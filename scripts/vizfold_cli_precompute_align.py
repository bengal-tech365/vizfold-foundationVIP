import argparse
import json
import logging
import os
import sys
import threading
import tempfile
from functools import partial
from multiprocessing import cpu_count
from shutil import copyfile

import openfold.data.mmcif_parsing as mmcif_parsing
from openfold.data.data_pipeline import AlignmentRunner
from openfold.data.parsers import parse_fasta
from openfold.data.tools import hhsearch, hmmsearch
from openfold.np import protein, residue_constants


logging.basicConfig(level=logging.WARNING)

def validate_args(args):
    errors = []

    if(not os.path.isdir(args.input_dir)):
        errors.append(f"--input_dir does not exist: {args.input_dir}")

    if(args.checkpoint is not None and not os.path.isfile(args.checkpoint)):
        errors.append(f"--checkpoint not found: {args.checkpoint}")

    if(not args.uniref90_database_path):
        errors.append("--uniref90_database_path is required")
    elif(not os.path.isfile(args.uniref90_database_path)):
        errors.append(f"--uniref90_database_path not found: {args.uniref90_database_path}")

    if(args.hmmsearch_binary_path and not args.pdb_seqres_database_path):
        errors.append("--hmmsearch_binary_path requires --pdb_seqres_database_path")

    if(args.hhsearch_binary_path and not args.pdb70_database_path):
        errors.append("--hhsearch_binary_path requires --pdb70_database_path")

    if(args.no_tasks < 1):
        errors.append(f"--no_tasks must be >= 1, got: {args.no_tasks}")

    if(errors):
        for e in errors:
            logging.warning(e)
        sys.exit(1)


def save_run_manifest(args, intermediates_dir):
    manifest = {
        "input_dir": args.input_dir,
        "output_dir": args.output_dir,
        "checkpoint": args.checkpoint,
        "capture_intermediates": args.capture_intermediates,
        "intermediates_layers": args.intermediates_layers if args.capture_intermediates else [],
        "intermediates_dir": intermediates_dir,
        "no_tasks": args.no_tasks,
        "cpus_per_task": args.cpus_per_task,
        "slurm_nodes": os.environ.get("SLURM_JOB_NUM_NODES", "N/A"),
        "slurm_node_id": os.environ.get("SLURM_NODEID", "N/A"),
    }
    manifest_path = os.path.join(args.output_dir, "vizfold_run_manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)


def run_seq_group_alignments(seq_groups, alignment_runner, args):
    dirs = set(os.listdir(args.output_dir))
    for seq, names in seq_groups:
        first_name = names[0]
        alignment_dir = os.path.join(args.output_dir, first_name)

        try:
            os.makedirs(alignment_dir)
        except Exception as e:
            logging.warning(f"Failed to create directory for {first_name} with exception {e}...")
            continue

        fd, fasta_path = tempfile.mkstemp(suffix=".fasta")
        with os.fdopen(fd, 'w') as fp:
            fp.write(f'>query\n{seq}')

        try:
            alignment_runner.run(fasta_path, alignment_dir)
        except Exception as e:
            logging.warning(e)
            logging.warning(f"Failed to run alignments for {first_name}. Skipping...")
            os.remove(fasta_path)
            os.rmdir(alignment_dir)
            continue

        os.remove(fasta_path)

        for name in names[1:]:
            if(name in dirs):
                logging.warning(f'{name} has already been processed. Skipping...')
                continue
            cp_dir = os.path.join(args.output_dir, name)
            os.makedirs(cp_dir, exist_ok=True)
            for f in os.listdir(alignment_dir):
                copyfile(os.path.join(alignment_dir, f), os.path.join(cp_dir, f))


def parse_and_align(files, alignment_runner, args):
    for f in files:
        path = os.path.join(args.input_dir, f)
        file_id = os.path.splitext(f)[0]
        seq_group_dict = {}

        if(f.endswith('.cif')):
            with open(path, 'r') as fp:
                mmcif_str = fp.read()
            mmcif = mmcif_parsing.parse(file_id=file_id, mmcif_string=mmcif_str)
            if(mmcif.mmcif_object is None):
                logging.warning(f'Failed to parse {f}...')
                if(args.raise_errors):
                    raise list(mmcif.errors.values())[0]
                else:
                    continue
            mmcif = mmcif.mmcif_object
            for chain_letter, seq in mmcif.chain_to_seqres.items():
                chain_id = '_'.join([file_id, chain_letter])
                l = seq_group_dict.setdefault(seq, [])
                l.append(chain_id)
        elif(f.endswith('.fasta') or f.endswith('.fa')):
            with open(path, 'r') as fp:
                fasta_str = fp.read()
            input_seqs, _ = parse_fasta(fasta_str)
            if(len(input_seqs) != 1):
                msg = f'More than one input_sequence found in {f}'
                if(args.raise_errors):
                    raise ValueError(msg)
                else:
                    logging.warning(msg)
            input_sequence = input_seqs[0]
            seq_group_dict[input_sequence] = [file_id]
        elif(f.endswith('.core')):
            with open(path, 'r') as fp:
                core_str = fp.read()
            core_prot = protein.from_proteinnet_string(core_str)
            aatype = core_prot.aatype
            seq = ''.join([
                residue_constants.restypes_with_x[aatype[i]]
                for i in range(len(aatype))
            ])
            seq_group_dict[seq] = [file_id]
        else:
            continue

        seq_group_tuples = [(k, v) for k, v in seq_group_dict.items()]
        run_seq_group_alignments(seq_group_tuples, alignment_runner, args)


def main(args):
    validate_args(args)

    os.makedirs(args.output_dir, exist_ok=True)

    intermediates_dir = None
    if(args.capture_intermediates):
        run_name = os.path.basename(args.input_dir.rstrip("/")) or "run"
        intermediates_dir = os.path.join(args.output_dir, f"{run_name}_intermediates")
        os.makedirs(intermediates_dir, exist_ok=True)

    save_run_manifest(args, intermediates_dir)

    if(args.hmmsearch_binary_path is not None and args.pdb_seqres_database_path is not None):
        template_searcher = hmmsearch.Hmmsearch(
            binary_path=args.hmmsearch_binary_path,
            hmmbuild_binary_path=args.hmmbuild_binary_path,
            database_path=args.pdb_seqres_database_path,
        )
    elif(args.hhsearch_binary_path is not None and args.pdb70_database_path is not None):
        template_searcher = hhsearch.HHSearch(
            binary_path=args.hhsearch_binary_path,
            databases=[args.pdb70_database_path],
        )
    else:
        template_searcher = None

    alignment_runner = AlignmentRunner(
        jackhmmer_binary_path=args.jackhmmer_binary_path,
        hhblits_binary_path=args.hhblits_binary_path,
        uniref90_database_path=args.uniref90_database_path,
        mgnify_database_path=args.mgnify_database_path,
        bfd_database_path=args.bfd_database_path,
        uniref30_database_path=args.uniref30_database_path,
        uniclust30_database_path=args.uniclust30_database_path,
        uniprot_database_path=args.uniprot_database_path,
        template_searcher=template_searcher,
        use_small_bfd=args.bfd_database_path is None,
        no_cpus=args.cpus_per_task,
    )

    files = list(os.listdir(args.input_dir))

    if(args.mmcif_cache is not None):
        with open(args.mmcif_cache, "r") as fp:
            cache = json.load(fp)
    else:
        cache = None

    dirs = []
    if(cache is not None and args.filter):
        dirs = set(os.listdir(args.output_dir))
        def prot_is_done(f):
            prot_id = os.path.splitext(f)[0]
            if(prot_id in cache):
                chain_ids = cache[prot_id]["chain_ids"]
                for c in chain_ids:
                    full_name = prot_id + "_" + c
                    if(not full_name in dirs):
                        return False
            else:
                return False
            return True
        files = [f for f in files if not prot_is_done(f)]

    def split_up_arglist(arglist):
        if(os.environ.get("SLURM_JOB_NUM_NODES", 0)):
            num_nodes = int(os.environ["SLURM_JOB_NUM_NODES"])
            if(num_nodes > 1):
                node_id = int(os.environ["SLURM_NODEID"])
                logging.warning(f"Num nodes: {num_nodes}")
                logging.warning(f"Node ID: {node_id}")
                arglist = arglist[node_id::num_nodes]

        t_arglist = []
        for i in range(args.no_tasks):
            t_arglist.append(arglist[i::args.no_tasks])

        return t_arglist

    if(cache is not None and "seqs" in next(iter(cache.values()))):
        seq_group_dict = {}
        for f in files:
            prot_id = os.path.splitext(f)[0]
            if(prot_id in cache):
                prot_cache = cache[prot_id]
                chains_seqs = zip(prot_cache["chain_ids"], prot_cache["seqs"])
                for chain, seq in chains_seqs:
                    chain_name = prot_id + "_" + chain
                    if(chain_name not in dirs):
                        l = seq_group_dict.setdefault(seq, [])
                        l.append(chain_name)

        func = partial(run_seq_group_alignments, alignment_runner=alignment_runner, args=args)
        seq_groups = [(k, v) for k, v in seq_group_dict.items()]
        seq_groups = sorted(seq_groups, key=lambda x: len(x[1]))
        task_arglist = [[a] for a in split_up_arglist(seq_groups)]
    else:
        func = partial(parse_and_align, alignment_runner=alignment_runner, args=args)
        task_arglist = [[a] for a in split_up_arglist(files)]

    threads = []
    for i, task_args in enumerate(task_arglist):
        print(f"Started thread {i}...")
        t = threading.Thread(target=func, args=task_args)
        threads.append(t)
        t.start()

    for t in threads:
        t.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True,
        help="Path to directory containing mmCIF, FASTA and/or ProteinNet .core files")
    parser.add_argument("--output_dir", type=str, required=True,
        help="Directory in which to output alignments")
    parser.add_argument("--checkpoint", type=str, default=None,
        help="Path to pretrained OpenFold model checkpoint")
    parser.add_argument("--capture_intermediates", action="store_true", default=False,
        help="Save intermediate attention maps and hidden states for VizFold visualization")
    parser.add_argument("--intermediates_layers", type=int, nargs="+", default=None,
        help="Which layers to capture intermediates from. Example: --intermediates_layers 4 8 12")
    parser.add_argument("--uniref90_database_path", type=str, default=None)
    parser.add_argument("--mgnify_database_path", type=str, default=None)
    parser.add_argument("--bfd_database_path", type=str, default=None)
    parser.add_argument("--uniref30_database_path", type=str, default=None)
    parser.add_argument("--uniclust30_database_path", type=str, default=None)
    parser.add_argument("--uniprot_database_path", type=str, default=None)
    parser.add_argument("--hhsearch_binary_path", type=str, default=None)
    parser.add_argument("--pdb70_database_path", type=str, default=None)
    parser.add_argument("--hmmsearch_binary_path", type=str, default=None)
    parser.add_argument("--hmmbuild_binary_path", type=str, default=None)
    parser.add_argument("--pdb_seqres_database_path", type=str, default=None)
    parser.add_argument("--jackhmmer_binary_path", type=str, default="jackhmmer")
    parser.add_argument("--hhblits_binary_path", type=str, default="hhblits")
    parser.add_argument("--raise_errors", action="store_true", default=False,
        help="Whether to crash on parsing errors")
    parser.add_argument("--cpus_per_task", type=int, default=cpu_count(),
        help="Number of CPUs to use")
    parser.add_argument("--mmcif_cache", type=str, default=None,
        help="Path to mmCIF cache. Used to filter files to be parsed")
    parser.add_argument("--no_tasks", type=int, default=1)
    parser.add_argument("--filter", type=bool, default=True)

    args = parser.parse_args()

    main(args)