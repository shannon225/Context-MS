import argparse
import sys
from pathlib import Path
from . import pipeline


def _add_input_args(p):
    g = p.add_argument_group("input")
    g.add_argument("--nontarget", type=Path, required=True,
                   help="Nontarget (background) feature TSV used to train Percolator.")
    g.add_argument("--target", type=Path, required=True,
                   help="Target panel feature TSV to score and calibrate.")
    g.add_argument("--prefix", required=True,
                   help="Output prefix for results files.")


def _add_runtime_args(p):
    g = p.add_argument_group("runtime")
    g.add_argument("--outdir", type=Path, default=Path("results"))
    g.add_argument("--container-cmd", default="podman",
                   help="podman or docker; container runtime used as a fallback when percolator or pyIsoPEP isn't available locally.")


def _resolve_inputs(args):
    if not args.nontarget.is_file():
        sys.exit(f"Not found: {args.nontarget}")
    if not args.target.is_file():
        sys.exit(f"Not found: {args.target}")
    return args.prefix, args.nontarget, args.target


def cmd_run(args):
    prefix, nt, tg = _resolve_inputs(args)
    pipeline.run(
        nt, tg,
        outdir=args.outdir, prefix=prefix, seed=args.seed,
        container_cmd=args.container_cmd,
    )
    return 0


def build_parser():
    p = argparse.ArgumentParser(
        prog="context",
        description="Context-aware confidence estimation for targeted proteomics",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="end-to-end run")
    _add_input_args(p_run)
    _add_runtime_args(p_run)
    p_run.add_argument("--seed", type=int, default=1,
                       help="default 1")
    p_run.set_defaults(func=cmd_run)

    return p


def main(argv):
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


def cli():
    sys.exit(main(sys.argv[1:]))


if __name__ == "__main__":
    cli()
