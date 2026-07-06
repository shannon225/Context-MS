import argparse
import sys
from pathlib import Path
from . import pipeline
from .io import INPUT_PROFILES


def _add_input_args(p):
    g = p.add_argument_group("input")
    g.add_argument("--background", type=Path, required=True,
                   help="Background feature TSV used to train the engine.")
    g.add_argument("--reference", type=Path, required=True,
                   help="Reference panel feature TSV to score and calibrate.")
    g.add_argument("--prefix", required=True,
                   help="Output prefix for results files.")


def _add_runtime_args(p):
    g = p.add_argument_group("runtime")
    g.add_argument("--outdir", type=Path, default=Path("results"),
                   help="Output directory for results files.")
    g.add_argument("--engine", choices=pipeline.ENGINES, default="percolator",
                   help="percolator or mprophet; engine used to train the "
                        "model on the background and score the reference "
                        "panel.")
    g.add_argument("--container-cmd", default="podman",
                   help="podman or docker; container runtime used as a "
                        "fallback when the engine or pyIsoPEP isn't available "
                        "locally.")
    g.add_argument("--psm-out", default=None,
                   help="File name (or path) for the PSM-level reference "
                        "output. If relative, written inside --outdir. "
                        "Default: <prefix>.psm.reference.txt")
    g.add_argument("--peptide-out", default=None,
                   help="File name (or path) for the peptide-level reference "
                        "output. If relative, written inside --outdir. "
                        "Default: <prefix>.peptide.reference.txt")
    g.add_argument("--rescored-out", default=None,
                   help="File name (or path) for the rescored-features TSV "
                        "(reference targets + decoys). If relative, written "
                        "inside --outdir. Default: "
                        "<prefix>.rescored_features.tsv")
    g.add_argument("--weights-out", default=None,
                   help="File name (or path) for the trained weights file. "
                        "If relative, written inside --outdir/weights. "
                        "Default: <prefix>.weights.txt")


def _add_mprophet_args(p):
    g = p.add_argument_group("mprophet engine")
    g.add_argument("--input-profile", choices=INPUT_PROFILES,
                   default=pipeline.DEFAULT_INPUT_PROFILE,
                   help="Feature-column selection profile for the mprophet "
                        "engine. 'auto' detects var_/main_var_ prefixed "
                        "inputs and falls back to 'encyclopedia' otherwise. "
                        "Has no effect on the percolator engine.")
    g.add_argument("--seed-coefficients",
                   default=pipeline.DEFAULT_SEED_COEFFICIENTS,
                   help="Seed-model coefficients for the mprophet engine. "
                        "Accepts a built-in name ('encyclopedia', 'none') or "
                        "a path to a JSON file mapping feature names to "
                        "coefficients. 'none' disables the seed model and "
                        "uses an auto-selected single-feature start for "
                        "inner iter 0.")


def _resolve_inputs(args):
    if not args.background.is_file():
        sys.exit(f"Not found: {args.background}")
    if not args.reference.is_file():
        sys.exit(f"Not found: {args.reference}")
    return args.prefix, args.background, args.reference


def cmd_run(args):
    prefix, bg, ref = _resolve_inputs(args)
    pipeline.run(
        bg, ref,
        outdir=args.outdir, prefix=prefix, seed=args.seed,
        engine=args.engine, container_cmd=args.container_cmd,
        psm_out=args.psm_out, peptide_out=args.peptide_out,
        rescored_out=args.rescored_out, weights_out=args.weights_out,
        input_profile=args.input_profile,
        seed_coefficients=args.seed_coefficients,
    )
    return 0


def build_parser():
    p = argparse.ArgumentParser(
        prog="context",
        description="Context-aware confidence estimation for targeted proteomics.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="end-to-end run")
    _add_input_args(p_run)
    _add_runtime_args(p_run)
    _add_mprophet_args(p_run)
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
