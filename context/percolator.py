from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np

PERCOLATOR_IMAGE = "ghcr.io/percolator/percolator:master"


def run_percolator(input_tsv, weights_out, seed, *, container_cmd="podman",
                   image=PERCOLATOR_IMAGE, extra_args=()):
    input_tsv = Path(input_tsv).resolve()
    weights_out = Path(weights_out).resolve()
    weights_out.parent.mkdir(parents=True, exist_ok=True)

    if shutil.which("percolator"):
        cmd = _local_cmd(input_tsv, weights_out, seed, extra_args)
    else:
        cmd = _container_cmd(input_tsv, weights_out, seed, container_cmd,
                             image, extra_args)

    console = weights_out.with_suffix(".log")
    with console.open("w") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        sys.exit(f"Percolator failed (rc={rc}). See {console}")


def _local_cmd(input_tsv, weights_out, seed, extra_args):
    return [
        "percolator",
        "--seed", str(seed),
        "--weights", str(weights_out),
        "--protein-report-duplicates",
        "--post-processing-tdc",
        "--override",
        *extra_args,
        str(input_tsv),
    ]


def _container_cmd(input_tsv, weights_out, seed, container_cmd, image, extra_args):
    in_dir = input_tsv.parent
    out_dir = weights_out.parent
    return [
        container_cmd, "run", "--rm",
        "--entrypoint", "percolator",
        "-v", f"{in_dir}:/in:ro",
        "-v", f"{out_dir}:/out",
        image,
        "--seed", str(seed),
        "--weights", f"/out/{weights_out.name}",
        "--protein-report-duplicates",
        "--post-processing-tdc",
        "--override",
        *extra_args,
        f"/in/{input_tsv.name}",
    ]


def parse_weights(path):
    lines = [
        ln.rstrip("\n")
        for ln in Path(path).read_text().splitlines()
        if ln.strip() and not ln.lstrip().startswith("#")
    ]
    if len(lines) % 3 != 0:
        raise ValueError(f"{path}: expected 3 lines per CV bin, got {len(lines)}")
    k = len(lines) // 3
    feature_names = None
    w_raw = []
    b_raw = []
    for i in range(k):
        header = [c.strip() for c in lines[3 * i].split("\t")][:-1]
        if feature_names is None:
            feature_names = header
        elif header != feature_names:
            raise ValueError(f"{path}: feature-name drift between CV bins")
        raw = np.array([float(v) for v in lines[3 * i + 2].split("\t")])
        w_raw.append(raw[:-1])
        b_raw.append(float(raw[-1]))
    return feature_names, np.stack(w_raw, axis=0), np.array(b_raw)
