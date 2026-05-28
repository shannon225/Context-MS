# Context workflow

`Context` trains a linear discriminant on the background, transfers the
learned weights to the target panel, and calculates q-values and PEPs with
pyIsoPEP. The discriminant is learned by one of two interchangeable engines.

For one (`nontarget`, `target`) pair, given a seed `K`, prefix `P`, and
engine `E ∈ {percolator, pyprophet}`:

```
0. Pairing check + feature pruning.
   The target file must declare exactly the same columns, in the same
   order, as the nontarget file. Errors out otherwise.
   Feature columns with σ < 1e-6 on the nontarget file are then dropped
   from BOTH files before either engine sees them.

1. Train discriminant on the pruned nontarget file.

   E = percolator:
     Run percolator on the nontarget file. The weights file written
     to <outdir>/weights/P.weights.txt is Percolator's native --weights
     output: a tab-separated text file with 3 rows per CV bin
     (header, normalised weights, raw weights).

   E = pyprophet:
     Compute μ, σ on the nontarget feature matrix.
     Standardize Z = (X − μ) / σ, rename columns to pyProphet's
     transition_group_id / decoy / main_var_* / var_* convention, write
     a pyProphet TSV, run `pyprophet score --classifier LDA`. pyProphet
     averages internal CV folds; we read the resulting <prefix>_weights.csv
     (standardized weights).
     De-standardize to raw-feature space: w_raw[i] = w_std[i] / σ[i].
     Bias term is dropped as pyProphet does.
     Persist <outdir>/weights/P.weights.txt as a tab-separated text file
     with columns: feature, weight, nontarget_mu, nontarget_sigma. The
     μ and σ columns record the standardisation that was applied to the
     nontarget feature matrix before training, so the raw-space weights
     can be reproduced or audited.

2. Combine weights into a single (w, b):
     Percolator: w = mean over CV bins, b = mean over CV bins.
     pyProphet:  w = w_raw, b = 0.

3. Score the target file: score = w · x + b.

4. PSM-level pyIsoPEP: run pyIsoPEP's q2pep with --calc-q-from-fdr on the
   rescored table, passing the resolved label column. pyIsoPEP
     (i) computes FDR + q from TDC,
     (ii) filters to targets,
     (iii) fits the I-Spline q -> PEP,
     (iv) returns a target-only frame with estimated q-values and PEPs

5. Peptide level: dedup the rescored target+decoy table by sequence
   (highest score wins), then repeat step 4.

6. Write outputs to <outdir>/ (default names shown; each can be
   overridden via --psm-out / --peptide-out / --rescored-out /
   --weights-out):
   <outdir>/P.psm.target.txt
   <outdir>/P.peptide.target.txt
```
