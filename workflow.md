# Context workflow

`Context` trains Percolator on the background, transfers the learned
weights to the target panel, and calculates q-values and PEPs with pyIsoPEP.


For one (`nontarget`, `target`) pair, given a seed `K` and prefix `P`:

```
0. Pairing check: the target file must declare exactly the same columns,
   in the same order, as the nontarget file. Errors out otherwise.

1. Run Percolator on the nontarget file.
   The weights file is written to results/weights/P.seedK.weights.txt
   (K=3 CV-bin raw-weight rows).

2. Parse weights.txt. Seed weight = mean across the K CV bins.

3. Score the target file: score = w . x + b.

4. PSM-level pyIsoPEP: invoke pyIsoPEP's q2pep with --calc-q-from-fdr on the
   rescored table, passing the resolved label column. pyIsoPEP
     (i) computes FDR + q from TDC,
     (ii) filters to targets,
     (iii) fits the I-Spline q -> PEP,
     (iv) returns a target-only frame with estimated q-values and PEPs

5. Peptide level: dedup the rescored target+decoy table by sequence
   (highest score wins), then repeat step 4.

6. Write output:
   results/P.seedK.psm.target.txt
   results/P.seedK.peptide.target.txt
```
