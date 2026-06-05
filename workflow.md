# Context workflow

`Context` trains a linear discriminant on the background, transfers the
learned weights to the target panel, and calculates q-values and PEPs with
pyIsoPEP. The discriminant is learned by one of two interchangeable engines.

For one (`nontarget`, `target`) pair, given a seed `K`, prefix `P`, and
engine `E ∈ {percolator, mprophet}`:

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

   E = mprophet:
     Operates directly on raw feature values. The training loop:

       a. Resolve seed coefficients from --seed-coefficients (built-in
          'encyclopedia', 'none', or a JSON path). Names not present in
          the input contribute 0; if everything collapses, the seed
          model is disabled.
       b. Resolve the training feature set from --input-profile:
            auto         — keep var_/main_var_-prefixed columns if any,
                           else fall back to 'encyclopedia'.
            pin          — keep all pin-derived feature columns.
            encyclopedia — drop Encyclopedia metadata columns (pepLength,
                           charge1..4, precursorMass, RTinMin, midTime,
                           numberOfMatchingPeaksAboveThreshold, primary,
                           TD).
       c. Outer loop: 50 iterations of random 2-fold splits over the
          target and decoy rows.
       d. Inner loop on the training fold: 10 semi-supervised iterations.
            i = 0: rank targets by the seed model (or the auto-selected
                   best single feature when the seed is empty); keep the
                   top 15% by score percentile.
            i = 1: rerank targets by the current LDA; keep targets with
                   Storey q < 0.02 (p-values from a Gaussian N(mu, sigma)
                   fit to the decoy scores, matching Encyclopedia; pi0
                   floored at 0.05).
            i ≥ 2: same scheme at q < 0.01.
            Decoys are capped at 10× the kept-positive count. The inner
            loop stops early when passing count stops growing.
       e. Score the held-out fold with the final inner-loop LDA. If the
          seed model beats it (more Storey-q<0.01 targets on the held-out
          fold), keep the seed model for this outer iteration.
       f. After 50 outer iterations, sort by held-out passing count, keep
          the top 25 models, and average their coefficients and bias to
          produce the final (w, b).

     The weights file written to <outdir>/weights/P.weights.txt is a
     two-column TSV with `feature` and `weight`, plus a final `__bias__`
     row carrying the LDA constant.

2. Combine weights into a single (w, b):
     Percolator: w = mean over CV bins, b = mean over CV bins.
     mprophet:   w, b as returned by the training routine.

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
