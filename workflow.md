# Context workflow

`Context` trains a linear discriminant on the background, transfers the
learned weights to the reference panel, and reports q-values and PEPs. The
discriminant is learned by one of two interchangeable engines, and each
engine comes with its own confidence-estimation back-end.

For one (`background`, `reference`) pair, given a seed `K`, prefix `P`, and
engine `E ∈ {percolator, mprophet}`:

```
0. Pairing check + feature pruning.
   The reference file must declare exactly the same columns, in the same
   order, as the background file. Errors out otherwise.
   Feature columns with σ < 1e-6 on the background file are then dropped
   from BOTH files before either engine sees them.

1. Train discriminant on the pruned background file.

   E = percolator:
     Run percolator on the background file. The weights file written
     to <outdir>/weights/P.weights.txt is Percolator's native --weights
     output: a tab-separated text file with 3 rows per CV bin
     (header, normalised weights, raw weights).

   E = mprophet:
     A Python re-implementation of EncyclopeDIA's PRM rescoring and statistical validation components.
     Operates directly on raw feature values. The training loop:

       a. Resolve seed coefficients from --seed-coefficients (built-in
          'encyclopedia', 'none', or a JSON path). Names not present in
          the input contribute 0; if everything collapses, the seed
          model is disabled.
       b. Resolve the training feature set from --input-profile:
            auto         — keep var_/main_var_-prefixed columns if any,
                           else fall back to 'encyclopedia'.
            pin          — keep all pin-derived feature columns.
            encyclopedia — drop EncyclopeDIA metadata columns (pepLength,
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
                   q < 0.02.
            i ≥ 2: same scheme at q < 0.01.
            Decoys are capped at 10× the kept-positive count. The inner
            loop stops early when passing count stops growing.
       e. Score the held-out fold with the final inner-loop LDA. If the
          seed model beats it (more q<0.01 references on the held-out
          fold), keep the seed model for this outer iteration.
       f. After 50 outer iterations, sort by held-out passing count, keep
          the top 25 models, and average their coefficients and bias to
          produce the final (w, b).

     The q-values are calculated using EncyclopeDIA's method exactly:
     fit a Gaussian N(µ_d, σ_d) to the decoy scores, compute complementary-CDF
     p-values for the targets, and take Benjamini-Hochberg-adjusted p-values as
     q-values.

     The weights file written to <outdir>/weights/P.weights.txt is a
     two-column TSV with `feature` and `weight`, plus a final `__bias__`
     row carrying the LDA constant.

2. Combine weights into a single (w, b):
     Percolator: w = mean over CV bins, b = mean over CV bins.
     mprophet:   w, b as returned by the training routine.

3. Score the reference file: score = w · x + b.

4. Confidence estimation on the target+decoy scored table.

   E = percolator (pyIsoPEP pathway):
     Run pyIsoPEP's q2pep with --calc-q-from-fdr on the rescored PSM
     table, passing the resolved label column. pyIsoPEP
       (i) computes FDR + q from TDC,
       (ii) filters to targets,
       (iii) fits the I-Spline q -> PEP,
       (iv) returns a target-only frame with estimated q-values and PEPs.

   E = mprophet (EncyclopeDIA's PRM rescoring and statistical validation components):
     Match EncyclopeDIA's method exactly:
       (i)   fit N(µ_d, σ_d) to the decoy scores;
       (ii)  p = complementary CDF of the reference scores under that null;
       (iii) q = Benjamini-Hochberg-adjusted p-values;
       (iv)  PEP = Storey qvalue::lfdr — KDE on the logit-transformed
             p-values (Silverman bandwidth), monotone non-decreasing in
             p, with π0 floored at 0.05.
     pyIsoPEP is not used on this pathway.

5. Peptide level: dedup the rescored target+decoy table by sequence
   (highest score wins), then repeat step 4 (percolator: pyIsoPEP,
   mprophet: EncyclopeDIA PRM).

6. Write outputs to <outdir>/ (default names shown; each can be
   overridden via --psm-out / --peptide-out / --rescored-out /
   --weights-out):
   <outdir>/P.psm.reference.txt
   <outdir>/P.peptide.reference.txt
```
