# A Numerical Walkthrough of `scanpy.pp.normalize_total()`

This post gives a concrete, numbers-first explanation of how Scanpy’s `pp.normalize_total()` rescales UMI counts across cells.

---

## What `normalize_total` does

For each cell \(i\), Scanpy computes a per-cell total \(T_i\) and rescales every gene count \(x_{gi}\) by a factor that makes totals comparable:

\[
x'_{gi} = x_{gi} \times \frac{\text{target\_sum}}{T_i}.
\]

- If `target_sum=None`, Scanpy uses the **median** of \(T_i\) across cells as the target.
- If `exclude_highly_expressed=True`, the totals \(T_i\) are computed **excluding** any “highly expressed” genes (defined by `max_fraction`).

---

## Example dataset (3 cells × 6 genes)

We’ll use one dominant gene **H** and five smaller genes **S1–S5**.

| Cell | H  | S1 | S2 | S3 | S4 | S5 | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| C1 | 60 | 3 | 3 | 3 | 3 | 3 | 75 |
| C2 |120 | 6 | 6 | 6 | 6 | 6 |150 |
| C3 | 30 | 1 | 1 | 1 | 1 | 1 | 35 |

Per-cell totals are \([75, 150, 35]\), so the **median is 75**.  
Therefore with `target_sum=None`, we use:

\[
\text{target\_sum} = 75.
\]

---

## Case A: `exclude_highly_expressed=False` (default)

Here the per-cell total is the full total \(T_i\).

Scaling factor:

\[
s_i = \frac{75}{T_i}.
\]

- C1: \(s_1 = 75/75 = 1\)
- C2: \(s_2 = 75/150 = 0.5\)
- C3: \(s_3 = 75/35 \approx 2.142857\)

Now rescale counts: \(x' = x \times s_i\).

### Example outputs
**Gene H**
- C1: \(60 \times 1 = 60\)
- C2: \(120 \times 0.5 = 60\)
- C3: \(30 \times 2.142857 \approx 64.2857\)

**Gene S1**
- C1: \(3 \times 1 = 3\)
- C2: \(6 \times 0.5 = 3\)
- C3: \(1 \times 2.142857 \approx 2.142857\)

✅ After this normalization, each cell’s **total** is ~75 (up to rounding), meaning differences in sequencing depth are reduced.

---

## Case B: `exclude_highly_expressed=True`, `max_fraction=0.05`

### Step 1 — Identify “highly expressed” genes

A gene is flagged if it exceeds `max_fraction` of a cell’s total in **any** cell.

For gene **H** in C1:

\[
\frac{60}{75} = 0.80 = 80\% > 5\%,
\]

so **H is considered highly expressed** and is excluded from the total used to compute normalization factors.

### Step 2 — Compute totals excluding H

Compute totals using only **S1–S5**:

- C1: \(3+3+3+3+3 = 15\)
- C2: \(6+6+6+6+6 = 30\)
- C3: \(1+1+1+1+1 = 5\)

Now compute scaling factors:

\[
s_i = \frac{75}{\text{(sum of non-excluded genes)}_i}.
\]

- C1: \(s_1 = 75/15 = 5\)
- C2: \(s_2 = 75/30 = 2.5\)
- C3: \(s_3 = 75/5 = 15\)

### Step 3 — Apply scaling to *all* genes (including H)

**Small genes become comparable across cells**

For S1:
- C1: \(3 \times 5 = 15\)
- C2: \(6 \times 2.5 = 15\)
- C3: \(1 \times 15 = 15\)

So S1–S5 now contribute a consistent total:

\[
S1+S2+S3+S4+S5 = 75 \quad \text{in every cell.}
\]

**But H becomes large**

- C1: \(60 \times 5 = 300\)
- C2: \(120 \times 2.5 = 300\)
- C3: \(30 \times 15 = 450\)

✅ This behavior is **intentional**: excluding a dominant gene prevents it from “soaking up” the normalization and pushing all other genes artificially low.

---

## Summary (what to remember)

- With `exclude_highly_expressed=False`, the **entire library size** drives scaling, and each cell’s total becomes `target_sum`.
- With `exclude_highly_expressed=True`, scaling is driven by the **non-dominant genes**, so they become comparable across cells—even if one gene dominates counts.
- After normalization, people typically apply `sc.pp.log1p(adata)` to stabilize variance and compress large values.

---

## (Optional) Code snippet (Scanpy)

```python
import scanpy as sc

# Default behavior: rescale each cell to median total (if target_sum=None)
sc.pp.normalize_total(adata, target_sum=None)

# Exclude highly expressed genes from the size-factor computation
sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=True, max_fraction=0.05)

# Common next step
sc.pp.log1p(adata)
```


---

## R equivalents

Below are simple R functions that reproduce the **core scaling logic** of `scanpy.pp.normalize_total()` for a raw counts matrix `counts` (genes × cells).

### 1) Default behavior (`target_sum = NULL`, `exclude_highly_expressed = FALSE`)

```r
# counts: genes x cells matrix (raw UMI counts)
normalize_total_r <- function(counts, target_sum = NULL) {
  stopifnot(is.matrix(counts) || inherits(counts, "Matrix"))

  totals <- colSums(counts)
  if (is.null(target_sum)) target_sum <- stats::median(totals)

  sf <- target_sum / totals                 # per-cell scaling factors
  norm <- t(t(counts) * sf)                 # scale each column (cell)

  list(X = norm, size_factors = sf, target_sum = target_sum)
}

# Example:
# out <- normalize_total_r(counts)
# norm_counts <- out$X
```

### 2) Excluding highly expressed genes (`exclude_highly_expressed = TRUE`, `max_fraction = 0.05`)

This matches Scanpy’s idea: identify genes that exceed `max_fraction` of a cell’s total **in any cell**, exclude them when computing per-cell totals (size factors), and then apply the resulting scaling to **all** genes.

```r
normalize_total_exclude_high_r <- function(counts, target_sum = NULL, max_fraction = 0.05) {
  stopifnot(is.matrix(counts) || inherits(counts, "Matrix"))

  totals_full <- colSums(counts)
  if (is.null(target_sum)) target_sum <- stats::median(totals_full)

  # Identify genes that are > max_fraction of a cell's total in ANY cell
  # (Compute fractions per cell; avoid division by zero)
  denom <- ifelse(totals_full == 0, 1, totals_full)
  frac <- t(t(counts) / denom)              # each entry: x_gi / total_i

  highly_expressed <- apply(frac, 1, function(v) any(v > max_fraction))
  counts_non_high <- counts[!highly_expressed, , drop = FALSE]

  totals_non_high <- colSums(counts_non_high)
  totals_non_high[totals_non_high == 0] <- 1  # guard against all-zero

  sf <- target_sum / totals_non_high
  norm <- t(t(counts) * sf)

  list(
    X = norm,
    size_factors = sf,
    target_sum = target_sum,
    highly_expressed_genes = rownames(counts)[highly_expressed]
  )
}

# Example:
# out2 <- normalize_total_exclude_high_r(counts, max_fraction = 0.05)
# norm_counts2 <- out2$X
# out2$highly_expressed_genes
```

### Optional next step: log1p (like `sc.pp.log1p`)

```r
log1p_transform <- function(norm_counts) {
  log1p(norm_counts)   # natural log(1 + x)
}

# logX <- log1p_transform(norm_counts)
```
