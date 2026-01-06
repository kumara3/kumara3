# Scanpy preprocessing and SCVI training in R (via reticulate)

This walkthrough shows how to:

- Configure **reticulate** to use a conda environment (`scvi-tools`)
- Load an example **AnnData** dataset from `scvi-tools`
- Inspect AnnData structure (`obs`, `var`, sparse `X`)
- Preprocess with Scanpy (filtering, normalization, log1p, HVGs)
- Train an **SCVI** model and plot the training curve (negative ELBO)
- Save the trained model to disk
- Reference: https://docs.scvi-tools.org/en/stable/tutorials/notebooks/r/api_overview_in_R.html; https://github.com/scverse/scanpy/issues/2578

---

## Setup: R + Python (conda) + imports

```r
library(reticulate)
library(anndata)
library(ggplot2)
# library(IRdisplay)  # optional

options(reticulate.conda_binary = "/home/miniforge3/condabin/conda")
conda_list()
use_condaenv("scvi-tools", required = TRUE)
py_config()

sc   <- import("scanpy", convert = FALSE)
scvi <- import("scvi",   convert = FALSE)
```

---

## Load example AnnData

```r
# Load a sample dataset from scvi-tools
adata <- scvi$data$heart_cell_atlas_subsampled()

class(adata)
# [1] "anndata._core.anndata.AnnData" "python.builtin.object"

obs <- adata$obs
var <- adata$var

class(adata$X)
# [1] "scipy.sparse._csr.csr_matrix"        "scipy.sparse._matrix.spmatrix"
# ...
# [11] "python.builtin.object"

adata$X$shape
```

---

## AnnData Summary

- **Dimensions:** `n_obs × n_vars = 18641 × 26662`

### `obs` (cell metadata)
- `NRP`
- `age_group`
- `cell_source`
- `cell_type`
- `donor`
- `gender`
- `n_counts`
- `n_genes`
- `percent_mito`
- `percent_ribo`
- `region`
- `sample`
- `scrublet_score`
- `source`
- `type`
- `version`
- `cell_states`
- `Used`

### `var` (gene / feature metadata)
- `gene_ids-Harvard-Nuclei`
- `feature_types-Harvard-Nuclei`
- `gene_ids-Sanger-Nuclei`
- `feature_types-Sanger-Nuclei`
- `gene_ids-Sanger-Cells`
- `feature_types-Sanger-Cells`
- `gene_ids-Sanger-CD45`
- `feature_types-Sanger-CD45`
- `n_counts`

### `uns` (unstructured annotations)
- `cell_type_colors`

---

## Preprocessing using Scanpy 

> **Storage convention used here**
>
> - Raw counts are preserved in: `adata$layers[["counts"]]`
> - Normalized + log1p values are stored in: `adata$X`
> - A snapshot of the current state is stored in: `adata$raw`

```r
# Filtering
sc$pp$filter_genes(adata, min_counts = 3L)
sc$pp$filter_cells(adata, min_genes  = 200L)

# Preserve raw counts
adata$layers["counts"] <- adata$X$copy()

# Normalize total counts per cell (target_sum = 1e4) into adata$X
sc$pp$normalize_total(adata, target_sum = 1e4)


# Log-transform normalized data in adata$X
sc$pp$log1p(adata)

# Store a snapshot for later
adata$raw <- adata
```

---

## Highly variable genes

```r
sc$pp$highly_variable_genes(
  adata,
  n_top_genes = r_to_py(2000),
  subset      = TRUE,
  layer       = "counts",
  flavor      = "seurat_v3",
  batch_key   = "cell_source"
)
```

---

## Create and train an SCVI model

```r
# Register AnnData with scvi-tools
scvi$model$SCVI$setup_anndata(
  adata,
  layer = "counts",
  categorical_covariate_keys = c("cell_source", "donor"),
  continuous_covariate_keys  = c("percent_mito", "percent_ribo")
)

# Initialize model
model <- scvi$model$SCVI(adata)
model

# Train model
model$train(max_epochs = 5L)
```

---

## Plot training curve (negative ELBO)

```r
elbo_train <- py_to_r(model$history[["elbo_train"]])
elbo_train$epoch <- as.numeric(row.names(elbo_train))
elbo_train$elbo_train <- as.numeric(unlist(elbo_train$elbo_train))

ggplot(elbo_train, aes(x = epoch, y = elbo_train)) +
  geom_line(color = "steelblue") +
  labs(title = "Negative ELBO over training epochs")
```

---

## Save trained model

```r
model_dir <- file.path(getwd(), "scvi_model")
model$save(model_dir, overwrite = TRUE)
```
