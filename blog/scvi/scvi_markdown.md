# scVI Probabilistic Model
## Single-cell Variational Inference - A Deep Learning Framework for Gene Expression

---

## Model Overview (https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html)

scVI (single-cell Variational Inference) is a probabilistic model that uses deep learning to analyze gene expression data from individual cells. It handles technical noise and batch effects while learning meaningful biological patterns.

### Real-World Analogy

Think of scVI like a sophisticated photo editor that:
- Removes camera-specific artifacts (batch effects)
- Handles missing pixels (dropout/zero-inflation)
- Compresses images while keeping important features (dimensionality reduction)
- Reconstructs the original with improvements

---

## Generative Process: How Data is Created

The model assumes each gene expression value is generated through this hierarchical process:

### Step 1: Latent Space
```
z_n ~ Normal(0, I)
```
**Low-dimensional representation of cell state**

↓

### Step 2: Library Size
```
ℓ_n ~ LogNormal(ℓ_mean,s, ℓ_var,s)
```
**Scaling factor correlating with sequencing depth**

↓

### Step 3: Neural Network Transformations
```
w_ng = f_w(z_n, s_n)^g
h_ng = f_h(z_n, s_n)^g
```
**Map latent space to gene expression patterns**

↓

### Step 4: Generate Expression
```
x_ng ~ ZINB(ℓ_n · w_ng, θ_g, h_ng)
```
**Observed gene expression count**

---

## Numerical Example

### For Cell n=1, Gene g=5 (e.g., CD4):

| Variable | Value | Interpretation |
|----------|-------|----------------|
| z₁ | [0.5, -0.3, 0.8, ..., 0.1] | 10-dimensional latent representation |
| ℓ₁ | 8.7 (log scale) | ≈ 6,000 total transcripts |
| w₁₅ | 0.002 | CD4 is 0.2% of total transcripts |
| h₁₅ | 0.1 | 10% probability of technical dropout |
| θ₅ | 5.0 | Inverse dispersion for CD4 gene |
| **x₁₅** | **12** | **Observed: 12 CD4 transcripts** |

**ZINB Mean:** ℓ₁ · w₁₅ = 6,000 × 0.002 = **12 expected counts**

---

## Neural Network Architecture

### Encoder Networks (Inference)

Map observed data to latent representations:

```
Input (x_n, s_n) → Hidden (128-256 nodes) → Hidden (128-256 nodes) → Output (μ, σ²)
```

**Encoder Output Example:**
- **Input:** Gene expression vector (20,000 genes) + batch ID
- **Output:** q(z|x) with mean = [0.5, -0.3, ...] and variance = [0.1, 0.2, ...]
- **Output:** q(ℓ|x) with log-mean = 8.7, log-variance = 0.5

### Decoder Networks (Generation)

Map latent space back to gene expression:

#### f_w: Expression Frequencies
- **Architecture:** 1-3 hidden layers, softmax output
- **Output:** Vector summing to 1 across all genes
- **Example:** [0.002, 0.015, 0.001, ..., 0.003]

#### f_h: Dropout Probabilities
- **Architecture:** 1-3 hidden layers, sigmoid output
- **Output:** Probability per gene
- **Example:** [0.1, 0.05, 0.2, ..., 0.15]

### Key Constraints:
```
• Sum over genes: Σ_g w_ng = 1 (via softmax)
• Dropout: 0 ≤ h_ng ≤ 1 (via sigmoid)
• Activation: ReLU between hidden layers
• Regularization: Dropout + Batch Normalization
```

---

## Variational Inference

Since exact posterior is intractable, scVI uses variational inference to approximate it.

### Objective Function (ELBO)

```
ELBO = 𝔼_q[log p(x|z,ℓ)] - KL(q(z|x) || p(z)) - KL(q(ℓ|x) || p(ℓ))
```

### Components:

#### Reconstruction Term
```
𝔼_q[log p(x|z,ℓ)]
```
- How well the model reconstructs the observed data
- **Example:** -1250.3 (higher is better)

#### Regularization Terms
```
KL(q || p)
```
- Keep latent distributions close to priors
- **Example:** KL_z = 5.2, KL_ℓ = 1.8

---

## Optimization Strategy

### Mini-batch Stochastic Optimization:

| Parameter | Value | Purpose |
|-----------|-------|---------|
| Batch Size | 128 cells | Memory efficiency, GPU optimization |
| Optimizer | Adam | Adaptive learning rate |
| Epochs | 120-250 | Until convergence |
| Learning Rate | Dataset-specific | Grid search on validation set |

### Reparameterization Trick

**Instead of:** `z ~ 𝒩(μ, σ²)` [non-differentiable]

**Use:** `z = μ + σ · ε, where ε ~ 𝒩(0, 1)` [differentiable!]

**Why This Matters:** Allows gradients to flow through random sampling, making end-to-end training possible with standard backpropagation.

---

## Zero-Inflated Negative Binomial (ZINB)

Handles both biological zeros and technical dropout:

```
P(x = 0) = h + (1-h)·NB(0|μ, θ)
P(x = k) = (1-h)·NB(k|μ, θ), for k > 0
```

### Numerical Example: Gene Expression Distribution

**Parameters:** μ = 12, θ = 5.0, h = 0.1

| Count (x) | Probability | Interpretation |
|-----------|-------------|----------------|
| 0 | 0.15 | Technical dropout (10%) + biological zeros |
| 10-14 | 0.45 | Most likely range (around mean) |
| 20+ | 0.08 | High expression events (heavy tail) |

### ZINB Parameters:

#### Mean (μ)
- Expected expression level
- μ = ℓ_n · w_ng
- Example: 6,000 × 0.002 = 12

#### Dispersion (θ)
- Controls variance and tail heaviness
- Higher θ = less overdispersion
- Gene-specific, learned during training

---

## Key Advantages & Features

### Scalability
- Mini-batch training (128 cells)
- No need to load entire dataset
- GPU acceleration
- Handles millions of cells

### Batch Correction
- Batch annotation (s_n) as input
- Conditional independence in latent space
- No explicit MMD penalty needed
- Preserves biological signal

### Flexibility
- Beyond generalized linear models
- Neural networks learn complex patterns
- 3 hyperparameters to tune
- Automatic hyperparameter selection

### Biological Accuracy
- ZINB handles technical dropout
- Gene-specific dispersion
- Library size normalization
- Preserves rare cell types

---

## Hyperparameters

| Parameter | Typical Range | Selection Method |
|-----------|---------------|------------------|
| Learning Rate | 1e-4 to 1e-3 | Grid search on validation set |
| Number of Layers | 1-3 | Maximize held-out log likelihood |
| Layer Width | 128 or 256 | Balance capacity vs. overfitting |
| Latent Dimensions | 10-20 | Fixed (typically 10) |
| Dropout Rate | 0.1-0.3 | Regularization (prevent overfitting) |

**Automatic Selection:** Hyperparameters are chosen via grid search that maximizes held-out log likelihood on a validation set—a standard practice for training deep generative models.

---

## Complete Training Workflow

### 1. Data Preparation
- **Input:** X (cells × genes), batch annotations
- **Preprocessing:** Log-normalize, compute library sizes

↓

### 2. Initialize Networks
- Encoder & decoder networks with random weights
- Set priors for z and ℓ based on data statistics

↓

### 3. Mini-batch Sampling
- Randomly sample M=128 cells from dataset
- Load corresponding gene expression + batch info

↓

### 4. Forward Pass
- **Encoder:** (x,s) → q(z|x), q(ℓ|x)
- **Sample:** z, ℓ via reparameterization
- **Decoder:** (z,s) → w, h → reconstruct x

↓

### 5. Compute ELBO
- Reconstruction + KL divergences
- ELBO = log p(x|z,ℓ) - KL(z) - KL(ℓ)

↓

### 6. Backpropagation
- Compute gradients via automatic differentiation
- Update weights using Adam optimizer

↓

### 7. Convergence
- Repeat steps 3-6 for 120-250 epochs
- Monitor validation loss for early stopping

---

## Training Time Example

**Dataset:** 100,000 cells × 20,000 genes

**Epochs:** 150

**Batch Size:** 128 cells

**Iterations per epoch:** 100,000 / 128 ≈ 780

**Total iterations:** 150 × 780 = 117,000

**Training time:** ~2-4 hours on GPU

---

## Key Mathematical Properties

### Mean-Field Variational Distribution
```
q(z, ℓ | x) = q(z | x) · q(ℓ | x)
```
- Assumes independence between latent variables
- Simplifies optimization
- Computationally efficient

### Conditional Independence
```
p(x | z, ℓ, s) = ∏_g p(x_g | z, ℓ, s)
```
- Genes are conditionally independent given latent variables
- Allows parallel computation
- Scalable to high-dimensional data

### Integration of Technical Variables
- **w_ng, h_ng, y_ng** are integrated out analytically
- Results in closed-form ZINB distribution
- Reduces computational burden

---

## Variable Glossary (https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html)

| Symbol | Meaning | Dimensions |
|--------|---------|------------|
| x_ng | Observed gene expression count | scalar |
| z_n | Latent representation of cell | d-dimensional vector (typically d=10) |
| ℓ_n | Library size (log scale) | scalar |
| s_n | Batch annotation | categorical |
| w_ng | Expected gene expression frequency | scalar ∈ [0,1] |
| h_ng | Dropout probability | scalar ∈ [0,1] |
| θ_g | Gene-specific inverse dispersion | scalar > 0 |
| f_w | Decoder network for frequencies | neural network |
| f_h | Decoder network for dropout | neural network |
| N | Number of cells | integer |
| G | Number of genes | integer (typically 20,000) |
| B | Number of batches | integer |

---

## The Why's

### Why Neural Networks?
Traditional generalized linear models (GLMs) assume a specific functional form relating covariates to outcomes. Neural networks allow scVI to:
- Learn arbitrary non-linear relationships
- Capture complex gene-gene interactions
- Adapt to different data distributions
- Scale to high-dimensional spaces

### Conditional Independence vs. MMD Penalty
- scVI relies on conditional independence: p(z | x, s) does not depend on s
- Alternative: Maximum Mean Discrepancy (MMD) penalty explicitly enforces batch mixing
- scVI's approach is simpler and often sufficient
- MMD can be added if stronger batch correction is needed

### Why ZINB over Other Distributions?
- **Poisson:** Cannot model overdispersion (variance > mean)
- **Negative Binomial:** Cannot model excess zeros from technical dropout
- **ZINB:** Handles both overdispersion AND zero-inflation
- Biologically motivated: separates technical zeros from biological zeros

---

## Practical Considerations

### When to Use scVI:
- Large-scale single-cell RNA-seq datasets
- Multiple batches/experiments to integrate
- Need for batch correction and dimensionality reduction
- Downstream tasks: clustering, differential expression, imputation

### Computational Requirements:
- **GPU:** Highly recommended for datasets > 10,000 cells
- **Memory:** Proportional to batch size, not dataset size
- **Storage:** Only mini-batches loaded at a time

### Model Selection:
- Use validation set to choose hyperparameters
- Monitor convergence via ELBO on validation data
- Early stopping prevents overfitting (12 epochs without improvement)

---

## Summary

scVI combines:
1. **Probabilistic modeling** (ZINB distribution)
2. **Deep learning** (neural network encoders/decoders)
3. **Variational inference** (approximate posterior)
4. **Stochastic optimization** (mini-batch training)

This creates a powerful, scalable framework for single-cell gene expression analysis that handles technical noise while preserving biological signal.

**Key Innovation:** Using neural networks within a probabilistic framework allows flexible modeling while maintaining interpretability through the generative process.

---