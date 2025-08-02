# Competing Risks Analysis: Educational Study

## 📚 Project Overview

This project provides an educational framework for understanding competing risks analysis, specifically comparing different implementations of Fine-Gray models and examining the "robust SE paradox" discussed in statistical literature.

## 🎯 Learning Objectives

1. **Compare 4 approaches** to competing risks analysis
2. **Observe robust SE phenomena** - when robust standard errors are actually smaller
3. **Understand tie handling** differences between methods
4. **Experience true value-based learning** through simulation

## 🚀 Quick Start

```r
# Run the main educational analysis
source("educational_cr_analysis.R")

# The analysis will:
# 1. Calculate true values from large dataset (30,000 observations)
# 2. Perform single analysis with sample data (400 observations)  
# 3. Run parallel simulation (200 iterations)
# 4. Provide comprehensive educational interpretation
```

## 📊 Four Methods Compared

| Method | Description | Key Features |
|--------|-------------|--------------|
| `cmprsk::crr()` | Direct Fine-Gray implementation | Breslow tie handling, direct subdistribution hazard |
| `finegray + coxph (Standard SE)` | Counting process transformation | Model-based standard errors |
| `finegray + coxph (Robust SE, no cluster)` | Sandwich estimator | Robust variance without clustering |
| `finegray + coxph (Robust SE + cluster)` | Full robust approach | Subject-level clustering for correct inference |

## 🔍 Key Phenomena to Observe

### 1. **Robust SE Paradox**
- **Expectation**: Robust SE should be larger (more conservative)
- **Reality**: In weighted regression, robust SE can be smaller
- **Reason**: Subject-level score aggregation reduces noise

### 2. **Tie Handling Impact**
- **Breslow** (CRR): Simple but less accurate with many ties
- **Efron** (Fine-Gray): More complex but handles ties better

### 3. **Clustering Importance**
- Fine-Gray creates multiple rows per subject
- `cluster(id)` essential for correct standard errors
- Without clustering: SE underestimation

## 📈 Performance Metrics

The analysis provides comprehensive performance evaluation:

- **Bias**: Distance from true values
- **RMSE**: Overall accuracy
- **Coverage**: 95% confidence interval performance  
- **SE Ratios**: Relative standard error comparisons

## ⚙️ System Requirements

- **R packages**: `survival`, `cmprsk`, `dplyr`, `ggplot2`, `tidyr`, `data.table`, `knitr`, `parallel`
- **Computing**: Automatically uses (total cores - 2) for parallel processing
- **Memory**: Handles large datasets (30,000+ observations) efficiently

## 📁 Project Structure

```
competing_risk_package_compare/
├── educational_cr_analysis.R    # Main educational analysis
├── perplexity_conv.md          # Reference discussion about SE phenomena
├── code_explanation.txt        # Legacy documentation
├── README.md                   # This file
└── archive/                    # Archived legacy files
    ├── survey_weight_demo.R    # IPTW demonstration
    ├── legacy_competing_comp.R # Original comparison code
    ├── legacy_package_comp.R   # Package comparison code
    ├── legacy_se_comparison.R  # SE comparison code
    └── legacy_se_estimate.R    # SE estimation code
```

## 🎓 Educational Value

This code is designed for learning and demonstrates:

1. **Methodological differences** between competing risks approaches
2. **Statistical phenomena** that contradict common intuitions
3. **Proper simulation practices** with true value benchmarks
4. **Parallel computing** for efficient statistical computation
5. **Comprehensive evaluation** of method performance

## 🔬 Based on Statistical Literature

The analysis incorporates insights from recent statistical literature on:
- Robust standard errors in weighted regression
- Fine-Gray model implementations
- Tie handling in survival analysis
- Clustering in counting process data

## 📞 Usage Notes

- Run time: ~5-10 minutes depending on system
- Output: Comprehensive tables and educational interpretation
- Parallelization: Automatically optimized for your system
- Results: Reproducible with set.seed() for consistency