# Statistical Inference and Modeling in R

## Analyses

| Analysis | Methods | Inference|
| --- | --- | --- |
| Blood pressure | Multiple linear regression; logistic regression | Confidence intervals for population means; coefficient tests; conditional effects and odds |
| Ecological count data | Poisson and negative binomial GLMs; custom maximum likelihood estimation | Distributional assumptions; overdispersion; likelihood optimization; AIC comparison |
| Competing risks | Cumulative incidence functions; Fine–Gray and cause-specific Cox models | Distinguishing cumulative event probabilities, subdistribution hazards, and cause-specific hazards |
| Treatment comparison | One-way ANOVA; estimated marginal means; planned contrasts | Omnibus hypothesis testing; mean differences; multiplicity-adjusted simultaneous confidence intervals |

## Statistical Approach

- **Maximum likelihood:** The negative binomial analysis defines a negative log-likelihood explicitly and estimates regression and dispersion parameters using `optim()` with BFGS.
- **Model assumptions:** Poisson regression assumes equality of conditional mean and variance; negative binomial regression relaxes that restriction.
- **Competing events:** Cumulative incidence accounts for death as a competing event when estimating recovery probability. Fine–Gray and cause-specific Cox regressions answer different questions.
- **Simultaneous inference:** Four treatment-versus-control contrasts use a multivariate-t adjustment to account for multiple comparisons.

## Repository structure

```text
statistical-inference-and-modeling/
├── statistical_inference.R
├── README.md
├── data/
│   ├── data 1.txt
│   ├── data 2.txt
│   ├── data 3.txt
│   └── data 4.txt
└── figure/                 # Created when the script runs
    └── competing_risks_cif.png
```

## Run

Install the R dependencies once:

```r
install.packages(c("MASS", "cmprsk", "survival", "dplyr", "ggplot2", "emmeans"))
```

Set your working directory to the repository root (the folder containing `statistical_inference.R`) and run:

```r
source("statistical_inference.R", echo = TRUE)
```

The competing-risks plot is written to `figure/competing_risks_cif.png` rather than relying on the size of the RStudio Plots pane. Other analyses print model summaries and/or display plots when run interactively.

## Notes

The project is intended to showcase the statistical workflow and code, rather than a complete diagnostic or validation study.