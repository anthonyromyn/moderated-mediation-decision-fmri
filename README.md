# Statistical Modeling of Decision-Making & fMRI

A quick repo for my trial-level decision-making modeling with fMRI BOLD signals. The pipeline models choices with hierarchical mixed-effects and explores a moderated-mediation structure via SEM. This work was completed during my time as a PhD student and includes a partially finished paper and accompanying R code.

---

- Trial-level behavior (**`Decision`**) is modeled from expected values (EVs) and neural predictors (e.g., dlPFC / vMPFC), with participants (`ptp`) and runs (`Run`) as grouping factors.
- A second model predicts **`NAcc`** from the same predictors.
- A **lavaan** SEM then evaluates (moderated) mediation paths using centered variables and interaction terms.
- The paper (PDF) documents the experimental design and rationale that accompany the code.

---

## Repository Contents

```
.
├── analysis.R
├── write_up.pdf
└── README.md
```

---

## Model Outline (informal)

- **Behavioral model (GLMM, binomial)**  
  `Decision ~ (selfEV + otherEV) × dlPFC_trial × (vMPFC_self + vMPFC_other) + (1 + slopes | ptp/Run)`

- **Mediator model (LMM, gaussian)**  
  `NAcc ~ (selfEV + otherEV) × dlPFC_trial × (vMPFC_self + vMPFC_other) + (1 + slopes | ptp/Run)`

- **SEM (lavaan)**  
  Moderated mediation with mean-centered main effects and **explicit product terms** (no latent variables). Bootstrap used for inference on indirect effects.

---

## Citing / Acknowledgements

If you use this code or text, please cite the accompanying paper and this repository.
Authored by **Anthony Romyn** & **Dr. Wil Cunningham**.

---

## Data & Paths

**Data are not included in the repo.** 

---

