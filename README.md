# BasicPLSPM_Prime

### Author:
Heungsun Hwang and Gyeongcheol Cho

## Description:
- The **BasicPLSPM_Prime** package enables users to estimate and evaluate basic PLSPM models.

## Features:
- Estimate GSCA model parameters and calculate their standard errors (SE) along with 95% confidence intervals (CI).
- Enable parallel computing for bootstrap sampling.
- Allow users to specify a sign-fixing indicator for each component.

## Installation:
To use this package in MATLAB:
1. Clone or download the repository:
   ```bash
   git clone https://github.com/GyeongcheolCho/BasicPLSPM_Prime.git
   ```
2. Add the package to your MATLAB path:
   ```matlab
    addpath(genpath('BasicPLSPM_Prime'))
   ```

## Usage:
- For examples on how to use the package, refer to the `Run_Example_BasicPLSPM.m` file. This file demonstrates the implementation of `BasicPLSPM()` using the ACSI dataset.

## Compatibility:
- Tested on MATLAB R2023b.
- Likely compatible with earlier MATLAB versions.

### Citation (APA):
- If you use **BasicPLSPM_Prime** in your research or publications, please cite it in APA format as follows:

```plaintext
Hwang, H. & Cho, G. (2024). BasicPLSPM_Prime: A package for basic partial least squares path modeling [Computer software]. GitHub. https://github.com/GyeongcheolCho/BasicPLSPM_Prime
```
