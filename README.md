# Heston-IG-asian-options-and-Variance-reduction
# Heston IG Asian Options and Variance Reduction

Replication of two quantitative finance research papers as part of a 
Laboratory on Quantitative Finance course. The project implements and 
validates Monte Carlo simulation methods for option pricing under the 
Heston stochastic volatility model in MATLAB.

## Report

A full written report (`report.pdf`) is included in the repository. 
It covers the theoretical background of both papers, the mathematical 
derivations, the implementation details of all algorithms, and a 
comparison of the numerical results obtained against the original papers.

## Papers Replicated

- **Tse & Wan (2014)** - Low-bias simulation scheme for the Heston model
  using Inverse Gaussian (IG) approximation (IPZ-IG scheme). Note: the 
  IPZ interpolation component for sampling V(t2)|V(t1) was not replicated;
  only the IG approximation for the integrated variance is implemented.
- **Dingeç, Sak & Hörmann (2014)** - Variance reduction for arithmetic 
  Asian options via Control Variate (CVN) and Conditional Monte Carlo (CMC)
  under the Heston model.

## Code Structure
```
src/
├── tse_wan_2014/
│   ├── IG.m                  % Inverse Gaussian sampler
│   └── simHestonIG.m         % European/Asian option pricing via IG scheme
│
└── dingec_et_al_2014/
    ├── Glasserman_pathV.m    % Variance path simulation (Glasserman 2004)
    ├── pathsim.m             % Alternative variance path via IG scheme
    ├── charfun.m             % Joint characteristic function (Heston model)
    ├── g.m                   % Function g(w) for control variate expectation
    ├── conditionalexp.m      % Conditional expectation for CMC
    ├── Varreduction.m        % Asian option pricing with CVN+CMC
    ├── Newton.m              % Newton's method for root finding
    └── greek_computation.m   % Delta and Gamma via Pathwise Derivative
```

## Requirements

- MATLAB

## References

- Tse, S., Wan, J. (2014). *Low bias simulation scheme for Heston model 
  by Inverse Gaussian approximation*.
- Dingeç, K.D., Sak, H., Hörmann, W. (2014). *Variance Reduction for 
  Asian Options under a General Model Framework*.
- Glasserman, P. (2004). *Monte Carlo Methods in Financial Engineering*. 
  Springer.
- Quarteroni, A., Saleri, F., Gervasio, P. *Matematica Numerica*. Springer.
- Broadie, M., Kaya, O. (2006). Exact simulation of stochastic volatility 
  and other affine jump diffusion processes.
- Andersen, L. (2007). Efficient simulation of the Heston stochastic 
  volatility model.
