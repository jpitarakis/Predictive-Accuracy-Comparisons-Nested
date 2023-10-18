This repository contains the codes for using the predictive accuracy comparison tests developed
in Pitarakis, J. (2023). "A novel approach to predictive accuracy testing in nested environments", Econometric Theory, 2023. 

Funding from the ESRC through grant ES/W000989/1 is gratefully acknowledged.

The folder "Full Packet of Files" aggregates all files used in the project. 

If you wish to use these tests on your own data you only need the following files to be installed in your working directory: recursive_hstep_fast.m, recursive_hstep_interceptonly.m, recursive_hstep_slow.m, Nested_Stats_S0.m, Nested_Stats_Sbar.m together with your data file, say dummy_data_2.xlsx

**Implementation Example**: 
The demo data file contains three time-series  $y_{t}$, $x_{1t}$ and $x_{2t}$. Suppose you wish to implement predictive accuracy comparisons between model 1: $y_{t}=\theta_{0}+\theta_{1}x_{1,t-1}+u_{t}$ and model 2: $y_{t}=\theta_{0}+\theta_{1}x_{1,t-1}+\theta_{2} x_{2,t-1}+u_{t}$

