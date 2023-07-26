# Volatility-Model-comparison

This Project is an assignment in FINA 5423(Financial Econometrics and Computational Methods II) at the Carlson School of Management, University of Minnesota.

Volatility has several facts. First, the volatility tends to stay at a high level if it is already high and stays at a low level if it is low. Second, the volatility reacts differently to good or bad market time. When modeling volatility, we can use ARCH and GARCH models. 

In the script, I used CRSP monthly CRSP- value-weighted portfolio return and constructed the GJR-GARCH model and it takes this form: 
$$
r_{t} & =\mu+\epsilon_{t} \\
\epsilon_{t} & =\sigma_{t} z_{t}, \quad z_{t} \stackrel{i i d}{\sim} \mathcal{N}(0,1) \\
\sigma_{t}^{2} & =\omega+\left(\alpha+\gamma \mathbf{1}\left\{r_{t-1}<\mu\right\}\right) \epsilon_{t-1}^{2}+\beta \sigma_{t-1}^{2}
$$
