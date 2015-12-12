# Monte Carlo Simulations

Option prices randomly evolve over time. 
They are usually modeled through Black and Scholes equations and depend on some parameters such as :
- the interest rate of the option.
- the contract's term.
- the option strike : the fixed price at which the owner of the option can buy, or sell the underlying asset.

For fixed values of these parameters, the option price can vary much and so does the estimator's variance.

Implementing Monte Carlo's techniques can enable to reduce the variance of the estimator. 

Here are the steps we followed :
- Simulate prices as a brownian motion.
- Estimating mean and variance with Standard Monte Carlo.
- Estimating mean and variance with Antithetic Variates.
- Estimating mean and variance with Control Variates.

This is group project of 3 engineers students at ENSAE ParisTech : 
- Capucine Columelli 
- Marie-Laura Luisada
- Charles Tanguy
