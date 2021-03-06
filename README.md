# Sensitivity Analysis
> Project in collaboration with Soroush Zamanian and Ge Liu

Failure in wastewater infrastructure systems is recognized as a serious, worldwide concern which can have irreversible impacts on health, environment, and the economy. Concrete sewer pipes are most commonly used in wastewater systems. Cracks form in them when the tensile stress on the pipes exceeds its tensile strain. The geometry of the pipe, materials used, properties of the soil in which the pipe is buried, etc. can affect the tensile stress on pipes. There are currently 11 known factors that can influence the tensile stress. In this project, we considered 5 of them. Our goal was to identify the dominant factors that influence the tensile stress on the pipe. We also estimated a safety region for the significant variables that would keep the tensile stress within permissible limits.


We used a simulator that provides the tensile stress as the output response. Due to the high computational complexity, the simulator takes about 2 hours for a single run. We fitted the simulator data to two models; Bayesian Gaussian Process (GP) and Bayesian Additive Regression Trees (BART). We then used the fitted models to perform Sensitivity Analysis (SA). We also selected one of the two models based on Mean Squared Prediction Error (MSPE) and it was used for some further analysis.


## Experiment design
We chose a sample size of 50 where 40 samples were used as training data to fit the two models (Bayesian GP and BART) and the remaining 10 samples were used as testing data and to measure the predictive performance of each model. To generate a randomized design of the experiment, we used a space-filling design method, the Latin Hypercube Sampling (LHS) on the domain of the five predictors. Based on the prior information we have on the five predictors, all of them are independent lognormal random variables with corresponding mean and standard deviation of the Gaussian component as given in the table below.

Simulator inputs


|<img src="http://latex.codecogs.com/gif.latex?X_i" title="X_i" />   |Predictor name|mean (<img src="http://latex.codecogs.com/gif.latex?\mu_i" title="\mu_i" />)|sd (<img src="http://latex.codecogs.com/gif.latex?\sigma_i" title="\sigma_i" />)|
|-----|:------------:|:-------:|------:|
|<img src="http://latex.codecogs.com/gif.latex?X_1" title="X_1" />   |compressive strength of concrete|8.597|0.184|
|<img src="http://latex.codecogs.com/gif.latex?X_2" title="X_2" />   |elastic modulus of bedding soil|8.153|0.198|
|<img src="http://latex.codecogs.com/gif.latex?X_3" title="X_3" />   |density of bedding soil|8.759|0.198|
|<img src="http://latex.codecogs.com/gif.latex?X_4" title="X_4" />   |elastic modulus of backfilling soil|−2.763|0.071|
|<img src="http://latex.codecogs.com/gif.latex?X_5" title="X_5" />   |density of backfilling soil|−2.595|0.071|

For i = 1,...,5, we generate the inputs using LHS such that <img src="http://latex.codecogs.com/gif.latex?X_i\in[a_i,b_i]" title="X_i\in[a_i,b_i]" />, where <img src="http://latex.codecogs.com/gif.latex?a_i" title="a_i" />  and <img src="http://latex.codecogs.com/gif.latex?b_i" title="b_i" /> are the 5th and the 95th quantile of the lognormal distribution of <img src="http://latex.codecogs.com/gif.latex?X_i" title="X_i" /> respectively. The response <img src="http://latex.codecogs.com/gif.latex?Y(x)" title="Y(x)" />, tensile stress, is recorded by running the simulator. Before fitting the dataset to any model, the samples for each Xi are normalized to the interval <img src="http://latex.codecogs.com/gif.latex?[0,1]" title="[0,1]" />.


## Bayesian Gaussian Process
The first model considered was the Bayesian Gaussian Process (GP) model with a non-zero trend function. Since the five material variables were correlated (the backfill shear velocity and the density of backfill were correlated and the bedding shear velocity and the density of bedding were correlated), we used the Bayesian GP model with trend function that captured the known structure in them. The output of a GP is given as

<img src="http://latex.codecogs.com/gif.latex?z(x)&space;\sim&space;GP(f(x)^T\beta,&space;c(\cdot;&space;\lambda-1,\rho))" title="z(x) \sim GP(f(x)^T\beta, c(\cdot; \lambda-1,\rho))" />, 

where <img src="http://latex.codecogs.com/gif.latex?f(x)&space;=&space;(1,&space;f_1(x),&space;f_2(x),&space;f_3(x),&space;f_4(x),&space;f_5(x))" title="f(x) = (1, f_1(x), f_2(x), f_3(x), f_4(x), f_5(x))" /> is the assumed trend function in the GP and <img src="http://latex.codecogs.com/gif.latex?\beta\in\mathbb{R}^6" title="\beta\in\mathbb{R}^6" /> are the unknown regression coefficients. We use the separable Gaussian covariance function,

<img src="http://latex.codecogs.com/gif.latex?c(x,x')&space;=&space;\lambda^{-1}\prod_{i=1}^5\rho_i^{|x_i-x'_i|}" title="c(x,x') = \lambda^{-1}\prod_{i=1}^5\rho_i^{|x_i-x'_i|}" />, 

where <img src="http://latex.codecogs.com/gif.latex?\{\rho_i\}" title="\{\rho_i\}" /> are the unknown correlation parameters and <img src="http://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /> is the precision. The trend functions are specified as,
<img src="http://latex.codecogs.com/gif.latex?\begin{align*}&space;f_1(x)&space;&=&space;\sqrt{x_1},&space;\quad&space;f_2(x)&space;=&space;x_4x_2,\\&space;f_3(x)&space;&=&space;x_5&space;x_2^3,&space;\quad&space;f_4(x)&space;=&space;x_4,&space;\quad&space;f_5(x)&space;=&space;x_5.&space;\end{align*}" title="\begin{align*} f_1(x) &= \sqrt{x_1}, \quad f_2(x) = x_4x_2,\\ f_3(x) &= x_5 x_2^3, \quad f_4(x) = x_4, \quad f_5(x) = x_5. \end{align*}" />

### Sampling
In the Bayesian setting, the unknown parameters in our model <img src="http://latex.codecogs.com/gif.latex?\rho_i" title="\rho_i" />, <img src="http://latex.codecogs.com/gif.latex?\beta_i" title="\beta_i" />, <img src="http://latex.codecogs.com/gif.latex?\lambda" title="\lambda" /> are treated as random variables with known prior distributions. The aim is to estimate the parameters from the samples of posterior distribution. The choice of priors are <img src="http://latex.codecogs.com/gif.latex?\lambda\sim&space;Gamma(2,0.5)" title="\lambda\sim Gamma(2,0.5)" />, and then the posterior is also a Gamma distribution. For <img src="http://latex.codecogs.com/gif.latex?\rho" title="\rho" />, the priors are <img src="http://latex.codecogs.com/gif.latex?\rho_1,\rho_3&space;\sim&space;Beta(5,1),&space;\rho_2&space;\sim&space;Beta(2,&space;5)&space;\text{&space;and&space;}&space;\rho_4,&space;\rho_5&space;\sim&space;Beta(2,&space;2)" title="\rho_1,\rho_3 \sim Beta(5,1), \rho_2 \sim Beta(2, 5) \text{ and } \rho_4, \rho_5 \sim Beta(2, 2)" />, the choice of the hyperparameters in the Beta distribution is based on the prior information of the strength of the correlation between each Xi and Y(x). 

