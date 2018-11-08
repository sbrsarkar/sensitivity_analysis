# Sensitivity Analysis
> Project in collaboration with Soroush Zamanian and Ge Liu

Failure in wastewater infrastructure systems is recognized as a serious, worldwide concern which can have irreversible impacts on health, environment, and the economy. Concrete sewer pipes are most commonly used in wastewater systems. Cracks form in them when the tensile stress on the pipes exceeds its tensile strain. The geometry of the pipe, materials used, properties of the soil in which the pipe is buried, etc. can affect the tensile stress on the pipe. There are currently 11 known factors that can influence the tensile stress. In this project, we considered 5 of them. Our goal was to identify the dominant factors that influence the tensile stress on the pipe. We also estimated a safety region for the significant variables that would keep the tensile stress within permissible limits.


We used a simulator that provides the tensile stress as the output response. Due to the high com- putational complexity, the simulator takes about 2 hours for a single run. We fitted the simulator data to two models that we studied in the course; Bayesian Gaussian Process (GP) and Bayesian Additive Regression Trees (BART). We then used the fitted models to perform Sensitivity Analysis (SA). We also selected one of the two models based on Mean Squared Prediction Error (MSPE) and it was used for some sensitivity_analysisfurther analysis.


## Experiment design
We chose a sample size of 50 where 40 samples were used as training data to fit the two models (Bayesian GP and BART) and the remaining 10 samples were used as the testing data for validation and to measure the predicting performance of each model. To generate a randomized design of the experiment, we implemented a space-filling design method, the Latin Hypercube Sampling (LHS) on the domain of the five predictors. Based on the prior information we have on the five predictors, all of them are distributed independently as the lognormal distribution with corresponding mean and standard deviation of the Gaussian component are given in Table ??.

Table 2.1: Simulator inputs


|<img src="http://latex.codecogs.com/gif.latex?X_i" title="X_i" />   |Predictor name|mean (<img src="http://latex.codecogs.com/gif.latex?\mu_i" title="\mu_i" />)|sd (<img src="http://latex.codecogs.com/gif.latex?\sigma_i" title="\sigma_i" />)|
|-----|:------------:|:-------:|------:|
|<img src="http://latex.codecogs.com/gif.latex?X_1" title="X_1" />   |compressive strength of concrete|8.597|0.184|
|<img src="http://latex.codecogs.com/gif.latex?X_2" title="X_2" />   |elastic modulus of bedding soil|8.153|0.198|
|<img src="http://latex.codecogs.com/gif.latex?X_3" title="X_3" />   |density of bedding soil|8.759|0.198|
|<img src="http://latex.codecogs.com/gif.latex?X_4" title="X_4" />   |elastic modulus of backfilling soil|−2.763|0.071|
|<img src="http://latex.codecogs.com/gif.latex?X_5" title="X_5" />   |density of backfilling soil|−2.595|0.071|

For i = 1,...,5, we generate the inputs using LHS such that <img src="http://latex.codecogs.com/gif.latex?X_i\in[a_i,b_i]" title="X_i\in[a_i,b_i]" />, where <img src="http://latex.codecogs.com/gif.latex?a_i" title="a_i" />  and <img src="http://latex.codecogs.com/gif.latex?b_i" title="b_i" /> are the 5th and the 95th quantile of the lognormal distribution of <img src="http://latex.codecogs.com/gif.latex?X_i" title="X_i" /> respectively. The response <img src="http://latex.codecogs.com/gif.latex?Y(x)" title="Y(x)" />, tensile stress, is recorded by running the simulator. Before fitting the dataset to any model, the samples for each Xi are normalized to the interval <img src="http://latex.codecogs.com/gif.latex?[0,1]" title="[0,1]" />.


## Bayesian Gaussian Process
The first model considered is the Bayesian Gaussian Process (GP) model with a non-zero trend function. Since the five material variables are correlated (the backfill shear velocity and the density of backfill are correlated and the bedding shear velocity and the density of bedding are correlated), we propose using the Bayesian GP model with specified linear trend as functions of the 5 variables capturing the correlated structure in them.
z(x) ∼ GP(f(x)T β, c(·; λ−1, ρ), (2.1)
where f(x) = (1, f1(x), f2(x), f3(x), f4(x), f5(x)) is the assumed trend function in the GP and
β = (β1, β2, . . . , β6) ∈ R6 are the unknown regression coefficients. We use the separable Gaussian
covariance function, c(x, x′) = λ−1 􏰀5 ρ|xi−x′i|, where {ρi} are the unknown correlation parameters i=1 i
and λ is the precision. The trend functions are specified as,
f1(x) = √x1, f2(x) = x4x2,
f3(x) = x5x23, f4(x) = x4, f5(x) = x5.
