# sensitivity_analysis
> Project in collaboration with Soroush Zamanian and Ge Liu

Failure in wastewater infrastructure systems is recognized as a serious, worldwide concern which can have irreversible impacts on health, environment, and the economy. Concrete sewer pipes are most commonly used in wastewater systems. Cracks form in them when the tensile stress on the pipes exceeds its tensile strain. The geometry of the pipe, materials used, properties of the soil in which the pipe is buried, etc. can affect the tensile stress on the pipe. There are currently 11 known factors that can influence the tensile stress. In this project, we considered 5 of them. Our goal was to identify the dominant factors that influence the tensile stress on the pipe. We also estimated a safety region for the significant variables that would keep the tensile stress within permissible limits.


We used a simulator that provides the tensile stress as the output response. Due to the high com- putational complexity, the simulator takes about 2 hours for a single run. We fitted the simulator data to two models that we studied in the course; Bayesian Gaussian Process (GP) and Bayesian Additive Regression Trees (BART). We then used the fitted models to perform Sensitivity Analysis (SA). We also selected one of the two models based on Mean Squared Prediction Error (MSPE) and it was used for some further analysis.


## Experiment design
We chose a sample size of 50 where 40 samples were used as training data to fit the two models (Bayesian GP and BART) and the remaining 10 samples were used as the testing data for validation and to measure the predicting performance of each model. To generate a randomized design of the experiment, we implemented a space-filling design method, the Latin Hypercube Sampling (LHS) on the domain of the five predictors. Based on the prior information we have on the five predictors, all of them are distributed independently as the lognormal distribution with corresponding mean and standard deviation of the Gaussian component are given in Table ??.

Table 2.1: Simulator inputs
|Xi   |Predictor name|mean (μi)|sd (σi)|
|----:|-------------:|--------:|------:|
|X1   |compressive strength of concrete|8.597|0.184|
|X2.  |elastic modulus of bedding soil|8.153|0.198|
|X3.  |density of bedding soil|8.759|0.198|
|X4.  |elastic modulus of backfilling soil|−2.763|0.071|
|X5.  |density of backfilling soil|−2.595|0.071|

For i = 1,...,5, we generate the inputs using LHS such that Xi ∈ [ai,bi], where ai and bi are the 5th and the 95th quantile of the lognormal distribution of Xi respectively. The response Y (x), tensile stress, is recorded by running the simulator. Before fitting the dataset to any model, the samples for each Xi are normalized to the interval [0, 1].
