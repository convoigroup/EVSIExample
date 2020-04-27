# EVSI Example
A repository for R code to compute EVSI using one of four methods, developed in the following papers:
1) Regression Based Method: Strong, M., Oakley, J.E., Brennan, A. and Breeze, P., 2015. Estimating the expected value of sample information using the probabilistic sensitivity analysis sample: a fast, nonparametric regression-based method. Medical Decision Making, 35(5), pp.570-583.
2) Importance Sampling Method: Menzies, N.A., 2016. An efficient estimator for the expected value of sample information. Medical Decision Making, 36(3), pp.308-320.
3) Gaussian Approximation Method: Jalal, H. and Alarid-Escudero, F., 2018. A Gaussian approximation approach for value of information analysis. Medical Decision Making, 38(2), pp.174-188.
4) Moment Matching Method: Heath, A., Manolopoulou, I. and Baio, G., 2019. Estimating the Expected Value of Sample Information across Different Sample Sizes Using Moment Matching and Nonlinear Regression. Medical Decision Making, 39(4), pp.347-359.

# Description:
The repository has four folders with R scripts that contain general functions for estimating EVSI using different methods. The three example decision models used for demonstration are a chemotherapy trial, a chronic pain study and a model proposed by Ades et al.

## "EVSI Estimation Methods" folder

This folder has five script files for computing EVSI estimates on any types of decision models.

**evsi_moment_matching_method.R** uses moment matching method to calculate EVSI proposed by Heath et al. This method requires estimating Expected Value of Partial Perfect Information (EVPPI) first to obtain the inner expectation of Incremental Net Benefit (INB) conditioned on the parameters that will be updated by the future data (focal parameters). The method then samples quantiles from simulations of the focal parameters. Next, the method simulates future datasets from the sampling distribution of the data and uses Bayesian methods to update the distribution of model parameters. This method requires reruns of the Probabilistic Sensitivity Analysis (PSA) to update the distribution of the net monetary benefit.

**evsi_Gaussian_approximation_method.R** uses a Gaussian approximation (GA) approach to calculate EVSI as proposed by Jalal et al. This approach uses a single PSA dataset and a linear meta-model step to compute the EVSI by using the GA method to compute the preposterior distribution of the parameters of interest.

**predict_ga.R** contains helper functions for the EVSI method by Jalal et al. method. The script will compute the preposterior for each of the basis functions of the Generalized Additive Model (GAM) model.

**evsi_importance_sampling_method.R** performs the method of estimating EVSI proposed by Menzies et al. This method requires EVPPI estimation in order to obtain the inner conditional expectation of INB on the focal parameters. Then the script will compute the likelihood of the simulated dataset conditional on all PSA simulations and calculate the weighted sum of samples. Each weighted sum is used to estimate EVSI.

**evsi_regression_based_method.R** will carry out a nonparametric regression-based method to estimate EVSI proposed by Strong et al. It requires only the PSA sample. After simulating and summarizing dataset, fitted values from regression model (eg. generalized additive model) estimate EVSI.

## "Decision Model Example - Chemotherapy" folder

**chemotherapy_decision_model.R** will run the chemotherapy decision model to calculate the INB of reducing chemotherapy treatment side effects and to generate PSA simulation data. The underlying model is a Markov model with four states (ambulatory care, hospital care, recovery and death) that evaluates two chemotherapy interventions. It assumes that the patients can transition between four different states with certain probabilities that depending on the state.

**chemotherapy_evsi_four_methods.R** will run the EVSI calculation steps using the four estimation methods for the chemotherapy decision model.

## "Decision Model Example - Chronic Pain" folder

**chronic_pain_decision_model.R** will run a Markov model to simulate a chronic pain treatment trial with 10 states, where each state has an associated Quality of Life (QoL) and cost. Patients have risk of experiencing adverse events and may withdraw from the initial treatment due to side effects or lack of efficacy. Following this, they can be offered an alternative therapy or withdraw completely. Similarly, patients can also experience adverse events from this second line of treatment and can withdraw either due to adverse events or lack of efficacy. If patients withdraw from the second line of treatment, they can receive further treatment or discontinue, as both states considered end points in the model. Initially, patients can receive morphine or a new treatment and if required, a second line treatment of oxycodone.

**chronic_pain_evsi_four_methods.R** will run the EVSI calculation steps using the four estimation methods for the chronic pain decision model.

## "Decision Model Example - Ades et al" folder

**Ades_et_al_decision_model.R** will run simulations to compare a hypothetical pair of standard and new treatments. The aim of the study is evaluate the reduction of critical events associated with the intervention. This critical event leads to decrease in Quality-Adjusted Life Years (QALYs) for the remainder of the patient’s life. The new treatment lowers the risk of critical event but the patient may also experience side effect, which gives a short-term reduction in QALYs along with direct cost of additional treatment. The model has 11 parameters, of which 4 are subject to uncertainty, which is modelled using 4 mutually independent distributions. For this case study, we consider 3 different data collection scenarios: EVSI for the Probability of Side Effects (Pse), EVSI for Quality of Life after Critical Event (Qe) and EVSI for the Treatment Effect Size in Odds Ratio (OR) measure.

**Ades_et_al_evsi_four_methods.R** will run the EVSI calculation steps using the four estimation methods for the Ades et al model.
