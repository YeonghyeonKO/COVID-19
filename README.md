# COVID-19
(20-1) Regression Analysis Project

- The Raw Data have been collected by John's Hopkins University, "https://data.world/covid-19-data-resource-hub/covid-19-case-counts/workspace"

1) Preprocess the confirmed cases data from the day the first case had appeared to May 31th in 2020.

2) Using nls2 and nls functions, estimate the initial values of parameters for Logistic, Bertalanffy, Gompertz models from South Korea, Brazil, Hungary and Paraguay.

3) Divide the preprocessed data for training sets(from the day the first case had appeared to April 19th in 2020) and testing sets(from April 20th to May 31th in 2020).

4) Predict the new daily confirmed cases by three models and evaluate each models based on MSE which is equivalent to R-Squared(Coefficient of Determination).

5) Since these three models analyzed above do not properly reflect new changes, build a new model with Segmented Poisson Regression model and compare all models based on MSE.
