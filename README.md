# COVID-19

- The Raw Data have been collected by John's Hopkins University, "https://data.world/covid-19-data-resource-hub/covid-19-case-counts/workspace"
- This is the repository for the project in [BIBS](http://bibs.snu.ac.kr) (BioInformatics & BioStatistics, Seoul National University)

<br>

 1) Preprocess the confirmed cases data from the day the first case had appeared to August 31th (2020) / December 14th (2020) / **March 20th (2021) (latest)**
    > Using **Peak Detection Algorithm**, we found breakpoints of each growth curve for each country.

 2) Using nls2 and nls functions, estimate the initial values of parameters for Logistic and Gompertz models

 3) Linear Regression with variable indicators such as population, GDP, Imports, etc.

