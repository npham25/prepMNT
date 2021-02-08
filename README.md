# prepMNT

#### This package preprocesses the MNTs data by removing MNTs which have few available values and imputing the zero values using Quantile Regression Imputation of Left-Censored data


```{r}
devtools::install_github('npham25/prepMNT')
library(prepMNT)

data_prepMNT <- prepMNT(data,pctZeroCut=75, freqCut = 19, uniqueCut = 20)
```
