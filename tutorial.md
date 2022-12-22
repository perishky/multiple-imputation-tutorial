
# Missing data

## Load the dataset


```r
library(mice)
```

```
## Loading required package: lattice
```

```
## 
## Attaching package: 'mice'
```

```
## The following objects are masked from 'package:base':
## 
##     cbind, rbind
```


```r
load("dataset/dataset.rda")

idx <- match(rownames(sva.smoking), colnames(meth))
meth <- meth[,idx]
samples <- samples[idx,]
sva.age <- sva.age[idx,]
```

Double-check that data objects are synchronized.

```r
stopifnot(identical(rownames(sva.smoking), colnames(meth)))
stopifnot(identical(rownames(sva.smoking), as.character(samples$gsm)))
stopifnot(identical(rownames(sva.smoking), rownames(sva.age)))
```

Make sure that character variables are factors.

```r
for (i in 1:ncol(samples))
    if (is.character(samples[[i]]))
        samples[[i]] <- as.factor(samples[[i]])
```

Remove the 'gsm' (individual id) and 'gse' (dataset) variables
since no individual has more than methylation profile.

```r
samples$gsm <- samples$gse <- NULL
```



## Exploring missing patterns

We'll be randomly removing data.
To make these analyses repeatable, we
set a random seed.

```r
set.seed(20180312)
```

Remove 20 values from smoking status and age.

```r
samples.mis <- samples
samples.mis$smoking[sample(1:nrow(samples.mis), 20)] <- NA
samples.mis$age[sample(1:nrow(samples.mis),20)] <- NA
```

Have a look at the missingness pattern.

```r
md.pattern(samples.mis)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.svg)

```
##     sex smoking age   
## 165   1       1   1  0
## 16    1       1   0  1
## 16    1       0   1  1
## 4     1       0   0  2
##       0      20  20 40
```

According to the `mice::md.pattern()` documentation:
> A matrix with ‘ncol(x)+1’ columns, in which each row corresponds
> to a missing data pattern (1=observed, 0=missing).  Rows and
> columns are sorted in increasing amounts of missing information.
> The last column and row contain row and column counts,
> respectively.


```r
sum(is.na(samples.mis$age) & is.na(samples.mis$smoking))
```

```
## [1] 4
```

```r
sum(is.na(samples.mis$age) & !is.na(samples.mis$smoking))
```

```
## [1] 16
```

```r
sum(!is.na(samples.mis$age) & is.na(samples.mis$smoking))
```

```
## [1] 16
```

```r
sum(!is.na(samples.mis$age) & !is.na(samples.mis$smoking))
```

```
## [1] 165
```

## Missing completely at random

To keep things simple for now,
we'll analyse data with only one variable (smoking status)
missing values.  Well remove 20% of the values.

```r
pct.missing <- 0.2
```

The way we remove the values, data should be 'missing completely at random'.

```r
samples.mcar <- samples
missing.idx <- sample(1:nrow(samples.mcar),
                      floor(pct.missing*nrow(samples.mcar)))
samples.mcar$smoking[missing.idx] <- NA
```

Have a look at the missingness patterns.

```r
table(samples.mcar$smoking, useNA="a")
```

```
## 
## current  former   never    <NA> 
##      20       0     141      40
```

```r
md.pattern(data.frame(samples.mcar))
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.svg)

```
##     age sex smoking   
## 161   1   1       1  0
## 40    1   1       0  1
##       0   0      40 40
```

Now we test if smoking status missingness is associated with
any other variables of interest.

```r
coef(summary(glm(is.na(smoking) ~ age + sex, data=samples.mcar, family="binomial")))
```

```
##                 Estimate Std. Error    z value  Pr(>|z|)
## (Intercept) -1.138172924  1.5286906 -0.7445411 0.4565492
## age         -0.009932406  0.0261631 -0.3796341 0.7042170
## sexmale      0.414839240  0.4079895  1.0167889 0.3092539
```

Assuming that age and sex are the only external variables of interest here,
we have shown that smoking status is missing completely at random.

## Missing at random

Now suppose we remove smoking status values
so that more are missing from females than males.

```r
samples.mar <- samples[,c("smoking","sex","age")]
missing.idx <- sample(1:nrow(samples.mar),
                      floor(pct.missing*nrow(samples.mar)),
                      prob=ifelse(samples.mar$sex == "male", 0.2, 0.8))
samples.mar$smoking[missing.idx] <- NA
```

Have a look at the patterns.

```r
table(samples.mar$smoking, useNA="a")
```

```
## 
## current  former   never    <NA> 
##      16       0     145      40
```

```r
md.pattern(samples.mar)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.svg)

```
##     sex age smoking   
## 161   1   1       1  0
## 40    1   1       0  1
##       0   0      40 40
```

Investigate missingness by sex.

```r
table(sex=samples.mar$sex, missing=is.na(samples.mar$smoking))
```

```
##         missing
## sex      FALSE TRUE
##   female    41   24
##   male     120   16
```

We should see a correlation of missingness with sex.

```r
coef(summary(glm(is.na(smoking) ~ age + sex, data=samples.mar, family="binomial")))
```

```
##                Estimate Std. Error    z value     Pr(>|z|)
## (Intercept)  1.64203786 1.66209832  0.9879306 3.231866e-01
## age         -0.03819851 0.02886948 -1.3231451 1.857871e-01
## sexmale     -1.60720872 0.38880536 -4.1337103 3.569534e-05
```

## Ad hoc solutions for testing associations

### Solution: Complete case (CC) analysis

Here we just remove samples with missing values.

Here are the full data association statistics with CpG site cg07178945.

```r
data <- cbind(samples, sva.smoking)
data$cpg <- meth["cg07178945",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##      Estimate    Std. Error       t value      Pr(>|t|) 
## -4.643091e-02  6.608516e-03 -7.025921e+00  4.984723e-11
```

In the complete case analysis of MCAR, the association with smoking status
is still there but much weaker due to reduced sample size.

```r
data <- cbind(samples.mcar, sva.smoking)
data$cpg <- meth["cg07178945",]
data <- na.omit(data)
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##      Estimate    Std. Error       t value      Pr(>|t|) 
## -4.254030e-02  7.164111e-03 -5.937973e+00  2.511050e-08
```

The full data analysis of cg14024579.

```r
data <- cbind(samples, sva.smoking)
data$cpg <- meth["cg14024579",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##    Estimate  Std. Error     t value    Pr(>|t|) 
## 0.009756944 0.005914970 1.649533908 0.100895782
```

This time, with complete case analysis of MAR, the association
becomes much stronger.

```r
data <- cbind(samples.mar, sva.smoking)
data$cpg <- meth["cg14024579",]
data <- na.omit(data)
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##    Estimate  Std. Error     t value    Pr(>|t|) 
## 0.021588732 0.007140281 3.023512913 0.003015148
```

### Solution: Impute mean/median/mode

To avoid losing power, we could replace missing values
with 'reasonable' values, e.g. mean, median or mode.


The full data analysis of cg07178945.

```r
data <- cbind(samples, sva.smoking)
data$cpg <- meth["cg07178945",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##      Estimate    Std. Error       t value      Pr(>|t|) 
## -4.643091e-02  6.608516e-03 -7.025921e+00  4.984723e-11
```

Complete case:

```r
data <- cbind(samples.mis, sva.smoking)
data$cpg <- meth["cg07178945",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##      Estimate    Std. Error       t value      Pr(>|t|) 
## -4.806207e-02  8.053398e-03 -5.967924e+00  2.055330e-08
```

By imputing with the mode (most people are never smokers),
we reclassify many smokers as never smokers and reduce power.

```r
data <- cbind(samples.mcar, sva.smoking)
data$smoking[which(is.na(data$smoking))] <- names(which.max(table(samples.mcar$smoking)))
data$cpg <- meth["cg07178945",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##      Estimate    Std. Error       t value      Pr(>|t|) 
## -4.384964e-02  7.212197e-03 -6.079928e+00  7.766277e-09
```

In the missing-at-random dataset, imputation with the mode
weakens the association compared to complete case analysis.

```r
data <- cbind(samples.mar, sva.smoking)
data$smoking[which(is.na(data$smoking))] <- names(which.max(table(samples.mar$smoking)))
data$cpg <- meth["cg07178945",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##      Estimate    Std. Error       t value      Pr(>|t|) 
## -4.517944e-02  8.141135e-03 -5.549526e+00  1.089015e-07
```

### Solution: Imputing to conserve statistical properties

We can make the imputation a bit more interesting
by imputing randomly according to properties of the existing data.
In the smoking example, we impute randomly so that the ratio
of smokers to never-smokers remains the same.

Not an improvement in the association of cg21161138 with smoking in the MCAR
dataset compared to imputing the mode.

Mode:

```r
data <- cbind(samples.mcar, sva.smoking)
data$smoking[which(is.na(data$smoking))] <- names(which.max(table(samples.mcar$smoking)))
data$cpg <- meth["cg21161138",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##     Estimate   Std. Error      t value     Pr(>|t|) 
## 9.185805e-02 5.556028e-03 1.653304e+01 3.803867e-37
```


```r
data <- cbind(samples.mcar, sva.smoking)
idx <- which(is.na(data$smoking))
p.smoker <- sum(data$smoking == "current", na.rm=T)/nrow(data)
data$smoking[idx] <- sample(c("current","never"),
                            size=length(idx),
                            prob=c(p.smoker, 1-p.smoker),
                            replace=T)
data$cpg <- meth["cg21161138",]
coef(summary(lm(cpg ~ ., data)))["smokingnever",]
```

```
##     Estimate   Std. Error      t value     Pr(>|t|) 
## 6.434755e-02 6.214994e-03 1.035360e+01 9.170254e-20
```

### Solution: Multiple imputation

In multiple imputation, we attempt to allow other variables
to help decide on the best values to impute for each sample
instead of assigning values blindly.

Multiple imputation begins by imputing 'reasonable' values
for each variable (e.g. mean) and then repeatedly using
the currently imputed dataset to predict values for one of the
variables with missing values.
This is repeated some maximum number of iterations or until
values don't change much.

To take into account the uncertainty in the imputed values,
multiple datasets are imputed.  Downstream association
testing is then applied to each dataset and the summary
statistics pooled.

#### Simple example

To make things interesting, we remove values for both smoking and age.

```r
samples.mis <- samples
samples.mis$smoking[sample(1:nrow(samples.mis), 30)] <- NA
samples.mis$age[sample(1:nrow(samples.mis),30)] <- NA
```

For future reference, we identify the positions of missing values.

```r
smoking.idx <- which(is.na(samples.mis$smoking))
age.idx <- which(is.na(samples.mis$age))
```

We begin multiple imputation by filling in missing values with reasonable
starting values, e.g. mean or mode.

```r
samples.mis$smoking[smoking.idx] <- names(which.max(table(samples.mis$smoking)))
samples.mis$age[age.idx] <- mean(samples.mis$age, na.rm=T)
```

For comparison, we show the relationship of the imputed dataset with the
original dataset.

```r
table(samples.mis$smoking, samples$smoking)
```

```
##          
##           current former never
##   current      20      0     0
##   former        0      0     0
##   never         2      0   179
```

```r
cor(samples.mis$age, samples$age)
```

```
## [1] 0.9383978
```

Next we predict age missing values using all other variables in the dataset.

```r
fit <- glm(age ~ ., samples.mis[-age.idx,], family="gaussian")
samples.mis$age[age.idx] <- predict(fit, newdata=samples.mis[age.idx,], type="response")
```

We do the same for smoking status.

```r
fit <- glm(sign(smoking == "current") ~ ., samples.mis[-smoking.idx,], family="binomial")
probs <- predict(fit, newdata=samples.mis[smoking.idx,], type="response")
samples.mis$smoking[smoking.idx] <- ifelse(probs > 0.5, "current", "never")
```

Unfortunately we observe no improvement.

```r
table(samples.mis$smoking, samples$smoking)
```

```
##          
##           current former never
##   current      20      0     0
##   former        0      0     0
##   never         2      0   179
```

```r
cor(samples.mis$age, samples$age)
```

```
## [1] 0.9358325
```

This actually makes sense because none of the variables
in the dataset are strongly associated with age or smoking status.

```r
coef(summary(glm(smoking ~ ., samples, family=binomial)))
```

```
##               Estimate Std. Error   z value  Pr(>|z|)
## (Intercept) 0.44535941 1.92931718 0.2308378 0.8174408
## age         0.02939856 0.03346428 0.8785056 0.3796694
## sexmale     0.06019524 0.49696980 0.1211245 0.9035924
```

```r
coef(summary(glm(age ~ ., samples, family=gaussian)))
```

```
##               Estimate Std. Error    t value     Pr(>|t|)
## (Intercept)  55.968180   1.586237 35.2836268 2.450239e-87
## smokingnever  1.328764   1.514053  0.8776207 3.812133e-01
## sexmale      -2.753330   1.010538 -2.7246172 7.013788e-03
```

#### Including DNA methylation to improve imputation

To improve imputation,
we include sites well-known to be predictive of smoking (cg05575921)
and age (cg21572722).

```r
samples.mis <- samples
samples.mis$cg05575921 <- meth["cg05575921",]
samples.mis$cg21572722 <- meth["cg21572722",]  
```

We also reset the dataset.

```r
samples.mis$age[age.idx] <- NA
samples.mis$smoking[smoking.idx] <- NA
samples.mis$smoking[smoking.idx] <- names(which.max(table(samples.mis$smoking)))
samples.mis$age[age.idx] <- mean(samples.mis$age, na.rm=T)
```

Once again, for reference we note comparison with the full dataset.

```r
table(samples.mis$smoking, samples$smoking)
```

```
##          
##           current former never
##   current      20      0     0
##   former        0      0     0
##   never         2      0   179
```

```r
cor(samples.mis$age, samples$age)
```

```
## [1] 0.9383978
```

Again, we predict missing values using the other variables.

```r
fit <- glm(age ~ ., samples.mis[-age.idx,], family="gaussian")
samples.mis$age[age.idx] <- predict(fit, newdata=samples.mis[age.idx,], type="response")

fit <- glm(sign(smoking == "current") ~ ., samples.mis[-smoking.idx,], family="binomial")
probs <- predict(fit, newdata=samples.mis[smoking.idx,], type="response")
samples.mis$smoking[smoking.idx] <- ifelse(probs > 0.5, "current", "never")
```

We observe greater agreement with the full dataset.

```r
table(samples.mis$smoking, samples$smoking)
```

```
##          
##           current former never
##   current      22      0     0
##   former        0      0     0
##   never         0      0   179
```

```r
cor(samples.mis$age, samples$age)
```

```
## [1] 0.9714208
```

We repeat prediction.

```r
fit <- glm(age ~ ., samples.mis[-age.idx,], family="gaussian")
samples.mis$age[age.idx] <- predict(fit, newdata=samples.mis[age.idx,], type="response")

fit <- glm(sign(smoking == "current") ~ ., samples.mis[-smoking.idx,], family="binomial")
probs <- predict(fit, newdata=samples.mis[smoking.idx,], type="response")
samples.mis$smoking[smoking.idx] <- ifelse(probs > 0.5, "current", "never")
```

This time we see no further improvement.

```r
table(samples.mis$smoking, samples$smoking)
```

```
##          
##           current former never
##   current      22      0     0
##   former        0      0     0
##   never         0      0   179
```

```r
cor(samples.mis$age, samples$age)
```

```
## [1] 0.972284
```

#### A more interesting example

To observe more interesting behaviour, we need
two variables that have missing values and
are associated with one another.

This time, we'll create missing values in the age variable
and the age-associated CpG site (cg21572722).

```r
set.seed(20180312)
samples.mis <- samples
samples.mis$cg21572722 <- meth["cg21572722",]
cg.idx <- sample(1:nrow(samples.mis), 20)
samples.mis$cg21572722[cg.idx] <- NA
samples.mis$age[age.idx] <- NA
samples.mis$cg21572722[cg.idx] <- mean(samples.mis$cg21572722, na.rm=T)
samples.mis$age[age.idx] <- mean(samples.mis$age, na.rm=T)
```

Ten times we'll predict values for each variable and compare to the full dataset.

```r
for (i in 1:8) {
    cat(i, cor(samples.mis$age, samples$age), cor(samples.mis$cg21572722, meth["cg21572722",]), "\n")
        
    fit <- glm(age ~ ., samples.mis[-age.idx,], family="gaussian")
    samples.mis$age[age.idx] <- predict(fit, newdata=samples.mis[age.idx,], type="response")

    fit <- glm(cg21572722 ~ ., samples.mis[-cg.idx,], family="gaussian")
    samples.mis$cg21572722[cg.idx] <- predict(fit, newdata=samples.mis[cg.idx,], type="response")
}
```

```
## 1 0.9383978 0.9354883 
## 2 0.9710155 0.9726933 
## 3 0.9723632 0.9728255 
## 4 0.972379 0.9728259 
## 5 0.9723809 0.9728258 
## 6 0.9723815 0.9728257 
## 7 0.9723817 0.9728257 
## 8 0.9723818 0.9728257
```

## Imputation using mice

Fortunately an R package is available to do most of this for us.

Once again, we create a dataset with missing age and CpG site values.
We include surrogate variables because we'll use `mice`
to illustrate the steps of an EWAS for a couple of CpG sites.

```r
set.seed(20180312)
samples.mis <- cbind(samples, sva=sva.age)
samples.mis$cg21572722 <- meth["cg21572722",]
cg.idx <- sample(1:nrow(samples.mis), 35)
samples.mis$cg21572722[cg.idx] <- NA
samples.mis$age[age.idx] <- NA
```

Imputation is completed with a single call to `mice()`.
Here ask to impute 10 datasets, each time predicting values
for each variable with missing values at most 10 times (`maxit`).

```r
library(mice)
imp <- mice(samples.mis, m=10, maxit=10, print=F, seed=20180312)
```

```
## Warning: Number of logged events: 200
```

The prediction matrix determines which variables (columns) are used to
impute values for other variables (rows).

```r
imp$pred
```

```
##            smoking age sex sva.1 sva.2 sva.3 sva.4 sva.5 sva.6 sva.7 sva.8
## smoking          0   1   1     1     1     1     1     1     1     1     1
## age              1   0   1     1     1     1     1     1     1     1     1
## sex              1   1   0     1     1     1     1     1     1     1     1
## sva.1            1   1   1     0     1     1     1     1     1     1     1
## sva.2            1   1   1     1     0     1     1     1     1     1     1
## sva.3            1   1   1     1     1     0     1     1     1     1     1
## sva.4            1   1   1     1     1     1     0     1     1     1     1
## sva.5            1   1   1     1     1     1     1     0     1     1     1
## sva.6            1   1   1     1     1     1     1     1     0     1     1
## sva.7            1   1   1     1     1     1     1     1     1     0     1
## sva.8            1   1   1     1     1     1     1     1     1     1     0
## sva.9            1   1   1     1     1     1     1     1     1     1     1
## sva.10           1   1   1     1     1     1     1     1     1     1     1
## sva.11           1   1   1     1     1     1     1     1     1     1     1
## sva.12           1   1   1     1     1     1     1     1     1     1     1
## sva.13           1   1   1     1     1     1     1     1     1     1     1
## sva.14           1   1   1     1     1     1     1     1     1     1     1
## sva.15           1   1   1     1     1     1     1     1     1     1     1
## sva.16           1   1   1     1     1     1     1     1     1     1     1
## sva.17           1   1   1     1     1     1     1     1     1     1     1
## sva.18           1   1   1     1     1     1     1     1     1     1     1
## sva.19           1   1   1     1     1     1     1     1     1     1     1
## sva.20           1   1   1     1     1     1     1     1     1     1     1
## sva.21           1   1   1     1     1     1     1     1     1     1     1
## sva.22           1   1   1     1     1     1     1     1     1     1     1
## sva.23           1   1   1     1     1     1     1     1     1     1     1
## sva.24           1   1   1     1     1     1     1     1     1     1     1
## sva.25           1   1   1     1     1     1     1     1     1     1     1
## sva.26           1   1   1     1     1     1     1     1     1     1     1
## sva.27           1   1   1     1     1     1     1     1     1     1     1
## sva.28           1   1   1     1     1     1     1     1     1     1     1
## sva.29           1   1   1     1     1     1     1     1     1     1     1
## sva.30           1   1   1     1     1     1     1     1     1     1     1
## sva.31           1   1   1     1     1     1     1     1     1     1     1
## sva.32           1   1   1     1     1     1     1     1     1     1     1
## sva.33           1   1   1     1     1     1     1     1     1     1     1
## sva.34           1   1   1     1     1     1     1     1     1     1     1
## sva.35           1   1   1     1     1     1     1     1     1     1     1
## sva.36           1   1   1     1     1     1     1     1     1     1     1
## sva.37           1   1   1     1     1     1     1     1     1     1     1
## sva.38           1   1   1     1     1     1     1     1     1     1     1
## sva.39           1   1   1     1     1     1     1     1     1     1     1
## sva.40           1   1   1     1     1     1     1     1     1     1     1
## sva.41           1   1   1     1     1     1     1     1     1     1     1
## sva.42           1   1   1     1     1     1     1     1     1     1     1
## sva.43           1   1   1     1     1     1     1     1     1     1     1
## sva.44           1   1   1     1     1     1     1     1     1     1     1
## sva.45           1   1   1     1     1     1     1     1     1     1     1
## sva.46           1   1   1     1     1     1     1     1     1     1     1
## sva.47           1   1   1     1     1     1     1     1     1     1     1
## sva.48           1   1   1     1     1     1     1     1     1     1     1
## sva.49           1   1   1     1     1     1     1     1     1     1     1
## sva.50           1   1   1     1     1     1     1     1     1     1     1
## sva.51           1   1   1     1     1     1     1     1     1     1     1
## sva.52           1   1   1     1     1     1     1     1     1     1     1
## sva.53           1   1   1     1     1     1     1     1     1     1     1
## sva.54           1   1   1     1     1     1     1     1     1     1     1
## sva.55           1   1   1     1     1     1     1     1     1     1     1
## sva.56           1   1   1     1     1     1     1     1     1     1     1
## sva.57           1   1   1     1     1     1     1     1     1     1     1
## sva.58           1   1   1     1     1     1     1     1     1     1     1
## sva.59           1   1   1     1     1     1     1     1     1     1     1
## sva.60           1   1   1     1     1     1     1     1     1     1     1
## sva.61           1   1   1     1     1     1     1     1     1     1     1
## sva.62           1   1   1     1     1     1     1     1     1     1     1
## sva.63           1   1   1     1     1     1     1     1     1     1     1
## cg21572722       1   1   1     1     1     1     1     1     1     1     1
##            sva.9 sva.10 sva.11 sva.12 sva.13 sva.14 sva.15 sva.16 sva.17
## smoking        1      1      1      1      1      1      1      1      1
## age            1      1      1      1      1      1      1      1      1
## sex            1      1      1      1      1      1      1      1      1
## sva.1          1      0      0      0      0      0      0      0      0
## sva.2          1      1      1      1      1      1      1      1      1
## sva.3          1      1      1      1      1      1      1      1      1
## sva.4          1      1      1      1      1      1      1      1      1
## sva.5          1      1      1      1      1      1      1      1      1
## sva.6          1      1      1      1      1      1      1      1      1
## sva.7          1      1      1      1      1      1      1      1      1
## sva.8          1      1      1      1      1      1      1      1      1
## sva.9          0      1      1      1      1      1      1      1      1
## sva.10         1      0      1      1      1      1      1      1      1
## sva.11         1      1      0      1      1      1      1      1      1
## sva.12         1      1      1      0      1      1      1      1      1
## sva.13         1      1      1      1      0      1      1      1      1
## sva.14         1      1      1      1      1      0      1      1      1
## sva.15         1      1      1      1      1      1      0      1      1
## sva.16         1      1      1      1      1      1      1      0      1
## sva.17         1      1      1      1      1      1      1      1      0
## sva.18         1      1      1      1      1      1      1      1      1
## sva.19         1      1      1      1      1      1      1      1      1
## sva.20         1      1      1      1      1      1      1      1      1
## sva.21         1      1      1      1      1      1      1      1      1
## sva.22         1      1      1      1      1      1      1      1      1
## sva.23         1      1      1      1      1      1      1      1      1
## sva.24         1      1      1      1      1      1      1      1      1
## sva.25         1      1      1      1      1      1      1      1      1
## sva.26         1      1      1      1      1      1      1      1      1
## sva.27         1      1      1      1      1      1      1      1      1
## sva.28         1      1      1      1      1      1      1      1      1
## sva.29         1      1      1      1      1      1      1      1      1
## sva.30         1      1      1      1      1      1      1      1      1
## sva.31         1      1      1      1      1      1      1      1      1
## sva.32         1      1      1      1      1      1      1      1      1
## sva.33         1      1      1      1      1      1      1      1      1
## sva.34         1      1      1      1      1      1      1      1      1
## sva.35         1      1      1      1      1      1      1      1      1
## sva.36         1      1      1      1      1      1      1      1      1
## sva.37         1      1      1      1      1      1      1      1      1
## sva.38         1      1      1      1      1      1      1      1      1
## sva.39         1      1      1      1      1      1      1      1      1
## sva.40         1      1      1      1      1      1      1      1      1
## sva.41         1      1      1      1      1      1      1      1      1
## sva.42         1      1      1      1      1      1      1      1      1
## sva.43         1      1      1      1      1      1      1      1      1
## sva.44         1      1      1      1      1      1      1      1      1
## sva.45         1      1      1      1      1      1      1      1      1
## sva.46         1      1      1      1      1      1      1      1      1
## sva.47         1      1      1      1      1      1      1      1      1
## sva.48         1      1      1      1      1      1      1      1      1
## sva.49         1      1      1      1      1      1      1      1      1
## sva.50         1      1      1      1      1      1      1      1      1
## sva.51         1      1      1      1      1      1      1      1      1
## sva.52         1      1      1      1      1      1      1      1      1
## sva.53         1      1      1      1      1      1      1      1      1
## sva.54         1      1      1      1      1      1      1      1      1
## sva.55         1      1      1      1      1      1      1      1      1
## sva.56         1      1      1      1      1      1      1      1      1
## sva.57         1      1      1      1      1      1      1      1      1
## sva.58         1      1      1      1      1      1      1      1      1
## sva.59         1      1      1      1      1      1      1      1      1
## sva.60         1      1      1      1      1      1      1      1      1
## sva.61         1      1      1      1      1      1      1      1      1
## sva.62         1      1      1      1      1      1      1      1      1
## sva.63         1      1      1      1      1      1      1      1      1
## cg21572722     1      1      1      1      1      1      1      1      1
##            sva.18 sva.19 sva.20 sva.21 sva.22 sva.23 sva.24 sva.25 sva.26
## smoking         1      1      1      1      1      1      1      1      1
## age             1      1      1      1      1      1      1      1      1
## sex             1      1      1      1      1      1      1      1      1
## sva.1           0      0      1      1      1      1      1      1      1
## sva.2           1      1      0      0      0      0      0      0      0
## sva.3           1      1      1      1      1      1      1      1      1
## sva.4           1      1      1      1      1      1      1      1      1
## sva.5           1      1      1      1      1      1      1      1      1
## sva.6           1      1      1      1      1      1      1      1      1
## sva.7           1      1      1      1      1      1      1      1      1
## sva.8           1      1      1      1      1      1      1      1      1
## sva.9           1      1      1      1      1      1      1      1      1
## sva.10          1      1      1      1      1      1      1      1      1
## sva.11          1      1      1      1      1      1      1      1      1
## sva.12          1      1      1      1      1      1      1      1      1
## sva.13          1      1      1      1      1      1      1      1      1
## sva.14          1      1      1      1      1      1      1      1      1
## sva.15          1      1      1      1      1      1      1      1      1
## sva.16          1      1      1      1      1      1      1      1      1
## sva.17          1      1      1      1      1      1      1      1      1
## sva.18          0      1      1      1      1      1      1      1      1
## sva.19          1      0      1      1      1      1      1      1      1
## sva.20          1      1      0      1      1      1      1      1      1
## sva.21          1      1      1      0      1      1      1      1      1
## sva.22          1      1      1      1      0      1      1      1      1
## sva.23          1      1      1      1      1      0      1      1      1
## sva.24          1      1      1      1      1      1      0      1      1
## sva.25          1      1      1      1      1      1      1      0      1
## sva.26          1      1      1      1      1      1      1      1      0
## sva.27          1      1      1      1      1      1      1      1      1
## sva.28          1      1      1      1      1      1      1      1      1
## sva.29          1      1      1      1      1      1      1      1      1
## sva.30          1      1      1      1      1      1      1      1      1
## sva.31          1      1      1      1      1      1      1      1      1
## sva.32          1      1      1      1      1      1      1      1      1
## sva.33          1      1      1      1      1      1      1      1      1
## sva.34          1      1      1      1      1      1      1      1      1
## sva.35          1      1      1      1      1      1      1      1      1
## sva.36          1      1      1      1      1      1      1      1      1
## sva.37          1      1      1      1      1      1      1      1      1
## sva.38          1      1      1      1      1      1      1      1      1
## sva.39          1      1      1      1      1      1      1      1      1
## sva.40          1      1      1      1      1      1      1      1      1
## sva.41          1      1      1      1      1      1      1      1      1
## sva.42          1      1      1      1      1      1      1      1      1
## sva.43          1      1      1      1      1      1      1      1      1
## sva.44          1      1      1      1      1      1      1      1      1
## sva.45          1      1      1      1      1      1      1      1      1
## sva.46          1      1      1      1      1      1      1      1      1
## sva.47          1      1      1      1      1      1      1      1      1
## sva.48          1      1      1      1      1      1      1      1      1
## sva.49          1      1      1      1      1      1      1      1      1
## sva.50          1      1      1      1      1      1      1      1      1
## sva.51          1      1      1      1      1      1      1      1      1
## sva.52          1      1      1      1      1      1      1      1      1
## sva.53          1      1      1      1      1      1      1      1      1
## sva.54          1      1      1      1      1      1      1      1      1
## sva.55          1      1      1      1      1      1      1      1      1
## sva.56          1      1      1      1      1      1      1      1      1
## sva.57          1      1      1      1      1      1      1      1      1
## sva.58          1      1      1      1      1      1      1      1      1
## sva.59          1      1      1      1      1      1      1      1      1
## sva.60          1      1      1      1      1      1      1      1      1
## sva.61          1      1      1      1      1      1      1      1      1
## sva.62          1      1      1      1      1      1      1      1      1
## sva.63          1      1      1      1      1      1      1      1      1
## cg21572722      1      1      1      1      1      1      1      1      1
##            sva.27 sva.28 sva.29 sva.30 sva.31 sva.32 sva.33 sva.34 sva.35
## smoking         1      1      1      1      1      1      1      1      1
## age             1      1      1      1      1      1      1      1      1
## sex             1      1      1      1      1      1      1      1      1
## sva.1           1      1      1      1      1      1      1      1      1
## sva.2           0      0      0      1      1      1      1      1      1
## sva.3           1      1      1      0      0      0      0      0      0
## sva.4           1      1      1      1      1      1      1      1      1
## sva.5           1      1      1      1      1      1      1      1      1
## sva.6           1      1      1      1      1      1      1      1      1
## sva.7           1      1      1      1      1      1      1      1      1
## sva.8           1      1      1      1      1      1      1      1      1
## sva.9           1      1      1      1      1      1      1      1      1
## sva.10          1      1      1      1      1      1      1      1      1
## sva.11          1      1      1      1      1      1      1      1      1
## sva.12          1      1      1      1      1      1      1      1      1
## sva.13          1      1      1      1      1      1      1      1      1
## sva.14          1      1      1      1      1      1      1      1      1
## sva.15          1      1      1      1      1      1      1      1      1
## sva.16          1      1      1      1      1      1      1      1      1
## sva.17          1      1      1      1      1      1      1      1      1
## sva.18          1      1      1      1      1      1      1      1      1
## sva.19          1      1      1      1      1      1      1      1      1
## sva.20          1      1      1      1      1      1      1      1      1
## sva.21          1      1      1      1      1      1      1      1      1
## sva.22          1      1      1      1      1      1      1      1      1
## sva.23          1      1      1      1      1      1      1      1      1
## sva.24          1      1      1      1      1      1      1      1      1
## sva.25          1      1      1      1      1      1      1      1      1
## sva.26          1      1      1      1      1      1      1      1      1
## sva.27          0      1      1      1      1      1      1      1      1
## sva.28          1      0      1      1      1      1      1      1      1
## sva.29          1      1      0      1      1      1      1      1      1
## sva.30          1      1      1      0      1      1      1      1      1
## sva.31          1      1      1      1      0      1      1      1      1
## sva.32          1      1      1      1      1      0      1      1      1
## sva.33          1      1      1      1      1      1      0      1      1
## sva.34          1      1      1      1      1      1      1      0      1
## sva.35          1      1      1      1      1      1      1      1      0
## sva.36          1      1      1      1      1      1      1      1      1
## sva.37          1      1      1      1      1      1      1      1      1
## sva.38          1      1      1      1      1      1      1      1      1
## sva.39          1      1      1      1      1      1      1      1      1
## sva.40          1      1      1      1      1      1      1      1      1
## sva.41          1      1      1      1      1      1      1      1      1
## sva.42          1      1      1      1      1      1      1      1      1
## sva.43          1      1      1      1      1      1      1      1      1
## sva.44          1      1      1      1      1      1      1      1      1
## sva.45          1      1      1      1      1      1      1      1      1
## sva.46          1      1      1      1      1      1      1      1      1
## sva.47          1      1      1      1      1      1      1      1      1
## sva.48          1      1      1      1      1      1      1      1      1
## sva.49          1      1      1      1      1      1      1      1      1
## sva.50          1      1      1      1      1      1      1      1      1
## sva.51          1      1      1      1      1      1      1      1      1
## sva.52          1      1      1      1      1      1      1      1      1
## sva.53          1      1      1      1      1      1      1      1      1
## sva.54          1      1      1      1      1      1      1      1      1
## sva.55          1      1      1      1      1      1      1      1      1
## sva.56          1      1      1      1      1      1      1      1      1
## sva.57          1      1      1      1      1      1      1      1      1
## sva.58          1      1      1      1      1      1      1      1      1
## sva.59          1      1      1      1      1      1      1      1      1
## sva.60          1      1      1      1      1      1      1      1      1
## sva.61          1      1      1      1      1      1      1      1      1
## sva.62          1      1      1      1      1      1      1      1      1
## sva.63          1      1      1      1      1      1      1      1      1
## cg21572722      1      1      1      1      1      1      1      1      1
##            sva.36 sva.37 sva.38 sva.39 sva.40 sva.41 sva.42 sva.43 sva.44
## smoking         1      1      1      1      1      1      1      1      1
## age             1      1      1      1      1      1      1      1      1
## sex             1      1      1      1      1      1      1      1      1
## sva.1           1      1      1      1      1      1      1      1      1
## sva.2           1      1      1      1      1      1      1      1      1
## sva.3           0      0      0      0      1      1      1      1      1
## sva.4           1      1      1      1      0      0      0      0      0
## sva.5           1      1      1      1      1      1      1      1      1
## sva.6           1      1      1      1      1      1      1      1      1
## sva.7           1      1      1      1      1      1      1      1      1
## sva.8           1      1      1      1      1      1      1      1      1
## sva.9           1      1      1      1      1      1      1      1      1
## sva.10          1      1      1      1      1      1      1      1      1
## sva.11          1      1      1      1      1      1      1      1      1
## sva.12          1      1      1      1      1      1      1      1      1
## sva.13          1      1      1      1      1      1      1      1      1
## sva.14          1      1      1      1      1      1      1      1      1
## sva.15          1      1      1      1      1      1      1      1      1
## sva.16          1      1      1      1      1      1      1      1      1
## sva.17          1      1      1      1      1      1      1      1      1
## sva.18          1      1      1      1      1      1      1      1      1
## sva.19          1      1      1      1      1      1      1      1      1
## sva.20          1      1      1      1      1      1      1      1      1
## sva.21          1      1      1      1      1      1      1      1      1
## sva.22          1      1      1      1      1      1      1      1      1
## sva.23          1      1      1      1      1      1      1      1      1
## sva.24          1      1      1      1      1      1      1      1      1
## sva.25          1      1      1      1      1      1      1      1      1
## sva.26          1      1      1      1      1      1      1      1      1
## sva.27          1      1      1      1      1      1      1      1      1
## sva.28          1      1      1      1      1      1      1      1      1
## sva.29          1      1      1      1      1      1      1      1      1
## sva.30          1      1      1      1      1      1      1      1      1
## sva.31          1      1      1      1      1      1      1      1      1
## sva.32          1      1      1      1      1      1      1      1      1
## sva.33          1      1      1      1      1      1      1      1      1
## sva.34          1      1      1      1      1      1      1      1      1
## sva.35          1      1      1      1      1      1      1      1      1
## sva.36          0      1      1      1      1      1      1      1      1
## sva.37          1      0      1      1      1      1      1      1      1
## sva.38          1      1      0      1      1      1      1      1      1
## sva.39          1      1      1      0      1      1      1      1      1
## sva.40          1      1      1      1      0      1      1      1      1
## sva.41          1      1      1      1      1      0      1      1      1
## sva.42          1      1      1      1      1      1      0      1      1
## sva.43          1      1      1      1      1      1      1      0      1
## sva.44          1      1      1      1      1      1      1      1      0
## sva.45          1      1      1      1      1      1      1      1      1
## sva.46          1      1      1      1      1      1      1      1      1
## sva.47          1      1      1      1      1      1      1      1      1
## sva.48          1      1      1      1      1      1      1      1      1
## sva.49          1      1      1      1      1      1      1      1      1
## sva.50          1      1      1      1      1      1      1      1      1
## sva.51          1      1      1      1      1      1      1      1      1
## sva.52          1      1      1      1      1      1      1      1      1
## sva.53          1      1      1      1      1      1      1      1      1
## sva.54          1      1      1      1      1      1      1      1      1
## sva.55          1      1      1      1      1      1      1      1      1
## sva.56          1      1      1      1      1      1      1      1      1
## sva.57          1      1      1      1      1      1      1      1      1
## sva.58          1      1      1      1      1      1      1      1      1
## sva.59          1      1      1      1      1      1      1      1      1
## sva.60          1      1      1      1      1      1      1      1      1
## sva.61          1      1      1      1      1      1      1      1      1
## sva.62          1      1      1      1      1      1      1      1      1
## sva.63          1      1      1      1      1      1      1      1      1
## cg21572722      1      1      1      1      1      1      1      1      1
##            sva.45 sva.46 sva.47 sva.48 sva.49 sva.50 sva.51 sva.52 sva.53
## smoking         1      1      1      1      1      1      1      1      1
## age             1      1      1      1      1      1      1      1      1
## sex             1      1      1      1      1      1      1      1      1
## sva.1           1      1      1      1      1      1      1      1      1
## sva.2           1      1      1      1      1      1      1      1      1
## sva.3           1      1      1      1      1      1      1      1      1
## sva.4           0      0      0      0      0      1      1      1      1
## sva.5           1      1      1      1      1      0      0      0      0
## sva.6           1      1      1      1      1      1      1      1      1
## sva.7           1      1      1      1      1      1      1      1      1
## sva.8           1      1      1      1      1      1      1      1      1
## sva.9           1      1      1      1      1      1      1      1      1
## sva.10          1      1      1      1      1      1      1      1      1
## sva.11          1      1      1      1      1      1      1      1      1
## sva.12          1      1      1      1      1      1      1      1      1
## sva.13          1      1      1      1      1      1      1      1      1
## sva.14          1      1      1      1      1      1      1      1      1
## sva.15          1      1      1      1      1      1      1      1      1
## sva.16          1      1      1      1      1      1      1      1      1
## sva.17          1      1      1      1      1      1      1      1      1
## sva.18          1      1      1      1      1      1      1      1      1
## sva.19          1      1      1      1      1      1      1      1      1
## sva.20          1      1      1      1      1      1      1      1      1
## sva.21          1      1      1      1      1      1      1      1      1
## sva.22          1      1      1      1      1      1      1      1      1
## sva.23          1      1      1      1      1      1      1      1      1
## sva.24          1      1      1      1      1      1      1      1      1
## sva.25          1      1      1      1      1      1      1      1      1
## sva.26          1      1      1      1      1      1      1      1      1
## sva.27          1      1      1      1      1      1      1      1      1
## sva.28          1      1      1      1      1      1      1      1      1
## sva.29          1      1      1      1      1      1      1      1      1
## sva.30          1      1      1      1      1      1      1      1      1
## sva.31          1      1      1      1      1      1      1      1      1
## sva.32          1      1      1      1      1      1      1      1      1
## sva.33          1      1      1      1      1      1      1      1      1
## sva.34          1      1      1      1      1      1      1      1      1
## sva.35          1      1      1      1      1      1      1      1      1
## sva.36          1      1      1      1      1      1      1      1      1
## sva.37          1      1      1      1      1      1      1      1      1
## sva.38          1      1      1      1      1      1      1      1      1
## sva.39          1      1      1      1      1      1      1      1      1
## sva.40          1      1      1      1      1      1      1      1      1
## sva.41          1      1      1      1      1      1      1      1      1
## sva.42          1      1      1      1      1      1      1      1      1
## sva.43          1      1      1      1      1      1      1      1      1
## sva.44          1      1      1      1      1      1      1      1      1
## sva.45          0      1      1      1      1      1      1      1      1
## sva.46          1      0      1      1      1      1      1      1      1
## sva.47          1      1      0      1      1      1      1      1      1
## sva.48          1      1      1      0      1      1      1      1      1
## sva.49          1      1      1      1      0      1      1      1      1
## sva.50          1      1      1      1      1      0      1      1      1
## sva.51          1      1      1      1      1      1      0      1      1
## sva.52          1      1      1      1      1      1      1      0      1
## sva.53          1      1      1      1      1      1      1      1      0
## sva.54          1      1      1      1      1      1      1      1      1
## sva.55          1      1      1      1      1      1      1      1      1
## sva.56          1      1      1      1      1      1      1      1      1
## sva.57          1      1      1      1      1      1      1      1      1
## sva.58          1      1      1      1      1      1      1      1      1
## sva.59          1      1      1      1      1      1      1      1      1
## sva.60          1      1      1      1      1      1      1      1      1
## sva.61          1      1      1      1      1      1      1      1      1
## sva.62          1      1      1      1      1      1      1      1      1
## sva.63          1      1      1      1      1      1      1      1      1
## cg21572722      1      1      1      1      1      1      1      1      1
##            sva.54 sva.55 sva.56 sva.57 sva.58 sva.59 sva.60 sva.61 sva.62
## smoking         1      1      1      1      1      1      1      1      1
## age             1      1      1      1      1      1      1      1      1
## sex             1      1      1      1      1      1      1      1      1
## sva.1           1      1      1      1      1      1      1      1      1
## sva.2           1      1      1      1      1      1      1      1      1
## sva.3           1      1      1      1      1      1      1      1      1
## sva.4           1      1      1      1      1      1      1      1      1
## sva.5           0      0      0      0      0      0      1      1      1
## sva.6           1      1      1      1      1      1      0      0      0
## sva.7           1      1      1      1      1      1      1      1      1
## sva.8           1      1      1      1      1      1      1      1      1
## sva.9           1      1      1      1      1      1      1      1      1
## sva.10          1      1      1      1      1      1      1      1      1
## sva.11          1      1      1      1      1      1      1      1      1
## sva.12          1      1      1      1      1      1      1      1      1
## sva.13          1      1      1      1      1      1      1      1      1
## sva.14          1      1      1      1      1      1      1      1      1
## sva.15          1      1      1      1      1      1      1      1      1
## sva.16          1      1      1      1      1      1      1      1      1
## sva.17          1      1      1      1      1      1      1      1      1
## sva.18          1      1      1      1      1      1      1      1      1
## sva.19          1      1      1      1      1      1      1      1      1
## sva.20          1      1      1      1      1      1      1      1      1
## sva.21          1      1      1      1      1      1      1      1      1
## sva.22          1      1      1      1      1      1      1      1      1
## sva.23          1      1      1      1      1      1      1      1      1
## sva.24          1      1      1      1      1      1      1      1      1
## sva.25          1      1      1      1      1      1      1      1      1
## sva.26          1      1      1      1      1      1      1      1      1
## sva.27          1      1      1      1      1      1      1      1      1
## sva.28          1      1      1      1      1      1      1      1      1
## sva.29          1      1      1      1      1      1      1      1      1
## sva.30          1      1      1      1      1      1      1      1      1
## sva.31          1      1      1      1      1      1      1      1      1
## sva.32          1      1      1      1      1      1      1      1      1
## sva.33          1      1      1      1      1      1      1      1      1
## sva.34          1      1      1      1      1      1      1      1      1
## sva.35          1      1      1      1      1      1      1      1      1
## sva.36          1      1      1      1      1      1      1      1      1
## sva.37          1      1      1      1      1      1      1      1      1
## sva.38          1      1      1      1      1      1      1      1      1
## sva.39          1      1      1      1      1      1      1      1      1
## sva.40          1      1      1      1      1      1      1      1      1
## sva.41          1      1      1      1      1      1      1      1      1
## sva.42          1      1      1      1      1      1      1      1      1
## sva.43          1      1      1      1      1      1      1      1      1
## sva.44          1      1      1      1      1      1      1      1      1
## sva.45          1      1      1      1      1      1      1      1      1
## sva.46          1      1      1      1      1      1      1      1      1
## sva.47          1      1      1      1      1      1      1      1      1
## sva.48          1      1      1      1      1      1      1      1      1
## sva.49          1      1      1      1      1      1      1      1      1
## sva.50          1      1      1      1      1      1      1      1      1
## sva.51          1      1      1      1      1      1      1      1      1
## sva.52          1      1      1      1      1      1      1      1      1
## sva.53          1      1      1      1      1      1      1      1      1
## sva.54          0      1      1      1      1      1      1      1      1
## sva.55          1      0      1      1      1      1      1      1      1
## sva.56          1      1      0      1      1      1      1      1      1
## sva.57          1      1      1      0      1      1      1      1      1
## sva.58          1      1      1      1      0      1      1      1      1
## sva.59          1      1      1      1      1      0      1      1      1
## sva.60          1      1      1      1      1      1      0      1      1
## sva.61          1      1      1      1      1      1      1      0      1
## sva.62          1      1      1      1      1      1      1      1      0
## sva.63          1      1      1      1      1      1      1      1      1
## cg21572722      1      1      1      1      1      1      1      1      1
##            sva.63 cg21572722
## smoking         1          1
## age             1          1
## sex             1          1
## sva.1           1          1
## sva.2           1          1
## sva.3           1          1
## sva.4           1          1
## sva.5           1          1
## sva.6           0          1
## sva.7           1          1
## sva.8           1          1
## sva.9           1          1
## sva.10          1          1
## sva.11          1          1
## sva.12          1          1
## sva.13          1          1
## sva.14          1          1
## sva.15          1          1
## sva.16          1          1
## sva.17          1          1
## sva.18          1          1
## sva.19          1          1
## sva.20          1          1
## sva.21          1          1
## sva.22          1          1
## sva.23          1          1
## sva.24          1          1
## sva.25          1          1
## sva.26          1          1
## sva.27          1          1
## sva.28          1          1
## sva.29          1          1
## sva.30          1          1
## sva.31          1          1
## sva.32          1          1
## sva.33          1          1
## sva.34          1          1
## sva.35          1          1
## sva.36          1          1
## sva.37          1          1
## sva.38          1          1
## sva.39          1          1
## sva.40          1          1
## sva.41          1          1
## sva.42          1          1
## sva.43          1          1
## sva.44          1          1
## sva.45          1          1
## sva.46          1          1
## sva.47          1          1
## sva.48          1          1
## sva.49          1          1
## sva.50          1          1
## sva.51          1          1
## sva.52          1          1
## sva.53          1          1
## sva.54          1          1
## sva.55          1          1
## sva.56          1          1
## sva.57          1          1
## sva.58          1          1
## sva.59          1          1
## sva.60          1          1
## sva.61          1          1
## sva.62          1          1
## sva.63          0          1
## cg21572722      1          0
```
A `1` indicates that the variable in the column is used to impute values for the
variable in the row.
Notice how there are 1's in the rows only of variables with missing values.

Within each imputed dataset, we can calculate the correlation of age with
the age variable in the full dataset.

```r
unlist(with(imp, cor(age, samples$age))$analyses)
```

```
##  [1] 0.9071685 0.9219443 0.9295520 0.9183630 0.9351234 0.9184417 0.9322300
##  [8] 0.9113708 0.9274307 0.9426592
```

```r
quantile(unlist(with(imp, cor(age, samples$age))$analyses))
```

```
##        0%       25%       50%       75%      100% 
## 0.9071685 0.9183826 0.9246875 0.9315605 0.9426592
```

Here we fit a linear model in each imputed dataset.

```r
fit <- with(imp, lm(cg21572722 ~ age + smoking + sex + sva.1 + sva.2 + sva.3 + sva.4))
```

Therefore we have 10 model fits. We show the summary statistics
for the first two.

```r
length(fit$analyses)
```

```
## [1] 10
```

```r
coef(summary(fit$analyses[[1]]))
```

```
##                  Estimate   Std. Error    t value     Pr(>|t|)
## (Intercept)   0.374550317 0.0416368941  8.9956354 2.205903e-16
## age           0.001637678 0.0001910164  8.5734920 3.223878e-15
## smokingnever -0.006255448 0.0042312569 -1.4783899 1.409336e-01
## sexmale      -0.053938452 0.0567221826 -0.9509234 3.428326e-01
## sva.1         0.560810000 0.5601839960  1.0011175 3.180235e-01
## sva.2         0.224786211 0.0315122320  7.1333002 1.922739e-11
## sva.3         0.005782113 0.0334017516  0.1731081 8.627479e-01
## sva.4         0.032992510 0.0296300678  1.1134808 2.668868e-01
```

```r
coef(summary(fit$analyses[[2]]))
```

```
##                  Estimate   Std. Error    t value     Pr(>|t|)
## (Intercept)   0.334689244 0.0432913974  7.7310797 5.729764e-13
## age           0.001725174 0.0001953258  8.8322896 6.265707e-16
## smokingnever -0.001885681 0.0043567044 -0.4328229 6.656267e-01
## sexmale      -0.008606018 0.0587175111 -0.1465665 8.836272e-01
## sva.1         0.097750579 0.5798431303  0.1685811 8.663027e-01
## sva.2         0.240467200 0.0325920399  7.3780960 4.643845e-12
## sva.3         0.015012083 0.0345397712  0.4346318 6.643153e-01
## sva.4         0.012306293 0.0306373860  0.4016757 6.883671e-01
```

To obtain summary statistics across all model fits,
we 'pool' the summary statistics of each.

```r
pool.fit <- pool(fit)
```

The pooled summary statistics for 'age':

```r
summary(pool.fit)["age",]
```

```
##        estimate    std.error statistic       df      p.value
## age 0.001594739 0.0002283402  6.984048 75.45443 7.469492e-11
```

The complete case association in this case is quite similar.

```r
fit.cc <- lm(cg21572722 ~ age + smoking + sex + sva.1 + sva.2 + sva.3 + sva.4, samples.mis)
coef(summary(fit.cc))["age",]
```

```
##     Estimate   Std. Error      t value     Pr(>|t|) 
## 1.767848e-03 2.322437e-04 7.612038e+00 4.865750e-12
```

What if we consider another CpG site known to be associated with age but
was not included in the imputation?

```r
fit <- with(imp, lm(meth["cg27193080",] ~ age + smoking + sex + sva.1 + sva.2 + sva.3 + sva.4))
summary(pool(fit))["age",]
```

```
##         estimate    std.error statistic       df     p.value
## age -0.001348686 0.0003573995 -3.773609 103.8462 0.000214735
```

The complete case is much stronger ...

```r
fit.cc <- lm(meth["cg27193080",] ~ age + smoking + sex + sva.1 + sva.2 + sva.3 + sva.4, samples.mis)
coef(summary(fit.cc))["age",]
```

```
##      Estimate    Std. Error       t value      Pr(>|t|) 
## -1.629459e-03  3.461588e-04 -4.707258e+00  5.334813e-06
```

## Imputation for EWAS

### Option 1: Separate CpGs

1. Identify available variables known to be associated with each covariate with missing values.
2. For each CpG site:
2a. Impute N datasets including the CpG site methylation levels, all additional variables and covariates.
2b. Within each dataset fit a model to test association of the CpG site with variable of interest.
2c. Pool model fits to obtain association summary statistics (e.g. estimate, se, p-value).

In other words, for each CpG site `cg`:

```r
data <- ... ## var, additional variables, covariates
for (cg in cgs) {
    data$cg <- meth[cg,]
    imp <- mice(data, ...)
    fit <- with(imp, lm(cg ~ var + score1 + score2 + ... + var1 + var2 + ... + cov1 + cov2 + ...))
    fit <- pool(fit)
    summary(fit)["var",] ## save summary statistics
}
```

Problem: Running time is very long.

### Option 2: Using associated CpGs – Naïve Method

1. To reduce running time, impute N datasets one time as above but
   excluding methylation levels of CpG sites (except for all those
   associated with the variable of interest).
2. For each CpG site, fit a model with each imputed dataset and pool model fits.


```r
data <- ... ## var, additional variables, covariates, all CpG sites associated with var
imp <- mice(data, ...)
for (cg in cgs) {
    data$cg <- meth[cg,]
    fit <- with(imp, lm(cg ~ var + score1 + score2 + ... + var1 + var2 + ... + cov1 + cov2 + ...))
    fit <- pool(fit)
    summary(fit)["var",] ## save summary statistics    
}
```

Problem: Models can be very large and association testing
for CpG sites not included in the imputation are biased to the null.

### Option 3: Using associated CpGs – Wu Method

As above, but instead of just including all associated CpG sites,
construct a methylation predictor of the variable of interest
and include only those in the imputation model.

Wu et al. generated a predictor using forward-stepwise selection
from the set of 100 CpG most strongly associated with the
variable of interest in the complete-case dataset.

Problem: Also tends to bias CpG sites not included in the imputation
to null but not as badly as the Naive method.

### Option 4: Random bins of fixed size always including associated CpGs – Wu bins

Construct bins of CpG sites of size B such that:
1. Each CpG site in the Wu predictor appears in every bin, and
2. Each other CpG site appears in exactly one bin.

For each bin:
1. Impute N datasets with all CpG sites and all variables and covariates.
2. For each CpG site in the bin (except We predictor sites), fit a
   model for each dataset, pool the model fits and save summary
   statistics.

Bin size B is selected so that there is a 3:1 or 10:1 ratio of samples to
variables (i.e. variables, covariates, CpG sites) in the imputation models.

### Option 5: Random bins of fixed size

As above except that no predictors are generated
so each CpG site is randomly assigned to a single bin.

## R package for imputation in EWAS

Coming soon to the `ewaff` R package.  Here is what an imputation EWAS will be executed.
We'll use the random bins method.

```r
imp <- ewaff.impute(methylation ~ age + smoking + sex, data, methylation, isv=TRUE)
```

Here we print out the Bonferroni-adjusted associations.

```r
imp$table[which(imp$table$p.adjust < 0.05),]
```

