---
title: "Lecture 10 - Forecasting and Fitting ARIMA Models"
author: "Ethan Shen"
date: "5/17/2020"
output: 
  html_document:
    keep_md: TRUE
---



# AR(1) Example

$\phi = 0.75$, $\delta=0.5$, and $\sigma_w^2=1$,


```r
ar1 = arima.sim(n=1000, model = list(order=c(1,0,0), ar=0.75), mean=0.5)
forecast::ggtsdisplay(ar1)
```

![](Lec10-Forecasting_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

### Using ARIMA function


```r
ar1_arima = forecast::Arima(ar1, order = c(1,0,0)) 
summary(ar1_arima)
```

```
## Series: ar1 
## ARIMA(1,0,0) with non-zero mean 
## 
## Coefficients:
##          ar1    mean
##       0.7514  1.7887
## s.e.  0.0209  0.1300
## 
## sigma^2 estimated as 1.054:  log likelihood=-1444.47
## AIC=2894.93   AICc=2894.95   BIC=2909.65
## 
## Training set error measures:
##                        ME     RMSE       MAE       MPE     MAPE      MASE
## Training set -0.002122881 1.025429 0.8274295 -4293.926 4519.452 0.9468485
##                      ACF1
## Training set -0.007473343
```

The "mean" reported by the `ARIMA` model is $E(y_t) = \frac{\delta}{1-\phi}$.

### Using lm() function


```r
d = data_frame(y = ar1, t=seq_along(ar1))
ar1_lm = lm(y~lag(as.vector(y)), data=d)
summary(ar1_lm)
```

```
## 
## Call:
## lm(formula = y ~ lag(as.vector(y)), data = d)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -3.11107 -0.69229  0.01323  0.69632  2.70814 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)    
## (Intercept)        0.44204    0.04972   8.891   <2e-16 ***
## lag(as.vector(y))  0.75097    0.02097  35.809   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.026 on 997 degrees of freedom
##   (1 observation deleted due to missingness)
## Multiple R-squared:  0.5626,	Adjusted R-squared:  0.5621 
## F-statistic:  1282 on 1 and 997 DF,  p-value: < 2.2e-16
```

We see the coefficients for `ar1` and `lag(y)` are very very similar. 


# Bayesian AR(1) Model 

Model follows decomposition of multivariate normal density into product of conditional densities.


```r
ar1_model = "model{
# likelihood
  y[1] ~ dnorm(delta/(1-phi), (sigma2w/(1-phi^2))^-1) # taking inverse of variance bc JAGS takes precision, not variance
   
   # posterior predictive
  y_hat[1] ~ dnorm(delta/(1-phi), (sigma2w/(1-phi^2))^-1)

 # conditional likelihoods
  for (t in 2:length(y)) {
    y[t] ~ dnorm(delta + phi*y[t-1], 1/sigma2w)
    y_hat[t] ~ dnorm(delta + phi*y[t-1], 1/sigma2w)
  }
  
  # expected value of y_t
  mu = delta/(1-phi)

# priors
  delta ~ dnorm(0,1/1000)
   # could be problematic bc model might not be stationary (in this case, we know the model is stationary, so not an issue)
  # however, if we don't know whether or not model is stationary, the mean and variance of y_t is not what we coded, because those are derived after assuming process is stationary
  
  phi ~ dnorm(0,1)
  tau ~ dgamma(0.001,0.001)
  sigma2w <- 1/tau
}"

m_ar1 = rjags::jags.model(textConnection(ar1_model),
                          data = list(y=ar1))
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 1000
##    Unobserved stochastic nodes: 1003
##    Total graph size: 4016
## 
## Initializing model
```

```r
update(m_ar1, n.iter = 1000)

samp_ar1 = rjags::coda.samples(m_ar1, n.iter = 5000,
                               variable.names = c("delta", "phi", "sigma2w", "mu", "y_hat"))
```

### MCMC Diagnostics


```r
tidybayes::gather_draws(samp_ar1, delta, phi, sigma2w) %>%
  ggplot(aes(x = .iteration, y = .value, color = as.factor(.variable))) + 
  geom_line(alpha = 0.8) + 
  facet_wrap(~.variable, ncol=1, scales = "free_y") + 
  labs(color = "Parameter")
```

![](Lec10-Forecasting_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

## Posterior Density & Comparison to ARIMA fit and lm() fit


```r
ar1_res = bind_rows(
  data_frame(
    model = "truth",
    term = c("delta", "phi", "sigma2_w"), 
    .value = c(0.5, 0.75, 1)
  ),
  data_frame(
    model = "lm",
    term = c("delta", "phi", "sigma2_w"), 
    .value = c(coef(ar1_lm), var(ar1_lm$residuals))
  ),
  data_frame(
    model = "ARIMA",
    term = c("delta", "phi", "sigma2_w"), 
    .value = c(ar1_arima$model$phi, (1-ar1_arima$model$phi)*ar1_arima$coef[2], ar1_arima$sigma2)
  )
)

tidybayes::gather_draws(samp_ar1, delta, phi, sigma2w) %>%
  ggplot(aes(x=.value)) +
  geom_density(fill="lightgrey") +
  geom_vline(data=ar1_res, aes(xintercept = .value, linetype=model, color=model), size=1.5, alpha=0.75) +
  facet_wrap(~term, ncol=3, scales = "free_x")
```

![](Lec10-Forecasting_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

We see that $\sigma^2$ is predicted pretty well, others not so much. 

## Predictions


```r
models = d %>% 
  modelr::add_predictions(ar1_lm,var = "lm") %>%
  mutate(
    ARIMA = ar1_arima %>% fitted(),
    bayes = tidybayes::gather_draws(samp_ar1, y_hat[i]) %>%  summarize(post_mean = mean(.value)) %>% pull(post_mean)
  ) %>%
  tidyr::gather(model, y, -t)

gridExtra::grid.arrange(
  models %>%
    filter(model=="y" | model=="lm") %>%
    ggplot(aes(x=t, y=y, color=forcats::as_factor(model))) +
    geom_line(alpha=0.75) +
    scale_color_manual(values=c("black", "red")) +
    labs(color="Model")
  ,
  models %>%
    filter(model=="y" | model=="bayes") %>%
    ggplot(aes(x=t, y=y, color=forcats::as_factor(model))) +
    geom_line(alpha=0.75) +
    scale_color_manual(values=c("black", "light blue")) +
    labs(color="Model")
  ,
  models %>%
    filter(model=="y" | model=="ARIMA") %>%
    ggplot(aes(x=t, y=y, color=forcats::as_factor(model))) +
    geom_line(alpha=0.75) +
    scale_color_manual(values=c("black", "light green")) +
    labs(color="Model")
)
```

![](Lec10-Forecasting_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

We see the predictions (ARIMA, lm, and Bayesian) are all very similar to the truth.

# ARMA(2,2)

$\phi = (1.3,-0.5)$, $\theta = (0.5,0.2)$, $\delta=0$, and $\sigma_w^2=1$ using the same models 


```r
y = arima.sim(n=500, model=list(ar=c(1.3,-0.5), ma=c(0.5,0.2))) 

forecast::ggtsdisplay(y, points = FALSE)
```

![](Lec10-Forecasting_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

PACF shows MA(2), ACF implies AR of some sort. 

## ARIMA


```r
forecast::Arima(y, order = c(2,0,2), include.mean = FALSE) %>% summary()
```

```
## Series: y 
## ARIMA(2,0,2) with zero mean 
## 
## Coefficients:
##          ar1      ar2     ma1     ma2
##       1.2791  -0.4555  0.4522  0.1921
## s.e.  0.0858   0.0785  0.0902  0.0736
## 
## sigma^2 estimated as 1:  log likelihood=-709.29
## AIC=1428.57   AICc=1428.69   BIC=1449.65
## 
## Training set error measures:
##                       ME      RMSE       MAE       MPE     MAPE      MASE
## Training set 0.004646019 0.9962157 0.7979356 -47.45931 124.1918 0.6629925
##                      ACF1
## Training set -0.005743304
```

$\theta$ generally harder to estimate, but all estimates are relatively close to true values. 

## AR only lm


```r
lm(y ~ lag(as.vector(y),1) + lag(as.vector(y),2)) %>% summary()
```

```
## 
## Call:
## lm(formula = y ~ lag(as.vector(y), 1) + lag(as.vector(y), 2))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -3.3314 -0.7340 -0.0014  0.6936  3.3618 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           0.01260    0.04597   0.274    0.784    
## lag(as.vector(y), 1)  1.56194    0.03187  49.008   <2e-16 ***
## lag(as.vector(y), 2) -0.70395    0.03178 -22.153   <2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.026 on 495 degrees of freedom
##   (2 observations deleted due to missingness)
## Multiple R-squared:  0.9193,	Adjusted R-squared:  0.9189 
## F-statistic:  2818 on 2 and 495 DF,  p-value: < 2.2e-16
```

These are estimates of $\phi_1$ and $\phi_2$. However, we can't do `lag(w_t)` in a simple `lm` model. 

# Hannan-Rissanen Algorithm 

## Steps 1 and 2


```r
ar_yw=ar.yw(y, order.max=10)
ar_yw
```

```
## 
## Call:
## ar.yw.default(x = y, order.max = 10)
## 
## Coefficients:
##       1        2        3  
##  1.5581  -0.7799   0.0867  
## 
## Order selected 3  sigma^2 estimated as  1.299
```

```r
ar_mle=ar.mle(y, order.max=10)
```

```
## Warning in arima0(x, order = c(i, 0L, 0L), include.mean = demean): possible
## convergence problem: optim gave code = 1
```

```r
ar_mle
```

```
## 
## Call:
## ar.mle(x = y, order.max = 10)
## 
## Coefficients:
##       1        2        3        4        5  
##  1.7326  -1.0692   0.1797   0.1247  -0.0886  
## 
## Order selected 5  sigma^2 estimated as  0.9886
```

## Step 3


```r
d = data_frame(
  y = y, 
  # create column of residuals 
  w_hat1 = ar_mle$resid 
)

# then take the lag of residuals
(lm1 = lm(y ~ lag(as.vector(y),1) + lag(as.vector(y),2) + lag(as.vector(w_hat1),1) + lag(as.vector(w_hat1),2), data=d)) %>%
  summary()
```

```
## 
## Call:
## lm(formula = y ~ lag(as.vector(y), 1) + lag(as.vector(y), 2) + 
##     lag(as.vector(w_hat1), 1) + lag(as.vector(w_hat1), 2), data = d)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -3.15879 -0.70605 -0.00185  0.66505  2.97365 
## 
## Coefficients:
##                           Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                0.01205    0.04531   0.266     0.79    
## lag(as.vector(y), 1)       1.31249    0.07155  18.345  < 2e-16 ***
## lag(as.vector(y), 2)      -0.48682    0.06272  -7.762 4.94e-14 ***
## lag(as.vector(w_hat1), 1)  0.40900    0.08561   4.777 2.35e-06 ***
## lag(as.vector(w_hat1), 2)  0.15529    0.07903   1.965     0.05 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.005 on 488 degrees of freedom
##   (7 observations deleted due to missingness)
## Multiple R-squared:  0.9226,	Adjusted R-squared:  0.922 
## F-statistic:  1454 on 4 and 488 DF,  p-value: < 2.2e-16
```

There is uncertainty on the residuals, so we refit by taking residuals of `lm1` and using those instead. 

## Step 4.1


```r
d = modelr::add_residuals(d,lm1,"w_hat2")

(lm2 = lm(y ~ lag(as.vector(y),1) + lag(as.vector(y),2) + lag(as.vector(w_hat2),1) + lag(as.vector(w_hat2),2), data=d)) %>%
  summary()
```

```
## 
## Call:
## lm(formula = y ~ lag(as.vector(y), 1) + lag(as.vector(y), 2) + 
##     lag(as.vector(w_hat2), 1) + lag(as.vector(w_hat2), 2), data = d)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -3.12173 -0.72661 -0.02084  0.66228  3.04655 
## 
## Coefficients:
##                           Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                0.02077    0.04544   0.457   0.6479    
## lag(as.vector(y), 1)       1.31474    0.07140  18.413  < 2e-16 ***
## lag(as.vector(y), 2)      -0.49038    0.06224  -7.879 2.18e-14 ***
## lag(as.vector(w_hat2), 1)  0.40787    0.08490   4.804 2.08e-06 ***
## lag(as.vector(w_hat2), 2)  0.14593    0.07848   1.859   0.0636 .  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.006 on 486 degrees of freedom
##   (9 observations deleted due to missingness)
## Multiple R-squared:  0.9226,	Adjusted R-squared:  0.922 
## F-statistic:  1448 on 4 and 486 DF,  p-value: < 2.2e-16
```

## Step 4.2


```r
d = modelr::add_residuals(d,lm2,"w_hat3")

(lm3 = lm(y ~ lag(as.vector(y),1) + lag(as.vector(y),2) + lag(as.vector(w_hat3),1) + lag(as.vector(w_hat3),2), data=d)) %>%
  summary()
```

```
## 
## Call:
## lm(formula = y ~ lag(as.vector(y), 1) + lag(as.vector(y), 2) + 
##     lag(as.vector(w_hat3), 1) + lag(as.vector(w_hat3), 2), data = d)
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -3.09934 -0.73640 -0.00399  0.66967  3.00808 
## 
## Coefficients:
##                           Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                0.01987    0.04556   0.436   0.6630    
## lag(as.vector(y), 1)       1.30362    0.07169  18.185  < 2e-16 ***
## lag(as.vector(y), 2)      -0.48103    0.06245  -7.703 7.58e-14 ***
## lag(as.vector(w_hat3), 1)  0.42080    0.08530   4.933 1.11e-06 ***
## lag(as.vector(w_hat3), 2)  0.16010    0.07897   2.027   0.0432 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 1.006 on 484 degrees of freedom
##   (11 observations deleted due to missingness)
## Multiple R-squared:  0.9228,	Adjusted R-squared:  0.9221 
## F-statistic:  1446 on 4 and 484 DF,  p-value: < 2.2e-16
```

## RMSEs


```r
modelr::rmse(lm1, data = d)
```

```
## [1] 0.9996785
```

```r
modelr::rmse(lm2, data = d)
```

```
## [1] 1.000582
```

```r
modelr::rmse(lm3, data = d)
```

```
## [1] 1.000943
```

Once coefficients get "close" to each other, we can say they've converged. We are done. 

# Bayesian ARMA(2,2)

This model looks like $y_t = \phi_1 y_{t-1} + \phi_2 y_{t-2} + w_t + \theta_1 w-{t-1} + \theta_2 w-{t-2} + e_t$, where $e_t \sim N(0, \sigma_e ^2)$. $e_t$ is NOT autoregressive, so this makes we are fitting a slightly different model compared to the `lm` model, but very similar to the Hannan Rissanen. The set of residuals that are left over (are not part of the $w_t$) are similar to the $e_t$. 


```r
arma22_model = "model{
# Likelihood
  for (t in 1:length(y)) {
    y[t] ~ dnorm(mu[t], 1/sigma2_e)
  }                    
  
  # estimating first two steps with latent errors 

  mu[1] = phi[1] * y_0  + phi[2] * y_n1 + w[1] + theta[1]*w_0  + theta[2]*w_n1
  mu[2] = phi[1] * y[1] + phi[2] * y_0  + w[2] + theta[1]*w[1] + theta[2]*w_0   
  for (t in 3:length(y)) { 
    mu[t] = phi[1] * y[t-1] + phi[2] * y[t-2] + w[t] + theta[1] * w[t-1] + theta[2] * w[t-2]
  }
  
# Priors

# similar to HR algorithm, we are estimating w_t
  for(t in 1:length(y)){
    w[t] ~ dnorm(0,1/sigma2_w)
  }

  sigma2_w = 1/tau_w; tau_w ~ dgamma(0.001, 0.001) 
  sigma2_e = 1/tau_e; tau_e ~ dgamma(0.001, 0.001) 
  for(i in 1:2) {
    phi[i] ~ dnorm(0,1)
    theta[i] ~ dnorm(0,1)
  }

# Latent errors and series values
  w_0  ~ dt(0,tau_w,2)
  w_n1 ~ dt(0,tau_w,2)
  y_0  ~ dnorm(0,1/1000)
  y_n1 ~ dnorm(0,1/1000)
}"
```

## Bayesian Fit


```r
if (!file.exists("arma22_coda.rds")) {
  m = rjags::jags.model(
    textConnection(arma22_model), 
    data = list(y = y), 
    quiet = TRUE
  ) 
  
  update(m, n.iter=50000, progress.bar="none")
  
  arma22_coda = rjags::coda.samples(
    m, variable.names=c("phi", "theta", "sigma2_w"), 
    n.iter=50000, progress.bar="none", thin = 50
  )
  
  saveRDS(arma22_coda, "arma22_coda.rds")
} else {
  arma22_coda = readRDS("arma22_coda.rds")
}
```


```r
arma22_res = bind_rows(
  data_frame(
    model = "True",
    term = c("phi[1]", "phi[2]", "theta[1]", "theta[2]"), 
    estimate = c(1.3, -0.5, 0.5, 0.2)
  ),
  data_frame(
    model = "Arima",
    term = c("phi[1]", "phi[2]", "theta[1]", "theta[2]"), 
    estimate = c(1.3171,  -0.5142,  0.4332,  0.1651)
  ),
  data_frame(
    model = "HR",
    term = c("phi[1]", "phi[2]", "theta[1]", "theta[2]"), 
    estimate = c(1.34414, -0.54288, 0.38065, 0.11431)
  )
)


arma22_params = tidybayes::gather_samples(arma22_coda, phi[i], theta[i]) %>%
  ungroup() %>%
  mutate(term = paste0(term,"[",i,"]"))
```


### MCMC Diagnostics 


```r
arma22_params %>%
  group_by(term) %>%
  slice(seq(1,n(),n()/200)) %>%
  ggplot(aes(x=.iteration, y=estimate, color=term)) +
    geom_line() +
    facet_grid(term~., scales = "free_y")
```

![](Lec10-Forecasting_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

These are pretty hard to sample with JAGS. 

### Posterior Densities


```r
arma22_params %>%
  ggplot(aes(x=estimate)) +
    geom_density(fill="lightgrey") +
    geom_vline(data=arma22_res, aes(xintercept = estimate, color=forcats::as_factor(model), linetype=forcats::as_factor(model)), size=1, alpha=0.5) +
    facet_wrap(~term, ncol=2, scales = "free")
```

![](Lec10-Forecasting_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

All these predictions (ARIMA vs. HR) are all pretty similar to each other. Only Bayesian posterior that is struggling is $\theta_1$. 
