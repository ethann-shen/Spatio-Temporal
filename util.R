normalize_weights = function(w) {
  diag(w) = 0
  rs = rowSums(w)
  rs[rs == 0] = 1
  w/rs
}

strip_attrs = function(obj)
{
  attributes(obj) = NULL
  obj
}

strip_class = function(obj)
{
  attr(obj,"class") = NULL
  obj
}

fixterm = function(term, i) {
  ifelse(is.na(i), term, paste0(term,"[",i,"]"))
}


inv_logit = function(x) exp(x) / (1+exp(x))
logit = function(p) log(p/(1-p))

get_coda_parameter = function(coda, pattern)
{
  w = coda[[1]] %>% colnames() %>% stringr::str_detect(pattern)
  coda[[1]][,w,drop=FALSE]
}

post_summary = function(m, ci_width=0.95) {
  d = data_frame(
    post_mean  = apply(m, 2, mean),
    post_med   = apply(m, 2, median),
    post_lower = apply(m, 2, quantile, probs=(1-ci_width)/2),
    post_upper = apply(m, 2, quantile, probs=1 - (1-ci_width)/2)
  )
  
  if (!is.null(colnames(m)))
    d = d %>% mutate(param = colnames(m)) %>% select(param,post_mean:post_upper)
  
  d
}

nugget_cov   = function(d, sigma2=1) { ifelse(d==0, sigma2, 0) }
exp_cov      = function(d, sigma2=1, l=1, sigma2_w=0) { sigma2 * exp(-abs(d)*l) + nugget_cov(d,sigma2_w) }
sq_exp_cov   = function(d, sigma2=1, l=1, sigma2_w=0) { sigma2 * exp(-(abs(d)*l)^2) + nugget_cov(d,sigma2_w) }
pow_exp_cov  = function(d, sigma2=1, l=1, p=1.5) { sigma2 * exp(-(abs(d)*l)^p) }
rquad_cov    = function(d, sigma2=1, l=1, a=1) { sigma2 * (1+d^2*l^2/a)^(-a) }
periodic_cov = function(d, sigma2=1, l=1, p=1) { sigma2 * exp(-2*l^2*sin(pi*d/p)^2) }
matern_cov   = function(d, sigma2=1, l=1, nu=1/2) { fields::Matern(d, alpha=l, nu=nu, phi=sigma2) }
sphere_cov   = function(d, sigma2=1, l=1) { ifelse(d > 1/l, 0, sigma2*(1 - 1.5*d*l + 0.5*(d*l)^3)) }

exp_sv     = function(d, sigma2=1, l=1, sigma2_w=0)      { sigma2 + sigma2_w - exp_cov(d,sigma2,l) - nugget_cov(d,sigma2_w) }
sq_exp_sv  = function(d, sigma2=1, l=1, sigma2_w=0)      { sigma2 + sigma2_w - sq_exp_cov(d,sigma2,l) - nugget_cov(d,sigma2_w) }
pow_exp_sv = function(d, sigma2=1, l=1, p=2, sigma2_w=0) { sigma2 + sigma2_w - pow_exp_cov(d,sigma2,l,p) - nugget_cov(d,sigma2_w) }

linear_cov   = function(x1, x2, sigma2_b=1, sigma2_v=1, c=0)
{
  if (!missing(x2))
  {
    sigma2_b + sigma2_v * (x1-c) * (x2-c)
  } else {
    expand.grid(t_i=x1, t_j=x1) %>%
      {linear_cov(.[[1]], .[[2]], sigma2_b=sigma2_b, sigma2_v=sigma2_v, c=c)} %>% 
      matrix(ncol=length(x1))
  }
}



# From https://github.com/hadley/dplyr/issues/271
bin = function(df, var, binwidth, origin = NULL) {
  n = nrow(df)
  
  var = as.character(substitute(var))
  x = df[[var]]
  
  if (is.null(origin)) {
    origin = min(x)
  }
  
  bin = (x - origin) %/% binwidth
  indices = unname(split(seq_len(n) - 1, bin))
  
  mid = origin + (0:(max(bin)+1)) * binwidth + binwidth/2
  
  df[["bin_mid"]] = mid[bin+1]
  
  attr(df, "indices") = indices
  attr(df, "drop") = FALSE
  attr(df, "group_sizes") = sapply(indices, length)
  attr(df, "biggest_group_size") = max(attr(df, "group_sizes"))
  attr(df, "labels") = data.frame(bin = seq_along(indices))
  attr(df, "vars") = list(quote(bin))
  class(df) = c("grouped_df", "tbl_df", "tbl", "data.frame")
  
  df
}



dist_long = function(d)
{
  d = as.matrix(d)
  d[upper.tri(d, diag = TRUE)] = NA
  
  data.frame(
    expand.grid(i=1:nrow(d), j=1:nrow(d)),
    c(d)
  ) %>%
    setNames(c("i","j","dist")) %>%
    filter(!is.na(dist))
}

emp_semivariogram = function(d, y, x, bin=FALSE, binwidth, range_max)
{
  y_col = as.character(substitute(y))
  x_col = as.character(substitute(x))
  
  d = d[[x_col]] %>%
    dist() %>% 
    dist_long() %>%
    mutate(y_i = d[[y_col]][i], y_j = d[[y_col]][j])
  
  
  if (bin)
  {
    d = d %>% bin(dist, binwidth = binwidth)
  } else {
    d = d %>% mutate(bin_mid = dist) %>% rowwise()
  }
  
  d = d %>%
    group_by(bin_mid) %>% 
    summarize(
      gamma = sum( (y_i - y_j)^2 / (2*n()) ),
      h = mean(bin_mid),
      n = n()
    )
  
  if (!missing(range_max))
    d = d %>% filter(h < range_max)
  
  d
}


rmvnorm = function(n, mu=rep(0, nrow(sigma)), sigma, diag_adjust = 1e-6)
{
  diag(sigma) = diag(sigma) + diag_adjust
  mu %*% matrix(1, ncol=n) + t(chol(sigma)) %*% matrix(rnorm(n*nrow(sigma)), ncol=n)
}

clean_spdynlm = function(m, start, end, thin=1) {
  m$p.beta.0.samples = window(m$p.beta.0.samples, start, end, thin)
  m$p.beta.samples = window(m$p.beta.samples, start, end, thin)
  m$p.sigma.eta.samples = window(m$p.sigma.eta.samples, start, end, thin)
  m$p.theta.samples = window(m$p.theta.samples, start, end, thin)
  m$p.y.samples = t(m$p.y.samples) %>% coda::as.mcmc() %>% window(start, end, thin) %>% t()
  m$p.u.samples = NULL
  
  m
}