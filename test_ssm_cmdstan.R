# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# library(cmdstanr)
# check_cmdstan_toolchain()
# install_cmdstan()

library(tidyverse)
library(cmdstanr)
library(bayesplot)


# generate random variables -----------------------------------------------

n_state = 2
#G_sea = diag(10) %>%
#  cbind(., matrix(rep(0, 10), nrow=10)) %>%
#  rbind(matrix(rep(-1, 11), nrow=1), .)
#G_sta = diag(n_state)
#G = cbind(G_sta, matrix(0, nrow=n_state, ncol=11)) %>%
#  rbind(., cbind(matrix(0, nrow=11, ncol=n_state), G_sea))
m_G = diag(n_state)

#F1 = matrix(c(0.3, -0.5, 1.0, rep(0, 10)), nrow=1)
#F2 = matrix(c(-0.4, 0.8, 1.0, rep(0, 10)), nrow=1)
#F3 = matrix(c(0.1, 0.3, 1.0, rep(0, 10)), nrow=1)
#F_array = array(dim=c(3, 1, n_state+11))
#F_array[1, , ] = F1
#F_array[2, , ] = F2
#F_array[3, , ] = F3
m_F = matrix(c(0.2, -0.6, 0.4, -0.4, 0.1, 0.9), ncol=n_state, byrow=T)
m_B = matrix(c(0.1, 0.2, 0.25, 0.3, 0.22, 0.15, 0.08, 0.02, -0.04, -0.02, -0.01,
               -0.4, -0.5, -0.4, -0.3, -0.2, 0.0, 0.1, 0.15, 0.23, 0.16, 0.1,
               0.2, 0.3, 0.35, 0.18, 0.12, 0.04, -0.04, -1.6, -2.3, -1.5, -0.4),
             nrow=3,
             byrow=T)
m_mu = matrix(c(2, -0.1, -3), nrow=3)

m_d = tibble(month=rep(1:12, 5), one=1) %>%
  mutate(id=row_number()) %>%
  pivot_wider(names_from="month", values_from="one", values_fill=0) %>%
  select(-c(id, `1`)) %>%
  as.matrix()

v_sd_x = c(5, 4)
v_sd_y = c(3, 1, 1.4)

Time = nrow(m_d)
x = matrix(nrow=Time, ncol=n_state)
y = matrix(nrow=Time, ncol=3)

set.seed(42)
for(t in 1:Time){
  if(t==1){
    x[1, ] = rnorm(n_state, sd=v_sd_x)
    y[1, ] = m_mu + m_F%*%x[1, ] + m_B%*%m_d[1, ] + rnorm(3, sd=v_sd_y)
  }else{
    x[t, ] = m_G%*%x[t-1, ] + rnorm(n_state, sd=v_sd_x)
    y[t, ] = m_mu + m_F%*%x[t, ] + m_B%*%m_d[t, ] + rnorm(3, sd=v_sd_y)
  }
}

df_x = as_tibble(x, .name_repair="unique") %>%
  rename_at(vars(contains("...")), ~str_replace(., "...", "x")) %>%
  mutate(time=row_number()) %>%
  relocate(time)

ggplot(df_x %>% pivot_longer(-time), aes(x=time, y=value, color=name)) +
  geom_line()

df_y = as_tibble(y, .name_repair="unique") %>%
  rename_at(vars(contains("...")), ~str_replace(., "...", "y")) %>%
  mutate(time=row_number()) %>%
  relocate(time)

ggplot(df_y %>% pivot_longer(-time), aes(x=time, y=value, color=name)) +
  geom_line()

# modeling ----------------------------------------------------------------

m = cmdstan_model("test.stan")

data = list(Time=Time,
            n_val=3,
            n_state=n_state,
            y=y,
            d_month=m_d)

init = list(eps_x=matrix(0, nrow=Time, ncol=n_state),
            v_lambda_others=rep(0, 3*n_state),
            mean_y=rep(0.4, 3),
            coef_month=m_d,
            sig_x=v_sd_x,
            sig_y=v_sd_y)

fit = m$sample(data=data,
               seed=123,
               iter_sampling=2000,
               iter_warmup=1000,
               chains=2,
               refresh=500,
               parallel_chains=2,
               init=rep(list(init), 2))

s = fit$summary()
draws_df = fit$draws(format="df")

mcmc_hist(fit$draws("lambda"))
