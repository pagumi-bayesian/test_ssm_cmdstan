data {
   int Time;
   int n_val;
   int n_state;
   matrix[Time, n_val] y;
   matrix[Time, 11] d_month;
}

transformed data {
  real sig_x;
  sig_x = 0.01;
}

parameters {
  matrix[Time-1, n_state] eps_x;
  // 因子負荷は識別性のため，行列の右上部分が全て0になるよう設計
  cholesky_factor_cov[n_val, n_state] lambda;
   
  vector[n_val*11] v_coef_month;
  vector[n_val] mean_y;
  vector<lower=0>[n_val] sig_y;
  //real<lower=0> sig_x;
}

transformed parameters {
  matrix[Time-1, n_state] d_x;
  matrix[Time-1, n_state] sum_d_x;
  matrix[Time-1, n_state] v_d_x;
  vector[n_state] sum_v_d_x;
  matrix[Time, n_state] x;
  matrix[Time, n_val] mu;
  matrix[n_val, 11] coef_month;
  
  coef_month = to_matrix(v_coef_month, n_val, 11);
  
  // 状態変数の期間トータルが0になるよう変換(コーホート分析を参考)
  d_x = sig_x * eps_x;
  
  for(t in 1:(Time-1)){
    v_d_x[t, ] = (Time - t) * d_x[t, ];
  }
  
  for(s in 1:n_state){
    sum_v_d_x[s] = sum(v_d_x[, s]);
    
    for(t in 2:Time){
      sum_d_x[t-1, s] = sum(d_x[1:(t-1), s]);
    }
  }
  
  for(t in 1:Time){
    if(t==1){
      x[1, ] = -(1.0/Time) * sum_v_d_x';
    }else{
      x[t, ] = x[1, ] + sum_d_x[t-1, ];
    }
    
    mu[t, ] = mean_y' + (lambda * x[t, ]' + coef_month * d_month[t, ]')';
    
  }
}

model {
  to_vector(eps_x) ~ normal(0, 1);
  
  for(t in 1:Time){
    y[t, ] ~ normal(mu[t, ], sig_y');
  }
  
  to_vector(lambda) ~ normal(0, 1);
  
  v_coef_month ~ normal(0, 1);
  mean_y ~ normal(0, 1);
  sig_y ~ cauchy(0, 0.01);
  //sig_x ~ exponential(1000);
}

generated quantities {
  matrix[Time, n_val] y_rng;
  
  for(t in 1:Time){
    for(n in 1:n_val){
      y_rng[t, n] = normal_rng(mu[t, n], sig_y[n]);
    }
  }
}
