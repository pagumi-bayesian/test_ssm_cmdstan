data {
   int Time;
   int n_val;
   int n_state;
   matrix[Time, n_val] y;
   matrix[Time, 11] d_month;
}

parameters {
   matrix[Time, n_state] eps_x;
   real<lower=0> lambda_first;
   vector[n_val*n_state-1] v_lambda_others;
   vector[n_val*11] v_coef_month;
   vector[n_val] mean_y;
   vector<lower=0>[n_state] sig_x;
   vector<lower=0>[n_val] sig_y;
}

transformed parameters {
   matrix[Time, n_state] x;
   matrix[Time, n_val] mu;
   matrix[Time, n_val] y_std;
   vector[n_val*n_state] v_lambda;
   matrix[n_val, n_state] lambda;
   matrix[n_val, 11] coef_month;
   
   v_lambda = append_row(lambda_first, v_lambda_others);
   lambda = to_matrix(v_lambda, n_val, n_state);
   coef_month = to_matrix(v_coef_month, n_val, 11);

   for (t in 1:Time) {
       if(t==1){
           x[1, ] = eps_x[1, ];
       }
       else if(t>=2){
           x[t, ] = x[t-1, ] + sig_x' .* eps_x[t, ];
       }

       mu[t, ] = mean_y' + (lambda * x[t, ]' + coef_month * d_month[t, ]')';
       y_std[t, ] = (y[t, ] - mu[t, ]) ./ sig_y';
   }

}

model {
   for(t in 1:Time){
       if(t==1){
           eps_x[1, ] ~ normal(0, 5);
       }
       if(t>=2){
           eps_x[t, ] ~ normal(0, 1);
       }
       y_std[t, ] ~ normal(0, 1);
   }

   v_lambda_others ~ double_exponential(0, 0.1);
   v_coef_month ~ normal(0, 1);
   mean_y ~ normal(0, 1);
   sig_x ~ normal(0, 0.1);
   sig_y ~ normal(0, 1);
}

generated quantities {
   matrix[Time, n_state] x_gen;

   for(t in 1:Time){

   }
}