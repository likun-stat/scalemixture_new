// -------------------------------------------------------------------------- //
// -------------------------     Matrix Algebra     ------------------------- //
// -------------------------------------------------------------------------- //
//  Created by ZhangLikun on 9/17/18.
//  Copyright © 2018 ZhangLikun. All rights reserved.
//




#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

/*
#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
*/

using namespace Rcpp;

// [[Rcpp::export()]]
double eig2logdet_c (arma::vec d){
    arma::vec log_d = log(d);
    return sum(log_d);
}

// [[Rcpp::export()]]
double eig2inv_quadform_vector (arma::mat V, arma::vec d_inv, arma::vec x) {
    arma::mat step_r = V.t()*x;
    arma::mat step_m = diagmat(d_inv)*step_r;
    return dot(step_r.t(), step_m);
}

// [[Rcpp::export()]]
arma::mat eig2inv_times_vector (arma::mat V, arma::vec d_inv, arma::vec x) {
    arma::mat step_r = V.t()*x;
    arma::mat step_m = diagmat(d_inv)*step_r;
    
    return V*step_m;
}

// [[Rcpp::export()]]
double dmvn_eig (arma::mat R, arma::mat V, arma::vec d_inv){
    int n_rep = R.n_cols;
    double left = -0.5*n_rep*eig2logdet_c(1/d_inv);
    
    arma::vec right(n_rep);
    for(int i=0; i<n_rep; i++){
        arma::vec R_i = R.col(i);
        right(i) = eig2inv_quadform_vector(V, d_inv, R_i);
    }
    
    return left-0.5* sum(right);
}

// [[Rcpp::export()]]
arma::vec trans_NVec (NumericVector x){
    arma::vec y = as<arma::vec>(x);
    return y;
}

// [[Rcpp::export()]]
double X_s_likelihood_conditional(arma::vec X_s, double R, arma::mat V, arma::vec d){
    NumericVector tmp = qnorm(as<NumericVector>(wrap(1-R/X_s)));
    arma::vec X_s_to_Z = trans_NVec(tmp);
    double part1 = -0.5*eig2inv_quadform_vector(V, 1/d, X_s_to_Z);
    double part2 = 0.5*sum(X_s_to_Z % X_s_to_Z)-2*sum(log(X_s));
    return part1+part2;
}



// [[Rcpp::export()]]
double X_s_likelihood_conditional_on_X (arma::vec X_s, arma::vec X, double R, arma::mat V, arma::vec d, double tau_sqd){
    
    // treat X_s and X as vec
    // arma::vec X_s = vectorise(X_s_m);
    // arma::vec X = vectorise(X_m);
    
    if(any(X_s<R)) return -std::numeric_limits<double>::infinity();
    else{
        NumericVector tmp = qnorm(as<NumericVector>(wrap(1-R/X_s)));
        arma::vec X_s_to_Z = trans_NVec(tmp);
        double part1 = -0.5*eig2inv_quadform_vector(V, 1/d, X_s_to_Z);
        double part2 = 0.5*sum(X_s_to_Z % X_s_to_Z)-2*sum(log(X_s));
        double part3 = - 0.5*sum((X_s-X) % (X_s-X))/tau_sqd;
        
        return part1+part2+part3;
    }
}



// [[Rcpp::export()]]
List var_at_a_time_update_X_s (arma::vec X, double R, arma::mat V, arma::vec d, double tau_sqd, double v_q=0.5, int n_chain=100){
    
    int n_s = X.n_elem;
    arma::vec X_s = X;
    arma::uvec tmp =find(X_s<R);
    
    if(tmp.n_elem>0){
        for(unsigned int i=0; i<tmp.n_elem; i++){
            int ind = tmp(i);
            arma::vec nugget = 0.5*arma::randn(1);
            X_s(ind) = R + fabs(as_scalar(nugget));
        }
    }
    
    arma::vec accept(n_s);
    accept.fill(0);
    
    arma::vec X_s_update(n_s);
    double log_num=0, log_denom=0, r=0;
    for(int i=0; i<n_chain; i++){
        for(int iter=0; iter<n_s; iter++){
            X_s_update = X_s;
            X_s_update(iter) = X_s[iter]+v_q*as_scalar(arma::randn(1));
            log_num = X_s_likelihood_conditional_on_X(X_s_update, X, R, V, d, tau_sqd);
            log_denom = X_s_likelihood_conditional_on_X(X_s, X, R, V, d, tau_sqd);
            
            r = exp(log_num - log_denom);
            if(as_scalar(arma::randu(1))<r){
                X_s = X_s_update;
                accept(iter) = accept(iter) + 1;
            }
        }
    }
    
    List result;
    result["X.s"] = X_s;
    result["accept"] = accept;
    
    return result;
}


// ------------------------------------------------------------------------------ //
// ------------- Apply OpenMP to parallelize over time replications ------------- //
// ------------------------------------------------------------------------------ //
/*
// [[Rcpp::export]]
double long_computation_omp(int nb, int threads=1) {
#ifdef _OPENMP
    if ( threads > 0 )
        omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif
    
    double sum = 0;
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nb; ++i) {
        double thread_sum = 0;
        for (int j = 0; j < nb; ++j) {
            thread_sum += R::dlnorm(i+j, 0.0, 1.0, 0);
        }
        sum += thread_sum;
    }
    return sum + nb;
}

// check OpenMP
//[[Rcpp::export]]
int tryfor (int iter){
    
    double A[iter];
    for (int i =0; i<iter; i++){
        A[i] = sin(i)*sin(i);
    }
    return A[0];
}

//[[Rcpp::export]]
int tryfor_omp (int iter, int threads = 1){
#ifdef _OPENMP
    if ( threads > 0 ) omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif
    
    double A[iter];
    
#pragma omp parallel for
    for (int i =0; i<iter; i++){
        A[i] = sin(i)*sin(i);
    }
    return A[0];
}

*/

/*
// [[Rcpp::export]]
List var_at_a_time_update_X_s_cols (arma::vec R,  arma::mat X, double tau_sqd, arma::mat V, arma::vec d, double v_q=0.5, int n_chain=100, int threads = 1) {
#ifdef _OPENMP
    if ( threads > 0 ) omp_set_num_threads( threads );
    REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif
    
    int n_s = X.n_rows;
    int n_t = X.n_cols;
    
    arma::vec log_rat(n_t);
    log_rat.fill(0);
    
    arma::mat prop_X_s(n_s, n_t);
    
    // Proposal
    // Treat the X process as the response, ignoring the marginal transformation
    // i.e. prop.X.s ~ X.s | X, R, other.params
    
#pragma omp parallel for schedule(static)
    for(int t=0; t<n_t; t++){
        
        arma::vec prop = X.col(t);
        arma::uvec tmp =find(prop<R(t));
        
        if(tmp.n_elem>0){
            for(unsigned int i=0; i<tmp.n_elem; i++){
                int ind = tmp(i);
                arma::vec nugget = 0.5*arma::randn(1);
                prop(ind) = R(t) + fabs(as_scalar(nugget));
            }
        }
        
        
        arma::vec prop_update(n_s);
        double log_num=0, log_denom=0, r=0;
        int i =0;
        while(i<2){
            for(int iter=0; iter<n_s; iter++){
                prop_update = prop;
                prop_update(iter) = prop(iter)+v_q*as_scalar(arma::randn(1));
                log_num = X_s_likelihood_conditional_on_X(prop_update, X.col(t), R(t), V, d, tau_sqd);
                log_denom = X_s_likelihood_conditional_on_X(prop, X.col(t), R(t), V, d, tau_sqd);
                
                r = exp(log_num - log_denom);
                if(as_scalar(arma::randu(1))<r){
                    prop = prop_update;
                }
            }
            i++;
        }
        
        
        prop_X_s.col(t) = prop;
        
        // partial M-H ratio
       
         log_rat[t] = X_s_likelihood_conditional_on_X(X_s.col(t), X.col(t), R(t), V, d, tau_sqd) +
         X_s_likelihood_conditional(prop_X_s.col(t), R(t), V, d)-
         X_s_likelihood_conditional_on_X(prop_X_s.col(t), X.col(t), R(t), V, d, tau_sqd)-
         X_s_likelihood_conditional(X_s.col(t), R(t), V, d);
        
    }
    
    List result;
    result["prop_X_s"] = prop_X_s;
    result["log_rat"] = 1;
    
    return result;
    
}
*/
