// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD 1
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#include <RcppArmadillo.h>

#include <stdexcept>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export()]]
arma::mat ridge_closed_form_cpp_econ(arma::mat x,arma::vec y,arma::vec lam,arma::vec w,arma::vec eta,arma::mat x_val,arma::vec y_val) 
{
    int length_l_2D = lam.n_elem;
    int length_e_2D = eta.n_elem;
    arma::mat mse(length_l_2D,length_e_2D); //mse matrix

    int n = x.n_rows;
    int p = x.n_cols;
    int n_val = x_val.n_rows; //used for MSE calculation

    arma::mat U=zeros(n,n);
    arma::vec s_value; //singular values
    arma::mat s=zeros(n,p); //corresponding singular value matrix (if p>n)
    arma::mat V=zeros(p,p);
    arma::mat I=eye(p,p);
    arma::mat cov=zeros(p,p);
    arma::vec beta_est=zeros(p);
    arma::vec predicted_Y=zeros(n_val);
    arma::vec residual=zeros(n_val);

    
    List ret_temp(length_e_2D*length_l_2D);

    if (n>=p){
        for (int ll=0; ll!=length_l_2D;++ll){
            for (int ee=0; ee!=length_e_2D;++ee){
                beta_est = inv_sympd(x.t()*x+n*lam(ll)*I)*(x.t()*y +n*eta(ee)*w);
                predicted_Y=x_val*beta_est;  
                residual=predicted_Y-y_val;
                mse(ll,ee)=dot(residual,residual)/n_val;
                ret_temp[ll*length_e_2D+ee]=join_rows(predicted_Y,y_val);
            }
        }
    }
    else
    {
        svd_econ(U, s_value, V, x);

        List ret_svd(3);
        ret_svd["U"] = U; //n by n
        ret_svd["s"] = s_value; 
        ret_svd["V"] = V ; //p by p
        
        //initialize the elements first. Do not repeatly declare elements
        int k=s_value.n_elem;
        arma::vec A=zeros(k);
        arma::vec B=zeros(k);
        arma::vec C=zeros(p);
        arma::vec Sigma1_vec=zeros(k);
        arma::vec Sigma2_vec=zeros(k);
        arma::mat Sigma1=zeros(k,k);
        arma::mat Sigma2=zeros(k,k);

        A=U.t()*y;
        B=V.t()*w;
        C=(I-V*V.t())*w;

        for (int l=0; l!=length_l_2D;++l){

            for(int jj=0; jj!=k;++jj){
                Sigma1_vec(jj)=s_value(jj)/(pow(s_value(jj),2)+n*lam(l));
            }

            for(int jj=0; jj!=k;++jj){
                Sigma2_vec(jj)=1/(pow(s_value(jj),2)+n*lam(l));
            }

            Sigma1=diagmat(Sigma1_vec);
            Sigma2=diagmat(Sigma2_vec);

            for (int e=0; e!=length_e_2D;++e){

                beta_est=V*Sigma1*A+n*eta(e)*V*Sigma2*B+eta(e)/lam(l)*C;
                predicted_Y=x_val*beta_est;  
                residual=predicted_Y-y_val;
                mse(l,e)=dot(residual,residual)/n_val; 
                ret_temp[l*length_l_2D+e]=join_rows(predicted_Y,y_val);        
            }
        }        
    }

    
    //find the penalty parameter with smallest MSE
    double smallest = mse(0, 0);  // Initialize with the first element
    int rowIdx = 0;
    int colIdx = 0;
    for (int l=0; l!=length_l_2D;++l) {
        for (int e=0; e!=length_e_2D;++e) {
            if (mse(l, e) < smallest) {
                smallest = mse(l, e);
                rowIdx = l;
                colIdx = e;
            }
        }
    }

    //organize the results
    arma::mat ret_fin=zeros(n_val,4);

    // Access the desired element of the list
    arma::mat mat = as<arma::mat>(ret_temp[rowIdx*length_l_2D+colIdx]);

    ret_fin.col(0)=mat.col(0);
    ret_fin.col(1)=mat.col(1);
    ret_fin.col(2)=lam(rowIdx)*ones(n_val);
    ret_fin.col(3)=eta(colIdx)*ones(n_val);

    return(ret_fin);

}

// [[Rcpp::export()]]
arma::vec ridge_closed_form_best_cpp_econ (arma::mat x, 
                                arma::vec y,
                                double lam, 
                                arma::vec w, 
                                double eta) 
{

    int n = x.n_rows;
    int p = x.n_cols;

    arma::mat U=zeros(n,n);
    arma::vec s_value; //singular values
    arma::mat s=zeros(n,p); //corresponding singular value matrix (if p>n)
    arma::mat V=zeros(p,p);
    arma::mat I=eye(p,p);
    arma::mat cov=zeros(p,p);
    arma::vec beta_est=zeros(p);

    if (n>=p){
        beta_est = inv_sympd(x.t()*x+n*lam*I)*(x.t()*y +n*eta*w);           
    }
    else
    {
        svd_econ(U, s_value, V, x);

        List ret_svd(3);
        ret_svd["U"] = U; //n by n
        ret_svd["s"] = s_value; 
        ret_svd["V"] = V ; //p by p
        
        //initialize the elements first. Do not repeatly declare elements
        int k=s_value.n_elem;
        arma::vec A=zeros(k);
        arma::vec B=zeros(k);
        arma::vec C=zeros(p);
        arma::vec Sigma1_vec=zeros(k);
        arma::vec Sigma2_vec=zeros(k);
        arma::mat Sigma1=zeros(k,k);
        arma::mat Sigma2=zeros(k,k);

        A=U.t()*y;
        B=V.t()*w;
        C=(I-V*V.t())*w;

        for(int jj=0; jj!=k;++jj){
            Sigma1_vec(jj)=s_value(jj)/(pow(s_value(jj),2)+n*lam);
        }

        for(int jj=0; jj!=k;++jj){
            Sigma2_vec(jj)=1/(pow(s_value(jj),2)+n*lam);
        }

        Sigma1=diagmat(Sigma1_vec);
        Sigma2=diagmat(Sigma2_vec);

        beta_est=V*Sigma1*A+n*eta*V*Sigma2*B+eta/lam*C;      
    
    }

    return(beta_est);
}