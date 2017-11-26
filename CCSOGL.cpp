#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// ########################################## ******* ###########################################
// ###                                        probability 
// [[Rcpp::export]]
mat beta0AddXB(int n, int k, rowvec beta0, mat XB)
{
	mat addi(n, k);
	for(int i = 0; i < n; i++)
	{
		addi.row(i) = beta0 + XB.row(i);
	}
	return addi;
}

// [[Rcpp::export]]
mat pCalc(int n, int k, mat beta0_add_XB)
{
   mat prob(n, k);
   mat exp_beta0_XB = exp(beta0_add_XB);
   colvec s = sum(exp_beta0_XB, 1) ;
   for(int i = 0; i < n; i++)
   {
	   prob.row(i) = exp_beta0_XB.row(i) / s(i);
   }
   return prob;
}

// ########################################## ******* ###########################################
// ###                                       x %*% beta
// [[Rcpp::export]]
cube cube_xTimesBeta(int n, int j, int k,
                     rowvec group_size, 
					 List dup_x, mat l_beta)
{
	cube XB_cube(n, k, j);
	int head = 0;
	int tail = -1;
	for(int i = 0; i < j; i++)
	{
		tail = tail + group_size(i);
		XB_cube.slice(i) = as<mat>(dup_x[i]) * l_beta.rows(head, tail);
		head = head + group_size(i);	
	}
	return XB_cube;
}

// [[Rcpp::export]]
mat xTimesBeta(int n, int j, int k, cube XB_cube)
{
	return sum(XB_cube, 2);
}


// ########################################## ******* ###########################################
// ###                                       Compute cost function
// ### Scaled Negative loglikelihood function
// [[Rcpp::export]]
double NegLoglikelihoodCalc(int n, mat prob, mat y)
{
	double logLik = 0;
	for(int i = 0; i < n; i++)
	{
		logLik = logLik + sum(log(prob.row(i)) % y.row(i));
	}
	return -logLik / n;
}

// ### l_1 norm (||beta||_1)
// [[Rcpp::export]]
double L1Norm(mat l_beta)
{
	return accu(abs(l_beta));
}

// ### group size scaled l_2 norm ( \sum sqrt(d_j) * ||beta_j||_2 )
// [[Rcpp::export]]
double L2Norm(int j, rowvec group_size, mat l_beta)
{
	double s = 0;
	int head = 0;
	int tail = -1;
	mat ma = pow(l_beta, 2);
	for(int i = 0; i < j; i++)
	{
		tail = tail + group_size(i);
		s = s + accu(sqrt(group_size(i) * sum(ma.rows(head, tail), 0)));
		head = head + group_size(i);
	}
	return s;
}

// ### cost function
// [[Rcpp::export]]
double costFunc(int n, int j, 
                mat prob, mat y, 
				rowvec group_size, mat l_beta, 
				double lambda1, double lambda2)
{
	return NegLoglikelihoodCalc(n, prob, y) + 
	       lambda1 * L1Norm(l_beta) +
		   lambda2 * L2Norm(j, group_size, l_beta);
}				
	
// ########################################## ******* ###########################################
// ###                                       beta0 solver
// [[Rcpp::export]]
double k_th_betaZeroSolver(int n, int k, int k_th, mat beta0_add_XB, 
                           rowvec beta0, mat XB, colvec y_k_th,
                           int iter = 1000, double  thresh = 0.000001)
{
	int count = 0;
  double num = 0;
  double denom = 0;
  double delta = 10;
	double beta0_k_th = beta0(k_th);
  mat pr;
    
	while(pow(delta, 2) > pow(thresh,2) && count < iter)
  {
    pr = pCalc(n, k, beta0_add_XB);
    num = 0;
    denom = 0;
  
    for(int i = 0; i < n; i++)
    {
      num = num + pr.col(k_th)(i) - y_k_th(i);
      denom = denom + pr.col(k_th)(i) * (1 - pr.col(k_th)(i));
    }

    delta = num / denom;
    beta0_k_th = beta0_k_th - delta;
		beta0_add_XB.col(k_th) = beta0_k_th + XB.col(k_th);
    count += 1;
  }
	//printf("#### beta_%d0 iteration: %d \n", k_th + 1, count);
  return beta0_k_th;
}

// ########################################## ******* ###########################################
// ###                                       block gradient
// [[Rcpp::export]]
colvec blockGradient(int n, int d_j, int k_th, 
                     colvec prob_k_th, mat x_block, colvec y_k_th)
{
	rowvec blockgrad(d_j); blockgrad.zeros();
	for(int i = 0; i < n; i++)
	{
		blockgrad = blockgrad + (prob_k_th(i) - y_k_th(i)) * x_block.row(i);
	}
	return blockgrad.t() / n;
}

// ########################################## ******* ###########################################
// ###                                       soft-thresholding
//[[Rcpp::export]]
colvec SoftThresh(colvec x, double lambda)
{
	int n = x.n_rows;
	colvec s(n);
	colvec abs_x = abs(x);
	colvec chopper(n);
	colvec xnew(n);
	
	for(int i = 0; i < n; i++)
  {
    if(x[i] > 0)
    {
      s[i] = 1;
    }
    if(x[i] == 0)
    {
      s[i] = 0;
    }
    if(x[i] < 0)
    {
      s[i] = -1;
    }
  }
  
  for(int i = 0; i < n; i++)
  {
    if(abs_x[i] >= lambda)
    {
      chopper[i] = abs_x[i] - lambda;
    }
    else 
    {
      chopper[i] = 0;
    }
  }
  
  for(int i = 0; i < n; i++)
  {
    xnew[i] = s[i] * chopper[i];
  }
  
  return xnew;
}

// ###                                       group minimizer in coordinate descent
// [[Rcpp::export]]
colvec groupMinimizer(double eta, colvec tau,
                      double lambda1, double lambda2)
{
	int n = tau.n_rows;
	double s = 0;  // norm of soft-thresholded tau
	colvec chopped_tau = SoftThresh(tau, lambda1 * eta);
	colvec minimizer(n); minimizer.zeros();
	s = sqrt(accu(pow(chopped_tau,2)));
	if(s > lambda2 * eta)
  {
    minimizer = (1 - lambda2 * eta / s) * chopped_tau;
  }
  return minimizer;
}

// ########################################## ******* ###########################################
// ###                               recover true beta from latent_beta
// [[Rcpp::export]]
mat betaRecover(int p, int j, int k, mat group, rowvec group_size, mat l_beta)
{
  int head = 0;
  int tail = -1;
  int counter;
  mat beta(p,k); beta.zeros();
  
  for(int h = 0; h < j; h++)
  {
    counter = 0;
    tail = tail + group_size(h);
    for(int l = 0; l < p; l++)
    {
      if(group(l, h) == 1)
      {
        counter += 1;
        beta.row(l) = beta.row(l) + (l_beta.rows(head,tail)).row(counter-1);
      }
    }
    head = head + group_size(h);
  }
  return beta;
}

// ########################################## ******* ###########################################
// ###                                        ultimate Solver
// [[Rcpp::export]]
List cppMultiSOGL(mat x, mat y, mat group, List dup_x,
                  double lambda1, double lambda2, 
                  double step = 0.1, int iter =3000, double thresh = 0.00002)
{
  int k = y.n_cols;      // # of classes
  int j = group.n_cols;  // # of groups
  int p = x.n_cols;      // # of features
  int n = x.n_rows;      // # of observations
  rowvec group_size = sum(group, 0); // a vector, recording size of each group
  
  mat l_beta(accu(group_size),k); l_beta.zeros();  // initialize latent coefficents
  rowvec beta0(k); beta0.zeros();                 // ####  initialize beta_0
  mat beta0_add_XB(n, k); beta0_add_XB.zeros();
  
  cube XB_cube(n, k, j); XB_cube.zeros();
  mat XB = sum(XB_cube, 2);
  mat prob = pCalc(n, k, beta0_add_XB);
  colvec gradient_block;
  
  vec cost(iter); cost.zeros(); // ###  recording cost for algorithm convergence
  cost(0) = costFunc(n, j, prob, y, group_size, l_beta, lambda1, lambda2);
  
  int head;
  int tail;
  int count = 1;
  while(count < iter)  // update in cyclic order
  {
    if(count%50 == 0){
      printf("Iteration: %d \n", count);
    }
    for(int k_th = 0; k_th < k; k_th++)
    {
      beta0(k_th) = k_th_betaZeroSolver(n, k, k_th, beta0_add_XB, beta0, XB, y.col(k_th));
      //printf("beta_%d0: %f \n", k_th + 1, beta0(k_th));
      beta0_add_XB.col(k_th) = beta0(k_th) + XB.col(k_th);
      prob = pCalc(n, k, beta0_add_XB);
      head = 0;
      tail = -1;
      for(int h = 0; h < j; h++)
      {
        tail = tail + group_size(h);
        gradient_block = blockGradient(n, group_size(h), k_th, prob.col(k_th), 
                                       as<mat>(dup_x[h]), y.col(k_th));
        (l_beta.col(k_th)).rows(head,tail) = 
          groupMinimizer(step, (l_beta.col(k_th)).rows(head,tail) - step * gradient_block,
                         lambda1, lambda2 * sqrt(group_size(h)));
        
        // update XB
        colvec old = XB_cube.slice(h).col(k_th);
        
        (XB_cube.slice(h)).col(k_th) = as<mat>(dup_x[h]) * (l_beta.col(k_th)).rows(head,tail);
        XB.col(k_th) = XB.col(k_th) - old + (XB_cube.slice(h)).col(k_th);
        beta0_add_XB.col(k_th) = beta0(k_th) + XB.col(k_th);
        prob = pCalc(n, k, beta0_add_XB);
        head = head + group_size(h);
      }
    }
    
    cost(count) = costFunc(n, j, prob, y, group_size, l_beta, lambda1, lambda2);
    //printf("Cost: %f \n", cost(count));
    
    if(cost(count) > cost(count-1) || count == iter-1)
    {
      printf("Algorithm does NOT converge! Stepsize is too large.\n");
      break;
    }
    if(std::abs(cost(count) - cost(count-1)) / std::abs(cost(count-1)) <= thresh)
    {
      break;
    }
    count = count + 1;
    
  }
  printf("# of iterations: %d \n", count);
  mat beta = betaRecover(p, j, k, group, group_size, l_beta);
  return List::create(Named("beta0") = beta0,
                      Named("beta") = beta,
                      Named("iteration") = count + 1,
                      Named("cost") = cost,
                      Named("latent_beta") = l_beta);  
}








