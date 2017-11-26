# CCSOGL
Multinomial classification with class-conditional feature group selection

1) This instruction provides an implementation of CCSOGL published in:
       http://www.sciencedirect.com/science/article/pii/S0167865517304087
       
   Li, Xiangrui, Dongxiao Zhu, and Ming Dong. "Multinomial classification with class-conditional overlapping sparse feature groups." Pattern Recognition Letters (2017).
   
2) Package dependencies:
   Rcpp, RcppArmadillo

3. To run the method, you first need to compile "CCSOGL.cpp" file using "Rcpp::sourceCpp" function from Rcpp package.

4. The main function cppMultiSOGL implements the CCSOGL algorithm. It has 9 arguments:
* x: the data matrix of dimension (n, p), where n is the number of samples and p is the number of features.
* y: the 0-1 label matrix of dimension (n, k) using one-hot encoding, where n is the number of samples, k is the number of classes.
* group: the 0-1 membership matrix of dimension (p, J), using multi-hot encoding, showing the feature membership for each feature group, where p is the number of features and J is the number of feature groups.
* dup_x: a list of length J, where each element is a submatrix corresponding to feature group.
* lambda1: the tuning parameter before L1-norm (i.e. lambda X alpha).
* lambda2: the tuning parameter before L2-norm (i.e. lambda X (1-alpha)).
* step: step size in the block coordinate descent.
* iter: maximal number of iteration.
* thresh: threshold for stopping algorithm (i.e. the difference of loss function values in two consecutive iterations is less than the threshold).

5. The values returned by cppMultiSOCGL:
* beta0: the intercepts
* beta: the coefficient matrix (calculated using latent_beta according to the coefficient decomposition)
* latent_beta: the latent coefficient matrix, corresponding to the decomposition of origninal coefficients beta (remember that we duplicate features that belong to multiple feature groups.)

6. An exmaple of CCSOGL along with the dataset is provided in r script. 
