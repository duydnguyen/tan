//Includes/namespaces
#include <Rcpp.h>
using namespace Rcpp;

//' Compute the adaptive Neyman test
//'
//' @param X A coverage vector for first condition
//' @param Y A coverage vector for second condition
//' @param na_rm TRUE to remove NA values
//' @param pool TRUE if variance is pooled. In this case, pre-computed poolVarX and poolVarY must be provided.
//'    If FALSE, pointwise variance is computed.
//' @param poolVarX A precomputed variance vector from the first condition. The pooled variance is implemented
//'    according to our paper.
//' @param poolVarY A precomputed variance vector from the second condition. The pooled variance is implemented
//'    according to our paper.
//' @export
// [[Rcpp::export]]
List AN_test(NumericMatrix X, NumericMatrix Y, LogicalVector na_rm, LogicalVector pool, NumericVector poolVarX, NumericVector poolVarY ) {
    LogicalVector naX;
    LogicalVector naY;
    /*  Remove NA values  */
    std::vector<bool> id(X.ncol());
    for(int i = 0; i < X.ncol(); i++) {
        id[i] = TRUE;
    }
    LogicalVector idxs(id.begin(), id.end());
    // Remove NA values from both X and Y
    if (na_rm) {
        for(int i = 0; i < X.nrow(); i++) {
            naX = is_na(X(i,_));
            naY = is_na(Y(i,_));
            idxs = idxs & (naX == FALSE & naY == FALSE);
        }
    }
    // Avoid var = 0 when two counts at a position t are the same
    std::vector<int> id1_same, id2_same, id_same(X.ncol());
    std::vector<int>::iterator it;
    for(int i = 0; i < X.ncol(); i++) {
        if (max(X(_,i)) == min(X(_,i)) ) {
            id1_same.push_back(i);
        }
        if (max(Y(_,i)) == min(Y(_,i)) ) {
            id2_same.push_back(i);
        }
    }
    it = std::set_union(id1_same.begin(), id1_same.end(), id2_same.begin(), id2_same.end(), id_same.begin() );
    // id_same contains index whose var is 0
    id_same.resize(it - id_same.begin());

    // Compute the size of non-NA values
    int count = 0;
    for(int i = 0; i < idxs.size(); i++) {
        if (idxs[i] == TRUE) {
            count += 1;
        }
    }
    count = count - id_same.size();
    NumericMatrix X_(X.nrow(), count);
    NumericMatrix Y_(Y.nrow(), count);
    int i_ = 0;
    if (na_rm(0) == 1) {
        for(int i = 0; i < idxs.size(); i++) {
            // check if id_same contains i
            std::vector<int>::iterator it1 = find(id_same.begin(), id_same.end(), i);
            if (it1 == id_same.end()) {
                if (idxs[i] == TRUE) {
                    X_(0, i_) = X(0, i);
                    X_(1, i_) = X(1, i);
                    Y_(0, i_) = Y(0, i);
                    Y_(1, i_) = Y(1, i);
                    i_ += 1;
                }
            }
        }
    }
    /*  Main Computation */
    NumericVector Var1(X_.ncol());
    NumericVector Var2(X_.ncol());
    NumericVector VarX(X_.ncol());
    NumericVector VarY(X_.ncol());
    NumericVector colMeanX(X_.ncol()), colMeanY(X_.ncol()) ;
    int p = X_.ncol(), n1 = X_.nrow(), n2 = Y_.nrow();
    int kstar = 0;
    double Tstar = -99999999.0; // -infinity
    if (pool(0) == 0) { // pool == FALSE
        NumericVector Xfoo(X_.ncol());
        // obtain poitwise variance
        for(int kk = 0; kk < X_.ncol(); kk++) {
            Var1[kk] = pow(sd(X_(_,kk)), 2 );
            Var2[kk] = pow(sd(Y_(_,kk)), 2);
            colMeanX[kk] = mean(X_(_,kk));
            colMeanY[kk] = mean(Y_(_,kk));
        }
        VarX[0] = Var1[0];
        VarY[0] = Var2[0];
        for(int kk = 1; kk < X_.ncol(); kk++) {
            NumericVector vec1 = Var1[seq_len(kk+1)-1];
            NumericVector vec2 = Var2[seq_len(kk+1)-1];
            VarX[kk] = sum(vec1) / (kk + 1);
            VarY[kk] = sum(vec2) / (kk + 1);
        }
        for(int kk = 1; kk < X_.ncol() + 1; kk++) {
            double tmp = 0;
            Xfoo = colMeanX[seq_len(kk)-1]  - colMeanY[seq_len(kk)-1] ;
            Xfoo = Xfoo / sqrt( VarX[kk-1]/n1 + VarY[kk-1]/n2 );
            tmp = sum(pow(Xfoo, 2) - 1 )/ sqrt(2 * (kk) );
            if ( (tmp > Tstar) && ( !std::isnan(tmp) )) {
                Tstar = tmp;
                kstar = kk;
            }
        }
        Tstar = sqrt(2 * log(log(p))) * Tstar - (2 * log(log(p)) + log(log(log(p)))/2 - log(4 * 3.141593)/2);
    } // end of if (pool_)
    else { // pool_ == TRUE
        NumericVector Xfoo(X_.ncol());
        // obtain poitwise variance
        for(int kk = 0; kk < X_.ncol(); kk++) {
            Var1[kk] = pow(sd(X_(_,kk)), 2 );
            Var2[kk] = pow(sd(Y_(_,kk)), 2);
            colMeanX[kk] = mean(X_(_,kk));
            colMeanY[kk] = mean(Y_(_,kk));
        }
        VarX[0] = Var1[0];
        VarY[0] = Var2[0];
        for(int kk = 1; kk < X_.ncol(); kk++) {
            NumericVector vec1 = Var1[seq_len(kk+1)-1];
            NumericVector vec2 = Var2[seq_len(kk+1)-1];
            VarX[kk] = sum(vec1) / (kk + 1);
            VarY[kk] = sum(vec2) / (kk + 1);
        }
        for(int kk = 1; kk < X_.ncol() + 1; kk++) {
            // Compute Cov of columns: Var1 is the pooled var of all vars of columns of X
            double tmp = 0;
            Xfoo = colMeanX[seq_len(kk)-1]  - colMeanY[seq_len(kk)-1] ;
            Xfoo = Xfoo / sqrt( poolVarX[kk-1]/n1 + poolVarY[kk-1]/n2 );
            tmp = sum(pow(Xfoo, 2) - 1 )/ sqrt(2 * (kk) );
            if ( (tmp > Tstar) && ( !std::isnan(tmp) )) {
                Tstar = tmp;
                kstar = kk;
            }
        }
        Tstar = sqrt(2 * log(log(p))) * Tstar - (2 * log(log(p)) + log(log(log(p)))/2 - log(4 * 3.141593)/2);
    } // end of else (pool_ == TRUE)

    //return wrap( sqrt( VarX[kk-1]/n1 + VarY[kk-1]/n2 )  ) ;
    //return wrap( Tstar );
    return List::create(Rcpp::Named("statistics") = wrap(Tstar),
                        Rcpp::Named("varX") = VarX,
                        Rcpp::Named("varY") = VarY);
}
