//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <R.h>
#include <Rcpp.h>

// To use functions in Rcpp library without the need of specifying "Rcpp::"
using namespace Rcpp;

//----------------------------------------------------------------------------------------
// C++ script
// 
// title:
// variational approximation for Bayesian SEM using
// multiple global shrinkage priors
//
// authors:
// Gwenael G.R. Leday, Mark A. van de Wiel
//
// date: 05/06/2018
//
//----------------------------------------------------------------------------------------





//----------------------------------------------------------------------------------------
// Function that carries out a single variational update for a single regression model
//----------------------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List VarOneIter(arma::colvec myy, arma::mat myX, arma::colvec myP, arma::colvec aRand, arma::colvec bRand, arma::colvec bRandStarInit, double dSigmaStarInit, double cSigma, double dSigma, arma::mat lincomb){
  
  // Initialization
  int then = myX.n_rows;
  int thep = myX.n_cols;

  // Get p_1, p_2, ..., p_K (nb parameters per random effect) 
  arma::colvec prvals = sort(unique(myP));
  int K = prvals.n_elem;
  arma::colvec ss = arma::zeros(K);
  for(int k=0; k<K; k++){
    ss(k) = sum(myP==prvals(k));
  }

  // Update posterior shape parameters
  arma::colvec aRandStar = aRand + 0.5*ss;
  double cSigmaStar = cSigma + 0.5*then + 0.5*thep;

  // Initialization variational parameters associated with random effects
  arma::colvec u = arma::zeros(thep);
  arma::colvec mybRandStarInit = arma::zeros(thep);

  for(int k=0; k<K; k++){
    u.elem(find(myP==prvals(k))) += aRandStar(k);
    if(all(bRandStarInit==0)){
      mybRandStarInit.elem(find(myP==prvals(k))) += bRand(k);
    }else{
      mybRandStarInit.elem(find(myP==prvals(k))) += bRandStarInit(k);
    }//end if
  }//end for

  // Expectations
  double expSig;
  arma::mat D;
  if(all(bRandStarInit==0)){
    D = diagmat(u/mybRandStarInit);
    expSig = cSigmaStar/dSigma;
  }else{
    D = diagmat(u/mybRandStarInit);
    expSig = cSigmaStar/dSigmaStarInit;
  }//end if

  //Cross-product
  arma::mat XTX = myX.t() * myX;

  //Update sigma
  arma::mat postSigma;
  if(then>=thep){
    postSigma = inv(expSig*(XTX + D));
  }else{
    arma::mat Dinv = diagmat(1/D.diag());
    postSigma = (1/expSig)*(Dinv - Dinv*myX.t()*inv(arma::eye(then,then)+myX*Dinv*myX.t())*myX*Dinv);
  }//end if
  arma::colvec v = postSigma.diag();

  // Calculate posterior mean of beta
  arma::colvec postMean = expSig*postSigma*trans(myX)*myy;

  // Update posterior rate for sigma^2
  arma::colvec myresid = myy-myX*postMean;
  arma::colvec tempval = v+square(postMean);
  double dSigmaStar = dSigma + 0.5*(as_scalar(trans(myresid)*(myresid))+trace(XTX*postSigma)) + 0.5*sum(D.diag()%tempval);

  // Update posterior rates for tau^2
  arma::colvec bRandStar = arma::zeros(K);
  for(int i=0; i<K; i++){
    bRandStar(i) = bRand(i) + 0.5*expSig*sum(tempval.elem(find(myP==prvals(i))));
  }//end for

  // Lower bound marginal likelihood
  double Lrand = 0;
  for(int i=0; i<K; i++){
    Lrand += aRand(i)*log(bRand(i))-aRandStar(i)*log(bRandStar(i))+lgamma(aRandStar(i))-lgamma(aRand(i));
  }//end for
  double Lsig = cSigma*log(dSigma)-cSigmaStar*log(dSigmaStar)+lgamma(cSigmaStar)-lgamma(cSigma);
  double theplus = 0.5*(cSigmaStar/dSigmaStar)*sum((u/mybRandStarInit)%tempval);
  double valLogDet, sign;
  log_det(valLogDet, sign, postSigma);
  double L;
  try{
    L = 0.5*thep - 0.5*then*log(2*arma::datum::pi) + 0.5*valLogDet + Lrand + Lsig + theplus;
  }
  catch(...){
    L = arma::datum::nan;
  }

  // ---  OUTPUT ---
  // priorRand matrix
  arma::mat priorRand(K,2);
  priorRand.zeros();
  priorRand.col(0) += aRand;
  priorRand.col(1) += bRand;
  // priorSig vector
  arma::colvec priorSig(2);
  priorSig.zeros();
  priorSig(0) = cSigma;
  priorSig(1) = dSigma;
  // postRand matrix
  arma::mat postRand(K,2);
  postRand.zeros();
  arma::colvec col0 = postRand.col(0);
  arma::colvec col1 = postRand.col(1);
  for(int k=0; k<K; k++){
    col0.elem(find(prvals==prvals(k))) += aRandStar(k);
    col1.elem(find(prvals==prvals(k))) += bRandStar(k);
  }
  postRand.col(0) += col0;
  postRand.col(1) += col1;
  // postSig vector
  arma::colvec postSig(2);
  postSig.zeros();
  postSig(0) = cSigmaStar;
  postSig(1) = dSigmaStar;
  // postBeta matrix
  arma::mat postBeta(thep,3);
  postBeta.zeros();
  postBeta.col(0) += postMean;
  postBeta.col(1) += sqrt(v);
  postBeta.col(2) += abs(postMean)/sqrt(v);
  // lincomb
  arma::colvec postMeanLincomb;
  arma::mat postVarLincomb;
  if(lincomb.n_rows>1){
    postMeanLincomb = lincomb*postMean;
    postVarLincomb = lincomb*postSigma*trans(lincomb);
  }

  return Rcpp::List::create(
				Rcpp::Named("L") = L,
				Rcpp::Named("priorRand") = priorRand,
				Rcpp::Named("priorSig") = priorSig,
				Rcpp::Named("postRand") = postRand,
				Rcpp::Named("postSig") = postSig,
				Rcpp::Named("postBeta") = postBeta,
				Rcpp::Named("postBeta") = postBeta,
				Rcpp::Named("postMeanLincomb") = postMeanLincomb,
				Rcpp::Named("postVarLincomb") = postVarLincomb);

}// end VarOneIter

//----------------------------------------------------------------------------------------
// Digamma function for vector
//----------------------------------------------------------------------------------------
arma::colvec mydigamma(arma::colvec vec){
  arma::colvec out(vec.n_elem);
  for(int k = 0; k<vec.n_elem; k++){
    out(k) = R::digamma(vec(k));
  }
  return out;
}

//----------------------------------------------------------------------------------------
// Fixed-point iteration for gamma prior estimation (VB empirical Bayes)
//----------------------------------------------------------------------------------------

Rcpp::NumericVector arma2vec(arma::colvec x) {
    return Rcpp::NumericVector(x.begin(), x.end());
}

// [[Rcpp::export]]
Rcpp::NumericVector fixedPointIterEB(arma::colvec initab, arma::colvec myallaRandStar, arma::colvec myallbRandStar, int mymaxiter, double myeps){
  arma::colvec myallaRandStar2 = nonzeros(myallaRandStar);
  arma::colvec myallbRandStar2 = nonzeros(myallbRandStar);
  arma::colvec a(mymaxiter+1);
  arma::colvec b(mymaxiter+1);
  a.zeros();
  b.zeros();
  a(0) = initab(0);
  b(0) = initab(1);
  int cpt = 1;
  bool mybool = true;
  double tp = (myallaRandStar2.n_elem/sum(myallaRandStar2/myallbRandStar2));
  double tp2 = mean(log(myallbRandStar2) - mydigamma(myallaRandStar2))-log(tp);
  while(mybool){
    a(cpt) = a(cpt-1) + 0.5*(1/(R::digamma(a(cpt-1))-log(a(cpt-1)))) + 0.5*(1/tp2);
    b(cpt) = a(cpt)*tp;
    if((abs(a(cpt)-a(cpt-1))<myeps) && (abs(b(cpt)-b(cpt-1))<myeps)){
      mybool = false;
    }else{
      if(cpt==mymaxiter){
        mybool = false;
      }else{
        cpt++;
      }//end if
    }//end if
  }//end while
  arma::colvec ab(2);
  ab(0) = a(cpt);
  ab(1) = b(cpt);
  
  return arma2vec(ab);
}




//----------------------------------------------------------------------------------------
// Function that carries out a single variational update for the BSEM model
//----------------------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List BSEMVarOneIter(Rcpp::List ylist, Rcpp::List Xlist, Rcpp::List Plist, Rcpp::List alist, Rcpp::List blist, Rcpp::List bstarlist, double cSigma, double dSigma, arma::colvec dstarvec, Rcpp::List lincomblist){
  
	// Initialize vector of variational lower bounds
	arma::colvec allmargs(Xlist.size());
	// Get number of global shrinkage priors (SAME NUMBER OF SHRINKAGE PRIORS PER REGRESSION)
	arma::colvec prvals = sort(unique(Rcpp::as<arma::colvec>(Plist[1])));
	// int K = prvals.n_elem; (unused variable 'K' warning)
	// Initialize temporary list
	Rcpp::List tplist;
	// Initialize output list
	Rcpp::List postSigList(Xlist.size());
	Rcpp::List postRandList(Xlist.size());
	Rcpp::List postBetaList(Xlist.size());
	Rcpp::List postMeanLincombList(Xlist.size());
	Rcpp::List postVarLincombList(Xlist.size());
	arma::mat emptyMat(1,1);

	// Fit all regression models
	for(int j=0; j<Xlist.size(); j++){

		if(lincomblist.size()>1){

			// update variational parameters once
			tplist = VarOneIter(Rcpp::as<arma::colvec>(ylist[j]), Rcpp::as<arma::mat>(Xlist[j]), Rcpp::as<arma::colvec>(Plist[j]), Rcpp::as<arma::colvec>(alist[j]), Rcpp::as<arma::colvec>(blist[j]), Rcpp::as<arma::colvec>(bstarlist[j]), dstarvec(j), cSigma, dSigma, Rcpp::as<arma::mat>(lincomblist[j]));

		}else{

			// update variational parameters once
			tplist = VarOneIter(Rcpp::as<arma::colvec>(ylist[j]), Rcpp::as<arma::mat>(Xlist[j]), Rcpp::as<arma::colvec>(Plist[j]), Rcpp::as<arma::colvec>(alist[j]), Rcpp::as<arma::colvec>(blist[j]), Rcpp::as<arma::colvec>(bstarlist[j]), dstarvec(j), cSigma, dSigma, emptyMat);

		}
		// get variational lower bound
		allmargs(j) = Rcpp::as<double>(tplist["L"]);

		// get postRand, postSig and postBeta
		postRandList[j] = Rcpp::as<arma::mat>(tplist["postRand"]);
		postSigList[j] = Rcpp::as<arma::colvec>(tplist["postSig"]);
		postBetaList[j] = Rcpp::as<arma::mat>(tplist["postBeta"]);
		postMeanLincombList[j] = Rcpp::as<arma::colvec>(tplist["postMeanLincomb"]);
		postVarLincombList[j] = Rcpp::as<arma::mat>(tplist["postVarLincomb"]);

	}//end for
  
  return Rcpp::List::create(
	Rcpp::Named("allmargs") = allmargs,
	Rcpp::Named("postRandList") = postRandList,
	Rcpp::Named("postSigList") = postSigList,
	Rcpp::Named("postBetaList") = postBetaList,
	Rcpp::Named("postMeanLincombList") = postMeanLincombList,
	Rcpp::Named("postVarLincombList") = postVarLincombList);

}//end varAlgoP


