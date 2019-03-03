#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

#define ALLOC(a,b)  R_alloc(a,b)
double **dmatrix(double *array, int ncol, int nrow){
   int i;
   double **pointer;
   pointer = (double **) ALLOC(nrow, sizeof(double *));
   for (i=0; i<nrow; i++) {
     pointer[i] = array;
     array += ncol;
 	}
  return(pointer);
}



 SEXP ScoreAFUNC(SEXP time2,   SEXP status2,  SEXP covar2,  SEXP betahat2,
                 SEXP offset2, SEXP weights2, SEXP strata2               ) {

  double **covar, **imat,   **cmat,  **cmat2;
  double *time,   *weights, *offset, *beta,  *newbeta, *a, *a2, *scale, *means, *u;
  double denom =0, denom2, dtime, deadwt, temp, temp2, wtave, zbeta, risk;
  SEXP   beta2, means2, u2, imat2, rlist, rlistnames;
  int    *status, *flag, *strata;
  int    nprotect, nvar, nused, i, j, k, person, nrisk, ndead;

  nvar  = ncols(covar2);
  nused = LENGTH(time2);
  nprotect = 0;

  PROTECT(beta2 = duplicate(betahat2));
  nprotect++;

  beta = REAL(beta2);

  time    = REAL(time2);
  weights = REAL(weights2);
  offset  = REAL(offset2);
  status  = INTEGER(status2);
  strata  = INTEGER(strata2);
  covar   = dmatrix(REAL(covar2), nused, nvar);


  PROTECT(means2 = allocVector(REALSXP, nvar));
  means = REAL(means2);
  nprotect++;

  PROTECT(imat2 = allocVector(REALSXP, nvar*nvar));
  imat = dmatrix(REAL(imat2),  nvar, nvar);
  nprotect++;

  PROTECT(u2 = allocVector(REALSXP, nvar));
  u = REAL(u2);
  nprotect++;

  a = (double *) R_alloc(2*nvar*nvar + 4*nvar, sizeof(double));
  newbeta = a + nvar;

  a2    = newbeta + nvar;
  scale = a2 + nvar;
  cmat  = dmatrix(scale + nvar,   nvar, nvar);
  cmat2 = dmatrix(scale + nvar +  nvar*nvar, nvar, nvar);

  strata[nused-1] =1;
  for (i=0; i<nvar; i++) {
	  u[i] =0;
	  a2[i] =0;
	  for (j=0; j<nvar; j++) {
	    imat[i][j] =0 ;
	    cmat2[i][j] =0;
    }
  }
  // Efron method to handle ties//
  for (person=nused-1; person>=0; ){
    if (strata[person] == 1) {
        nrisk =0 ;
        denom = 0;
        for (i=0; i<nvar; i++) {
            a[i] = 0;
            for (j=0; j<nvar; j++) cmat[i][j] = 0;

        }
     }
     dtime = time[person];
     ndead  = 0;
     deadwt = 0;
     denom2 = 0;
     while (person >=0 && time[person]==dtime) {
       nrisk++;
       zbeta = offset[person];
       for (i=0; i<nvar; i++)
           zbeta += beta[i]*covar[i][person];
       risk = exp(zbeta) * weights[person];
       if (status[person] ==0) {
         denom += risk;
         for (i=0; i<nvar; i++) {
             a[i] += risk*covar[i][person];
             for (j=0; j<=i; j++)
               cmat[i][j] += risk*covar[i][person]*covar[j][person];
         }
       }
       else {
         ndead++;
         deadwt += weights[person];
         denom2 += risk;

         for (i=0; i<nvar; i++) {
           u[i]  +=  weights[person]*covar[i][person];
           a2[i] +=  risk*covar[i][person];
           for (j=0; j<=i; j++)
               cmat2[i][j] += risk*covar[i][person]*covar[j][person];
         }
      }
      person--;
      if (person>=0 && strata[person]==1) break;
    }

    if (ndead >0) {
        wtave = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          for (i=0; i<nvar; i++) {
            a[i] += a2[i]/ndead;
            temp2 = a[i]/denom;
            u[i] -= wtave *temp2;
            for (j=0; j<=i; j++) {
              cmat[i][j] += cmat2[i][j]/ndead;
              imat[j][i] += wtave*(cmat[i][j] - temp2*a[j])/denom;
            }
          }
        }
        for (i=0; i<nvar; i++) {
          a2[i]=0;
          for (j=0; j<nvar; j++) cmat2[i][j]=0;
        }
      }
  }

  for (i=0;i<nvar;i++)
    for (j=0;j<i;j++)
      imat[i][j] = imat[j][i];

  //for (i=0;i<nvar;i++){
  //  u[i] /= nused;
  //  for (j=0;j<nvar;j++)
  //    imat[i][j] /= -nused;
  //}

  PROTECT(rlist= allocVector(VECSXP, 2));
    SET_VECTOR_ELT(rlist, 0, u2);
    SET_VECTOR_ELT(rlist, 1, imat2);

  PROTECT(rlistnames= allocVector(STRSXP, 2));
    SET_STRING_ELT(rlistnames, 0, mkChar("score"));
    SET_STRING_ELT(rlistnames, 1, mkChar("neg_info_mat"));

  UNPROTECT(nprotect+2);
  return rlist;
}
