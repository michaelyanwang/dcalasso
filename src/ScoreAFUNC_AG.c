#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

#define ALLOC(a,b)  R_alloc(a,b)

double **dmatrix2(double *array, int ncol, int nrow){
   int i;
   double **pointer;
   pointer = (double **) ALLOC(nrow, sizeof(double *));
   for (i=0; i<nrow; i++) {
     pointer[i] = array;
     array += ncol;
 	}
  return(pointer);
}

// surv2: Surv(start,stop,status); covar2: covariates; method2: 1 if Efron (always)
// strata2: (not consider strata = n); weights2, offset2 as default is nil
// betahat2: last iteration's beta
SEXP ScoreAFUNC_AG(SEXP surv2,      SEXP covar2,   SEXP strata2,
                   SEXP weights2,   SEXP offset2,
                   SEXP sort12,     SEXP sort22,   SEXP betahat2) {

    int i,j,k, person;
    int indx1, istrat, p, p1;
    int nrisk;
    int nused, nvar;
    int rank, rank2, fail;

    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *oldbeta, *weights, *offset;
    double *scale;
    double *a2, **cmat2;
    double *eta;
    double  denom, zbeta, risk;
    double  dtime;
    double  temp, temp2;
    double  newlk =0;
    int     halving;    /*are we doing step halving at the moment? */
    double  tol_chol;
    double  meanwt;
    int deaths;
    double denom2, etasum;
    int *keep;               /* marker for useless obs */

    /* inputs */
    double *start, *tstop, *event;
    int *sort1, *sort2;
    int *strata, nstrat;
    double method;  /* saving this as double forces some double arithmetic */
    int doscale;

    /* returned objects */
    SEXP imat2, beta2, u2, loglik2;
    double *beta, *u, *loglik;
    SEXP sctest2, flag2, iter2;
    double *sctest;
    int *flag, *iter;
    SEXP rlist, rlistnames;

    int nprotect;  /* number of protect calls I have issued */

    /* get sizes and constants */
    nused = nrows(covar2);
    nvar  = ncols(covar2);
    method= 1;
    nstrat = LENGTH(strata2);

    /* input arguments */
    start = REAL(surv2);
    tstop  = start + nused;
    event = tstop + nused;
    weights = REAL(weights2);
    offset = REAL(offset2);
    sort1  = INTEGER(sort12);
    sort2  = INTEGER(sort22);
    strata = INTEGER(strata2);


    /*
    ** scratch space
    **  nvar: a, a2, oldbeta, scale
    **  nvar*nvar: cmat, cmat2
    **  nused:  eta, keep
    */
    eta = (double *) R_alloc(nused + 4*nvar + 2*nvar*nvar, sizeof(double));
    a = eta + nused;
    a2= a + nvar;
    scale  = a2 + nvar;
    oldbeta = scale + nvar;
    keep = (int *) R_alloc(nused, sizeof(int));


    /*
    **  Set up the ragged arrays
    **  covar2 might not need to be duplicated, even though
    **  we are going to modify it, due to the way this routine was
    **  was called.  In this case NAMED(covar2) will =0
    */
    PROTECT(imat2 = allocVector(REALSXP, nvar*nvar));
    nprotect =1;
    PROTECT(u2 = allocVector(REALSXP, nvar));
    u = REAL(u2);
    nprotect++;

    PROTECT(beta2 = duplicate(betahat2));
    nprotect++;
    if (NAMED(covar2)>0) {
        PROTECT(covar2 = duplicate(covar2));
        nprotect++;
    }
    covar= dmatrix2(REAL(covar2), nused, nvar);
    imat = dmatrix2(REAL(imat2),  nvar, nvar);
    cmat = dmatrix2(oldbeta+ nvar,   nvar, nvar);
    cmat2= dmatrix2(oldbeta+ nvar + nvar*nvar, nvar, nvar);
    beta = REAL(beta2);


    indx1 =0;
    person =0;
    for (k=0; k<nused; k++) keep[k] =1;
    for (istrat=0; istrat<nstrat; istrat++) {
       while(person < strata[istrat]) {
           /* find the next death */
           for (k=person; k< strata[istrat]; k++) {
               p = sort2[k];
               if (event[p] ==1) {
                   dtime = tstop[p];
                   break;
               }
           }
           if (k== strata[istrat]) {
               /* no more deaths in this strata */
               person = k;
               indx1 =k;  /* we can move on */
           }

           for (; indx1 < strata[istrat]; indx1++) {
               p1 = sort1[indx1];
               if (start[p1] < dtime) break;
               keep[p1]--;
           }
           for (; person < strata[istrat]; person++) {
               p = sort2[person];
               if (tstop[p] < dtime) break;
               if (keep[p] ==1) keep[p] =2;
           }
       }
    }

    /* First iteration, which has different ending criteria */
    for (person=0; person<nused; person++) {
        zbeta = 0;      /* form the term beta*z   (vector mult) */
        for (i=0; i<nvar; i++)
            zbeta += beta[i]*covar[i][person];
        eta[person] = zbeta + offset[person];
    }

    /*
    **  'person' walks through the the data from 1 to n,
    **     sort1[0] points to the largest stop time, sort1[1] the next, ...
    **  'dtime' is a scratch variable holding the time of current interest
    **  'indx1' walks through the start times.
    */
    newlk =0;
    for (i=0; i<nvar; i++) {
        u[i] =0;
        for (j=0; j<nvar; j++) imat[i][j] =0;
    }
    person =0;
    indx1 =0;
    istrat =0;

    /* this next set is rezeroed at the start of each stratum */
    denom=0;
    nrisk=0;
    etasum =0;
    for (i=0; i<nvar; i++) {
        a[i] =0;
        for (j=0; j<nvar; j++) cmat[i][j] =0;
    }
    /* end of the per-stratum set */

    while (person < nused) {
        /* find the next death time */
        for (k=person; k< nused; k++) {
            if (k == strata[istrat]) {
                /* hit a new stratum; reset temporary sums */
                istrat++;
                denom = 0;
                nrisk = 0;
                etasum =0;
                for (i=0; i<nvar; i++) {
                    a[i] =0;
                    for (j=0; j<nvar; j++) cmat[i][j] =0;
                }
                person =k;  /* skip to end of stratum */
                indx1  =k;
            }
            p = sort2[k];
            if (event[p] == 1) {
                dtime = tstop[p];
                break;
            }
        }
        if (k == nused) person =k;  /* no more deaths to be processed */
        else {
            /* remove any subjects no longer at risk */
            /*
            ** subtract out the subjects whose start time is to the right
            ** If everyone is removed reset the totals to zero.  (This happens when
            ** the survSplit function is used, so it is worth checking).
            */
            for (; indx1<strata[istrat]; indx1++) {
                p1 = sort1[indx1];
                if (start[p1] < dtime) break;
                if (keep[p1] == 0) continue;  /* skip any never-at-risk rows */
                nrisk--;
                if (nrisk ==0) {
                    etasum =0;
                    denom =0;
                    for (i=0; i<nvar; i++) {
                        a[i] =0;
                        for (j=0; j<=i; j++) cmat[i][j] =0;
                    }
                }
                else {
                    etasum -= eta[p1];
                    risk = exp(eta[p1]) * weights[p1];
                    denom -= risk;
                    for (i=0; i<nvar; i++) {
                        a[i] -= risk*covar[i][p1];
                        for (j=0; j<=i; j++)
                            cmat[i][j] -= risk*covar[i][p1]*covar[j][p1];
                    }
                }
                /*
                ** We must avoid overflow in the exp function (~750 on Intel)
                ** and want to act well before that, but not take action very often.
                ** One of the case-cohort papers suggests an offset of -100 meaning
                ** that etas of 50-100 can occur in "ok" data, so make it larger
                ** than this.
                ** If the range of eta is more then log(1e16) = 37 then the data is
                **  hopeless: some observations will have effectively 0 weight.  Keeping
                **  the mean sensible suffices to keep the max in check for all other
                *   data sets.
                */
                if (fabs(etasum/nrisk) > 200) {
                    flag[1]++;  /* a count, for debugging/profiling purposes */
                    temp = etasum/nrisk;
                    for (i=0; i<nused; i++) eta[i] -= temp;
                    temp = exp(-temp);
                    denom *= temp;
                    for (i=0; i<nvar; i++) {
                        a[i] *= temp;
                        for (j=0; j<nvar; j++) {
                            cmat[i][j]*= temp;
                        }
                    }
                    etasum =0;
                }
            }

            /*
            ** add any new subjects who are at risk
            ** denom2, a2, cmat2, meanwt and deaths count only the deaths
            */
            denom2= 0;
            meanwt =0;
            deaths=0;
            for (i=0; i<nvar; i++) {
                a2[i]=0;
                for (j=0; j<nvar; j++) {
                    cmat2[i][j]=0;
                }
            }

            for (; person<strata[istrat]; person++) {
                p = sort2[person];
                if (tstop[p] < dtime) break; /* no more to add */
                risk = exp(eta[p]) * weights[p];

                if (event[p] ==1 ){
                    nrisk++;
                    etasum += eta[p];
                    deaths++;
                    denom2 += risk*event[p];
                    meanwt += weights[p];
                    newlk += weights[p]* eta[p];
                    for (i=0; i<nvar; i++) {
                        u[i] += weights[p] * covar[i][p];
                        a2[i]+= risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat2[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                }
                else if (keep[p] >0) {
                    nrisk++;
                    etasum += eta[p];
                    denom += risk;
                    for (i=0; i<nvar; i++) {
                        a[i] += risk*covar[i][p];
                        for (j=0; j<=i; j++)
                            cmat[i][j] += risk*covar[i][p]*covar[j][p];
                    }
                }
            }
            /*
            ** Add results into u and imat for all events at this time point
            */
            if (method==0 || deaths ==1) { /*Breslow */
                denom += denom2;
                newlk -= meanwt*log(denom);  /* sum of death weights*/
                for (i=0; i<nvar; i++) {
                    a[i] += a2[i];
                    temp = a[i]/denom;   /*mean covariate at this time */
                    u[i] -= meanwt*temp;
                    for (j=0; j<=i; j++) {
                        cmat[i][j] += cmat2[i][j];
                        imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
                    }
                }
            }
            else {
                meanwt /= deaths;
                for (k=0; k<deaths; k++) {
                    denom += denom2/deaths;
                    newlk -= meanwt*log(denom);
                    for (i=0; i<nvar; i++) {
                        a[i] += a2[i]/deaths;
                        temp = a[i]/denom;
                        u[i] -= meanwt*temp;
                        for (j=0; j<=i; j++) {
                            cmat[i][j] += cmat2[i][j]/deaths;
                            imat[j][i] += meanwt*((cmat[i][j]- temp*a[j])/denom);
                        }
                        }
                }
            }
            /*
            ** We must avoid overflow in the exp function (~750 on Intel)
            ** and want to act well before that, but not take action very often.
            ** One of the case-cohort papers suggests an offset of -100 meaning
            ** that etas of 50-100 can occur in "ok" data, so make it larger
            ** than this.
            ** If the range of eta is more then log(1e16) = 37 then the data is
            **  hopeless: some observations will have effectively 0 weight.  Keeping
            **  the mean sensible suffices to keep the max in check for all other
            *   data sets.
            */
            if (fabs(etasum/nrisk) > 200) {
                flag[1]++;  /* a count, for debugging/profiling purposes */
                temp = etasum/nrisk;
                for (i=0; i<nused; i++) eta[i] -= temp;
                temp = exp(-temp);
                denom *= temp;
                for (i=0; i<nvar; i++) {
                    a[i] *= temp;
                    for (j=0; j<nvar; j++) {
                        cmat[i][j]*= temp;
                    }
                }
                etasum =0;
            }
        }
    }   /* end  of accumulation loop */


    for (i=0;i<nvar;i++)
        for (j=0;j<i;j++)
            imat[i][j] = imat[j][i];

//    for (i=0;i<nvar;i++){
//      u[i] /= nused;
//      for (j=0;j<nvar;j++)
//        imat[i][j] /= -nused;
//    }

    PROTECT(rlist= allocVector(VECSXP, 2));
    SET_VECTOR_ELT(rlist, 0, u2);
    SET_VECTOR_ELT(rlist, 1, imat2);

    PROTECT(rlistnames= allocVector(STRSXP, 2));
    SET_STRING_ELT(rlistnames, 0, mkChar("score"));
    SET_STRING_ELT(rlistnames, 1, mkChar("neg_info_mat"));

    UNPROTECT(nprotect+2);
    return(rlist);
}
