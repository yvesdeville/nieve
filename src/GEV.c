#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#define NODEBUG

/* ===========================================================================
   AUTHOR Yves Deville <deville.yves@alpestat.com>

   NOTE This program is part of the 'nieve' R package. See the LICENSE
   of the package.

   GOAL Compute the probability functions related to the Generalised
   Extreme Value (GEV) distribution, possibly including the gradient
   and the Hessian (for the density). The implementation in C is
   faster than a pure R implementation. Unlike some other
   implementations in existing R packages, NAs are returned when the
   shape is negative, which can be a desirable behaviour for when
   unconstrained optimisation is used to maximise the log-likelihood.
   =========================================================================== */  


/* ==========================================================================
 * Density function
 * ========================================================================== */

SEXP Call_dGEV(SEXP x,             /*  double                          */
	       SEXP loc,           /*  double                          */
	       SEXP scale,         /*  double                          */
	       SEXP shape,         /*  double                          */
	       SEXP logFlag,       /*  integer                         */
	       SEXP derivFlag,     /*  integer                         */ 
	       SEXP hessianFlag) { /*  integer                         */
 
  int n, nx, nloc, nscale, nshape, i, ix, iloc, iscale, ishape,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0];
  
  double eps = 2e-4, z, emz, V, xi, xiz, dlogfdz;
  
  SEXP val;
  
  PROTECT(x = coerceVector(x, REALSXP));
  PROTECT(loc = coerceVector(loc, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));
  PROTECT(shape = coerceVector(shape, REALSXP));

  double *rx = REAL(x), *rloc = REAL(loc), *rscale = REAL(scale), 
    *rshape = REAL(shape);

  nx = LENGTH(x);						
  nloc = LENGTH(loc);						
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);
  
  if ((nx == 0) || (nloc == 0) || (nscale == 0) || (nshape == 0)) {
    UNPROTECT(4);
    return(allocVector(REALSXP, 0));
  }
  
  n = nx;							
  if (n < nloc) n = nloc;						
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);
  
  if (deriv) {

    double U, W;
    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n * 3));
    double *rgrad = REAL(grad);

    /* =====================================================================
       NOTE: although some variables that are used only when the
       hessian is required and hence could be declared only in this
       case, it seems that this does not work. It seems that the
       compilator does not see declarations inside two nested 'if'.
       ===================================================================== */

    // if (hessian) { 

    SEXP hess;
    
    PROTECT(hess = allocVector(REALSXP, n * 3 * 3));
    double *rhess = REAL(hess);
    
    // specific auxiliary variables for hessian
    int j, k;
    double L, zVm1, zVm2, L2, sigma, sigma2, sigmaV;
    
    // }
    
    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

    for (i = ix = iloc = iscale = ishape = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	       
	 iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
	
	rgrad[i] = NA_REAL;
	rgrad[i + n] = NA_REAL;
	rgrad[i + 2 * n] = NA_REAL; 

	if (hessian) {
	  // row 'mu'
	  rhess[i] = NA_REAL;
	  rhess[i + n] = NA_REAL;
	  rhess[i + 3 * n] = NA_REAL;
	  rhess[i + 2 * n] = NA_REAL;
	  rhess[i + 6 * n] = NA_REAL;
	  // row 'sigma'
	  rhess[i + 4 * n] = NA_REAL;
	  rhess[i + 5 * n] = NA_REAL;
	  rhess[i + 7 * n] = NA_REAL;
	  // row 'xi'
	  rhess[i + 8 * n] = NA_REAL;
	}
	
      } else {
	
	z = (rx[ix] - rloc[iloc]) / rscale[iscale];
	xi = rshape[ishape];
	sigma = rscale[iscale]; 
	sigma2 = sigma * sigma;

	if (fabs(xi) < eps) {

	  if (z < -30.0) {
	    
	    rgrad[i] = 0.0;
	    rgrad[i + n] =  0.0;
	    rgrad[i + 2 * n] = 0.0;
	    
	    if (hessian) {
	      
	      // row 'mu'
	      rhess[i] =  0.0;
	      rhess[i + n] = 0.0;
	      rhess[i + 3 * n] = 0.0;
	      rhess[i + 2 * n] = 0.0;
	      rhess[i + 6 * n] = 0.0;
	      // row 'sigma'
	      rhess[i + 4 * n] = 0.0;
	      rhess[i + 5 * n] = 0.0;
	      rhess[i + 7 * n] = 0.0;
	      // row 'xi'
	      rhess[i + 8 * n] = 0.0;
	    
	    }

	  } else {
	    emz = exp(-z);

	    // improved approximation based on Taylor expansion
	    xiz = xi * z;
	    rval[i] = -log(rscale[iscale]) - z - emz +
	      xiz * (0.5 * ((1.0 - emz) * z - 2.0) -
		     xiz * (-12.0 + 8.0 * (1.0 - emz) * z + 3.0 * z * z * emz) / 24.0);
	    dlogfdz = -(1.0 - emz) - 0.5 * (2.0 - 2.0 * (1.0 - emz) * z - emz * z * z) * xi; 
	    rgrad[i] = - dlogfdz / sigma;
	    rgrad[i + n] =  (-1.0 - z * dlogfdz) / sigma;
	    rgrad[i + 2 * n] = z * z * (1.0 - emz) / 2.0 - z +
	      (12.0 - 8.0 * (1.0 - emz) * z - 3.0 * emz * z * z) * xiz * z / 12.0;
	    
	    if (hessian) {
	      
	      // row 'mu'
	      rhess[i] =  - emz / sigma2;
	      rhess[i + n] = - (1.0 + (z - 1.0) * emz) / sigma2;
	      rhess[i + 3 * n] = rhess[i + n];
	      rhess[i + 2 * n] =  (1 - z +  z * (1 - z / 2.0) * emz) / sigma;
	      rhess[i + 6 * n] = rhess[i + 2 * n];
	      // row 'sigma'
	      rhess[i + 4 * n] =  (((2.0 - z) * emz - 2.0) * z + 1.0) / sigma2;
	      rhess[i + 5 * n] = (z * (1.0 - z) +  z * z * (1.0 - z / 2.0) * emz) / sigma; 
	      rhess[i + 7 * n] = rhess[i + 5 * n];
	      // row 'xi'
	      rhess[i + 8 * n] =  z * z * 
		(1.0 - 2.0 * z / 3 + emz * z * (2.0 / 3 + - z / 4.0));
	    
	    }
	  }
	  
	} else {
	  
	  V = 1.0 + xi * z;
	  sigmaV = sigma * V;

	  // Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, V);
	  
	  if (V > 0.0) {
	
	    rval[i] = -log(rscale[iscale]) - pow(V, - 1.0 / xi) -
	      (1.0 / xi + 1.0) * log(V);	  
	    
	    L = log(V) / xi / xi;
	    W = pow(V, -1.0 / xi);
	    U = (1.0 + xi - W) / V / sigma;
	    zVm1 = z / V;
	    zVm2 = zVm1 * zVm1;
	    L2 = L - zVm1 / xi;

	    rgrad[i] = U;
	    rgrad[i + n] = -1.0 / sigma + z * U;
	    rgrad[i + 2 * n] = (1.0 - W) * L - z * U * sigma / xi;
	    
	    if (hessian) {
	      // row 'mu'
	      rhess[i] =  - (xi + 1.0) * (W - xi) / sigmaV / sigmaV;
	      rhess[i + n] = ((1.0 - (xi + 1.0) * zVm1) * (W - xi) - 1.0) / V / sigma2;
	      rhess[i + 3 * n] = rhess[i + n];
	      rhess[i + 2 * n] = (W * (- L2 + zVm1) + 1.0 - (xi + 1.0) * zVm1) / sigmaV;
	      rhess[i + 6 * n] = rhess[i + 2 * n];
	      // row 'sigma'
	      rhess[i + 4 * n] = (((W - xi) * (2.0 - (xi + 1.0) * zVm1) - 2.0) *  zVm1 + 1.0) / sigma2;
	      rhess[i + 5 * n] = (z * W * (-L2 + zVm1) / V  + zVm1 - (xi + 1.0) * zVm2) / sigma;
	      rhess[i + 7 * n] = rhess[i + 5 * n];
	      // row 'xi'
	      rhess[i + 8 * n] = - (L2 * L2  + (-2.0 * L2 + zVm2) / xi) * W +
		(- 2.0 * L2 + (xi + 1.0) * zVm2) / xi;

	    }
	    
	  } else {

	    rval[i] = R_NegInf;
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;
	    rgrad[i + 2 * n] = 0.0;

	    if (hessian) {
	      // row 'mu'
	      rhess[i] = 0.0;
	      rhess[i + n] = 0.0;
	      rhess[i + 3 * n] = 0.0;
	      rhess[i + 2 * n] = 0.0;
	      rhess[i + 6 * n] = 0.0;
	      // row 'sigma'
	      rhess[i + 4 * n] = 0.0;
	      rhess[i + 5 * n] = 0.0;
	      rhess[i + 7 * n] = 0.0;
	      // row 'xi'
	      rhess[i + 8 * n] = 0.0;
	    }

	  }

	} /* non-Gumbel case */
	
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	  rgrad[i] *= rval[i];
	  rgrad[i + n] *= rval[i];
	  rgrad[i + 2 * n] *= rval[i];

	  if (hessian) {
	    // multiply the hessian by the density and add the
	    // 'tcrossprod' of the gradient
	    for (j = 0; j < 3; j++) {
	      for (k = 0; k < 3; k++) {
		rhess[i + (j + k * 3) * n] *= rval[i];
		rhess[i + (j + k * 3) * n] += rgrad[i + j * n] * rgrad[i + k * n];
	      }
	    }
	  }
	  
	}
	
      }   /* non-NA case     */
      
    }    /* loop */

    SET_ATTR(val, attrNm, grad);

    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    UNPROTECT(8);

    return(val);
    
  } else {

    for (i = ix = iloc = iscale = ishape = 0;  i < n; 
	 ix = (++ix == nx) ? 0 : ix, 	       
	 iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rx[ix]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
		
      } else {
	
	z = (rx[ix] - rloc[iloc]) / rscale[iscale];
	xi = rshape[ishape];

	if (fabs(xi) < eps) {
	  
	  emz = exp(-z);
	  // rval[i] = -log(rscale[iscale]) - z - emz;	  
	  // improved approximation base on Taylor expansion
	  xiz = xi * z;
	  rval[i] = -log(rscale[iscale]) - z - emz +
	    xiz * (0.5 * ((1.0 - emz) * z - 2.0) -
		   xiz * (-12.0 + 8.0 * (1.0 - emz) * z + 3.0 * z * z * emz) / 24.0);
	} else {
	  
	  V = 1.0 + xi * z;
	  // Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, V);

	  if (V > 0.0) {
	    rval[i] = -log(rscale[iscale]) - pow(V, - 1.0 / xi) -
	      (1.0 / xi + 1.0) * log(V);
	  } else {
	    rval[i] = R_NegInf;
	  }
	  
	  
	} /* non-Gumbel case */
	
	if (!INTEGER(logFlag)[0]) {
	  rval[i] = exp(rval[i]);
	}
	
      }   /* non-NA case     */
      
    }
  
    UNPROTECT(5);
    return(val);
    
  }

}

/* ==========================================================================
 * Distribution function
 * ========================================================================== */

SEXP Call_pGEV(SEXP q,               /*  double                          */
	       SEXP loc,             /*  double                          */
	       SEXP scale,           /*  double                          */
	       SEXP shape,           /*  double                          */
	       SEXP lowerTailFlag,   /*  integer                         */
	       SEXP derivFlag) {     /*  integer                         */
 
  int n, nq, nloc, nscale, nshape, i, iq, iloc, iscale, ishape,
    deriv = INTEGER(derivFlag)[0];
  
  double eps = 2e-4, z, emz, V, Z, W, xi, ee, eemz = 1.0, dFdz = 0.0;
  
  SEXP val;
  
  PROTECT(q = coerceVector(q, REALSXP));
  PROTECT(loc = coerceVector(loc, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));
  PROTECT(shape = coerceVector(shape, REALSXP));

  double *rq = REAL(q), *rloc = REAL(loc), *rscale = REAL(scale), 
    *rshape = REAL(shape);

  nq = LENGTH(q);						
  nloc = LENGTH(loc);						
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);
  
  if ((nq == 0) || (nloc == 0) || (nscale == 0) || (nshape == 0)) {
    UNPROTECT(4);
    return(allocVector(REALSXP, 0));
  }
  
  n = nq;							
  if (n < nloc) n = nloc;						
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);
 

  if (deriv) {

    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n * 3));
    double *rgrad = REAL(grad);

    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

    for (i = iq = iloc = iscale = ishape = 0;  i < n; 
	 iq = (++iq == nq) ? 0 : iq, 	       
	 iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rq[iq]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
	
	rgrad[i] = NA_REAL;
	rgrad[i + n] = NA_REAL;
	rgrad[i + 2 * n] = NA_REAL; 
	

      } else {
	
	z = (rq[iq] - rloc[iloc]) / rscale[iscale];
	xi = rshape[ishape];

	if (fabs(xi) < eps) {
   
	  if (z < -30.0) {
	    rval[i] = 0.0;
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;
	    rgrad[i + 2 * n] = 0.0;
	  } else {
	    
	    // improved approximation
	    emz = exp(-z);
	    ee = exp(-emz);
	    rval[i] = ee * (1.0 + emz * z * z * xi *
			    (-0.5 + (8.0 - 3.0 * (1.0 - emz) * z) * xi * z / 24.0));
	    eemz = ee * emz;
	    dFdz = eemz * (1.0 +  0.5 * z * ( -2.0  + (1.0 - emz) * z) * xi);
	    rgrad[i] = - dFdz / rscale[iscale];
	    rgrad[i + n] =  -z * dFdz / rscale[iscale];
	    rgrad[i + 2 * n] = eemz * (z * z * (-0.5  + z * (8.0 - 3.0 * (1.0 - emz) * z) * xi / 12.0));
	  }

	  // Rprintf("%4d, %4d, %6.3f, %7.2f, %6.3f, %6.3f", i, ishape, xi, z, Z, rscale[iscale]);
	  // Rprintf(" grad: %6.3f, %6.3f, %6.3f\n", rgrad[i], rgrad[i + n], rgrad[i + 2 * n]);

	} else {
	  
	  V = 1.0 + xi * z;
	  // Rprintf("%d, %d, %6.3f, %6.3f\n", i, ishape, xi, V);

	  if (V > 0.0) {

	    W = pow(V, - 1.0 / xi);
	    Z = exp(-W);

	    rval[i] = Z;
	    Z *= W;
	   
	    rgrad[i] = -Z / V / rscale[iscale];
	    rgrad[i + n] = z * rgrad[i];
	    rgrad[i + 2 * n] = -Z * (log(V) / xi - z / V) / xi;
	  
	  } else {
	    if (xi > 0.0) {
	      rval[i] = 0.0;
	    } else {
	      rval[i] = 1.0;
	    }
	    rgrad[i] = 0.0;
	    rgrad[i + n] = 0.0;
	    rgrad[i + 2 * n] = 0.0;

	  }

	} /* non-Gumbel case */
	
	if (!INTEGER(lowerTailFlag)[0]) {
	  rval[i] = 1.0 - rval[i];
	  rgrad[i] = -rval[i];
	  rgrad[i + n] = -rgrad[i + n];
	  rgrad[i + 2 * n] = -rgrad[i + 2 * n];
	}
	
      }   /* non-NA case     */
      
    }
  
    SET_ATTR(val, attrNm, grad);
    UNPROTECT(7);
    return(val);
    
  } else {

    for (i = iq = iloc = iscale = ishape = 0;  i < n; 
	 iq = (++iq == nq) ? 0 : iq, 	       
	 iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rq[iq]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
		
      } else {
	
	z = (rq[iq ] - rloc[iloc]) / rscale[iscale];
	xi = rshape[ishape];

	if (fabs(xi) < eps) {
	  
	  // improved approximation
	  emz = exp(-z);
	  ee = exp(-emz);
	  rval[i] = ee * (1.0 + emz * z * z * xi *
			  (-0.5 + (8.0 - 3.0 * (1.0 - emz) * z) * xi * z / 24.0));
	} else {
	  
	  V = 1.0 + xi * z;
	  
	  if (V > 0.0) { 
	    W = pow(V, - 1.0 / xi);
	    rval[i] = exp(-W);
	  } else {
	    if (xi > 0.0) {
	      rval[i] = 0.0;
	    } else {
	      rval[i] = 1.0;
	    }
	  }
	  
	  
	} /* non-Gumbel case */
	
	if (!INTEGER(lowerTailFlag)[0]) {
	  rval[i] = 1.0 - rval[i];
	}
	
      }   /* non-NA case     */
      
    }
  
    UNPROTECT(5);
    return(val);
    
  }

}


/* ==========================================================================
 * Quantile function
 * ========================================================================== */

SEXP Call_qGEV(SEXP p,               /*  double                          */
	       SEXP loc,             /*  double                          */
	       SEXP scale,           /*  double                          */
	       SEXP shape,           /*  double                          */
	       SEXP lowerTailFlag,   /*  integer                         */
	       SEXP derivFlag,        /* integer                          */
               SEXP hessianFlag) {   /*  integer                         */
 
  int n, np, nloc, nscale, nshape, i, ip, iloc, iscale, ishape, j,
    deriv = INTEGER(derivFlag)[0], hessian = INTEGER(hessianFlag)[0],
    lowerTail = INTEGER(lowerTailFlag)[0];
    
  double eps = 2e-4, A, logA, V, xi, prob, mu, sigma;
  
  SEXP val;
  
  PROTECT(p = coerceVector(p, REALSXP));
  PROTECT(loc = coerceVector(loc, REALSXP));
  PROTECT(scale = coerceVector(scale, REALSXP));
  PROTECT(shape = coerceVector(shape, REALSXP));

  double *rp = REAL(p), *rloc = REAL(loc), *rscale = REAL(scale), 
    *rshape = REAL(shape);

  np = LENGTH(p);						
  nloc = LENGTH(loc);						
  nscale = LENGTH(scale);		
  nshape = LENGTH(shape);
  
  if ((np == 0) || (nloc == 0) || (nscale == 0) || (nshape == 0)) {
    UNPROTECT(4);
    return(allocVector(REALSXP, 0));				
  }
  
  n = np;							
  if (n < nloc) n = nloc;						
  if (n < nscale) n = nscale;
  if (n < nshape) n = nshape;
  
  PROTECT(val = allocVector(REALSXP, n));
  double *rval = REAL(val);

  if (deriv) {

    SEXP grad, attrNm;

    PROTECT(grad = allocVector(REALSXP, n * 3));
    double *rgrad = REAL(grad);

    PROTECT(attrNm = NEW_CHARACTER(1)); 
    SET_STRING_ELT(attrNm, 0, mkChar("gradient"));

 // if (hessian) { 

    SEXP hess;
    
    PROTECT(hess = allocVector(REALSXP, n * 3 * 3));
    double *rhess = REAL(hess);
    double dV, d2V;
    
    // }
   
    for (i = ip = iloc = iscale = ishape = 0;  i < n; 
	 ip = (++ip == np) ? 0 : ip, 	       
	 iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      
      if (ISNA(rp[ip]) || (rscale[iscale] <= 0.0)) {
	// Rprintf("NA case\n");

	rval[i] = NA_REAL;
	
	for (j = 0; j < 3; j++) {
	  rgrad[i + j * n] = NA_REAL;
	}
	  
	if (hessian) {
	  for (j = 0; j < 9; j++) {
	    rhess[i + j * n] = NA_REAL;
	  }
	}
	
      } else if (((rp[ip] == 0.0) && lowerTail) ||
		 ((rp[ip] == 1.0) && !lowerTail)) {

	xi = rshape[ishape];
	
	if (xi > 0.0) {

	  sigma = rscale[iscale];
	  mu = rloc[iloc];
	  
	  rval[i] = mu - sigma / xi;
	
	  rgrad[i] = 1.0;
	  rgrad[i + n] = - 1.0 / xi;
	  rgrad[i + 2 * n] = sigma / xi / xi;
	  
	  if (hessian) {
	    rhess[i] = 0.0;
	    rhess[i + n] = 0.0;
	    rhess[i + 2 * n] = 0.0;
	    rhess[i + 3 * n] = 0.0;
	    rhess[i + 4 * n] = 0.0;
	    rhess[i + 6 * n] = 0.0;
	    // row 'sigma'
	    rhess[i + 5 * n] = 1.0 / xi / xi;
	    rhess[i + 7 * n] = rhess[i + 5 * n];
	    // row 'xi'
	    rhess[i + 8 * n] = - 2.0 * sigma / xi / xi / xi;
	  }
	  
	} else {
	  
	  rval[i] = R_NegInf;
	  
	  for (j = 0; j < 3; j++) {
	    rgrad[i + j * n] = NA_REAL;
	  }
	  
	  if (hessian) {
	    for (j = 0; j < 9; j++) {
	      rhess[i + j * n] = NA_REAL;
	    }
	  }

	}
	
      } else if (((rp[ip] == 1.0) && lowerTail) ||
		 ((rp[ip] == 0.0) && !lowerTail)) {
	
	xi = rshape[ishape];

	if (xi >= 0.0) {
	  
	  rval[i] = R_PosInf;
	  
	  for (j = 0; j < 3; j++) {
	    rgrad[i + j * n] = NA_REAL;
	  }
	  
	  if (hessian) {
	    for (j = 0; j < 9; j++) {
	      rhess[i + j * n] = NA_REAL;
	    }
	  }
	  
	} else {

	  xi = rshape[ishape];
	  sigma = rscale[iscale];
	  mu = rloc[iloc];
	  
	  rval[i] =  mu - sigma / xi;
	  
	  rgrad[i] =  1.0;
	  rgrad[i + n] =  -1.0 / xi;
	  rgrad[i + 2 * n] =  sigma / xi / xi;
	  
	  if (hessian) {
	    rhess[i] = 0.0;
	    rhess[i + n] = 0.0;
	    rhess[i + 2 * n] = 0.0;
	    rhess[i + 3 * n] = 0.0;
	    rhess[i + 4 * n] = 0.0;
	    rhess[i + 6 * n] = 0.0;
	    // row sigma
	    rhess[i + 5 * n] = 1.0 / xi / xi;
	    rhess[i + 7 * n] = rhess[i + 5 * n];
	    // row 'xi'
	    rhess[i + 8 * n] = - 2.0 * sigma / xi / xi / xi;
	  } 
	    
	}
	
      } else {

	if (!lowerTail) {
	  prob = 1.0 - rp[ip];
	} else {
	  prob = rp[ip];
	}

	xi = rshape[ishape];
	A = -log(prob);
	logA = log(A);

	if (fabs(xi) < eps) {

	  // improved approximation base on Taylor expansion
	  
	  rval[i] = rloc[iloc] + rscale[iscale] * logA *
	    (-1.0 + logA * xi * (0.5 - logA * xi / 6.0)) ;

	  rgrad[i] = 1.0;
	  rgrad[i + n] = logA * (-1.0  + 0.5 * logA * xi);
	  rgrad[i + 2 * n] = rscale[iscale] * logA * logA * (0.5 - logA * xi / 3.0);

	  if (hessian) {
	    rhess[i] = 0.0;
	    rhess[i + n] = 0.0;
	    rhess[i + 2 * n] = 0.0;
	    rhess[i + 3 * n] = 0.0;
	    rhess[i + 4 * n] = 0.0;
	    rhess[i + 6 * n] = 0.0;
	    // row 'sigma'
	    rhess[i + 5 * n] = logA * logA / 2.0;
	    rhess[i + 7 * n] = rhess[i + 5 * n];
	    // row 'xi'
	    rhess[i + 8 * n] = - logA * logA * logA * rscale[iscale] / 3.0;
	  }

	} else {
	  
	  V = (1.0 - pow(A, - xi)) / xi;
	  rval[i] = rloc[iloc] - rscale[iscale] * V;  
	  
	  rgrad[i] = 1.0;
	  rgrad[i + n] = -V;
	  rgrad[i + 2 * n] = rscale[iscale] * (V - logA * (-xi * V + 1.0)) / xi;

	  if (hessian) {
	    dV = - (V - logA) / xi - V * logA;
            d2V =  (V - logA) / xi / xi - dV * (logA + 1.0 / xi);
	    rhess[i] = 0.0;
	    rhess[i + n] = 0.0;
	    rhess[i + 2 * n] = 0.0;
	    rhess[i + 3 * n] = 0.0;
	    rhess[i + 4 * n] = 0.0;
	    rhess[i + 6 * n] = 0.0;
	    // row 'sigma'
	    rhess[i + 5 * n] = -dV;
	    rhess[i + 7 * n] = rhess[i + 5 * n];
	    // row 'xi'
	    rhess[i + 8 * n] = - rscale[iscale] * d2V;
	  }


	} /* non-Gumbel case */

      }   /* non-NA case     */
      
    }     /* loop            */

    SET_ATTR(val, attrNm, grad);
    
    if (hessian) {
      SET_STRING_ELT(attrNm, 0, mkChar("hessian"));
      SET_ATTR(val, attrNm, hess);
    } 

    UNPROTECT(8);
    return(val);
    
  } else {

    for (i = ip = iloc = iscale = ishape = 0;  i < n; 
	 ip = (++ip == np) ? 0 : ip, 	       
	   iloc = (++iloc == nloc) ? 0 : iloc, 
	   iscale = (++iscale == nscale) ? 0 : iscale,
	   ishape = (++ishape == nshape) ? 0 : ishape,
	   ++i) {
      
      if (ISNA(rp[ip]) || (rscale[iscale] <= 0.0)) {
	
	rval[i] = NA_REAL;
	
      } else if (((rp[ip] == 0.0) && lowerTail) ||
		 ((rp[ip] == 1.0) && !lowerTail)) {

	xi = rshape[ishape];
	
	if (xi > 0.0) {
	  rval[i] = rloc[iloc] - rscale[iscale] / rshape[ishape];
	} else {
	  rval[i] = R_NegInf;
	}
	
      } else if (((rp[ip] == 1.0) && lowerTail) ||
		 ((rp[ip] == 0.0 && !lowerTail))) {
	
	xi = rshape[ishape];

	if (xi >= 0.0) {	  
	  rval[i] = R_PosInf;
	} else {
	  rval[i] = rloc[iloc] - rscale[iscale] / rshape[ishape];
	}
	
      } else {
	
	if (!lowerTail) {
	  prob = 1.0 - rp[ip];
	} else {
	  prob = rp[ip];
	}
	
	xi = rshape[ishape];
	A = -log(prob);

	if (fabs(xi) < eps) {

	  logA = log(A);

	  // rval[i] = rloc[iloc] - rscale[iscale] * logA;
	  
	  // improved approximation
	  rval[i] = rloc[iloc] + rscale[iscale] * logA *
	    (-1.0 + logA * xi * (0.5 - logA * xi / 6.0));
	} else {
	  V = (1.0 - pow(A, - xi)) / xi;
	  rval[i] = rloc[iloc] - rscale[iscale] * V;  
       	  	  
	} /* non-Gumbel case */
	
      }   /* non-NA case     */
      
    }     /* loop            */
  
    UNPROTECT(5);
    return(val);
    
  }
}
