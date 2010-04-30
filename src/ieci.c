/**************************************************************************** 
 *
 * Particle Learning of Gaussian Processes
 * Copyright (C) 2010, University of Cambridge
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
 * 02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
 *
 ****************************************************************************/


#include "matrix.h"
#include "linalg.h"
#include "rhelp.h"
#include <stdlib.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>

/* defined in covar.c */
void covar(const int col, double **X1, const int n1, double **X2,
	   const int n2, double d, double g, double **K);
void covar_sep(const int col, double **X1, const int n1, double **X2,
	       const int n2, double *d, double g, double **K);

/*
 * EI:
 *
 * calculates the expected improvement following
 * Williams et al by integrating over the parameters
 * to the GP predictive
 */

double EI(const double m, const double s2, const int df, 
	  const double fmin)
{
  double diff, sd, diffs, scale, ei;

  diff = fmin - m;
  sd = sqrt(s2);
  diffs = diff/sd;
  scale = (df*sd + sq(diff)/sd)/(df-1.0);
  ei = diff*pt(diffs, (double) df, 1, 0);
  ei += scale*dt(diffs, (double) df, 0);

  return(ei);
}


/*
 * calc_ecis:
 *
 * function that iterates over the m Xref locations, and the
 * statis calculated by previous calc_* function in order to 
 * use the EI function to calculate the ECI relative to the
 * reference locations -- writes over ktKik input
 */

void calc_ecis(const int m, double *ktKik, double *s2p, const double phi, 
	       const double g, double *badj, double *tm, const double tdf, 
	       const double fmin)
{
  int i;
  double zphi, ts2;

  for(i=0; i<m; i++) {
    zphi = (s2p[1] + phi)*(1.0 + g - ktKik[i]);
    ts2 = badj[i] * zphi / (s2p[0] + tdf);
    ktKik[i] = EI(tm[i], ts2, tdf, fmin);
  }
}


/*
 * calc_ieci:
 *
 * function that iterates over the m Xref locations, and the
 * statis calculated by previous calc_* function in order to 
 * use the EI function to calculate the ECI relative to the
 * reference locations, and then averaves to create the IECI
 */

double calc_ieci(const int m, double *ktKik, double *s2p, const double phi, 
		 const double g, double *badj, double *tm, const double tdf, 
		 const double fmin, double *w)
{
  int i;
  double zphi, ts2, eci, ieci;

  ieci = 0.0;
  for(i=0; i<m; i++) {
    zphi = (s2p[1] + phi)*(1.0 + g - ktKik[i]);
    ts2 = badj[i] * zphi / (s2p[0] + tdf);
    eci = EI(tm[i], ts2, tdf, fmin);
    if(w) ieci += w[i]*eci;
    else ieci += eci;
  }

  return (ieci/m);
}


/*
 * calc_g_mui_kxy:
 *
 * function for calculating the g vector, mui scalar, and
 * kxy vector for the IECI calculation
 */

void calc_g_mui_kxy(const int col, double *x, double **X, 
		    const int n, double **Ki, double **Xref, 
		    const int m, double *d, const int dlen, 
		    const double g, double *gvec, double *mui, 
		    double *kxy)
{
  double *kx;
  double mu_neg;
  int i;

  /* allocate kx */
  kx = new_vector(n);

  if(dlen == 1) { /* isotropic */
    /* kx <- drop(covar(X1=pall$X, X2=x, d=Zt$d, g=Zt$g)) */
    covar(col, &x, 1, X, n, *d, g, &kx);
    /* kxy <- drop(covar(X1=x, X2=Xref, d=Zt$d, g=0)) */
    covar(col, &x, 1, Xref, m, *d, 0.0, &kxy);
  } else {
    assert(dlen == col);
    /* kx <- drop(covar.sep(X1=pall$X, X2=x, d=Zt$d, g=Zt$g)) */
    covar_sep(col, &x, 1, X, n, d, g, &kx);
    /* kxy <- drop(covar.sep(X1=x, X2=Xref, d=Zt$d, g=0)) */
    covar_sep(col, &x, 1, Xref, m, d, 0.0, &kxy);
  }

  /* Kikx <- drop(util$Ki %*% kx) stored in gvex */
  linalg_dsymv(n,1.0,Ki,n,kx,1,0.0,gvec,1);

  /* mui <- drop(1 + Zt$g - t(kx) %*% Kikx) */
  *mui = 1.0 + g - linalg_ddot(n, kx, 1, gvec, 1);
  
  /* gvec <- - Kikx/mui */
  mu_neg = 0.0 - 1.0/(*mui);
  for(i=0; i<n; i++) gvec[i] *= mu_neg;

  /* clean up */
  free(kx);
}


/*
 * calc_ktKikx:
 *
 * function for calculating the ktKikx vector used in the
 * IECI calculation -- writes over the KtKik input
 */

void calc_ktKikx(double *ktKik, const int m, double **k, const int n,
		 double *g, const double mui, double *kxy, double *ktKikx)
{
  int i;
  double **Gmui;
  double *ktGmui;

  /* first calculate Gmui = g %*% t(g)/mu */
  Gmui = new_zero_matrix(n, n);
  linalg_dgemm(CblasNoTrans,CblasTrans,n,n,1,
               mui,&g,n,&g,n,0.0,Gmui,n);

  /* used in the for loop below */
  ktGmui = new_vector(n);

  /* loop over all of the m candidates */
  for(i=0; i<m; i++) {

    /* ktGmui = drop(t(k) %*% Gmui) */
    /* zerov(ktGmui, n); */
    linalg_dsymv(n,1.0,Gmui,n,k[i],1,0.0,ktGmui,1);

    /* ktKik += diag(t(k) %*% (g %*% t(g) * mui) %*% k) */
    ktKikx[i] = ktKik[i] + linalg_ddot(n, ktGmui, 1, k[i], 1);

    /* ktKik.x += + 2*diag(kxy %*% t(g) %*% k) */
    ktKikx[i] += 2.0*linalg_ddot(n, k[i], 1, g, 1)*kxy[i];

    /* ktKik.x + kxy^2/mui */
    ktKikx[i] += sq(kxy[i])/mui;
  }

  /* clean up */
  free(ktGmui);
  delete_matrix(Gmui);
}


/*
 * calc_ktKikx_R:
 *
 * function for calculating the ktKikx matrix used in the
 * IECI calculation -- R interface
 */
 
void calc_ktKikx_R(double *ktKik_inout, int *m_in, double *k_in, int *n_in,
		   double *g_in, double *mui_in, double *kxy_in)
{
  int m, n;
  double  **k;
  
  /* copy integers */
  m = *m_in;
  n = *n_in;

  /* make matrix bones */
  k = new_matrix_bones(k_in, m, n);

  /* do the work */
  calc_ktKikx(ktKik_inout, m, k, n, g_in, *mui_in, kxy_in, ktKik_inout);

  /* clean up */
  free(k);
}


/*
 * calc2_ktKikx_R:
 *
 * function for calculating the ktKikx matrix used in the
 * IECI calculation -- R interface
 */
 
void calc2_ktKikx_R(double *ktKik_inout, int *m_in, double *k_in, int *n_in,
		    double *x_in, int *col_in, double *X_in, 
		    double *Ki_in, double *Xref_in, double *d_in, 
		    int *dlen_in, double *g_in)
{
  int m, n, col, dlen;
  double **X, **Xref, **k, **Ki;
  double *gvec, *kxy;
  double mui = 0;
  
  /* copy integers */
  m = *m_in;
  n = *n_in;
  col = *col_in;
  dlen = *dlen_in;

  /* make matrix bones */
  X = new_matrix_bones(X_in, n, col);
  Ki = new_matrix_bones(Ki_in, n, n);
  Xref = new_matrix_bones(Xref_in, m, col);
  k = new_matrix_bones(k_in, m, n);  

  /* allocate g and kxy vector */
  gvec = new_vector(n);
  kxy = new_vector(m);

  /* do the work */
  calc_g_mui_kxy(col, x_in, X, n, Ki, Xref, m, d_in,
		 dlen, *g_in, gvec, &mui, kxy);

  /* do the work */
  calc_ktKikx(ktKik_inout, m, k, n, gvec, mui, kxy, ktKik_inout);

  /* clean up */
  free(gvec);
  free(kxy);
  free(X);
  free(Ki);
  free(Xref);
  free(k);
}


/*
 * calc_ecis_R:
 *
 * function for calculating the ECIs used in the
 * IECI calculation -- R interface; stores the output (EICs) in 
 * ktKik_inout, writing over the input there
 */
 
void calc_ecis_R(double *ktKik_inout, int *m_in, double *k_in, int *n_in,
		 double *x_in, int *col_in, double *X_in, double *Ki_in, 
		 double *Xref_in, double *d_in, int *dlen_in, double *g_in, 
		 double *s2p_in, double *phi_in, double *badj_in, 
		 double *tm_in, int *tdf_in, double *fmin_in)
{
  int m, n, col, dlen;
  double **X, **Xref, **k, **Ki;
  double *gvec, *kxy;
  double mui;

  mui = 0;
  
  /* copy integers */
  m = *m_in;
  n = *n_in;
  col = *col_in;
  dlen = *dlen_in;

  /* make matrix bones */
  X = new_matrix_bones(X_in, n, col);
  Ki = new_matrix_bones(Ki_in, n, n);
  Xref = new_matrix_bones(Xref_in, m, col);
  k = new_matrix_bones(k_in, m, n);  

  /* allocate g and kxy vector */
  gvec = new_vector(n);
  kxy = new_vector(m);

  /* calculate the g vector, mui, and kxy */
  calc_g_mui_kxy(col, x_in, X, n, Ki, Xref, m, d_in, 
		 dlen, *g_in, gvec, &mui, kxy);

  /* use g, mu, and kxy to calculate ktKik.x */
  calc_ktKikx(ktKik_inout, m, k, n, gvec, mui, kxy, ktKik_inout);

  calc_ecis(m, ktKik_inout, s2p_in, *phi_in, *g_in, badj_in, 
	    tm_in, *tdf_in, *fmin_in);

  /* clean up */
  free(gvec);
  free(kxy);
  free(X);
  free(Ki);
  free(Xref);
  free(k);
}


/*
 * calc_ieci_R:
 *
 * function for calculating the IECI for a particular input
 * x, with reference to loactions in Xref -- R interface; 
 * stores the output in ieci_out, but ktKik_in is also written
 * over
 */

void calc_ieci_R(double *ktKik_in, int *m_in, double *k_in, int *n_in,
		 double *x_in, int *col_in, double *X_in, double *Ki_in, 
		 double *Xref_in, double *d_in, int *dlen_in, double *g_in, 
		 double *s2p_in, double *phi_in, double *badj_in, double *tm_in, 
		 int *tdf_in, double *fmin_in, double *w_in, double *ieci_out)
{
  int m, n, col, dlen;
  double **X, **Xref, **k, **Ki;
  double *gvec, *kxy;
  double mui;

  mui = 0;
  
  /* copy integers */
  m = *m_in;
  n = *n_in;
  col = *col_in;
  dlen = *dlen_in;

  /* make matrix bones */
  X = new_matrix_bones(X_in, n, col);
  Ki = new_matrix_bones(Ki_in, n, n);
  Xref = new_matrix_bones(Xref_in, m, col);
  k = new_matrix_bones(k_in, m, n);  

  /* allocate g and kxy vector */
  gvec = new_vector(n);
  kxy = new_vector(m);

  /* calculate the g vector, mui, and kxy */
  calc_g_mui_kxy(col, x_in, X, n, Ki, Xref, m, d_in, 
		 dlen, *g_in, gvec, &mui, kxy);

  /* use g, mu, and kxy to calculate ktKik.x */
  calc_ktKikx(ktKik_in, m, k, n, gvec, mui, kxy, ktKik_in);

  /* calculate the IECI */
  *ieci_out = calc_ieci(m, ktKik_in, s2p_in, *phi_in, *g_in, badj_in, 
			tm_in, *tdf_in, *fmin_in, w_in);

  /* clean up */
  free(gvec);
  free(kxy);
  free(X);
  free(Ki);
  free(Xref);
  free(k);
}


/*
 * calc_iecis_R:
 *
 * function for calculating the IECIs at all candidate locations
 * Xcand, with reference to loactions in Xref -- R interface; 
 * stores the output in vector ieci_out
 */


void calc_iecis_R(double *ktKik_in, int *m_in, double *k_in, int *n_in,
		  double *Xcand_in, int *I_in, int *col_in, double *X_in, 
		  double *Ki_in, double *Xref_in, double *d_in, int *dlen_in, 
		  double *g_in, double *s2p_in, double *phi_in, double *badj_in, 
		  double *tm_in, int *tdf_in, double *fmin_in, double *w_in, 
		  int *verb_in, double *ieci_out)
{
  int m, n, col, dlen, I, i;
  double **X, **Xcand, **Xref, **k, **Ki;
  double *gvec, *kxy, *ktKikx;
  double mui;

  mui = 0;
  
  /* copy integers */
  m = *m_in;
  n = *n_in;
  col = *col_in;
  dlen = *dlen_in;
  I = *I_in;

  /* make matrix bones */
  X = new_matrix_bones(X_in, n, col);
  Ki = new_matrix_bones(Ki_in, n, n);
  Xcand = new_matrix_bones(Xcand_in, I, col);
  Xref = new_matrix_bones(Xref_in, m, col);
  k = new_matrix_bones(k_in, m, n);  

  /* allocate g, kxy, and ktKikx vectors */
  gvec = new_vector(n);
  kxy = new_vector(m);
  ktKikx = new_vector(m);

  /* calculate the IECI for each candidate */
  for(i=0; i<I; i++) {

    /* progress meter */
    if(*verb_in > 1)
      myprintf(stdout, "calculating ECI for point %d of %d\n", i, I);
    
    /* calculate the g vector, mui, and kxy */
    calc_g_mui_kxy(col, Xcand[i], X, n, Ki, Xref, m, d_in, 
		   dlen, *g_in, gvec, &mui, kxy);

    /* skip if numerical problems */
    if(mui <= sqrt(DOUBLE_EPS)) {
      ieci_out[i] = 1e300 * 1e300;
      continue;
    }

    /* use g, mu, and kxy to calculate ktKik.x */
    calc_ktKikx(ktKik_in, m, k, n, gvec, mui, kxy, ktKikx);
    
    /* calculate the IECI */
    ieci_out[i] = calc_ieci(m, ktKikx, s2p_in, *phi_in, *g_in, badj_in, 
			    tm_in, *tdf_in, *fmin_in, w_in);
  }

  /* clean up */
  free(ktKikx);
  free(gvec);
  free(kxy);
  free(X);
  free(Xcand);
  free(Ki);
  free(Xref);
  free(k);
}


/*
 * calc_eis_R:
 *
 * R interface function for calculating EI over many sets of
 * predictive (t) distributions in tmat
 */

void calc_eis_R(double *tmat_in, int *n_in, double *fmin_in, double *eis_out)
{
  int n, i;
  double **tmat;
  double fmin;

  /* copy scalars */
  fmin = *fmin_in;
  n = *n_in;

  /* make matrix bones */
  tmat = new_matrix_bones(tmat_in, n, 3);

  for(i=0; i<n; i++) 
    eis_out[i] = EI(tmat[i][0], tmat[i][1], tmat[i][2], fmin);

  /* clean up */
  free(tmat);
}
