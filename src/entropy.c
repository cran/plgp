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
#include <stdlib.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>


/*
 * entropu:
 *
 * calculates the entropy of the probabilities in
 * pvec; assumes 0 <= pvec <= 1
 */

double entropy(double *pvec, const int nc)
{
  int i;
  double ent;
  
  ent = 0.0;
  for(i=0; i<nc; i++) ent -= pvec[i] * log(pvec[i]);

  return(ent);
}


/*
 * calc_ents_R:
 *
 * R interface function for calculating entropies 
 * over many sets of predictive (class) distributions
 * in pmat
 */

void calc_ents_R(double *pmat_in, int *n_in, int *nc_in, double *ents_out)
{
  int n, nc, i;
  double **pmat;

  /* copy scalars */
  n = *n_in;
  nc = *nc_in;

  /* make matrix bones */
  pmat = new_matrix_bones(pmat_in, n, nc);

  for(i=0; i<n; i++) ents_out[i] = entropy(pmat[i], nc);

  /* clean up */
  free(pmat);
}

