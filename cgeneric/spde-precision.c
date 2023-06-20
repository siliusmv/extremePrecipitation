#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "cgeneric.h"
#include "spde-precision.h"
#include "smat-operations.h"

inla_cgeneric_smat_tp spde_precision(double log_rho,
				     double log_sigma,
				     inla_cgeneric_mat_tp ** B,
				     inla_cgeneric_smat_tp ** M) {

  // Compute phi0, phi1 and phi2. Then create diagonal matrices D0, D1 and D2 that
  // contain elements of these vectors. This is just as seen in the
  // R function INLA::inla.spde2.precision()
  double * phi0 = theta2phi(log_rho, log_sigma, B[0], 1);
  double * phi1 = theta2phi(log_rho, log_sigma, B[1], 1);
  double * phi2 = theta2phi(log_rho, log_sigma, B[2], 0);
  diag_mat_tp D0 = {B[0]->nrow, phi0};
  diag_mat_tp D1 = {B[1]->nrow, phi1};
  diag_mat_tp D12 = {B[2]->nrow, phi2};
  for (int i = 0; i < D12.dim; i++) {
    D12.x[i] *= phi1[i];
  }

  // Create the matrix Q, which we will build up to be the precision matrix of the
  // SPDE, together with a temporary matrix needed for some of the computations
  inla_cgeneric_smat_tp Q = copy_smat(M[0]);
  inla_cgeneric_smat_tp tmp = copy_smat(M[1]);

  // Compute Q = D1 * M0 * D1 in place
  diag_smat_diag_mult(&Q, &D1);

  // Compute tmp = D12 * M1 in place
  diag_smat_mult(&tmp, &D12);

  // Set Q = Q + tmp
  smat_smat_add_in_place(&Q, &tmp);

  // Transpose tmp from D12 * M1 to M1' * D12, and add it to Q
  transpose_smat(&tmp);
  smat_smat_add_in_place(&Q, &tmp);

  // Set Q = Q + M2
  smat_smat_add_in_place(&Q, M[2]);

  // Set Q = D0 * Q * D0
  diag_smat_diag_mult(&Q, &D0);

  // Only keep the upper diagonal of the precision matrix
  upper_diag(&Q);

  // Sort the entries of the smat so they correspond to the R-INLA requirements
  sort_smat(&Q);

  // Free up the allocated memory to avoid memory leaks
  free(phi0);
  free(phi1);
  free(phi2);
  free_smat(&tmp);

  return Q;
}

inla_cgeneric_smat_tp * spde_precision_multimesh(double log_rho,
						 double log_sigma,
						 inla_cgeneric_mat_tp ** B,
						 inla_cgeneric_smat_tp ** M,
						 int n_mesh) {
  inla_cgeneric_smat_tp * res = Malloc(n_mesh, inla_cgeneric_smat_tp);
  for (int i = 0; i < n_mesh; ++i) {
    res[i] = spde_precision(log_rho, log_sigma, B + i * 3, M + i * 3);
  }
  return res;
}

double * theta2phi(double log_rho,
		   double log_sigma,
		   inla_cgeneric_mat_tp const * B,
		   int exp_transform) {
  assert(B->ncol == 3);
  double * res = Malloc(B->nrow, double);
  for (int i = 0; i < B->nrow; ++i) {
    res[i] = B->x[i * B->ncol] +
      B->x[i * B->ncol + 1] * log_rho +
      B->x[i * B->ncol + 2] * log_sigma;
  }
  if (exp_transform) {
    for (int i = 0; i < B->nrow; ++i) {
      res[i] = exp(res[i]);
    }
  }
  return res;
}

void block_diag_smat(inla_cgeneric_smat_tp * A, int n) {
  // Reallocate all the elements of A to make room for more data
  A->i = (int *) realloc(A->i, A->n * n * sizeof(int));
  A->j = (int *) realloc(A->j, A->n * n * sizeof(int));
  A->x = (double *) realloc(A->x, A->n * n * sizeof(double));

  // Copy all the elements of x n times, but increase the values
  // of i and j for each time we add a new block matrix
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < A->n; ++j) {
      A->i[j + A->n * i] = A->i[j] + A->nrow * i;
      A->j[j + A->n * i] = A->j[j] + A->ncol * i;
      A->x[j + A->n * i] = A->x[j];
    }
  }

  // Update the number of rows, columns and elements in A
  A->nrow *= n;
  A->ncol *= n;
  A->n *= n;
}

void upper_diag(inla_cgeneric_smat_tp * A) {
  int count = 0;
  for (int i = 0; i < A->n; ++i) {
    // If we are in the upper diagonal, keep the element. If not, ovewrite the element
    if (A->i[i] <= A->j[i]) {
      if (count < i) {
	A->i[count] = A->i[i];
	A->j[count] = A->j[i];
	A->x[count] = A->x[i];
      }
      count += 1;
    }
  }

  // Update the number of elements available in A
  A->n = count;

  // Reallocate all the elements of A to free unused memory
  A->i = (int *) realloc(A->i, A->n * sizeof(int));
  A->j = (int *) realloc(A->j, A->n * sizeof(int));
  A->x = (double *) realloc(A->x, A->n * sizeof(double));
}

void sort_smat(inla_cgeneric_smat_tp * A) {
  // Create one smat_element for each element in A
  smat_element * elements = Malloc(A->n, smat_element);
  for (int i = 0; i < A->n; ++i) {
    elements[i].i = A->i[i];
    elements[i].j = A->j[i];
    elements[i].x = A->x[i];
  }
  // Sort all the smat_elements
  qsort(elements, A->n, sizeof(smat_element), smat_cmp_func_i_then_j);
  // Overwrite A with the values in the sorted smat_elements
  for (int i = 0; i < A->n; ++i) {
    A->i[i] = elements[i].i;
    A->j[i] = elements[i].j;
    A->x[i] = elements[i].x;
  }
  // Free all the smat_elements
  free(elements);
}

int smat_cmp_func_i_then_j(void const * a, void const * b) {
  if (((smat_element*)a)->i > ((smat_element*)b)->i) {
    return 1; // a > b if a.i > b.i
  } else if (((smat_element*)a)->i < ((smat_element*)b)->i) {
    return -1; // a < b if a.i < b.i
  } else if (((smat_element*)a)->j > ((smat_element*)b)->j) {
    return 1; // a > b if a.i == b.i and a.j > b.j
  } else if (((smat_element*)a)->j < ((smat_element*)b)->j) {
    return -1; // a < b if a.i == b.i and a.j < b.j
  } else {
    return 0; // a == b if a.i == b.i and a.j == b.j
  }
}
