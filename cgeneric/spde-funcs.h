#pragma once

#include <assert.h>
#if !defined(__FreeBSD__)
#endif
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <stdio.h>

#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_spblas.h>

#include "cgeneric.h"

#define Calloc(n, type)  (type *) calloc(n, sizeof(type))
#define Malloc(n, type)  (type *) malloc((n) * sizeof(type))

// This is a struct that describes a diagonal matrix of dimension (dim x dim),
// and with diagonal elements x
typedef struct {
  int dim;
  double * x;
} diag_mat_tp;

// This is a helper structure for keeping track of sparse matrices.
// Given an array that contains `num_ints` different integer values, where replicates are allowed,
// we can keep track of the number of the locations of integer nr. j (locs[j]), and
// the number of times integer nr. j is repeated (n[j])
typedef struct {
  int ** locs;
  int * n;
  int num_ints;
} int_loc;

// This is a helper structure for keeping track of sparse matrices.
// Within one row or column of a sparse matrix, we can keep track of all the
// nonzero elements. We kan know their locations in the row/column (index),
// and their values (val), and we can know how many nonzero elements there are (n)
typedef struct {
  int * index;
  double * val;
  int n;
} nonzero_index;


// Compute the max of two integers
int int_max(const int x, const int y) {
  return y > x ? y : x;
}

double double_max(const double x, const double y) {
  return y > x ? y : x;
}

// Compute the min of two integers
int int_min(const int x, const int y) {
  return y < x ? y : x;
}

int_loc locate_ints_in_arr(const int *arr, const int num_ints, const int arr_length) {
  int_loc res;
  res.locs = (int **) calloc(num_ints, sizeof(int *));
  res.n = Calloc(num_ints, int);
  res.num_ints = num_ints;
  int *max_size = Calloc(num_ints, int);
  for (int i = 0; i < num_ints; i++) {
    res.locs[i] = Calloc(5, int);
    max_size[i] = 5;
    res.n[i] = 0;
  }
  for (int i = 0; i < arr_length; i++) {
    if ((res.n[arr[i]] + 1) == max_size[arr[i]]) {
      max_size[arr[i]] *= 2;
      res.locs[arr[i]] = realloc(res.locs[arr[i]], max_size[arr[i]] * sizeof(int));
    }
    res.locs[arr[i]][res.n[arr[i]]] = i;
    res.n[arr[i]] += 1;
  }
  free(max_size);
  return res;
}

int_loc locate_nonzero_smat_rows(const inla_cgeneric_smat_tp *A) {
  return locate_ints_in_arr(A->i, A->nrow, A->n);
}

inla_cgeneric_smat_tp copy_smat(const inla_cgeneric_smat_tp *A) {
  inla_cgeneric_smat_tp res;
  res.name = NULL;
  res.nrow = A->nrow;
  res.ncol = A->ncol;
  res.n = A->n;
  res.x = Calloc(A->n, double);
  res.i = Calloc(A->n, int);
  res.j = Calloc(A->n, int);
  for (int k = 0; k < A->n; k++) {
    res.x[k] = A->x[k];
    res.i[k] = A->i[k];
    res.j[k] = A->j[k];
  }
  return res;
}

inla_cgeneric_mat_tp copy_mat(const inla_cgeneric_mat_tp *A) {
  inla_cgeneric_mat_tp res;
  res.name = NULL;
  res.nrow = A->nrow;
  res.ncol = A->ncol;
  res.x = Calloc(A->nrow * A->ncol, double);
  for (int i = 0; i < A->nrow * A->ncol; i++) {
    res.x[i] = A->x[i];
  }
  return res;
}


void free_smat(inla_cgeneric_smat_tp *A) {
  free(A->x);
  free(A->i);
  free(A->j);
  free(A->name);
}

void free_mat(inla_cgeneric_mat_tp *A) {
  free(A->x);
  free(A->name);
}

void diag_smat_mult(inla_cgeneric_smat_tp *A, const diag_mat_tp *D) {
  assert(D->dim == A->nrow);
  for (int i = 0; i < A->n; i++) {
    A->x[i] *= D->x[A->i[i]];
  }
}

void smat_diag_mult(inla_cgeneric_smat_tp *A, const diag_mat_tp *D) {
  assert(D->dim == A->ncol);
  for (int i = 0; i < A->n; i++) {
    A->x[i] *= D->x[A->j[i]];
  }
}

void transpose_smat(inla_cgeneric_smat_tp *A) {
  int *tmp_pnt = A->i;
  A->i = A->j;
  A->j = tmp_pnt;
  int tmp = A->nrow;
  A->nrow = A->ncol;
  A->ncol = tmp;
}

nonzero_index sum_of_nonzero_indices(const nonzero_index *x, const nonzero_index *y) {
  // We assume that there are no duplicate elements inside each of the index-structs
  nonzero_index res;
  res.index = Calloc(x->n + y->n, int);
  res.val = Calloc(x->n + y->n, double);
  int any_equal;
  int n = 0;
  for (int i = 0; i < x->n; i++) {
    res.index[n] = x->index[i];
    res.val[n] = x->val[i];
    for (int j = 0; j < y->n; j++) {
      if (x->index[i] == y->index[j]) {
	res.val[n] += y->val[j];
	break;
      }
    }
    n += 1;
  }
  for (int i = 0; i < y->n; i++) {
    any_equal = 0;
    for (int j = 0; j < x->n; j++) {
      if (y->index[i] == x->index[j]) {
	any_equal = 1;
	break;
      }
    }
    if (any_equal == 0) {
      res.index[n] = y->index[i];
      res.val[n] = y->val[i];
      n += 1;
    }
  }
  res.n = n;
  return res;
}

inla_cgeneric_smat_tp smat_smat_add(const inla_cgeneric_smat_tp *A, const inla_cgeneric_smat_tp *B) {
  assert(A->ncol == B->ncol);
  assert(A->nrow == B->nrow);

  inla_cgeneric_smat_tp res;
  res.name = NULL;
  res.ncol = A->ncol;
  res.nrow = A->nrow;

  int max_alloc_size = int_max(A->n, B->n);
  res.x = Calloc(max_alloc_size, double);
  res.i = Calloc(max_alloc_size, int);
  res.j = Calloc(max_alloc_size, int);
  int n = 0;
  
  // Where can I find the location of row 0, 1, 2, ...?
  int_loc A_row_locs = locate_smat_rows(A);
  int_loc B_row_locs = locate_smat_rows(B);

  nonzero_index A_nonzero_cols;
  nonzero_index B_nonzero_cols;
  nonzero_index A_B_sum;
  A_nonzero_cols.index = Calloc(A->ncol, int);
  A_nonzero_cols.val = Calloc(A->ncol, double);
  B_nonzero_cols.index = Calloc(B->ncol, int);
  B_nonzero_cols.val = Calloc(B->ncol, double);

  for (int i = 0; i < A->nrow; i++) {
    if (A_row_locs.n[i] == 0 && B_row_locs.n[i] > 0) {
      for (int j = 0; j < B_row_locs.n[i]; j++) {
	if (n == max_alloc_size) {
	  max_alloc_size *= 2;
	  res.x = realloc(res.x, max_alloc_size * sizeof(double));
	  res.i = realloc(res.i, max_alloc_size * sizeof(int));
	  res.j = realloc(res.j, max_alloc_size * sizeof(int));
	}
	res.x[n] = B->x[B_row_locs.locs[i][j]];
	res.j[n] = B->j[B_row_locs.locs[i][j]];
	res.i[n] = i;
	n += 1;
      }
    } else if (A_row_locs.n[i] > 0 && B_row_locs.n[i] == 0) {
      for (int j = 0; j < A_row_locs.n[i]; j++) {
	if (n == max_alloc_size) {
	  max_alloc_size *= 2;
	  res.x = realloc(res.x, max_alloc_size * sizeof(double));
	  res.i = realloc(res.i, max_alloc_size * sizeof(int));
	  res.j = realloc(res.j, max_alloc_size * sizeof(int));
	}
	res.x[n] = A->x[A_row_locs.locs[i][j]];
	res.j[n] = A->j[A_row_locs.locs[i][j]];
	res.i[n] = i;
	n += 1;
      }
    } else if (A_row_locs.n[i] + B_row_locs.n[i] > 0) {
      for (int j = 0; j < A_row_locs.n[i]; j++) {
	A_nonzero_cols.index[j] = A->j[A_row_locs.locs[i][j]];
	A_nonzero_cols.val[j] = A->x[A_row_locs.locs[i][j]];
      }
      A_nonzero_cols.n = A_row_locs.n[i];
      for (int j = 0; j < B_row_locs.n[i]; j++) {
	B_nonzero_cols.index[j] = B->j[B_row_locs.locs[i][j]];
	B_nonzero_cols.val[j] = B->x[B_row_locs.locs[i][j]];
      }
      B_nonzero_cols.n = B_row_locs.n[i];
      A_B_sum = sum_of_nonzero_indices(&A_nonzero_cols, &B_nonzero_cols);
      for (int j = 0; j < A_B_sum.n; j++) {
	if (n == max_alloc_size) {
	  max_alloc_size *= 2;
	  res.x = realloc(res.x, max_alloc_size * sizeof(double));
	  res.i = realloc(res.i, max_alloc_size * sizeof(int));
	  res.j = realloc(res.j, max_alloc_size * sizeof(int));
	}
	res.x[n] = A_B_sum.val[j];
	res.j[n] = A_B_sum.index[j];
	res.i[n] = i;
	n += 1;
      }
    }
  }
  res.n = n;
  return res;
}

void smat_smat_add_in_place(inla_cgeneric_smat_tp *A, const inla_cgeneric_smat_tp *B) {
  inla_cgeneric_smat_tp res = smat_smat_add(A, B);
  free_smat(A);
  A->name = res.name;
  A->i = res.i;
  A->j = res.j;
  A->x = res.x;
  A->n = res.n;
}

double *theta2phi2(const double log_rho,
		   const double log_sigma,
		   const inla_cgeneric_mat_tp * B,
		   const int exp_transform) {
  assert(B->ncol == 3);
  double * res = Calloc(B->nrow, double);
  for (int i = 0; i < B->nrow; i++) {
    res[i] = B->x[i * B->ncol] +
      B->x[i * B->ncol + 1] * log_rho +
      B->x[i * B->ncol + 2] * log_sigma;
  }
  if (exp_transform) {
    for (int i = 0; i < B->nrow; i++) {
      res[i] = exp(res[i]);
    }
  }
  return res;
}

inla_cgeneric_smat_tp spde_precision2(double log_rho,
				      double log_sigma,
				      inla_cgeneric_mat_tp **B,
				      inla_cgeneric_smat_tp **M) {
  double *phi0 = theta2phi2(log_rho, log_sigma, B[0], 1);
  double *phi1 = theta2phi2(log_rho, log_sigma, B[1], 1);
  double *phi2 = theta2phi2(log_rho, log_sigma, B[2], 0);
  diag_mat_tp D0 = {B[0]->nrow, phi0};
  diag_mat_tp D1 = {B[1]->nrow, phi1};
  diag_mat_tp D12 = {B[2]->nrow, phi2};
  for (int i = 0; i < D12.dim; i++) {
    D12.x[i] *= phi1[i];
  }
  inla_cgeneric_smat_tp res = copy_smat(M[0]);
  inla_cgeneric_smat_tp tmp = copy_smat(M[1]);
  smat_diag_mult(&res, &D1);
  diag_smat_mult(&res, &D1);
  diag_smat_mult(&tmp, &D12);
  smat_smat_add_in_place(&res, &tmp);
  transpose_smat(&tmp);
  smat_smat_add_in_place(&res, &tmp);
  smat_smat_add_in_place(&res, M[2]);
  diag_smat_mult(&res, &D0);
  smat_diag_mult(&res, &D0);
  free(phi0);
  free(phi1);
  free(phi2);
  free_smat(&tmp);
  return res;
}

void upper_diag(inla_cgeneric_smat_tp *A) {
  int count = 0;
  for (int i = 0; i < A->n; i++) {
    if (A->i[i] <= A->j[i]) {
      if (count < i) {
	A->i[count] = A->i[i];
	A->j[count] = A->j[i];
	A->x[count] = A->x[i];
      }
      count += 1;
    }
  }
  A->n = count;
}

typedef struct {
  int i;
  int j;
  double x;
} smat_element;

int smat_cmp_func_i_then_j(const void *a, const void *b) {
  if (((smat_element*)a)->i > ((smat_element*)b)->i) {
    return 1;
  } else if (((smat_element*)a)->i < ((smat_element*)b)->i) {
    return -1;
  } else if (((smat_element*)a)->j > ((smat_element*)b)->j) {
    return 1;
  } else if (((smat_element*)a)->j < ((smat_element*)b)->j) {
    return -1;
  } else {
    return 0;
  }
}

void sort_smat(inla_cgeneric_smat_tp *A) {
  smat_element *elements = Calloc(A->n, smat_element);
  for (int i = 0; i < A->n; i++) {
    elements[i].i = A->i[i];
    elements[i].j = A->j[i];
    elements[i].x = A->x[i];
  }
  qsort(elements, A->n, sizeof(smat_element), smat_cmp_func_i_then_j);
  for (int i = 0; i < A->n; i++) {
    A->i[i] = elements[i].i;
    A->j[i] = elements[i].j;
    A->x[i] = elements[i].x;
  }
  free(elements);
}

void block_diag_smat(inla_cgeneric_smat_tp *A, const int n) {
  A->i = (int *) realloc(A->i, A->n * n * sizeof(int));
  A->j = (int *) realloc(A->j, A->n * n * sizeof(int));
  A->x = (double *) realloc(A->x, A->n * n * sizeof(double));
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < A->n; j++) {
      A->i[j + A->n * i] = A->i[j] + A->nrow * i;
      A->j[j + A->n * i] = A->j[j] + A->ncol * i;
      A->x[j + A->n * i] = A->x[j];
    }
  }
  A->nrow *= n;
  A->ncol *= n;
  A->n *= n;
}

gsl_spmatrix *inla_smat_to_gsl_spmatrix(const inla_cgeneric_smat_tp *A) {
  gsl_spmatrix * tmp = gsl_spmatrix_alloc(A->nrow, A->ncol);
  for (int i = 0; i < A->n; i++) {
    gsl_spmatrix_set(tmp, A->i[i], A->j[i], A->x[i]);
  }
  gsl_spmatrix * res = gsl_spmatrix_ccs(tmp);
  gsl_spmatrix_free(tmp);
  return res;
}

gsl_spmatrix * csc_to_coo(gsl_spmatrix const * A) {
  gsl_spmatrix *res = gsl_spmatrix_alloc(A->size1, A->size2);
  int counter = 0, nums_in_col;
  for (int i = 0; i < A->size2; i++) {
    if (i < (A->size2 - 1)) {
      nums_in_col = A->p[i + 1] - A->p[i];
    } else {
      nums_in_col = A->nz - A->p[i];
    }
    for (int j = 0; j < nums_in_col; j++) {
      gsl_spmatrix_set(res, A->i[counter], i, A->data[counter]);
      counter += 1;
    }
  }
  return res;
}

inla_cgeneric_smat_tp gsl_spmatrix_to_inla_smat(gsl_spmatrix * A) {
  gsl_spmatrix * tmp;
  if (A->sptype == GSL_SPMATRIX_CSC) {
    tmp = csc_to_coo(A);
  } else {
    tmp = A;
  }
  assert(tmp->sptype == GSL_SPMATRIX_COO);
  inla_cgeneric_smat_tp res;
  res.name = NULL;
  res.nrow = tmp->size1;
  res.ncol = tmp->size2;
  res.n = tmp->nz;
  res.x = Calloc(res.n, double);
  res.i = Calloc(res.n, int);
  res.j = Calloc(res.n, int);
  for (int k = 0; k < tmp->nz; k++) {
    res.x[k] = tmp->data[k];
    res.i[k] = tmp->i[k];
    res.j[k] = tmp->p[k];
  }
  if (A->sptype != GSL_SPMATRIX_COO) {
    gsl_spmatrix_free(tmp);
  }
  return res;
}

gsl_vector *theta2phi(double const log_rho,
		       double const log_sigma,
		       inla_cgeneric_mat_tp const * B,
		       int const exp_transform) {
  assert(B->ncol == 3);
  gsl_vector * res = gsl_vector_alloc(B->nrow);
  double tmp;
  for (int i = 0; i < B->nrow; i++) {
    tmp = B->x[i * B->ncol] +
      B->x[i * B->ncol + 1] * log_rho +
      B->x[i * B->ncol + 2] * log_sigma;
    if (exp_transform) {
      tmp = exp(tmp);
    }
    gsl_vector_set(res, i, tmp);
  }
  return res;
}

void gsl_spmatrix_add_in_place(gsl_spmatrix * left, gsl_spmatrix const * right) {
  assert(left->size1 == right->size1);
  assert(left->size2 == right->size2);
  assert(left->sptype == right->sptype);
  gsl_spmatrix * tmp = gsl_spmatrix_alloc(left->size1, left->size2);
  tmp->sptype = left->sptype;
  gsl_spmatrix_memcpy(tmp, left);
  gsl_spmatrix_add(left, tmp, right);
  gsl_spmatrix_free(tmp);
}

inla_cgeneric_smat_tp spde_precision(double log_rho,
				     double log_sigma,
				     inla_cgeneric_mat_tp** B,
				     inla_cgeneric_smat_tp** M) {
  gsl_vector * phi0 = theta2phi(log_rho, log_sigma, B[0], 1);
  gsl_vector * phi1 = theta2phi(log_rho, log_sigma, B[1], 1);
  gsl_vector * phi12 = theta2phi(log_rho, log_sigma, B[2], 0);
  gsl_vector_mul(phi12, phi1);

  gsl_spmatrix * tmp1 = inla_smat_to_gsl_spmatrix(M[0]);
  gsl_spmatrix_scale_columns(tmp1, phi1);
  gsl_spmatrix_scale_rows(tmp1, phi1);

  gsl_spmatrix * tmp2 = inla_smat_to_gsl_spmatrix(M[1]);
  gsl_spmatrix_scale_columns(tmp2, phi12);

  gsl_spmatrix * Q = gsl_spmatrix_alloc(M[0]->nrow, M[0]->ncol);
  Q->sptype = GSL_SPMATRIX_CSC;
  gsl_spmatrix_add(Q, tmp1, tmp2);
  gsl_spmatrix_free(tmp1);

  tmp1 = gsl_spmatrix_alloc(tmp2->size1, tmp2->size2);
  tmp1->sptype = tmp2->sptype;
  gsl_spmatrix_transpose_memcpy(tmp1, tmp2);
  gsl_spmatrix_free(tmp2);
  gsl_spmatrix_add_in_place(Q, tmp1);
  gsl_spmatrix_free(tmp1);

  tmp1 = inla_smat_to_gsl_spmatrix(M[2]);
  gsl_spmatrix_add_in_place(Q, tmp1);
  gsl_spmatrix_free(tmp1);

  gsl_spmatrix_scale_columns(Q, phi0);
  gsl_spmatrix_scale_rows(Q, phi0);

  gsl_vector_free(phi0);
  gsl_vector_free(phi1);
  gsl_vector_free(phi12);

  inla_cgeneric_smat_tp res = gsl_spmatrix_to_inla_smat(Q);
  gsl_spmatrix_free(Q);

  return res;
}

