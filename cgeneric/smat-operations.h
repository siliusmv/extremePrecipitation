#pragma once

#include "cgeneric.h"

#define Calloc(n, type)  (type *) calloc(n, sizeof(type))
#define Malloc(n, type)  (type *) malloc((n) * sizeof(type))

// Create a sparse matrix that equals the sum A + B of the two sparse
// matrices A and B
inla_cgeneric_smat_tp smat_smat_add(inla_cgeneric_smat_tp const * A,
				    inla_cgeneric_smat_tp const * B);

// Add the sparse matrix B to the sparse matrix A
void smat_smat_add_in_place(inla_cgeneric_smat_tp * A,
			    inla_cgeneric_smat_tp const * B);

// Compute A * x and place the results in b
void mat_vec_mult(double * b, inla_cgeneric_mat_tp const * A, double const * x);

// Compute A * x and place the results in x
void mat_vec_mult_in_place(inla_cgeneric_mat_tp const * A, double * x);

// Compute A * x
void smat_vec_mult(double * out, inla_cgeneric_smat_tp const * A, double const * x);
double * smat_vec_mult_unsafe(inla_cgeneric_smat_tp const * A, double const * x);

// Remove column nr. n and row nr. n from the sparse matrix A.
// Do nothing if n < 0
void remove_col_and_row_from_smat(inla_cgeneric_smat_tp * A, int n);

void remove_cols_and_rows_from_smat(inla_cgeneric_smat_tp * A, int * n, int length);

// This is a struct that describes a diagonal matrix of dimension (dim x dim),
// and with diagonal elements x
typedef struct {
  int dim;
  double * x;
} diag_mat_tp;

// Create an (n x n)-dimensional diagonal matrix, containing x at each diagonal entry
diag_mat_tp diag(double const x, int const n);

// Given a sparse matrix A and a diagonal matrix D,
// Set A = D * A
void diag_smat_mult(inla_cgeneric_smat_tp * A, diag_mat_tp const * D);

// Given a sparse matrix A and a diagonal matrix D,
// Set A = A * D
void smat_diag_mult(inla_cgeneric_smat_tp * A, diag_mat_tp const * D);

// Given a sparse matrix A and a diagonal matrix D,
// Set A = D * A * D
void diag_smat_diag_mult(inla_cgeneric_smat_tp * A, diag_mat_tp const * D);

// Turn the sparse matrix A into the sparse matrix [A 0; 0 D], where D is a diagonal matrix
void smat_diag_merge(inla_cgeneric_smat_tp * A, diag_mat_tp const * D);

// Turn the sparse matrix A into the sparse matrix [A 0; 0 B]
void smat_smat_merge(inla_cgeneric_smat_tp * A, inla_cgeneric_smat_tp const * B);

// Transpose the sparse matrix A in place
void transpose_smat(inla_cgeneric_smat_tp * A);

// Create a deep copy of the sparse matrix A
inla_cgeneric_smat_tp copy_smat(inla_cgeneric_smat_tp const * A);

// Create a deep copy of the matrix A
inla_cgeneric_mat_tp copy_mat(inla_cgeneric_mat_tp const * A);

// Free allocated memory from the sparse matrix A
void free_smat(inla_cgeneric_smat_tp * A);

// Free allocated memory from the matrix A
void free_mat(inla_cgeneric_mat_tp * A);

// Free allocated memory from the diagonal matrix A
void free_diag(diag_mat_tp * A);

// This is a helper structure for keeping track of sparse matrices.
// Within one row or column of a sparse matrix, we can keep track of all the
// nonzero elements. We kan know their locations in the row/column (index),
// and their values (val), and we can know how many nonzero elements there are (n)
typedef struct {
  int * index;
  double * val;
  int n;
} nonzero_index;

// Given two nonzero_index structures x and y, representing two sparse
// vectors with equal length, compute the sum of these two vectors and return
// the result as a nonzero_index structure
nonzero_index sum_of_nonzero_indices(nonzero_index const * x,
				     nonzero_index const * y);

// Free allocated memory from the nonzero index A
void free_nonzero_index(nonzero_index A);

// This is a helper structure for keeping track of sparse matrices.
// Given an array that contains `num_ints` different integer values, where replicates are allowed,
// we can keep track of the number of the locations of integer nr. j (locs[j]), and
// the number of times integer nr. j is repeated (n[j])
typedef struct {
  int ** locs;
  int * n;
  int num_ints;
} int_loc;

// Create an int_loc struct for the sparse matrix A, that represents how many
// rows of A contain nonzero elements, how many nonzero elements are in each row,
// and where in A we can find all the elements for each row.
int_loc locate_nonzero_smat_rows(inla_cgeneric_smat_tp const * A);


// The array arr is an array of integers between 0 and num_ints - 1,
// of length arr_length.
// This function creates an int_loc struct for arr, that describes how many
// times each of the numbers from 0 to num_ints - 1 are represented in arr,
// and where in arr we can find each of the different ints.
int_loc locate_ints_in_arr(int const * arr,
			   int num_ints,
			   int arr_length);

// Compute the minimum of two integers
int int_min(int x, int y);
