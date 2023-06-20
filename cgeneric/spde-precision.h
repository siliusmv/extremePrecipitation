#pragma once

#include "cgeneric.h"

#define Calloc(n, type)  (type *) calloc(n, sizeof(type))
#define Malloc(n, type)  (type *) malloc((n) * sizeof(type))

// Compute the precision matrix of the SPDE approximation
// This is a modified version of the R function INLA::inla.spde2.precision()
inla_cgeneric_smat_tp spde_precision(double log_rho,
				     double log_sigma,
				     inla_cgeneric_mat_tp ** B,
				     inla_cgeneric_smat_tp ** M);

// Compute n_mesh different precision matrixes of the SPDE approximation.
inla_cgeneric_smat_tp * spde_precision_multimesh(double log_rho,
						 double log_sigma,
						 inla_cgeneric_mat_tp ** B,
						 inla_cgeneric_smat_tp ** M,
						 int n_mesh);

// Copy of the R functions INLA::theta2phiX(), where X = 0, 1, 2.
// Necessary for computing the precision matrix of the SPDE approximation.
double * theta2phi(double log_rho,
		   double log_sigma,
		   inla_cgeneric_mat_tp const * B,
		   int exp_transform);

// Take the inla_cgeneric_smat_tp A with dimension (m x k) and turn it into
// a block diagonal matrix of dimension ((n * m) x (n * k)) that contains
// n copies of the original A matrix.
void block_diag_smat(inla_cgeneric_smat_tp * A, int n);

// Remove all the lower diagonal elements of an inla_cgeneric_smat_tp.
// This is necessary because the cgeneric format only wants the upper diagonal
// parts of the precision matrix
void upper_diag(inla_cgeneric_smat_tp * A);

// Sort the elements of an inla_cgeneric_smat_tp, such that i is strictly increasing,
// while j is strictly increasing within each value of i. This is the
// format required by the cgeneric model for the precision matrix and graph.
void sort_smat(inla_cgeneric_smat_tp * A);

// This is a helper function for the sort_smat() function, used as input to the
// c function qsort(), which we use inside the sort_smat() function.
int smat_cmp_func_i_then_j(void const * a, void const * b);

// This is a helper struct used in the sort_smat() function and in the
// smat_cmp_func_i_then_j() function. The smat_element struct represents one
// element of an inla_cgeneric_smat_tp.
typedef struct {
  int i;
  int j;
  double x;
} smat_element;
