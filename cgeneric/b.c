#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cgeneric.h"
#include "smat-operations.c"
#include "spde-precision.c"

double * spde_b_model(inla_cgeneric_cmd_tp cmd,
		      double * theta,
		      inla_cgeneric_data_tp * data) {

  double * ret = NULL;

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================

  // Size of the latent field
  assert(!strcmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // Which s0 do the different replications come from?
  assert(!strcmp(data->ints[2]->name, "s0_index"));
  int n_repl = data->ints[2]->len;
  int * s0_index = data->ints[2]->ints;

  // Which of the parameters should be estimated, and which are fixed?
  assert(!strcmp(data->ints[3]->name, "is_fixed"));
  int * is_fixed = data->ints[3]->ints;
  int n_theta = data->ints[3]->len;
  int n_free = 0;
  for (int i = 0; i < n_theta; ++i) {
    if (!is_fixed[i]) ++n_free;
  }

  // Look at the matrices required for building the meshes
  int n_mesh = data->n_smats / 3;
  assert(data->n_smats == n_mesh * 3); // Ensure that we did not get rounding when dividing
  assert(data->n_smats == data->n_mats); // We need the same number of mats and smats

  // Where are each of the s0's located in their respective meshes?
  assert(!strcmp(data->ints[4]->name, "s0_location"));
  assert(data->ints[4]->len == n_mesh);
  int * s0_locs = data->ints[4]->ints;

  // Initial values
  assert(!strcmp(data->doubles[0]->name, "init"));
  double * init = data->doubles[0]->doubles;
  assert(data->doubles[0]->len == n_theta);

  // y0
  assert(!strcmp(data->doubles[1]->name, "y0"));
  assert(data->doubles[1]->len == n_repl);
  double * y0 = data->doubles[1]->doubles;

  int prior_start = 2; // This is the first element of data->doubles that contains a prior
  int dist_start = prior_start + n_free; // This is the first element of data->doubles that contains dist_to_s0
  assert(data->n_doubles == n_mesh + dist_start);

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, beta0, lambda, kappa;
  assert(n_theta == 5);
  if (theta) {
    double theta_full[5];
    int count = 0;
    for (int i = 0; i < n_theta; ++i) {
      if (is_fixed[i]) {
	theta_full[i] = init[i];
      } else {
	theta_full[i] = theta[count];
	++count;
      }
    }
    assert(count == n_free);
    log_rho = theta_full[0];
    log_sigma = theta_full[1];
    beta0 = exp(theta_full[2]);
    lambda = exp(theta_full[3]);
    kappa = exp(theta_full[4]);
  } else {
    log_rho = log_sigma = beta0 = lambda = kappa = NAN;
  }

  // ==========================================================
  // This switch statement is the required method for implementing
  // cgeneric models in R-INLA
  // ==========================================================
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }

  case INLA_CGENERIC_GRAPH:
    {
      // Compute all the different precision matrices
      inla_cgeneric_smat_tp * all_precisions =
	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);

      for (int i = 0; i < n_mesh; ++i) {
	remove_col_and_row_from_smat(all_precisions + i, s0_locs[i]);
      }

      // Merge them together
      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + s0_index[0]);
      for (int i = 1; i < n_repl; ++i) {
	smat_smat_merge(&precision, all_precisions + s0_index[i]);
      }

      // Extract all the necessary i and j values from the block matrix
      ret = Malloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  // i
	ret[2 + precision.n + i] = precision.j[i];                    // j
      }

      // Avoid memory leaks
      for (int i = 0; i < n_mesh; ++i) {
	free_smat(all_precisions + i);
      }
      free(all_precisions);
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {

      // Compute all the different precision matrices
      inla_cgeneric_smat_tp * all_precisions =
	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);

      for (int i = 0; i < n_mesh; ++i) {
	remove_col_and_row_from_smat(all_precisions + i, s0_locs[i]);
      }

      // Merge them together
      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + s0_index[0]);
      for (int i = 1; i < n_repl; ++i) {
	smat_smat_merge(&precision, all_precisions + s0_index[i]);
      }

      double ** beta_vals = Malloc(n_mesh, double *);
      for (int i = 0; i < n_mesh; ++i) {
	beta_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	for (int j = 0; j < data->doubles[dist_start + i]->len; ++j) {
	  beta_vals[i][j] = beta0 * exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda, kappa));
	}
      }
      // Create a diagonal matrix containing the values of 1 / b, that can be used
      // for rescaling the precision matrix with the correct b values
      diag_mat_tp b_inv;
      b_inv.x = Malloc(precision.nrow, double);
      b_inv.dim = precision.nrow;
      int count = 0;
      double log_y0;
      for (int i = 0; i < n_repl; ++i) {
	log_y0 = log(y0[i]);
	for (int j = 0; j < data->doubles[dist_start + s0_index[i]]->len; ++j) {
	  b_inv.x[count] = exp(-log_y0 * beta_vals[s0_index[i]][j]);
	  ++count;
	}
      }
      assert(count == precision.nrow);

      // Multiply the diagonal matrix with the precision matrix: B_inv * Q * B_inv
      diag_smat_diag_mult(&precision, &b_inv);

      // optimized format
      // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
      // where M is the length of Qij
      ret = Malloc(2 + precision.n, double);
      ret[0] = -1;                                   // code for optimized output
      ret[1] = precision.n;                          // number of (i <= j)
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.x[i];
      }

      // Avoid memory leaks
      for (int i = 0; i < n_mesh; ++i) {
	free(beta_vals[i]);
	free_smat(all_precisions + i);
      }
      free(all_precisions);
      free(beta_vals);
      free(b_inv.x);
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0

      // Z_b has zero mean
      ret = Calloc(1, double);
      break;
    }

  case INLA_CGENERIC_INITIAL:
    {
      // return c(M, initials)
      // where M is the number of hyperparameters

      // The initial values depend on how many parameters that are fixed or estimated.
      ret = Malloc(n_free + 1, double);
      ret[0] = n_free;
      int count = 0;
      for (int i = 0; i < n_theta; ++i) {
	if (!is_fixed[i]) {
	  ret[count + 1] = init[i];
	  ++count;
	}
      }
      assert(count == n_free);
      break;
    }

  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      // return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself
      ret = NULL;
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)

      // Only add terms for the parameters that are not fixed
      ret = Calloc(1, double);

      // This is an approximation for -0.5 * log(2 * pi)
      double gaussian_const = -0.91893853320467;

      int count = 0;
      double * prior;
      if (!is_fixed[0]) { // rho
	prior = data->doubles[prior_start + count]->doubles;
	double lambda0 = -log(prior[1]) * prior[0];
	ret[0] += log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
	++count;
      }
      if (!is_fixed[1]) { // sigma
	prior = data->doubles[prior_start + count]->doubles;
	double lambda1 = -log(prior[1]) / prior[0];
	ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
	++count;
      }
      if (!is_fixed[2]) { // beta0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(beta0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // lambda
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // kappa
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      assert(count == n_free);
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}

double * spde_model(inla_cgeneric_cmd_tp cmd,
		    double * theta,
		    inla_cgeneric_data_tp * data) {

  double * ret = NULL;

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================

  // Size of the latent field
  assert(!strcmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // Which mesh do the different realisations come from?
  assert(!strcmp(data->ints[2]->name, "mesh_index"));
  int n_realisations = data->ints[2]->len;
  int * mesh_index = data->ints[2]->ints;

  // Which of the parameters should be estimated, and which are fixed?
  assert(!strcmp(data->ints[3]->name, "is_fixed"));
  int * is_fixed = data->ints[3]->ints;
  int n_theta = data->ints[3]->len;
  int n_free = 0;
  for (int i = 0; i < n_theta; ++i) {
    if (!is_fixed[i]) ++n_free;
  }

  // Look at the matrices required for building the meshes
  int n_mesh = data->n_smats / 3;
  assert(data->n_smats == n_mesh * 3); // Ensure that we did not get rounding when dividing
  assert(data->n_smats == data->n_mats); // We need the same number of mats and smats

  // What is the first element of data->ints that contains constraining indices?
  int constr_start = 4;
  // There should either be 0 or n_mesh constraining index vectors
  assert(data->n_ints == constr_start ||
	 data->n_ints == constr_start + n_mesh);
  int any_constraints;
  if (data->n_ints == constr_start) {
    any_constraints = 0;
  } else {
    any_constraints = 1;
  }

  // Initial values
  assert(!strcmp(data->doubles[0]->name, "init"));
  double * init = data->doubles[0]->doubles;
  assert(data->doubles[0]->len == n_theta);

  int prior_start = 1; // This is the first element of data->doubles that contains a prior

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma;
  assert(n_theta == 2);
  if (theta) {
    double theta_full[2];
    int count = 0;
    for (int i = 0; i < n_theta; ++i) {
      if (is_fixed[i]) {
	theta_full[i] = init[i];
      } else {
	theta_full[i] = theta[count];
	++count;
      }
    }
    assert(count == n_free);
    log_rho = theta_full[0];
    log_sigma = theta_full[1];
  } else {
    log_rho = log_sigma = NAN;
  }

  // ==========================================================
  // This switch statement is the required method for implementing
  // cgeneric models in R-INLA
  // ==========================================================
  switch (cmd) {
  case INLA_CGENERIC_VOID:
    {
      assert(!(cmd == INLA_CGENERIC_VOID));
      break;
    }

  case INLA_CGENERIC_GRAPH:
    {
      // Compute all the different precision matrices
      inla_cgeneric_smat_tp * all_precisions =
	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);

      if (any_constraints) {
	for (int i = 0; i < n_mesh; ++i) {
	  if (data->ints[constr_start + i]->len == 1 && data->ints[constr_start + i]->ints[0] < 0) continue;
	  remove_cols_and_rows_from_smat(all_precisions + i,
					 data->ints[constr_start + i]->ints,
					 data->ints[constr_start + i]->len);
	}
      }

      // Merge them together
      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + mesh_index[0]);
      for (int i = 1; i < n_realisations; ++i) {
	smat_smat_merge(&precision, all_precisions + mesh_index[i]);
      }

      // Extract all the necessary i and j values from the block matrix
      ret = Malloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  // i
	ret[2 + precision.n + i] = precision.j[i];                    // j
      }

      // Avoid memory leaks
      for (int i = 0; i < n_mesh; ++i) {
	free_smat(all_precisions + i);
      }
      free(all_precisions);
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {

      // Compute all the different precision matrices
      inla_cgeneric_smat_tp * all_precisions =
	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);

      if (any_constraints) {
	for (int i = 0; i < n_mesh; ++i) {
	  if (data->ints[constr_start + i]->len == 1 && data->ints[constr_start + i]->ints[0] < 0) continue;
	  remove_cols_and_rows_from_smat(all_precisions + i,
					 data->ints[constr_start + i]->ints,
					 data->ints[constr_start + i]->len);
	}
      }

      // Merge them together
      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + mesh_index[0]);
      for (int i = 1; i < n_realisations; ++i) {
	smat_smat_merge(&precision, all_precisions + mesh_index[i]);
      }

      // optimized format
      // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
      // where M is the length of Qij
      ret = Malloc(2 + precision.n, double);
      ret[0] = -1;                                   // code for optimized output
      ret[1] = precision.n;                          // number of (i <= j)
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.x[i];
      }

      // Avoid memory leaks
      for (int i = 0; i < n_mesh; ++i) {
	free_smat(all_precisions + i);
      }
      free(all_precisions);
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0

      // Z_b has zero mean
      ret = Calloc(1, double);
      break;
    }

  case INLA_CGENERIC_INITIAL:
    {
      // return c(M, initials)
      // where M is the number of hyperparameters

      // The initial values depend on how many parameters that are fixed or estimated.
      ret = Malloc(n_free + 1, double);
      ret[0] = n_free;
      int count = 0;
      for (int i = 0; i < n_theta; ++i) {
	if (!is_fixed[i]) {
	  ret[count + 1] = init[i];
	  ++count;
	}
      }
      assert(count == n_free);
      break;
    }

  case INLA_CGENERIC_LOG_NORM_CONST:
    {
      // return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself
      ret = NULL;
      break;
    }

  case INLA_CGENERIC_LOG_PRIOR:
    {
      // return c(LOG_PRIOR)

      // Only add terms for the parameters that are not fixed
      ret = Calloc(1, double);

      int count = 0;
      double * prior;
      if (!is_fixed[0]) { // rho
	prior = data->doubles[prior_start + count]->doubles;
	double lambda0 = -log(prior[1]) * prior[0];
	ret[0] += log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
	++count;
      }
      if (!is_fixed[1]) { // sigma
	prior = data->doubles[prior_start + count]->doubles;
	double lambda1 = -log(prior[1]) / prior[0];
	ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
	++count;
      }
      assert(count == n_free);
      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  return (ret);
}
