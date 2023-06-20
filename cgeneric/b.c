#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "cgeneric.h"
#include "smat-operations.c"
#include "spde-precision.c"

double * probit_model(inla_cgeneric_cmd_tp cmd,
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
  for (int i = 0; i < n_theta; ++i) if (!is_fixed[i]) ++n_free;

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

  assert(!strcmp(data->doubles[1]->name, "dist_to_s0"));
  double * dist_to_s0 = data->doubles[1]->doubles;
  int n_dist = data->doubles[1]->len;

  int prior_start = 2; // This is the first element of data->doubles that contains a prior
  assert(data->n_doubles == prior_start + n_free);

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, mu0, mu1;
  assert(n_theta == 4);
  if (theta) {
    double theta_full[4];
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
    mu0 = exp(theta_full[2]);
    mu1 = exp(theta_full[3]);
  } else {
    log_rho = log_sigma = mu0 = mu1 = NAN;
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

      diag_mat_tp mean_precision = diag(1, n_dist);
      smat_diag_merge(&precision, &mean_precision);
      free_diag(&mean_precision);

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

      double high_prec = exp(15);

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

      diag_mat_tp mean_precision = diag(high_prec, n_dist);
      smat_diag_merge(&precision, &mean_precision);
      free_diag(&mean_precision);

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
      ret = Calloc(n + 1, double);
      ret[0] = n;

      double sigma = exp(log_sigma);
      int count_start = n - n_dist + 1;
      for (int i = 0; i < n_dist; ++i) {
	ret[count_start + i] = sigma * (mu0 - mu1 * dist_to_s0[i]);
      }

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

      // This is an approximation for -0.5 * log(2 * pi)
      double gaussian_const = -0.91893853320467;

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
      if (!is_fixed[2]) { // mu0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(mu0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // mu1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(mu1) - prior[0], 2) / (2 * pow(prior[1], 2));
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



double * probit_model_old(inla_cgeneric_cmd_tp cmd,
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
  for (int i = 0; i < n_theta; ++i) if (!is_fixed[i]) ++n_free;

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
  int dist_start = prior_start + n_free; // This is the first element of data->doubles that contains dist_to_s0
  assert(data->n_doubles == n_mesh + dist_start);

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, mu0, mu1;
  assert(n_theta == 4);
  if (theta) {
    double theta_full[4];
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
    mu0 = exp(theta_full[2]);
    mu1 = exp(theta_full[3]);
  } else {
    log_rho = log_sigma = mu0 = mu1 = NAN;
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
      double sigma = exp(log_sigma);

      double ** mean_vals = Malloc(n_mesh, double *);
      for (int i = 0; i < n_mesh; ++i) {
	mean_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	for (int j = 0; j < data->doubles[dist_start + i]->len; ++j) {
	  mean_vals[i][j] = sigma * (mu0 - mu1 * data->doubles[dist_start + i]->doubles[j]);
	}
      }

      /*
      if (any_constraints) {
	for (int i = 0; i < n_mesh; ++i) {
	  remove_elements_from_vec(mean_vals[i],
				   data->ints[constr_start + i]->ints,
				   data->doubles[dist_start + i]->len,
				   data->ints[constr_start + i]->len);
	}
      }
      */

      ret = Malloc(n + 1, double);
      ret[0] = n;
      int count = 0;
      for (int i = 0; i < n_realisations; ++i) {
	for (int j = 0; j < data->doubles[dist_start + mesh_index[i]]->len; ++j) {
	  ++count;
	  ret[count] = mean_vals[mesh_index[i]][j];
	}
      }

      for (int i = 0; i < n_mesh; ++i) {
	free(mean_vals[i]);
      }
      free(mean_vals);

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

      // This is an approximation for -0.5 * log(2 * pi)
      double gaussian_const = -0.91893853320467;

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
      if (!is_fixed[2]) { // mu0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(mu0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // mu1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(mu1) - prior[0], 2) / (2 * pow(prior[1], 2));
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




double * spde_b_model_final(inla_cgeneric_cmd_tp cmd,
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



double softplus(double x, double k) {
  double res = x * k;
  if (res < 10) res = log(1 + exp(res));
  res /= k;
  return res;
}

double * spde_scale_model(inla_cgeneric_cmd_tp cmd,
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
  for (int i = 0; i < n_theta; ++i) if (!is_fixed[i]) ++n_free;

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
  int dist_start = prior_start + n_free; // This is the first element of data->doubles that contains dist_to_s0
  assert(data->n_doubles == n_mesh + dist_start);

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, s, lambda_s, kappa_s;
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
    s = exp(theta_full[2]);
    lambda_s = exp(theta_full[3]);
    kappa_s = exp(theta_full[4]);
  } else {
    log_rho = log_sigma = s = lambda_s = kappa_s = NAN;
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

      // Rescale all the precision matrices
      diag_mat_tp inv_scale;
      int scale_size = all_precisions[0].nrow;
      for (int i = 1; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      inv_scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	inv_scale.dim = all_precisions[i].nrow;
	for (int j = 0; j < inv_scale.dim; ++j) {
	  inv_scale.x[j] = 1 / (1 + s * exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda_s, kappa_s)));
	}
	diag_smat_diag_mult(all_precisions + i, &inv_scale);
      }
      free(inv_scale.x);

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

      // This is an approximation for -0.5 * log(2 * pi)
      double gaussian_const = -0.91893853320467;

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
      if (!is_fixed[2]) { // s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // lambda_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda_s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // kappa_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa_s) - prior[0], 2) / (2 * pow(prior[1], 2));
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

double * spde_scale_model2(inla_cgeneric_cmd_tp cmd,
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
  for (int i = 0; i < n_theta; ++i) if (!is_fixed[i]) ++n_free;

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
  int dist_start = prior_start + n_free; // This is the first element of data->doubles that contains dist_to_s0
  assert(data->n_doubles == n_mesh + dist_start);

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, s, lambda_s, kappa_s;
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
    s = exp(theta_full[2]);
    lambda_s = exp(theta_full[3]);
    kappa_s = exp(theta_full[4]);
  } else {
    log_rho = log_sigma = s = lambda_s = kappa_s = NAN;
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

      // Rescale all the precision matrices
      diag_mat_tp inv_scale;
      int scale_size = all_precisions[0].nrow;
      for (int i = 1; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      inv_scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	inv_scale.dim = all_precisions[i].nrow;
	for (int j = 0; j < inv_scale.dim; ++j) {
	  inv_scale.x[j] = 1 / (1 + s * (1 - exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda_s, kappa_s))));
	}
	diag_smat_diag_mult(all_precisions + i, &inv_scale);
      }
      free(inv_scale.x);

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

      // This is an approximation for -0.5 * log(2 * pi)
      double gaussian_const = -0.91893853320467;

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
      if (!is_fixed[2]) { // s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // lambda_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda_s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // kappa_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa_s) - prior[0], 2) / (2 * pow(prior[1], 2));
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

double * spde_ab_model(inla_cgeneric_cmd_tp cmd,
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
  double log_rho, log_sigma, beta0, lambda_b, kappa_b, lambda0, lambda1, kappa0, kappa1, scale0, scale1, lambda_s, kappa_s;
  assert(n_theta == 13);
  if (theta) {
    double theta_full[13];
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
    lambda_b = exp(theta_full[3]);
    kappa_b = exp(theta_full[4]);
    lambda0 = exp(theta_full[5]);
    lambda1 = exp(theta_full[6]);
    kappa0 = exp(theta_full[7]);
    kappa1 = exp(theta_full[8]);
    scale0 = exp(theta_full[9]);
    scale1 = exp(theta_full[10]);
    lambda_s = exp(theta_full[11]);
    kappa_s = exp(theta_full[12]);

  } else {
    log_rho = log_sigma = beta0 = lambda_b = kappa_b = lambda0 = lambda1 = kappa0 = kappa1 = scale0 = scale1 = lambda_s = kappa_s = NAN;
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

      // Rescale all the precision matrices
      diag_mat_tp inv_scale;
      int scale_size = all_precisions[0].nrow;
      for (int i = 1; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      inv_scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	inv_scale.dim = all_precisions[i].nrow;
	for (int j = 0; j < inv_scale.dim; ++j) {
	  inv_scale.x[j] = 1 /
	    (scale0 + scale1 * (1 - exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda_s, kappa_s))));
	}
	diag_smat_diag_mult(all_precisions + i, &inv_scale);
      }
      free(inv_scale.x);

      // Constrain all the precision matrices
      for (int i = 0; i < n_mesh; ++i) {
	remove_col_and_row_from_smat(all_precisions + i, s0_locs[i]);
      }

      // Merge them together
      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + s0_index[0]);
      for (int i = 1; i < n_repl; ++i) {
	smat_smat_merge(&precision, all_precisions + s0_index[i]);
      }
      assert(precision.nrow == n);

      // Compute all unique value of beta
      double ** beta_vals = Malloc(n_mesh, double *);
      for (int i = 0; i < n_mesh; ++i) {
	beta_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	for (int j = 0; j < data->doubles[dist_start + i]->len; ++j) {
	  beta_vals[i][j] = beta0 * exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda_b, kappa_b));
	}
      }

      // Compute 1 / b
      diag_mat_tp b_inv;
      b_inv.x = Malloc(n, double);
      b_inv.dim = n;
      int count = 0;
      double lambda, kappa, a;
      for (int i = 0; i < n_repl; ++i) {
	lambda = softplus(lambda0 - y0[i] * lambda1, 10);
	kappa = softplus(kappa0 - y0[i] * kappa1, 10);
	for (int j = 0; j < data->doubles[dist_start + s0_index[i]]->len; ++j) {
	  a = y0[i] * exp(-pow(data->doubles[dist_start + s0_index[i]]->doubles[j] / lambda, kappa));
	  b_inv.x[count] = 1 / (1 + pow(a, beta_vals[s0_index[i]][j]));
	  ++count;
	}
      }
      assert(count == n);

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
	free_smat(all_precisions + i);
	free(beta_vals[i]);
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

      ret = Malloc(n + 1, double);
      ret[0] = n;

      int count = 1;
      double lambda, kappa;
      for (int i = 0; i < n_repl; ++i) {
	lambda = softplus(lambda0 - y0[i] * lambda1, 10);
	kappa = softplus(kappa0 - y0[i] * kappa1, 10);
	for (int j = 0; j < data->doubles[dist_start + s0_index[i]]->len; ++j) {
	  ret[count] = y0[i] * exp(-pow(data->doubles[dist_start + s0_index[i]]->doubles[j] / lambda, kappa));
	  ++count;
	}
      }
      assert(count == (n + 1));
      
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
      if (!is_fixed[3]) { // lambda_b
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda_b) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // kappa_b
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa_b) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }

      if (!is_fixed[5]) { // lambda0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[6]) { // lambda1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda1) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[7]) { // kappa0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[8]) { // kappa1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa1) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[9]) { // scale0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(scale0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[10]) { // scale1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(scale1) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[11]) { // lambda_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda_s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[12]) { // kappa_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa_s) - prior[0], 2) / (2 * pow(prior[1], 2));
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




double * spde_ab_model_old(inla_cgeneric_cmd_tp cmd,
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
  double log_rho, log_sigma, beta0, lambda_b, kappa_b, lambda0, lambda1, kappa0, kappa1;
  assert(n_theta == 9);
  if (theta) {
    double theta_full[9];
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
    lambda_b = exp(theta_full[3]);
    kappa_b = exp(theta_full[4]);
    lambda0 = exp(theta_full[5]);
    lambda1 = exp(theta_full[6]);
    kappa0 = exp(theta_full[7]);
    kappa1 = exp(theta_full[8]);

  } else {
    log_rho = log_sigma = beta0 = lambda_b = kappa_b = lambda0 = lambda1 = kappa0 = kappa1 = NAN;
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
      assert(precision.nrow == n);

      // Compute all unique value of beta
      double ** beta_vals = Malloc(n_mesh, double *);
      for (int i = 0; i < n_mesh; ++i) {
	beta_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	for (int j = 0; j < data->doubles[dist_start + i]->len; ++j) {
	  beta_vals[i][j] = beta0 * exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda_b, kappa_b));
	}
      }

      // Compute 1 / b
      diag_mat_tp b_inv;
      b_inv.x = Malloc(n, double);
      b_inv.dim = n;
      int count = 0;
      double lambda, kappa, a;
      for (int i = 0; i < n_repl; ++i) {
	lambda = softplus(lambda0 - y0[i] * lambda1, 10);
	kappa = softplus(kappa0 - y0[i] * kappa1, 10);
	for (int j = 0; j < data->doubles[dist_start + s0_index[i]]->len; ++j) {
	  a = y0[i] * exp(-pow(data->doubles[dist_start + s0_index[i]]->doubles[j] / lambda, kappa));
	  b_inv.x[count] = 1 / (1 + pow(a, beta_vals[s0_index[i]][j]));
	  ++count;
	}
      }
      assert(count == n);

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
	free_smat(all_precisions + i);
	free(beta_vals[i]);
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

      ret = Malloc(n + 1, double);
      ret[0] = n;

      int count = 1;
      double lambda, kappa;
      for (int i = 0; i < n_repl; ++i) {
	lambda = softplus(lambda0 - y0[i] * lambda1, 10);
	kappa = softplus(kappa0 - y0[i] * kappa1, 10);
	for (int j = 0; j < data->doubles[dist_start + s0_index[i]]->len; ++j) {
	  ret[count] = y0[i] * exp(-pow(data->doubles[dist_start + s0_index[i]]->doubles[j] / lambda, kappa));
	  ++count;
	}
      }
      assert(count == (n + 1));
      
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
      if (!is_fixed[3]) { // lambda_b
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda_b) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // kappa_b
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa_b) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }

      if (!is_fixed[5]) { // lambda0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[6]) { // lambda1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda1) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[7]) { // kappa0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[8]) { // kappa1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa1) - prior[0], 2) / (2 * pow(prior[1], 2));
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
  double log_rho, log_sigma, beta0, lambda_b, kappa_b, s, lambda_s, kappa_s;
  assert(n_theta == 8);
  if (theta) {
    double theta_full[8];
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
    lambda_b = exp(theta_full[3]);
    kappa_b = exp(theta_full[4]);
    s = exp(theta_full[5]);
    lambda_s = exp(theta_full[6]);
    kappa_s = exp(theta_full[7]);
  } else {
    log_rho = log_sigma = beta0 = lambda_b = kappa_b = s = lambda_s = kappa_s = NAN;
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

      // Rescale all the precision matrices
      diag_mat_tp inv_scale;
      int scale_size = all_precisions[0].nrow;
      for (int i = 1; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      inv_scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	inv_scale.dim = all_precisions[i].nrow;
	for (int j = 0; j < inv_scale.dim; ++j) {
	  inv_scale.x[j] = 1 /
	    (1 + s * (1 - exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda_s, kappa_s))));
	}
	diag_smat_diag_mult(all_precisions + i, &inv_scale);
      }
      free(inv_scale.x);

      // Constrain all the precision matrices
      for (int i = 0; i < n_mesh; ++i) {
	remove_col_and_row_from_smat(all_precisions + i, s0_locs[i]);
      }

      // Merge them together
      inla_cgeneric_smat_tp precision = copy_smat(all_precisions + s0_index[0]);
      for (int i = 1; i < n_repl; ++i) {
	smat_smat_merge(&precision, all_precisions + s0_index[i]);
      }

      // Compute the parts of beta that depends on distance from a conditioning site and not y0 
      double ** beta_vals = Malloc(n_mesh, double *);
      for (int i = 0; i < n_mesh; ++i) {
	beta_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	for (int j = 0; j < data->doubles[dist_start + i]->len; ++j) {
	  beta_vals[i][j] = beta0 * exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda_b, kappa_b));
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
	//pow(log(beta0) - log(1 - beta0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // lambda_b
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda_b) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // kappa_b
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa_b) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[5]) { // s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[6]) { // lambda_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda_s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[7]) { // kappa_s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa_s) - prior[0], 2) / (2 * pow(prior[1], 2));
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

double * spde_b_model_old(inla_cgeneric_cmd_tp cmd,
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
  double log_rho, log_sigma, beta0, beta1, lambda, kappa;
  assert(n_theta == 6);
  if (theta) {
    double theta_full[6];
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
    //beta0 = beta0 / (1 + beta0);
    beta1 = exp(theta_full[3]);
    lambda = exp(theta_full[4]);
    kappa = exp(theta_full[5]);
  } else {
    log_rho = log_sigma = beta0 = beta1 = lambda = kappa = NAN;
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

      // Compute the parts of beta that depends on distance from a conditioning site and not y0 
      double ** beta_dist_vals = Malloc(n_mesh, double *);
      for (int i = 0; i < n_mesh; ++i) {
	beta_dist_vals[i] = Malloc(data->doubles[dist_start + i]->len, double);
	for (int j = 0; j < data->doubles[dist_start + i]->len; ++j) {
	  beta_dist_vals[i][j] = exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda, kappa));
	}
      }
      // Compute the parts of beta that depends on y0 and not the distance
      double * beta_y0_vals = Malloc(n_repl, double);
      for (int i = 0; i < n_repl; ++i) {
	beta_y0_vals[i] = softplus(beta0 - beta1 * y0[i], 10);
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
	  b_inv.x[count] = exp(-log_y0 * beta_y0_vals[i] * beta_dist_vals[s0_index[i]][j]);
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
	free(beta_dist_vals[i]);
	free_smat(all_precisions + i);
      }
      free(all_precisions);
      free(beta_dist_vals);
      free(beta_y0_vals);
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
	//pow(log(beta0) - log(1 - beta0) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // beta1
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(beta1) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // lambda
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[5]) { // kappa
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

double * spde_model3(inla_cgeneric_cmd_tp cmd,
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

  assert(!strcmp(data->ints[2]->name, "n_spline"));
  int n_spline = data->ints[2]->ints[0];
  int n_theta = n_spline + 1;

  // Which mesh do the different realisations come from?
  assert(!strcmp(data->ints[3]->name, "mesh_index"));
  int n_realisations = data->ints[3]->len;
  int * mesh_index = data->ints[3]->ints;

  // Look at the matrices required for building the meshes
  int n_mesh = data->n_smats / 4;
  assert(data->n_smats == n_mesh * 4); // Ensure that we did not get rounding when dividing
  assert(data->n_mats == n_mesh * 3);

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

  assert(!strcmp(data->doubles[1]->name, "rho_prior"));
  double * rho_prior = data->doubles[1]->doubles;
  assert(!strcmp(data->doubles[2]->name, "spline_prior"));
  double * spline_prior = data->doubles[2]->doubles;

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho;
  double * spline_params = Malloc(n_spline, double);
  if (theta) {
    log_rho = theta[0];
    for (int i = 0; i < n_spline; ++i) {
      spline_params[i] = theta[1 + i];
    }
  } else {
    log_rho = NAN;
    for (int i = 0; i < n_spline; ++i) {
      spline_params[i] = NAN;
    }
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
	spde_precision_multimesh(log_rho, 0.0, data->mats, data->smats + n_mesh, n_mesh);

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
	spde_precision_multimesh(log_rho, 0.0, data->mats, data->smats + n_mesh, n_mesh);

      diag_mat_tp scale;
      int scale_size = all_precisions[0].nrow;
      for (int i = 1; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	scale.dim = all_precisions[i].nrow;
	smat_vec_mult(scale.x, data->smats[i], spline_params);
	for (int j = 0; j < scale.dim; ++j) {
	  scale.x[j] = exp(-scale.x[j]);
	}
	diag_smat_diag_mult(all_precisions + i, &scale);
      }
      free(scale.x);

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
      ret = Malloc(n_theta + 1, double);
      ret[0] = n_theta;
      for (int i = 0; i < n_theta; ++i) {
	ret[i + 1] = init[i];
      }
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

      ret[0] += gaussian_const - log(rho_prior[1]) -
	pow(log_rho - rho_prior[0], 2) / (2 * rho_prior[1] * rho_prior[1]);

      for (int i = 0; i < n_spline; ++i) {
	ret[0] += gaussian_const - log(spline_prior[1]) -
	  pow(spline_params[i] - spline_prior[0], 2) / (2 * spline_prior[1] * spline_prior[1]);
      }

      break;
    }

  case INLA_CGENERIC_QUIT:
  default:
    break;
  }

  free(spline_params);
  return (ret);
}

/*
double * spde_model3(inla_cgeneric_cmd_tp cmd,
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
  int dist_start = prior_start + n_free;

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, a, b;
  assert(n_theta == 4);
  if (theta) {
    double theta_full[4];
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
    a = exp(theta_full[2]);
    b = exp(theta_full[3]);
  } else {
    log_rho = log_sigma = a = b = NAN;
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

      diag_mat_tp scale;
      int scale_size = 0;
      for (int i = 0; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	scale.dim = all_precisions[i].nrow;
	for (int j = 0; j < scale.dim; ++j) {
	  scale.x[j] = 1 / sqrt(a + pow(data->doubles[dist_start + i]->doubles[j], -b));
	}
	diag_smat_diag_mult(all_precisions + i, &scale);
      }
      free(scale.x);

      if (any_constraints) {
	for (int i = 0; i < n_mesh; ++i) {
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
      if (!is_fixed[2]) { // a
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(a) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // slope
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(b) - prior[0], 2) / (2 * pow(prior[1], 2));
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
*/


// double * spde_model3(inla_cgeneric_cmd_tp cmd,
// 		     double * theta,
// 		     inla_cgeneric_data_tp * data) {
// 
//   double * ret = NULL;
// 
//   // ==========================================================
//   // Assert that all the input variables look as they shold
//   // ==========================================================
// 
//   // Size of the latent field
//   assert(!strcmp(data->ints[0]->name, "n"));
//   int n = data->ints[0]->ints[0];
//   assert(n > 0);
// 
//   // Which mesh do the different realisations come from?
//   assert(!strcmp(data->ints[2]->name, "mesh_index"));
//   int n_realisations = data->ints[2]->len;
//   int * mesh_index = data->ints[2]->ints;
// 
//   // Which of the parameters should be estimated, and which are fixed?
//   assert(!strcmp(data->ints[3]->name, "is_fixed"));
//   int * is_fixed = data->ints[3]->ints;
//   int n_theta = data->ints[3]->len;
//   int n_free = 0;
//   for (int i = 0; i < n_theta; ++i) {
//     if (!is_fixed[i]) ++n_free;
//   }
// 
//   // Look at the matrices required for building the meshes
//   int n_mesh = data->n_smats / 3;
//   assert(data->n_smats == n_mesh * 3); // Ensure that we did not get rounding when dividing
//   assert(data->n_smats == data->n_mats); // We need the same number of mats and smats
// 
//   // What is the first element of data->ints that contains constraining indices?
//   int constr_start = 4;
//   // There should either be 0 or n_mesh constraining index vectors
//   assert(data->n_ints == constr_start ||
// 	 data->n_ints == constr_start + n_mesh);
//   int any_constraints;
//   if (data->n_ints == constr_start) {
//     any_constraints = 0;
//   } else {
//     any_constraints = 1;
//   }
// 
//   // Initial values
//   assert(!strcmp(data->doubles[0]->name, "init"));
//   double * init = data->doubles[0]->doubles;
//   assert(data->doubles[0]->len == n_theta);
// 
//   int prior_start = 1; // This is the first element of data->doubles that contains a prior
//   int dist_start = prior_start + n_free;
// 
//   // ==========================================================
//   // Initiate the correct parameter values, depending on whether
//   // the parameters are fixed or estimated.
//   // ==========================================================
//   double log_rho, log_sigma, slope;
//   assert(n_theta == 3);
//   if (theta) {
//     double theta_full[3];
//     int count = 0;
//     for (int i = 0; i < n_theta; ++i) {
//       if (is_fixed[i]) {
// 	theta_full[i] = init[i];
//       } else {
// 	theta_full[i] = theta[count];
// 	++count;
//       }
//     }
//     assert(count == n_free);
//     log_rho = theta_full[0];
//     log_sigma = theta_full[1];
//     slope = exp(theta_full[2]);
//   } else {
//     log_rho = log_sigma = slope = NAN;
//   }
// 
//   // ==========================================================
//   // This switch statement is the required method for implementing
//   // cgeneric models in R-INLA
//   // ==========================================================
//   switch (cmd) {
//   case INLA_CGENERIC_VOID:
//     {
//       assert(!(cmd == INLA_CGENERIC_VOID));
//       break;
//     }
// 
//   case INLA_CGENERIC_GRAPH:
//     {
//       // Compute all the different precision matrices
//       inla_cgeneric_smat_tp * all_precisions =
// 	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);
// 
//       if (any_constraints) {
// 	for (int i = 0; i < n_mesh; ++i) {
// 	  remove_cols_and_rows_from_smat(all_precisions + i,
// 					 data->ints[constr_start + i]->ints,
// 					 data->ints[constr_start + i]->len);
// 	}
//       }
// 
//       // Merge them together
//       inla_cgeneric_smat_tp precision = copy_smat(all_precisions + mesh_index[0]);
//       for (int i = 1; i < n_realisations; ++i) {
// 	smat_smat_merge(&precision, all_precisions + mesh_index[i]);
//       }
// 
//       // Extract all the necessary i and j values from the block matrix
//       ret = Malloc(2 + 2 * precision.n, double);
//       ret[0] = n;
//       ret[1] = precision.n;
//       for (int i = 0; i < precision.n; i++) {
// 	ret[2 + i] = precision.i[i];                                  // i
// 	ret[2 + precision.n + i] = precision.j[i];                    // j
//       }
// 
//       // Avoid memory leaks
//       for (int i = 0; i < n_mesh; ++i) {
// 	free_smat(all_precisions + i);
//       }
//       free_smat(&precision);
// 
//       break;
//     }
// 
//   case INLA_CGENERIC_Q:
//     {
// 
//       // Compute all the different precision matrices
//       inla_cgeneric_smat_tp * all_precisions =
// 	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);
// 
//       diag_mat_tp scale;
//       int scale_size = 0;
//       for (int i = 0; i < n_mesh; ++i) {
// 	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
//       }
//       scale.x = Malloc(scale_size, double);
//       double tmp;
//       for (int i = 0; i < n_mesh; ++i) {
// 	scale.dim = all_precisions[i].nrow;
// 	for (int j = 0; j < scale.dim; ++j) {
// 	  tmp = 1 - slope * data->doubles[dist_start + i]->doubles[j];
// 	  if (tmp < 0.00001) tmp = 0.00001;
// 	  scale.x[j] = 1 / sqrt(tmp);
// 	}
// 	diag_smat_diag_mult(all_precisions + i, &scale);
//       }
//       free(scale.x);
// 
//       if (any_constraints) {
// 	for (int i = 0; i < n_mesh; ++i) {
// 	  remove_cols_and_rows_from_smat(all_precisions + i,
// 					 data->ints[constr_start + i]->ints,
// 					 data->ints[constr_start + i]->len);
// 	}
//       }
// 
//       // Merge them together
//       inla_cgeneric_smat_tp precision = copy_smat(all_precisions + mesh_index[0]);
//       for (int i = 1; i < n_realisations; ++i) {
// 	smat_smat_merge(&precision, all_precisions + mesh_index[i]);
//       }
// 
//       // optimized format
//       // return c(n, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
//       // where M is the length of Qij
//       ret = Malloc(2 + precision.n, double);
//       ret[0] = -1;                                   // code for optimized output
//       ret[1] = precision.n;                          // number of (i <= j)
//       for (int i = 0; i < precision.n; i++) {
// 	ret[2 + i] = precision.x[i];
//       }
// 
//       // Avoid memory leaks
//       for (int i = 0; i < n_mesh; ++i) {
// 	free_smat(all_precisions + i);
//       }
//       free_smat(&precision);
// 
//       break;
//     }
// 
//   case INLA_CGENERIC_MU:
//     {
//       // return (N, mu)
//       // if N==0 then mu is not needed as its taken to be mu[]==0
// 
//       // Z_b has zero mean
//       ret = Calloc(1, double);
//       break;
//     }
// 
//   case INLA_CGENERIC_INITIAL:
//     {
//       // return c(M, initials)
//       // where M is the number of hyperparameters
// 
//       // The initial values depend on how many parameters that are fixed or estimated.
//       ret = Malloc(n_free + 1, double);
//       ret[0] = n_free;
//       int count = 0;
//       for (int i = 0; i < n_theta; ++i) {
// 	if (!is_fixed[i]) {
// 	  ret[count + 1] = init[i];
// 	  ++count;
// 	}
//       }
//       assert(count == n_free);
//       break;
//     }
// 
//   case INLA_CGENERIC_LOG_NORM_CONST:
//     {
//       // return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself
//       ret = NULL;
//       break;
//     }
// 
//   case INLA_CGENERIC_LOG_PRIOR:
//     {
//       // return c(LOG_PRIOR)
// 
//       // Only add terms for the parameters that are not fixed
//       ret = Calloc(1, double);
// 
//       // This is an approximation for -0.5 * log(2 * pi)
//       double gaussian_const = -0.91893853320467;
// 
//       int count = 0;
//       double * prior;
//       if (!is_fixed[0]) { // rho
// 	prior = data->doubles[prior_start + count]->doubles;
// 	double lambda0 = -log(prior[1]) * prior[0];
// 	ret[0] += log(lambda0) - lambda0 * exp(-log_rho) - log_rho;
// 	++count;
//       }
//       if (!is_fixed[1]) { // sigma
// 	prior = data->doubles[prior_start + count]->doubles;
// 	double lambda1 = -log(prior[1]) / prior[0];
// 	ret[0] += log(lambda1) - lambda1 * exp(log_sigma) + log_sigma;
// 	++count;
//       }
//       if (!is_fixed[2]) { // slope
// 	prior = data->doubles[prior_start + count]->doubles;
// 	ret[0] += gaussian_const - log(prior[1]) -
// 	  pow(log(slope) - prior[0], 2) / (2 * pow(prior[1], 2));
// 	++count;
//       }
//       assert(count == n_free);
//       break;
//     }
// 
//   case INLA_CGENERIC_QUIT:
//   default:
//     break;
//   }
// 
//   return (ret);
// }



/*
double * spde_model3(inla_cgeneric_cmd_tp cmd,
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
  int dist_start = prior_start + n_free;

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, lambda, kappa, s0, delta;
  assert(n_theta == 6);
  if (theta) {
    double theta_full[6];
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
    lambda = exp(theta_full[2]);
    kappa = exp(theta_full[3]);
    s0 = exp(theta_full[4]);
    delta = exp(theta_full[5]);
  } else {
    log_rho = log_sigma = lambda = kappa = s0 = delta = NAN;
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
      free_smat(&precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {

      // Compute all the different precision matrices
      inla_cgeneric_smat_tp * all_precisions =
	spde_precision_multimesh(log_rho, log_sigma, data->mats, data->smats, n_mesh);

      diag_mat_tp scale;
      int scale_size = 0;
      for (int i = 0; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	scale.dim = all_precisions[i].nrow;
	for (int j = 0; j < scale.dim; ++j) {
	  if (data->doubles[dist_start + i]->doubles[j] < delta) {
	    scale.x[j] = 1 / sqrt(s0 + 1);
	  } else {
	    scale.x[j] = 1 / sqrt(s0 + exp(-pow((data->doubles[dist_start + i]->doubles[j] - delta) / lambda, kappa)));
	  }
	}
	diag_smat_diag_mult(all_precisions + i, &scale);
      }
      free(scale.x);

      if (any_constraints) {
	for (int i = 0; i < n_mesh; ++i) {
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
      if (!is_fixed[2]) { // lambda
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // kappa
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // s0
      	prior = data->doubles[prior_start + count]->doubles;
      	ret[0] += gaussian_const - log(prior[1]) -
      	  pow(log(s0) - prior[0], 2) / (2 * pow(prior[1], 2));
      	++count;
      }
      if (!is_fixed[5]) { // delta
      	prior = data->doubles[prior_start + count]->doubles;
      	ret[0] += gaussian_const - log(prior[1]) -
      	  pow(log(delta) - prior[0], 2) / (2 * pow(prior[1], 2));
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
*/



double * spde_model2(inla_cgeneric_cmd_tp cmd,
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
  int dist_start = prior_start + n_free;

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, lambda, kappa, s0, s1;
  assert(n_theta == 6);
  if (theta) {
    double theta_full[6];
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
    lambda = exp(theta_full[2]);
    kappa = exp(theta_full[3]);
    s0 = exp(theta_full[4]);
    s1 = exp(theta_full[5]);
  } else {
    log_rho = log_sigma = lambda = kappa = s0 = s1 = NAN;
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

      diag_mat_tp scale;
      int scale_size = 0;
      for (int i = 0; i < n_mesh; ++i) {
	if (all_precisions[i].nrow > scale_size) scale_size = all_precisions[i].nrow;
      }
      scale.x = Malloc(scale_size, double);
      for (int i = 0; i < n_mesh; ++i) {
	scale.dim = all_precisions[i].nrow;
	for (int j = 0; j < scale.dim; ++j) {
	  scale.x[j] = 1 / sqrt(s0 + s1 * exp(-pow(data->doubles[dist_start + i]->doubles[j] / lambda, kappa)));
	}
	diag_smat_diag_mult(all_precisions + i, &scale);
      }
      free(scale.x);

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
      if (!is_fixed[2]) { // lambda
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // kappa
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(kappa) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // s0
      	prior = data->doubles[prior_start + count]->doubles;
      	ret[0] += gaussian_const - log(prior[1]) -
      	  pow(log(s0) - prior[0], 2) / (2 * pow(prior[1], 2));
      	++count;
      }
      if (!is_fixed[5]) { // s1
      	prior = data->doubles[prior_start + count]->doubles;
      	ret[0] += gaussian_const - log(prior[1]) -
      	  pow(log(s1) - prior[0], 2) / (2 * pow(prior[1], 2));
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




double * spde_probit_model(inla_cgeneric_cmd_tp cmd,
			   double * theta,
			   inla_cgeneric_data_tp * data) {

  double * ret = NULL;

  // The fixed precision used in the precision matrix of the mean structure
  double high_prec = exp(15);

  // ==========================================================
  // Assert that all the input variables look as they shold
  // ==========================================================

  // Size of the latent field
  assert(!strcmp(data->ints[0]->name, "n"));
  int n = data->ints[0]->ints[0];
  assert(n > 0);

  // How many time points do we have?
  assert(!strcmp(data->ints[2]->name, "n_time"));
  int n_time = data->ints[2]->ints[0];

  // Which of the parameters should be estimated, and which are fixed?
  assert(!strcmp(data->ints[3]->name, "is_fixed"));
  int * is_fixed = data->ints[3]->ints;
  int n_theta = data->ints[3]->len;
  int n_free = 0;
  for (int i = 0; i < n_theta; ++i) {
    if (!is_fixed[i]) ++n_free;
  }

  // Initial values
  assert(!strcmp(data->doubles[0]->name, "init"));
  double * init = data->doubles[0]->doubles;
  assert(data->doubles[0]->len == n_theta);

  // dist_to_big
  assert(!strcmp(data->doubles[1]->name, "dist_to_big"));
  double * dist_to_big = data->doubles[1]->doubles;
  int n_y = data->doubles[1]->len;

  int prior_start = 2; // This is the first element of data->doubles that contains a prior

  // B0, B1 and B2, required for building the SPDE precision
  assert(!strcasecmp(data->mats[0]->name, "B0"));
  int n_spde = data->mats[0]->nrow;
  assert(n == n_y + n_time * n_spde);
  assert(data->mats[0]->ncol == 3);
  assert(!strcasecmp(data->mats[1]->name, "B1"));
  assert(data->mats[1]->nrow == n_spde);
  assert(data->mats[1]->ncol == 3);
  assert(!strcasecmp(data->mats[2]->name, "B2"));
  assert(data->mats[2]->nrow == n_spde);
  assert(data->mats[2]->ncol == 3);

  // M0, M1 and M2, required for building the SPDE precision
  assert(!strcasecmp(data->smats[0]->name, "M0"));
  assert(data->smats[0]->nrow == n_spde);
  assert(data->smats[0]->ncol == n_spde);
  assert(!strcasecmp(data->smats[1]->name, "M1"));
  assert(data->smats[1]->nrow == n_spde);
  assert(data->smats[1]->ncol == n_spde);
  assert(!strcasecmp(data->smats[2]->name, "M2"));
  assert(data->smats[2]->nrow == n_spde);
  assert(data->smats[2]->ncol == n_spde);

  // ==========================================================
  // Initiate the correct parameter values, depending on whether
  // the parameters are fixed or estimated.
  // ==========================================================
  double log_rho, log_sigma, mu0, s, lambda, kappa;
  assert(n_theta == 6);
  if (theta) {
    double theta_full[6];
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
    mu0 = theta_full[2];
    s = exp(theta_full[3]);
    lambda = exp(theta_full[4]);
    kappa = exp(theta_full[5]);
  } else {
    log_rho = log_sigma = mu0 = s = lambda = kappa = NAN;
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
      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);

      block_diag_smat(&precision, n_time);

      diag_mat_tp mean_precision = diag(high_prec, n_y);

      smat_diag_merge(&precision, &mean_precision);

      // Extract all the necessary i and j values from the block matrix
      ret = Malloc(2 + 2 * precision.n, double);
      ret[0] = n;
      ret[1] = precision.n;
      for (int i = 0; i < precision.n; i++) {
	ret[2 + i] = precision.i[i];                                  // i
	ret[2 + precision.n + i] = precision.j[i];                    // j
      }

      // Avoid memory leaks
      free_smat(&precision);
      free_diag(&mean_precision);

      break;
    }

  case INLA_CGENERIC_Q:
    {

      inla_cgeneric_smat_tp precision = spde_precision(log_rho, log_sigma, data->mats, data->smats);

      block_diag_smat(&precision, n_time);

      diag_mat_tp mean_precision = diag(high_prec, n_y);

      smat_diag_merge(&precision, &mean_precision);

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
      free_smat(&precision);
      free_diag(&mean_precision);

      break;
    }

  case INLA_CGENERIC_MU:
    {
      // return (N, mu)
      // if N==0 then mu is not needed as its taken to be mu[]==0

      // Z_b has zero mean
      ret = Calloc(n + 1, double);
      ret[0] = n;

      int start = 1 + n_spde * n_time;
      for (int i = 0; i < n_y; ++i) {
	ret[start + i] = mu0 + s * exp(-pow(dist_to_big[i] / lambda, kappa));
      }

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
      if (!is_fixed[2]) { // mu0
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(mu0 - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[3]) { // s
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(s) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[4]) { // lambda
	prior = data->doubles[prior_start + count]->doubles;
	ret[0] += gaussian_const - log(prior[1]) -
	  pow(log(lambda) - prior[0], 2) / (2 * pow(prior[1], 2));
	++count;
      }
      if (!is_fixed[5]) { // kappa
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
