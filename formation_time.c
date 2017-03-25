#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_erf.h>

#include "var.h"
#include "proto.h"
#include "define.h"


double bar_alpha_eff(double z, double M)
{
  double sigma_cdm_square = sigma_square(z, M);
  double alpha_eff;
  
  gsl_function F;
  double result, abserr;
  
  F.function = &sigma_diff_over_m;
  F.params = &z;
	
  gsl_deriv_central(&F, M, 1e-6 * M, &result, &abserr);
	
  alpha_eff = -1.0 * M / sigma_cdm_square * result;
  
  return alpha_eff * (0.6268 + 0.3058 * alpha_eff);
}                                         /* end bar_alpha_eff */

double sigma_diff_over_m(double M, void * params)
{
  double z = *(double *)params;
  
  return sigma_square(z, M);	
	 
}                                         /* end sigma_diff_over_m */

double form_prob(double M, double z_f, double z)   /* halo of mass M, observed at z formed at z_f */
{ 
  double d_z_f = growth_factor(z_f);
  double delta_c_z_f = DELTA_C / d_z_f;
  double delta_c_z = DELTA_C / growth_factor(z);
  double sigma_dimo = sqrt(sigma_square(z_i, M /2.0) - sigma_square(z_i, M));
  double w_f = (delta_c_z_f - delta_c_z) / sigma_dimo;
  double w_f_diff_over_z_f;
  
  gsl_function F;
  double result, abserr;
  
  F.function = &growth_factor_diff_over_z;
	
  gsl_deriv_central(&F, z_f, 1e-6 * z_f, &result, &abserr);
	
  w_f_diff_over_z_f = -1.0 * DELTA_C * result / d_z_f / d_z_f / sigma_dimo;
	
	
  double alpha = bar_alpha_eff(z_i, M);
  double c_alpha = 1.0 - (1.0 - alpha) / 25.0;
  double a_alpha = sqrt(8.0 / Pi) * (1.0 - alpha) * ( 0.0107 + 0.0163 * alpha);
  double b_alpha = 2.0 / a_alpha * (c_alpha - (pow(2, alpha) - 1.0 ) / alpha);
  
  double p_diff_over_w_f = a_alpha / (1.0 + b_alpha * w_f) * exp(-5.0 * w_f * w_f);
  p_diff_over_w_f += 2.0 * c_alpha * w_f * gsl_sf_erfc(w_f / sqrt(2));
  
	
  return p_diff_over_w_f * w_f_diff_over_z_f;
	
}                                            /* end form_prob */

double growth_factor_diff_over_z(double z, void * params)
{
  return growth_factor(z);
	
}                                           /* end growth_factor_diff_over_z */
