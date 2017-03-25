#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "var.h"
#include "proto.h"
#include "define.h"

double un_norm_p_nu(double nu_1, double nu_2, double sigma_m, double d_f)
{
#ifdef LEE	
  double deta_cf = DELTA_C / d_f;
  double lambda_cf = LAMBDA_C / d_f;	
#endif

#ifdef DL
  double deta_cf = DELTA_C / d_f;
  double lambda_cf = LAMBDA_C / d_f;	
#endif	

  double lambda_1, lambda_2;
  double p_lambda;
  double lambda_nu_deform;

#ifdef LEE  
  
  lambda_1 = (1.0 + (deta_cf * d_f - 2.0) * nu_2 * nu_2 + nu_1 * nu_1) /
             (d_f * (1.0 + nu_1 * nu_1 + nu_2 * nu_2));
  
  lambda_2 = (1.0 + (deta_cf * d_f - 2.0) * nu_1 * nu_1 + nu_2 * nu_2) /
             (d_f * (1.0 + nu_1 * nu_1 + nu_2 * nu_2));
                          
  lambda_nu_deform = 4.0 * pow(deta_cf * d_f - 3.0, 2) * nu_1 * nu_2 /
                     (d_f * d_f * pow(1.0 + nu_1 * nu_1 + nu_2 * nu_2 , 3.0));
                     
#endif

#ifdef DL
  lambda_1 = (1.0 + (deta_cf * d_f - 2.0) * nu_2 + nu_1) /
             (d_f * (1.0 + nu_1 + nu_2));
  
  lambda_2 = (1.0 + (deta_cf * d_f - 2.0) * nu_1 + nu_2) /
             (d_f * (1.0 + nu_1 + nu_2));

  lambda_nu_deform = pow(deta_cf * d_f - 3.0 , 2) / (d_f * d_f * pow(1.0 + nu_1 + nu_2, 3.0));

#endif                     
                     
   if(lambda_1 > (1.0 / d_f))
     return 0;
     
   if((lambda_1 + lambda_2) > (deta_cf - lambda_cf))
     return 0;		
	
   p_lambda = -2.5 * deta_cf * deta_cf / pow(sigma_m, 2);
   p_lambda += 7.5 * deta_cf * (lambda_1 + lambda_2) / pow(sigma_m , 2);
   p_lambda += -7.5 * (lambda_1 * lambda_1 + lambda_1 * lambda_2 + lambda_2 * lambda_2) / pow(sigma_m , 2);
   p_lambda = exp(p_lambda);
   p_lambda *= ( 2.0 * lambda_1 + lambda_2 - deta_cf);
   p_lambda *= (lambda_1 - lambda_2);
   p_lambda *= (lambda_1 + 2.0 * lambda_2 - deta_cf);
   p_lambda *= 3375.0 / 4.0 / sqrt(10 * Pi) / pow(sigma_m, 5);  /* coefficient differnt from Lee et al. 2005,
                                                                   it's wrong there. */
   
   return p_lambda * lambda_nu_deform;
	
}                       /* end un_norm_p_nu */

double un_norm_intg_kernel_over_nu_1(double nu_1, void * params)
{
  double * param_list = (double *) params;
  double nu_2 = param_list[0];
  double sigma_m = param_list[1];
  double d_f = param_list[2];
  
  return un_norm_p_nu(nu_1, nu_2, sigma_m, d_f);
		
}                      /* end un_norm_intg_kernel_over_nu_1 */

double un_norm_intg_over_nu_1(double nu_2, void * params)
{
 int WORKSIZE = 100000;

 double * param_list_1 = (double *) params;  
 double param_list_2[3];
 
 param_list_2[0] = nu_2;
 param_list_2[1] = param_list_1[0];  /* sigma_m */
 param_list_2[2] = param_list_1[1];  /* d_f */
  
 gsl_function F;
 gsl_integration_workspace *workspace;
 double result, abserr;

 workspace = gsl_integration_workspace_alloc(WORKSIZE);
 F.function = &un_norm_intg_kernel_over_nu_1;
 F.params = &param_list_2;
 
 gsl_integration_qag(&F, nu_2, 1, 0, TOL_ERROR, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  
 gsl_integration_workspace_free(workspace);

 return result;	
}                     /* end un_norm_intg_over_nu_1 */

double un_norm_intg_kernel_over_nu_2(double nu_2, void *params)
{
  double * param_list = (double *)params;	
  double nu_1 = param_list[0];
  double sigma_m = param_list[1];
  double d_f = param_list[2];
  
  return un_norm_p_nu(nu_1, nu_2, sigma_m, d_f);
	
}                      /* end un_norm_intg_kernel_over_nu_2 */

double un_norm_intg_over_nu_2(double nu_1, void * params)
{
 int WORKSIZE = 100000;

 double * param_list_1 = (double *) params;  
 double param_list_2[3];
 
 param_list_2[0] = nu_1;
 param_list_2[1] = param_list_1[0];  /* sigma_m */
 param_list_2[2] = param_list_1[1];  /* d_f */
  
 gsl_function F;
 gsl_integration_workspace *workspace;
 double result, abserr;

 workspace = gsl_integration_workspace_alloc(WORKSIZE);
 F.function = &un_norm_intg_kernel_over_nu_2;
 F.params = &param_list_2;
 
 gsl_integration_qag(&F, 0, nu_1, 0, TOL_ERROR, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  
 gsl_integration_workspace_free(workspace);

 return result;	
}                     /* end un_norm_intg_over_nu_1 */


void init_prob_norm(double z_f, double M)    /* for each (z, M) do the normalization */
{
  int WORKSIZE = 100000;
  double param_list[2];
  
#ifndef LOCAL_NG  
  double sigma_m = sqrt(sigma_square(z_i, M)); /* sigma is calculated at z_i instead of z_f */
#else
  double sigma_m = sqrt(local_ng_sigma_square(z_i, M)); /* sigma is calculated at z_i instead of z_f */
#endif  

  double d_f = growth_factor(z_f);
	
  param_list[0] = sigma_m;
  param_list[1] = d_f;
  
  gsl_function F;
  gsl_integration_workspace *workspace;
  double result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &un_norm_intg_over_nu_1;
  F.params = &param_list;
 
  gsl_integration_qag(&F, 0, 1, 0, TOL_ERROR, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  
  gsl_integration_workspace_free(workspace);

  norm_prob_A = 1.0 / result;
	
}                   /* end int_prob_norm */

double p_nu_2_fixed_z(double nu_2, void *params)  /* nu_2 = a/c, params[2] is a flag, if >0, renormalize  */
{
  double * param_list_1 = (double *) params; /* recieve the parameter list */
  double param_list_2[2];                    /* make a new parameter list */
  
  double z_f = param_list_1[0];
  double M = param_list_1[1];
  
  if(param_list_1[2] > 0)
    init_prob_norm(z_f, M);
#ifndef LOCAL_NG    
  param_list_2[0] = sqrt(sigma_square(z_i, M));
#else
  param_list_2[0] = sqrt(local_ng_sigma_square(z_i, M));
#endif
  param_list_2[1] = growth_factor(z_f);
  
  return norm_prob_A * un_norm_intg_over_nu_1(nu_2, param_list_2);
	
}                   /* end p_nu_2_fixed_z */


double p_nu_1_fixed_z(double nu_1, void *params)  /* nu_1 = b/c, params[2] is a flag, if >0, renormalize  */
{
  double * param_list_1 = (double *) params; /* recieve the parameter list */
  double param_list_2[2];                    /* make a new parameter list */
  
  double z_f = param_list_1[0];
  double M = param_list_1[1];
  
  if(param_list_1[2] > 0)
    init_prob_norm(z_f, M);
#ifndef LOCAL_NG     
  param_list_2[0] = sqrt(sigma_square(z_i, M));
#else
  param_list_2[0] = sqrt(local_ng_sigma_square(z_i, M));
#endif
  param_list_2[1] = growth_factor(z_f);
  
  return norm_prob_A * un_norm_intg_over_nu_2(nu_1, param_list_2);
	
}                   /* end p_nu_1_fixed_z */
