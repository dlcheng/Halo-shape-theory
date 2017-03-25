#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "var.h"
#include "proto.h"
#include "define.h"

double p_a_over_c(double nu_2, double M, double z)  /* nu_2 = a/c */
{
  int WORKSIZE = 10000;
  double a = 1.0 / (1.0 + z);
  double param_list[3];
  
  param_list[0] = nu_2;
  param_list[1] = z;
  param_list[2] = M;
  
  gsl_function F;
  gsl_integration_workspace *workspace;
  double result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &combine_intg_a_c;
  F.params = &param_list;
 
  gsl_integration_qag(&F, 0, a, 0, TOL_ERROR, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  
  gsl_integration_workspace_free(workspace);
	
  return result;
	
}                                                 /* end p_a_over_c */


double combine_intg_a_c(double a_f, void * params)
{
  double *param_list = (double *)params;
  double nu_2 = param_list[0];
  double z = param_list[1];
  double M = param_list[2];
  double z_f = 1.0 / a_f - 1.0;
  double param_list_1[3];
  double result = 1.0;

#ifdef LEE
  param_list_1[0] = z_f;
  param_list_1[1] = M * 2.0;
  param_list_1[2] = 1.0;                            /* renomalize the distribution */
   
  result = form_prob(2.0 * M, z_f, z) * p_nu_2_fixed_z(nu_2, param_list_1) / a_f / a_f;
#endif

#ifdef DL  
  param_list_1[0] = z_f;
  param_list_1[1] = M / 2.0;
  param_list_1[2] = 1.0;                            /* renomalize the distribution */
   
  result = form_prob(M, z_f, z) * p_nu_2_fixed_z(nu_2, param_list_1) / a_f / a_f;
#endif	

  return result;
}                                                   /* end combine_intg_a_c */


double p_b_over_c(double nu_1, double M, double z)  /* nu_1 = b/c */
{
  int WORKSIZE = 10000;
  double a = 1.0 / (1.0 + z);
  double param_list[3];
  
  param_list[0] = nu_1;
  param_list[1] = z;
  param_list[2] = M;
  
  gsl_function F;
  gsl_integration_workspace *workspace;
  double result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &combine_intg_b_c;
  F.params = &param_list;
 
  gsl_integration_qag(&F, 0, a, 0, TOL_ERROR, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  
  gsl_integration_workspace_free(workspace);
	
  return result;
	
}                                                    /* end p_b_over_c */


double combine_intg_b_c(double a_f, void * params)
{
  double *param_list = (double *)params;
  double nu_1 = param_list[0];
  double z = param_list[1];
  double M = param_list[2];
  double z_f = 1.0 / a_f - 1.0;
  double param_list_1[3];
  double result = 1.0;

#ifdef LEE
  param_list_1[0] = z_f;
  param_list_1[1] = M * 2.0;
  param_list_1[2] = 1.0;                        /* renomalize the distribution */
   
  result = form_prob(2.0 * M, z_f, z) * p_nu_1_fixed_z(nu_1, param_list_1) / a_f / a_f;
#endif

#ifdef DL  
  param_list_1[0] = z_f;
  param_list_1[1] = M / 2.0;
  param_list_1[2] = 1.0;                         /* renomalize the distribution */
   
  result = form_prob(M, z_f, z) * p_nu_1_fixed_z(nu_1, param_list_1) / a_f / a_f;
#endif	

  return result;
}                                                 /* end combine_intg_b_c */
