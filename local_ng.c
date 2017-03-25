#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "var.h"
#include "proto.h"
#include "define.h"

#ifdef LOCAL_NG
double local_ng_f_k(double k, double R)              /* defined as T_k^2  * k^(2+2*ns-1) * W(k*R)^2 */
{
 double tk = T_k(k);
 double k_ns = pow(k, 2*ns+1);
 double wk,k_R;

 k_R = R * k;        /* k*R */
 tk = tk * tk;

 if(k_R < 1e-3)
   wk = 1.0 / 3.0 * (1.0 - 0.1 * k_R * k_R);
 else 
   wk = sin(k_R) / k_R / k_R / k_R - cos(k_R) / k_R / k_R;

 wk = wk * wk;

 return tk * k_ns * wk;
}                                  /* end local_ng_f_k */

double local_ng_f_q(double q, void *params)
{
 double R = *(double *)params;
 return local_ng_f_k(1.0 / q - 1.0, R) / q / q;
}                                  /* end local_ng_f_q */

double local_ng_un_norm_sigma_square(double z, double R)
{
 int WORKSIZE = 100000;
 double d_growth = growth_factor(z);
 d_growth *= d_growth;

 double c = 4.5 * d_growth / Pi_square;

 gsl_function F;
 gsl_integration_workspace *workspace;
 double result,abserr;

 workspace = gsl_integration_workspace_alloc(WORKSIZE);
 F.function = &local_ng_f_q;
 F.params = &R;
  
 gsl_integration_qag(&F, 0, 1, 0, 1.0e-6, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  
 gsl_integration_workspace_free(workspace);

 result *= c;
 
 return result;
}                                   /* end local_ng_un_norm_sigma_square */


double local_ng_sigma_square(double z, double M)
{
  double R = pow(3.0 * M / 4.0 / Pi / mean_rho, 1.0 / 3.0);	
	
  double C_temp = 0.5 / Pi_square * F_NL * F_NL * norm_pow_A * norm_pow_A;
         C_temp = C_temp * pow(1.5 * Omega_m / 9e6, 2) * ps_local_factor(ns);
  double sigma_g_square = sigma_square(z, M);
  double sigma_ng_square = C_temp * local_ng_un_norm_sigma_square(z, R);
  
//  printf("C_temp/norm_pow_A = %.5e\n", C_temp/norm_pow_A);
  
  return sigma_g_square + sigma_ng_square;	  
}                                   /* end local_ng_sigma_square */

double ps_local_factor_ing_kernel(double t, void *params)
{
	double n = *(double *) params;
	double result, temp;
	
	if(t == 0)
	  result = 0;
	  
	if(t == 1)
	  result = -2.0;
	else
	  {
	  temp = pow(1+t, n-2) - pow(fabs(1-t), n-2);	
	  result = 1.0 / (n - 2.0) * pow(t, n-3) * temp - 2.0 * pow(t, n-2) - 1.0 / (n - 2.0) * t * temp;		
      }
    
    return result;
}                                    /* end ps_local_factor_ing_kernel */

double ps_local_factor(double n)
{
	int WORKSIZE = 1000000;
	gsl_function F;
	gsl_integration_workspace *workspace;
	double result,abserr;

	workspace = gsl_integration_workspace_alloc(WORKSIZE);
	F.function = &ps_local_factor_ing_kernel;
	F.params = &n;
 
 
	gsl_integration_qagiu(&F, 0, 0, 1.0e-7, WORKSIZE, workspace, &result, &abserr);
  
	gsl_integration_workspace_free(workspace);
 
	return result;		
}                                     /* end ps_local_factor */

#endif
