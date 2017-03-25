#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "var.h"
#include "proto.h"
#include "define.h"

void init_global()
{
  init_input_format();
  
#ifdef EH_T_K
  init_eh_t_k();  
#endif 
 
#ifdef WD_ST 
  WD_alpha = 0.189 * pow(M_WD, -0.858) * pow(Omega_m/0.26, -0.136) * pow(H0/0.7, 0.692);
#endif  
  
  init_constants();
  
  init_power_norm();
	
}                    /* end int_global */

void init_input_format()
{
  z_i = Z_NORM;
  Omega_m = OMEGA_M;
  Omega_b = OMEGA_B;
  Omega_v = OMEGA_V;
  ns = N_S;
  h0 = H0;
  Sigma_8 = SIGMA_8;
	
}                    /* end int_input_format */

void init_constants()
{
  mean_rho = 2.78e11 * Omega_m / h0;   /* in unit of h0^3*M_sun/Mpc^3 for flat universe */
		
}                    /* end int_constants */
