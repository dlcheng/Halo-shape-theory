#ifndef VAR_H
 #include "var.h"
#endif

void init_global();
void init_input_format();
void init_constants();

double growth_factor_anyz(double z);
double g_square_factor(double z);
double growth_factor(double z);
double T_k(double k);
double f_k(double k, double R);
double f_q(double q, void *params);
double un_norm_sigma_square(double z, double R);
void init_power_norm();
double sigma_square(double z, double M);

#ifdef LOCAL_NG
double local_ng_f_k(double k, double R);
double local_ng_f_q(double q, void *params);
double local_ng_un_norm_sigma_square(double z, double R);
double local_ng_sigma_square(double z, double M);
double ps_local_factor_ing_kernel(double t, void *params);
double ps_local_factor(double n);
#endif

double un_norm_p_nu(double nu_1, double nu_2, double sigma_m, double d_f);
double un_norm_intg_kernel_over_nu_1(double nu_1, void * params);
double un_norm_intg_kernel_over_nu_2(double nu_2, void *params);
double un_norm_intg_over_nu_1(double nu_2, void * params);
double un_norm_intg_over_nu_2(double nu_1, void * params);
void init_prob_norm(double z_f, double M);
double p_nu_2_fixed_z(double nu_2, void *params);
double p_nu_1_fixed_z(double nu_1, void *params);

double bar_alpha_eff(double z, double M);
double sigma_diff_over_m(double M, void * params);
double form_prob(double M, double z_f, double z);
double growth_factor_diff_over_z(double z, void *params);

double p_a_over_c(double nu_2, double M, double z);
double combine_intg_a_c(double a_f, void * params);
double p_b_over_c(double nu_1, double M, double z);
double combine_intg_b_c(double a_f, void * params);

#ifdef EH_T_K
void init_eh_t_k();
double eh_t_k(double k_0);
#endif
