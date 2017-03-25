/* Internal unit:
 *  Mass unit: M_sun, but input M_sun /h
 *  Length unit : Mpc/h
 */
 
#define Z_NORM    100                                                  /* arbitrary normalization redshift */      
#define H0        0.7                                                  /* hubbule constant = 100H0 Km/s/Mpc */
#define OMEGA_M   0.3                                                  /* matter fraction today */
#define OMEGA_V   (1.0 - OMEGA_M)                                      /* flat universe here */
#define OMEGA_B   0.024 / H0 / H0                                      /* baryon fraction */
#define N_S       0.96                                                 /* index of primordial P.S. */
#define SIGMA_8   0.8                                                  /* sigma 8 */          

#ifdef WD_ST
#define M_WD      1                                                    /* mass of sterile neutrino, 
                                                                        * fitting valid for range [0.3,15] 
                                                                        * in unit of keV */
#endif

#ifdef LOCAL_NG
#define F_NL      0                                                    /* local type f_nl */
#endif

#define DELTA_C   1.686                                                /* critical overdensity by linear theory */
#define LAMBDA_C  0.37                                                 /* critical value of smallest lambda */

#define TOL_ERROR 1e-2                                                 /*  better not too small */

#define OUT_PUT   "/home/dalong/Projects/Code/Code-to-be-refined/Halo-shape-theory/Result/Gaussian-plot/"   //"./Result/WDM-ST-1kev/"


#define Pi        3.1415926535897932384626433832795028842              /* pi */
#define Pi_square 9.86960440108935861883449099987615113531             /* pi^2 */
