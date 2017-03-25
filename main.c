#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "var.h"
#include "proto.h"
#include "define.h"

int main()
{
   int i;
   double z, M;
   FILE * fp;
   char filename[300];
   double p_ac, p_bc;      
   double p_fixed_time;
   double sigma_1, sigma_2;
   double param_list[3];   
   int op;
   int op_file;
   
   init_global();
   
   printf("Enter the mass in unit of M_sun/h, M=");
   if(scanf("%lf", &M) == 0)
      printf("Invalid value of M.\n");
   
   printf("Enter the redshift, z=");
   if(scanf("%lf", &z) == 0)
      printf("Invalid value of z.\n");
      
   printf("Save the result to file (1/0)?");
   if(scanf("%d", &op_file) == 0)
      printf("Invalid value of op_file.\n");
      
      
   printf("What to output?\n");
   printf("1. Gaussian case p(a/c) and p(b/c).\n");
   printf("2. Local type non-gaussian p(a/c) at formation time.\n");
   printf("3. Local type non-gaussian sigma at z_i=%.1f, eleven mass bins compares with the gaussian case.\n", z_i);
   printf("Your choice:");
   if(scanf("%d", &op) == 0)
      printf("Invalid value of op.\n");
      
   if(op_file != 0)
     {        
#ifdef LEE   
     sprintf(filename,"%s%s%.3e%s%.1f%s%.1f%s%d%s", 
             OUT_PUT, "Lee_M_", M, "_z_", z, "_sigma8_", SIGMA_8, "_op_", op, ".txt");
     fp = fopen(filename, "w+");
#endif   

#ifdef DL   
     sprintf(filename,"%s%s%.3e%s%.1f%s%.1f%s%d%s", 
             OUT_PUT, "Dl_M_", M, "_z_", z, "_sigma8_", SIGMA_8, "_op_", op, ".txt");
     fp = fopen(filename, "w+");
#endif  
     }
     
     
   M = M / H0;                 /* internal mass unit is M_sun */
   
   param_list[0] = z;
   param_list[1] = M;
   param_list[2] = 1.0;

   if(op == 1)
    {   
     for(i=1; i<100; i++)
       {
	   p_ac = p_a_over_c(0.01 * i, M, z);
       p_bc = p_b_over_c(0.01 * i, M, z);
       if(op_file != 0)
        fprintf(fp, "%.5e	%.5e	%.5e\n", 0.01 * i, p_ac, p_bc);
       printf("%.5e	%.5e	%.5e\n", 0.01 * i, p_ac, p_bc); 
       } 
    }
    
    if(op == 2)
      {
	   for(i=1; i<100; i++)	  
		 {
 	      p_fixed_time = p_nu_2_fixed_z(0.01 * i, param_list);  		  
          printf("%.5e	%.5e\n", 0.01 * i, p_fixed_time);
          if(op_file !=0)			 
		   fprintf(fp, "%.5e	%.5e\n", 0.01 * i, p_fixed_time);
		 }
	  }
    
    if(op == 3)
       {
	    for(i=0; i<11; i++)
	      {
		   sigma_1 = sqrt(sigma_square(z, M/pow(2, i-5)));
		   sigma_2 = sqrt(local_ng_sigma_square(z, M/pow(2, i-5)));	
		   printf("%.5e	%.7e	%.7e	%.7e\n", M * pow(2, i-5) * H0, sigma_2, sigma_1, fabs(sigma_2/sigma_1 - 1));		   
		   if(op_file != 0)
		     printf("%.5e	%.5e	%.5e\n", M * pow(2, i-5) * H0, sigma_2, sigma_1);	   
	      }
	   }
    
   if(op_file != 0)     
     fclose(fp);
      
   return 0;
}  /* end main */
