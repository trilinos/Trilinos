#include <stdio.h>
#include <stdlib.h>
/* Program to read in rigid body information and to write  */
/* out a series of files suitable for ml_readex to process */
/*                                                         */
/* The form of the input file is                           */
/*     Nrigid_body_modes  vector_length                    */
/*     rigid1[0]                                           */
/*         .                                               */
/*         .                                               */
/*     rigid1[vector_length-1]                             */
/*     rigid2[0]                                           */
/*         .                                               */
/*         .                                               */
/*     rigid2[vector_length-1]                             */
/*         .                                               */
/*         .                                               */
/*     rigidNrigid_body_modes[vector_length-1]             */

int main(int argc, char *argv[]) {

int Nrigid, length, i, j;
char input_file[80], rigid_file[80];
FILE *fp_in, *fp_out;
double dvalue;

   system("rm -f rigid_body_mode*");
   printf("Enter file name where rigid body modes are kept\n");
   scanf("%s",&input_file);
   fp_in = fopen(input_file,"r");
   if (fp_in == NULL) {
      printf("could not open file: %s\n",input_file);
      exit(1);
   }
   if (fscanf(fp_in,"%d",&Nrigid)!= 1){printf("Num rigid not found\n");exit(1);}
   if (fscanf(fp_in,"%d",&length)!= 1){printf("length not found\n");exit(1);}
   if (Nrigid > 1000) {
      printf("too many rigid body modes: %d %d\n",Nrigid,length);
      exit(1);
   }

   for (i = 0; i < Nrigid; i++) {
      sprintf(rigid_file,"rigid_body_mode%d",i+1);
      fp_out = fopen(rigid_file,"w");
      if (fp_out == NULL) {
         printf("could not open file: %s\n",rigid_file);
         exit(1);
      }
      fprintf(fp_out,"%d\n",length);
      for (j = 0; j < length; j++) {
         if (fscanf(fp_in,"%lf",&dvalue) != 1) 
            printf("Data missing for %dth mode %dth component\n",i+1,j);
         fprintf(fp_out,"%d %20.13e\n-1\n",j,dvalue);
      }
      fclose(fp_out);
   }
}
