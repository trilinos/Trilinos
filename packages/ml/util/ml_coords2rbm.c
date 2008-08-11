#include <stdlib.h>
#include <stdio.h>
#include "ml_rbm.h"
/* Program to read in coordinate information and the      */
/* create rigid body information suitable for ml_rigid_spread */
/*                                                         */

int main(int argc, char *argv[]) {

int Nnodes, Nelements, Nmaterials, Ndim,Ndof, maxpernode, Nrigid = 6;
int i, j, k;
char input_file[80], output_file[80], garbage[80];
FILE *fp_in, *fp_out;
double x,y,z, *rbm;

   rbm = (double *) malloc(36*sizeof(double));
   printf("Enter file name where coordinates are kept\n");
   scanf("%s",input_file);
   printf("Enter file name for rigid body mode file\n");
   scanf("%s",output_file);
   fp_out = fopen(output_file,"w");
   if (fp_out == NULL) {
      printf("could not open file: %s\n",output_file);
      exit(1);
   }
   i = 0;
  
   while ( i < Nrigid) {
      fp_in = fopen(input_file,"r");
      if (fp_in == NULL) {
         printf("could not open file: %s\n",input_file);
         exit(1);
      }
      if (fscanf(fp_in,"%d",&Nnodes)!= 1)
         {printf("Num nodes not found\n");exit(1);}
      if (fscanf(fp_in,"%d",&Nelements)!= 1)
         {printf("Num elements not found\n");exit(1);}
      if (fscanf(fp_in,"%d",&Nmaterials)!= 1)
         {printf("Num material not found\n");exit(1);}
      if (fscanf(fp_in,"%d",&Ndim)!= 1)
         {printf("Dimension not found\n");exit(1);}
      if (fscanf(fp_in,"%d",&Ndof)!= 1)
         {printf("Ndof not found\n");exit(1);}
      if (Ndof == 1) Nrigid = 1;
      if (fscanf(fp_in,"%d",&maxpernode)!= 1)
         {printf("maxpernode not found\n");exit(1);}
      if (fscanf(fp_in,"%s",garbage)!= 1)
         {printf("nopr not found\n");exit(1);}
      if (fscanf(fp_in,"%s",garbage)!= 1)
         {printf("nopa not found\n");exit(1);}
      if (fscanf(fp_in,"%s",garbage)!= 1)
         {printf("noco not found\n");exit(1);}
      if (fscanf(fp_in,"%s",garbage)!= 1)
         {printf("coordinates not found\n");exit(1);}

      if (i == 0) fprintf(fp_out,"%d %d\n",6, Nnodes*Ndof);

      for (j = 0; j < Nnodes; j++) {
         if (fscanf(fp_in,"%d",&k) != 1) 
            printf("Data missing before %dth x coordinate\n",j);
         if (fscanf(fp_in,"%d",&k) != 1) 
            printf("Data Missing before %dth x coordinate\n",j);
         if (fscanf(fp_in,"%lf",&x) != 1) 
            printf("Data missing for %dth x coordinate\n",j);
         if (fscanf(fp_in,"%lf",&y) != 1) 
            printf("Data missing for %dth y coordinate\n",j);
         if (fscanf(fp_in,"%lf",&z) != 1) 
            printf("Data missing for %dth z coordinate\n",j);

         ML_Coord2RBM(1, &x, &y, &z, rbm, Ndof);
         for (k = 0; k < Ndof ; k++)
            fprintf(fp_out,"%e\n",rbm[k+Ndof*i]);
      }
      i++;
      fclose(fp_in);
   }
   fclose(fp_out);
   free(rbm);
   return(1);
}
