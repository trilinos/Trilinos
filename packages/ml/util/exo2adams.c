/*
    I don't get this ....
*/
#include <stdio.h>
extern int read_up_to(char *pattern, int length);
main() {
   int index, Nnodes, i, status, Ndof = 3, j,k;
   double *x,*y,*z;
   char   comma;
   

#ifndef adams
   read_up_to("num_nodes =", 11);
#else
   read_up_to("Kasper", 6);
#endif
   
   if ( scanf("%d",&Nnodes) == 0) {
      printf("Nodes not found\n");
      exit(1);
   }
   x = (double *) malloc(sizeof(double)*Nnodes);
   y = (double *) malloc(sizeof(double)*Nnodes);
   z = (double *) malloc(sizeof(double)*Nnodes);
   if (z == NULL) {
      printf("not enough space for coordinates\n");
      exit(1);
   }

#ifndef adams
   read_up_to("coord =", 7);
#else
   read_up_to("coordinates", 11);
#endif

   /* this is adams header */

   printf("%d    0   0 3    %d   0\n",Nnodes, Ndof);
   printf("\nnopr\nnopa\nnoco\n\ncoordinates\n");

   for (i = 0; i < Nnodes; i++) {
#ifndef adams
      if (i != 0) { while ( (comma = getchar()) != ',') ; }
#else
      scanf("%d",&j);
      scanf("%d",&k);
#endif
      if ( scanf("%lf",&(x[i])) == 0) {
         printf("x(%d) not found\n",i); 
         exit(1);
      }
#ifndef adams
   }
   for (i = 0; i < Nnodes; i++) {
      while ( (comma = getchar()) != ',') ; 
#endif
      if ( scanf("%lf",&(y[i])) == 0) {
         printf("y(%d) not found\n",i); 
         exit(1);
      }
#ifndef adams
   }
   for (i = 0; i < Nnodes; i++) {
      while ( (comma = getchar()) != ',') ; 
#endif
      if ( scanf("%lf",&(z[i])) == 0) {
         printf("z(%d) not found\n",i); 
         exit(1);
      }
/*
      printf("%d 0 %20.13e %20.13e %20.13e\n",i+1,x,y,z);
*/
   }
   for (i = 0; i < Nnodes; i++) {
      printf("%d 0 %20.13e %20.13e %20.13e\n",i+1,x[i],y[i],z[i]);
   }   
}

int read_up_to(char *pattern, int length) {
   int index = 0;
   char current;

   while ( index < length ) {
      current = getchar();
      if (current == EOF) {
         printf("End of file reached\n");
         exit(1);
      }
      if (current == pattern[index]) index++;
      else index = 0;
   }
   return(1);
}
