#include <stdlib.h>
#include <stdio.h>

void write_vec(const char *filename, int n_equations,double *x)
 
/*  write ASCII data file:
*/

{
  FILE *out_file ;
  int i;

 
  if ( (out_file = fopen( filename, "w")) == NULL ) {
    fprintf(stderr,"Error: Cannot open file: %s\n",filename);
    return;
  }
 
  for (i=0; i< n_equations; i++)
  fprintf(out_file, "%20.15e\n", x[i]) ;

  fclose(out_file);

/* end write_vec */
}
