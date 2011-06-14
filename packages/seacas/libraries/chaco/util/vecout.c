/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	<math.h>
#include	"defs.h"

/* Print vertically range of double vector. */
void      vecout(double *vec, int beg, int end, char *tag, char *file_name)
{
    FILE     *file;
    int       i;
    int       print_indices;

    print_indices = FALSE;

    if (file_name != NULL)
	file = fopen(file_name, "w");
    else
	file = stdout;


    if (print_indices) {
	fprintf(file, "%s:\n", tag);
	for (i = beg; i <= end; i++) {
	    if (fabs(vec[i]) >= 1.0e-16)
		fprintf(file, "%2d.   %24.16f\n", i, vec[i]);
	    else
		fprintf(file, "%2d.         %g \n", i, vec[i]);
	}
    }

    else {
	fprintf(file, "%s:\n", tag);
	for (i = beg; i <= end; i++) {
	    if (fabs(vec[i]) >= 1.0e-16)
		fprintf(file, "%2d.   %24.16f\n", i, vec[i]);
	    else
		fprintf(file, "%2d.         %g \n", i, vec[i]);
	}
    }

    if (file_name != NULL)
	fclose(file);
}

/* Scale the eigenvector such that the first element is non-negative.
 * This is an attempt to reduce machine-to-machine variance which can
 * result from calculated eigenvectors being mirror-images of each
 * other due to small roundoff..
 */
void vecnorm(double *vec, int beg, int end)
{
  int i;
  if (vec[beg] < 0.0) {
    for (i = beg; i <= end; i++) {
      vec[i] *= -1.0;
    }
  }
}
