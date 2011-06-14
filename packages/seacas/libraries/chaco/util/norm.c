/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <math.h>

/* Returns 2-norm of a double n-vector over range. */
double    norm(double *vec, int beg, int end)
{
    double    dot(double *vec1, int beg, int end, double *vec2);

    return (sqrt(dot(vec, beg, end, vec)));
}

/* Returns 2-norm of a float n-vector over range. */
double    norm_float(float *vec, int beg, int end)
{
    double    temp;
    double    dot_float(float *vec1, int beg, int end, float *vec2);

    temp = sqrt(dot_float(vec, beg, end, vec));
    return (temp);
}
