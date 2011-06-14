/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Scale - fills vec1 with alpha*vec2 over range, double version */
void      vecscale(double *vec1, int beg, int end, double alpha, double *vec2)
{
    int       i;

    vec1 += beg;
    vec2 += beg;
    for (i = end - beg + 1; i; i--) {
	(*vec1++) = alpha * (*vec2++);
    }
}

/* Scale - fills vec1 with alpha*vec2 over range, float version */
void      vecscale_float(float *vec1, int beg, int end, float alpha, float *vec2)
{
    int       i;

    vec1 += beg;
    vec2 += beg;
    for (i = end - beg + 1; i; i--) {
	(*vec1++) = alpha * (*vec2++);
    }
}
