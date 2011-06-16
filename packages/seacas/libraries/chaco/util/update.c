/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* update - fills double vec1 with vec2 + alpha*vec3 over range*/
void      update(double *vec1, int beg, int end, double *vec2, double fac, double *vec3)
{
    int       i;

    vec1 += beg;
    vec2 += beg;
    vec3 += beg;
    for (i = end - beg + 1; i; i--) {
	(*vec1++) = (*vec2++) + fac * (*vec3++);
    }
}

/* update - fills float vec1 with vec2 + alpha*vec3 over range*/
void      update_float(float *vec1, int beg, int end, float *vec2, float fac, float *vec3)
{
    int       i;

    vec1 += beg;
    vec2 += beg;
    vec3 += beg;
    for (i = end - beg + 1; i; i--) {
	(*vec1++) = (*vec2++) + fac * (*vec3++);
    }
}
