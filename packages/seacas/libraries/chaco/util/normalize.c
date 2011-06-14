/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

/* Normalizes a double n-vector over range. */
double    normalize(double *vec, int beg, int end)
{
    int       i;
    double    scale;
    double    norm(double *vec, int beg, int end);

    scale = norm(vec, beg, end);
    vec = vec + beg;
    for (i = end - beg + 1; i; i--) {
	*vec = *vec / scale;
	vec++;
    }
    return (scale);
}

/* Normalizes such that element k is positive */
double    sign_normalize(double *vec, int beg, int end, int k)
{
    int       i;
    double    scale, scale2;
    double    norm(double *vec, int beg, int end);

    scale = norm(vec, beg, end);
    if (vec[k] < 0) {
	scale2 = -scale;
    }
    else {
	scale2 = scale;
    }
    vec = vec + beg;
    for (i = end - beg + 1; i; i--) {
	*vec = *vec / scale2;
	vec++;
    }
    return (scale);
}

/* Normalizes a float n-vector over range. */
double    normalize_float(float *vec, int beg, int end)
{
    int       i;
    float     scale;
    double    norm_float(float *vec, int beg, int end);

    scale = norm_float(vec, beg, end);
    vec = vec + beg;
    for (i = end - beg + 1; i; i--) {
	*vec = *vec / scale;
	vec++;
    }
    return ((double) scale);
}
