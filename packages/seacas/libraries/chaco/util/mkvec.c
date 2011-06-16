/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include "smalloc.h"

/* Allocates a double vector with range [nl..nh]. Dies. */
double   *mkvec(int nl, int nh)
{
    double   *v;

    v = smalloc((nh - nl + 1) * sizeof(double));
    return (v - nl);
}

/* Allocates a double vector with range [nl..nh]. Returns error code. */
double   *mkvec_ret(int nl, int nh)
{
    double   *v;

    v = smalloc_ret((nh - nl + 1) * sizeof(double));
    if (v == NULL) 
	return(NULL);
    else 
        return (v - nl);
}

/* Free a double vector with range [nl..nh]. */
void      frvec(double *v, int nl)
{


    sfree((v + nl));
    v = NULL;
}

/* Allocates a float vector with range [nl..nh]. Dies. */
float   *mkvec_float(int nl, int nh)
{
    float   *v;

    v = smalloc((nh - nl + 1) * sizeof(float));
    return (v - nl);
}

/* Allocates a float vector with range [nl..nh]. Returns error code. */
float   *mkvec_ret_float(int nl, int nh)
{
    float   *v;

    v = smalloc_ret((nh - nl + 1) * sizeof(float));
    if (v == NULL) 
	return(NULL);
    else 
        return (v - nl);
}

/* Free a float vector with range [nl..nh]. */
void      frvec_float(float *v, int nl)
{


    sfree((v + nl));
    v = NULL;
}
