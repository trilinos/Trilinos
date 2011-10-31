/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include "f2c.h"

double    pow_dd(doublereal *ap, doublereal *bp)
{
    double    pow();

    return (pow(*ap, *bp));
}
