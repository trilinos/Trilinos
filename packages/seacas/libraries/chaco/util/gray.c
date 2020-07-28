/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/* Compute the binary reflected Gray code of a value. */
int gray(int i) { return ((i >> 1) ^ i); }

/* Compute the inverse of the binary reflected Gray code of a value. */
/*
int       invgray(i)
int       i;
{
    int       k;

    k = i;
    while (k) {
        k >>= 1;
        i ^= k;
    }
    return (i);
}
*/
