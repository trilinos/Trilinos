#ifndef _IFP_ARCH_H_
#define _IFP_ARCH_H_

#define ifp_name2(a,b) a ## b

#if defined(CRAY) || defined(matched) || \
    defined(TRILINOS_HAVE_NO_FORTRAN_UNDERSCORE)
#define F77NAME(x) x
#else
#define F77NAME(x) ifp_name2(x,_)
#endif

typedef int      integer;
typedef long int integer8;    /* should be same size as a pointer */
typedef double   doublereal;

/* Note that C and Fortran integer sizes should match, since the storage
 * for Fortran arrays (for Harwell-Boeing matrices) is allocated by C. 
 */

#endif // _IFP_ARCH_H_
