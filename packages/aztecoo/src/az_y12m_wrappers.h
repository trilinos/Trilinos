#ifndef _AZ_Y12M_WRAPPERS_H_
#define _AZ_Y12M_WRAPPERS_H_

#include "az_f77func.h"

#   define Y12MBF_F77                 F77_FUNC(y12mbf,Y12MBF)
#   define Y12MCF_F77                 F77_FUNC(y12mcf,Y12MCF)
#   define Y12MDF_F77                 F77_FUNC(y12mdf,Y12MDF)

#ifdef __cplusplus
extern "C" {
#include <stdio.h>
#endif

void Y12MBF_F77(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, int ha[], int *iha, double aflag[],
                      int iflag[], int *ifail);

void Y12MCF_F77(int *n, int *z, double val[], int snr[], int *nn,
                      int rnr[], int *nn1, double pivot[], double b[], int ha[],
                      int *iha, double aflag[], int iflag[], int *ifail);

void Y12MDF_F77(int *n, double val[], int *nn, double b[], double pivot[],
                      int snr[], int ha[], int *iha, int iflag[], int *ifail);

#ifdef __cplusplus
}
#endif

#endif /* _AZ_Y12M_WRAPPERS_H_ */
