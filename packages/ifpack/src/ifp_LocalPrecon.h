#ifndef _IFP_LOCALPRECON_H_
#define _IFP_LOCALPRECON_H_

// Make the values explicit, so we can match them with FORTRAN parameters.

enum LocalPreconName
{                          // arguments
    LP_LU           =  1,  //
    LP_INVERSE      =  2,  //
    LP_SVD          =  3,  // rthresh, athresh
    LP_DIAG         = 10,  //
    LP_SOR          = 12,  // omega, iterations
    LP_SSOR         = 13,  // omega, iterations
    LP_DIAGDOM      = 15,  // 
    LP_GERSH        = 16   // alpha
};

class ifp_LocalPrecon
{
public:
    LocalPreconName name;
    int iarg1;
    int iarg2;
    double darg1;
    double darg2;
};

#endif // _IFP_LOCALPRECON_H_
