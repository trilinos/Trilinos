#ifndef _IFP_BILUK_H_
#define _IFP_BILUK_H_

#include "ifp_GlobalPrecon.h"

class ifp_LocalMat;
class ifp_BlockMat;

class ifp_biluk : public ifp_GlobalPrecon
{
private:
    const ifp_BlockMat *Ap;

    ifp_LocalMat **diag;  // inverse or factors of diagonal blocks

    ifp_LocalMat **al;    // lower triangular factor (strict lower part stored)
    int      *jal;
    int      *ial;
    ifp_LocalMat **au;    // upper triangular factor
    int      *jau;
    int      *iau;

    int NumEntries_;
    int NumNonzeros_;

public:
    ifp_biluk();
   ~ifp_biluk();

   int NumEntries() {return(NumEntries_);};
   int NumNonzeros() {return(NumNonzeros_);};

   ifp_BlockMat * A() {return((ifp_BlockMat *)Ap);};
 
    void setup(const ifp_BlockMat& A, int levfill);
    void apply (int, int, const double *, int, double *, int);
    void applyr(int, int, const double *, int, double *, int);
    void applyl(int, int, const double *, int, double *, int);

    void multiply(int, int, const double *, int, double *, int);
    double condest();

    static int growth;
};

#endif // _IFP_BILUK_H_
