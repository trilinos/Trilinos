#ifndef _IFP_BLOCKMAT_H_
#define _IFP_BLOCKMAT_H_

#include <iostream.h>
#include "ifp_Matrix.h"

class ifp_LocalMat;

class ifp_BlockMat : public ifp_Matrix
{
private:
    ifp_LocalMat **a;
    int *ja;
    int *ia;

    int nrow;
    int ncol;
    int nnz;
    int nnzs;

    int *kvstr;
    int *kvstc;

public:
    ifp_BlockMat(double *, int *, int *, int *, int *,
		   int *, int, int, int, int);
   ~ifp_BlockMat();

    const ifp_LocalMat& val(unsigned int i) const {return *a[i];}
    ifp_LocalMat& val(unsigned int i) {return *a[i];}
    const int& row_ptr(unsigned int i) const {return ia[i];}
    const int& col_ind(unsigned int i) const {return ja[i];}

    int numrow() const {return nrow;}
    int numcol() const {return ncol;}
    int numnz()  const {return nnz;}
    int NumEntries()  const {return nnz;}
    int NumNonzeros()  const {return nnzs;}

    int dimrow() const {return kvst_row(nrow);} // scalar dimension
    int dimcol() const {return kvst_col(ncol);}

    const int& kvst_row(unsigned int i) const {return kvstr[i];}
    const int& kvst_col(unsigned int i) const {return kvstc[i];}

    void mult(int, int, const double *, int, double *, int) const;
    void trans_mult(int, int, const double *, int, double *, int) const;
};

ostream& operator << (ostream& os, const ifp_BlockMat& mat);

#endif // _IFP_BLOCKMAT_H_
