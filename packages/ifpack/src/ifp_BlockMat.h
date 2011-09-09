/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef _IFP_BLOCKMAT_H_
#define _IFP_BLOCKMAT_H_

#include "Ifpack_config.h"

#include <iostream>
#include "ifp_Matrix.h"

class ifp_LocalMat;

class IFPACK_DEPRECATED ifp_BlockMat : public ifp_Matrix
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

std::ostream& operator << (std::ostream& os, const ifp_BlockMat& mat);

#endif // _IFP_BLOCKMAT_H_
