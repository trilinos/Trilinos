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

#ifndef _IFP_BILUK_H_
#define _IFP_BILUK_H_

#include "Ifpack_config.h"

#include "ifp_GlobalPrecon.h"

class ifp_LocalMat;
class ifp_BlockMat;

class IFPACK_DEPRECATED ifp_biluk : public ifp_GlobalPrecon
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
