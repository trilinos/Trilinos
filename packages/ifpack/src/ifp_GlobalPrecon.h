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

#ifndef _IFP_GLOBALPRECON_H_
#define _IFP_GLOBALPRECON_H_

#include "Ifpack_config.h"

#include "ifp_Precon.h"
#include "ifp_LocalPrecon.h"

class IFPACK_DEPRECATED ifp_GlobalPrecon : public ifp_Precon
{
protected:
    ifp_LocalPrecon local_precon;

public:
    ifp_GlobalPrecon() {local_precon.name = (LocalPreconName) 0;}
    virtual ~ifp_GlobalPrecon() {}

    // set the method for solving or inverting the blocks

    void localprecon(LocalPreconName b) 
        {local_precon.name = b;}
    void localprecon(LocalPreconName b, int i1) 
        {local_precon.name = b; local_precon.iarg1 = i1;}
    void localprecon(LocalPreconName b, double d1) 
        {local_precon.name = b; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, int i1, double d1) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, double d1, int i1) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.darg1 = d1;}
    void localprecon(LocalPreconName b, int i1, int i2) 
      {local_precon.name = b; local_precon.iarg1 = i1; local_precon.iarg2 = i2;}
    void localprecon(LocalPreconName b, double d1, double d2) 
      {local_precon.name = b; local_precon.darg1 = d1; local_precon.darg2 = d2;}
};

#endif // _IFP_GLOBALPRECON_H_
