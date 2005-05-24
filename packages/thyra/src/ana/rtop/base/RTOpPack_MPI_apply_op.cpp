// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#include "RTOpPack_MPI_apply_op.hpp"

extern "C" {

void RTOpPack_MPI_apply_op_reduction_op(
  void              *invec
  ,void             *inoutvec
  ,int              *len
  ,RTOp_Datatype    *datatype
  )
{
  RTOpPack::get_reduct_op()->reduce_reduct_objs(invec,inoutvec,len,datatype);
}

} // extern "C"

//
// This implementation of handling the reduction operator
// will work just fine in a single threaded program.
// For multi-threaded program this has to be reworked!
//

namespace {
Teuchos::RefCountPtr<const RTOpPack::MPIReductionOpBase> the_reduct_op = Teuchos::null;
}

void RTOpPack::set_reduct_op( const Teuchos::RefCountPtr<const MPIReductionOpBase>& reduct_op )
{
  the_reduct_op = reduct_op;
}

Teuchos::RefCountPtr<const RTOpPack::MPIReductionOpBase> RTOpPack::get_reduct_op()
{
  return the_reduct_op;
}
