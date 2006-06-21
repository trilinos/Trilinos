// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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

#include "Teuchos_MpiReductionOpSetter.hpp"
#include "Teuchos_OpaqueWrapper.hpp"

//
// This implementation of handling the reduction operator
// will work just fine in a single threaded program.
// For multi-threaded program this has to be reworked!
//

//
// The callback that an MPI implementation will actually call
//

extern "C" {

void Teuchos_MPI_reduction_op(
  void              *invec
  ,void             *inoutvec
  ,int              *len
  ,MPI_Datatype     *datatype
  );

} // extern "C"

//
// Manage the reduction operation as static data.  I have used access
// functions here to allow more insulation for other implementations other
// other than just single threaded programs.
//

namespace {

Teuchos::RefCountPtr<const Teuchos::MpiReductionOpBase> the_reduct_op = Teuchos::null;

Teuchos::RefCountPtr<const Teuchos::OpaqueWrapper<MPI_Op> > the_mpi_op = Teuchos::null;

Teuchos::RefCountPtr<const Teuchos::MpiReductionOpBase> get_reduct_op()
{
  return the_reduct_op;
}

void set_reduct_op( const Teuchos::RefCountPtr<const Teuchos::MpiReductionOpBase>& reduct_op )
{
  using Teuchos::null;
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( get_reduct_op() != null && reduct_op != null  );
#endif
  if(!the_mpi_op.get()) {
    MPI_Op mpi_op = MPI_OP_NULL;
    TEST_FOR_EXCEPT(
      0!=MPI_Op_create( &Teuchos_MPI_reduction_op ,1 ,&mpi_op ) // Assume op is commutative?
      );
    the_mpi_op = Teuchos::opaqueWrapper(mpi_op,MPI_Op_free);
  }
  the_reduct_op = reduct_op;
}

} // namespace

extern "C" {

void Teuchos_MPI_reduction_op(
  void              *invec
  ,void             *inoutvec
  ,int              *len
  ,MPI_Datatype     *datatype
  )
{
  get_reduct_op()->reduce(invec,inoutvec,len,datatype);
}

} // extern "C"

namespace Teuchos {

MpiReductionOpSetter::MpiReductionOpSetter(
  const Teuchos::RefCountPtr<const MpiReductionOpBase>& reduct_op
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(!reduct_op.get())
#endif
  set_reduct_op(reduct_op);
}

MpiReductionOpSetter::~MpiReductionOpSetter()
{
  set_reduct_op( null );
}

MPI_Op MpiReductionOpSetter::mpi_op() const
{
  return (*the_mpi_op)();
}

} // namespace Teuchos
