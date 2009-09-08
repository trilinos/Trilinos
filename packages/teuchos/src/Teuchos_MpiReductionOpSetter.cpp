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
#include "Teuchos_GlobalMPISession.hpp"

#ifdef MPIAPI
#define CALL_API MPIAPI
#else
#define CALL_API 
#endif

//
// This implementation of handling the reduction operator
// will work just fine in a single threaded program.
// For multi-threaded program this has to be reworked!
//

//
// The callback that an MPI implementation will actually call
//

extern "C" {
void CALL_API Teuchos_MPI_reduction_op(
  void              *invec
  ,void             *inoutvec
  ,int              *len
  ,MPI_Datatype     *datatype
  );

void Teuchos_MPI_Op_free( MPI_Op *op )
{
  //if(Teuchos::GlobalMPISession::mpiIsInitialized())
  //  MPI_Op_free(op);
  //else
  //  *op = MPI_OP_NULL;
  //
  // RAB: I have commented this out because this is getting called after
  // MPI_Finalize() is called when the Teuchos::GlobalMPISession class is not
  // being used..  On some systems, like MPICH, this was not problem.
  // However, there are some systems that complain when you do this.
  // Therefore, since I don't really know how to fix this problem, I am just
  // going to punt and just not delete this MPI_Op object.  I suspect that
  // many people do not clean up their MPI objects correctly so I would guess
  // that almost every MPI implementation allows you to not free objects and
  // end just fine.
}

} // extern "C"

//
// Manage the reduction operation as static data.  I have used access
// functions here to allow more insulation for other implementations other
// other than just single threaded programs.
//

namespace {

Teuchos::RCP<const Teuchos::MpiReductionOpBase> the_reduct_op = Teuchos::null;

Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Op> > the_mpi_op = Teuchos::null;

Teuchos::RCP<const Teuchos::MpiReductionOpBase> get_reduct_op()
{
  return the_reduct_op;
}

void set_reduct_op( const Teuchos::RCP<const Teuchos::MpiReductionOpBase>& reduct_op )
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
    the_mpi_op = Teuchos::opaqueWrapper(mpi_op,Teuchos_MPI_Op_free);
  }
  the_reduct_op = reduct_op;
}

} // namespace

extern "C" {
void CALL_API Teuchos_MPI_reduction_op(
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
  const Teuchos::RCP<const MpiReductionOpBase>& reduct_op
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
