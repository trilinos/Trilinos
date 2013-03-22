// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
  TEUCHOS_TEST_FOR_EXCEPT( get_reduct_op() != null && reduct_op != null  );
#endif
  if(!the_mpi_op.get()) {
    MPI_Op mpi_op = MPI_OP_NULL;
    TEUCHOS_TEST_FOR_EXCEPT(
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
  TEUCHOS_TEST_FOR_EXCEPT(!reduct_op.get())
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
