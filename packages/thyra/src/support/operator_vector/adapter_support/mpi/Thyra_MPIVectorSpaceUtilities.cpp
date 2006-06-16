// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#include "Thyra_MPIVectorSpaceUtilities.hpp"
#include "Teuchos_RawMPITraits.hpp"

namespace Thyra {

Index MPIVectorSpaceUtilities::computeMapCode(
  MPI_Comm mpiComm, const Index localSubDim
  )
{
  //
  // Here we will make a map code out of just the local
  // sub-dimension on each processor.  If each processor
  // has the same number of local elements, then the maps
  // will be the same and this is all you need for
  // RTOp compatibility unless the operations are not
  // coordinate invariant.  I will work on this issue
  // if it becomes a problem.
  //
  typedef Teuchos::RawMPITraits<Index> IRMT;
  int procRank = -1;
  MPI_Comm_rank(mpiComm,&procRank);
  Index mapCode = -1;
  Index localCode = localSubDim % (procRank+1) + localSubDim;
  MPI_Allreduce(
    &localCode                  // sendbuf
    ,&mapCode                   // recvbuf
    ,IRMT::adjustCount(1)       // count
    ,IRMT::type()               // datatype
    ,IRMT::sumOp()              // op
    ,mpiComm                    // comm
    );
  return mapCode;
}

Index MPIVectorSpaceUtilities::computeLocalOffset(
  MPI_Comm mpiComm, const Index localSubDim
  )
{
  typedef Teuchos::RawMPITraits<Index> IRMT;
  Index localOffset;
  Index _localOffset = localSubDim;
  MPI_Scan(
    &_localOffset               // sendbuf
    ,&localOffset               // recvbuf
    ,IRMT::adjustCount(1)       // count
    ,IRMT::type()               // datatype
    ,IRMT::sumOp()              // op
    ,mpiComm                    // comm
    );
  localOffset -= localSubDim;
  return localOffset;
}

Index MPIVectorSpaceUtilities::computeGlobalDim(
  MPI_Comm mpiComm, const Index localSubDim
  )
{
  typedef Teuchos::RawMPITraits<Index> IRMT;
  Index globalDim = -1;
  MPI_Allreduce(
    (void*)&localSubDim         // sendbuf
    ,&globalDim                 // recvbuf
    ,IRMT::adjustCount(1)       // count
    ,IRMT::type()               // datatype
    ,IRMT::sumOp()              // op
    ,mpiComm                    // comm
    );
  return globalDim;
}

void MPIVectorSpaceUtilities::broadcast(
  MPI_Comm mpiComm, const int rootRank, Index* value
  )
{
  typedef Teuchos::RawMPITraits<Index> IRMT;
  MPI_Bcast(
    (void*)value                // buffer [in/out]
    ,IRMT::adjustCount(1)       // count
    ,IRMT::type()               // datatype
    ,rootRank                   // root
    ,mpiComm                    // comm
    );
}

} // namespace Thyra
