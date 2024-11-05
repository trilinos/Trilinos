/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK_GREEDYPARTITIONER_H
#define IFPACK_GREEDYPARTITIONER_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_OverlappingPartitioner.h"
#include "Teuchos_ParameterList.hpp"
class Epetra_Comm;
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Ifpack_GreedyPartitioner: A class to decompose Ifpack_Graph's using a simple greedy algorithm.

class Ifpack_GreedyPartitioner : public Ifpack_OverlappingPartitioner {

public:

  //! Constructor.
  Ifpack_GreedyPartitioner(const Ifpack_Graph* Graph) :
    Ifpack_OverlappingPartitioner(Graph),
    RootNode_(0)
  {}

  //! Destructor.
  virtual ~Ifpack_GreedyPartitioner() {};

  //! Sets all the parameters for the partitioner (root node).
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    RootNode_ = List.get("partitioner: root node", RootNode_);

    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

private:

  int RootNode_;

}; // class Ifpack_GreedyPartitioner

#endif // IFPACK_GREEDYPARTITIONER_H
