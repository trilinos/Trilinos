/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef IFPACK2_GREEDYPARTITIONER_HPP
#define IFPACK2_GREEDYPARTITIONER_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Partitioner.hpp"
#include "Ifpack2_OverlappingPartitioner.hpp"
#include "Teuchos_ParameterList.hpp"
class Tpetra_Comm;
class Ifpack2_Graph;
class Tpetra_Map;
class Tpetra_BlockMap;
class Tpetra_Import;

//! Ifpack2_GreedyPartitioner: A class to decompose Ifpack2_Graph's using a simple greedy algorithm.

class Ifpack2_GreedyPartitioner : public Ifpack2_OverlappingPartitioner {

public:

  //! Constructor.
  Ifpack2_GreedyPartitioner(const Ifpack2_Graph* Graph) :
    Ifpack2_OverlappingPartitioner(Graph),
    RootNode_(0)
  {}

  //! Destructor.
  virtual ~Ifpack2_GreedyPartitioner() {};

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

}; // class Ifpack2_GreedyPartitioner

#endif // IFPACK2_GREEDYPARTITIONER_HPP
