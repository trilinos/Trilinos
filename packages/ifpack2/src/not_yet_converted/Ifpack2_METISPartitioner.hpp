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

#ifndef IFPACK2_METISPARTITIONER_HPP
#define IFPACK2_METISPARTITIONER_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Partitioner.hpp"
#include "Ifpack2_OverlappingPartitioner.hpp"
#include "Teuchos_ParameterList.hpp"
class Ifpack2_Graph;
class Tpetra_Map;
class Tpetra_BlockMap;
class Tpetra_Import;

//! Ifpack2_METISPartitioner: A class to decompose Ifpack2_Graph's using METIS.
/*!
Class Ifpack2_METISPartitioner enables the decomposition of the local
Ifpack2_Graph's using METIS. In order to work properly, this class requires
TIFPACK to be configured with option \c --enable-ifpack-metis.
Otherwise, this class will always create one partition.

\date Last modified: Oct-04.
*/

class Ifpack2_METISPartitioner : public Ifpack2_OverlappingPartitioner {

public:

  //! Constructor.
  Ifpack2_METISPartitioner(const Ifpack2_Graph* Graph) :
    Ifpack2_OverlappingPartitioner(Graph),
    UseSymmetricGraph_(true)
  {}

  //! Destructor.
  virtual ~Ifpack2_METISPartitioner() {};

  //! Sets all the parameters for the partitioner (none at moment).
  int SetPartitionParameters(Teuchos::ParameterList& List)
  {
    // `true' is the safest option, as singletons may still be
    // in the system, even after Ifpack2_SingletonFilter (think
    // for example to an upper triangular matrix, which has one
    // singleton, and still has one after the elimination of 
    // the first one...)
    UseSymmetricGraph_ = List.get("partitioner: use symmetric graph", 
				  UseSymmetricGraph_);

    return(0);
  }

  //! Computes the partitions. Returns 0 if successful.
  int ComputePartitions();

private:
  bool UseSymmetricGraph_;

}; // class Ifpack2_METISPartitioner

#endif // IFPACK2_METISPARTITIONER_HPP
