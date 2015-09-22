//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_ISORROPIA_CRSGRAPH_H
#define EPETRAEXT_ISORROPIA_CRSGRAPH_H

#include <EpetraExt_Transform.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Generates an Epetra_CrsGraph based on the repartitioning algorithms of Isorropia
 */
class Isorropia_CrsGraph : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  // Parameter list to specify partitioning algorithm for Isorropia.
  Teuchos::ParameterList partitionList_;

  // New graph for the repartitioned matrix.
  Teuchos::RCP<Epetra_CrsGraph> NewGraph_;

 public:

  ///
  /* Destructor
   */
  ~Isorropia_CrsGraph();

  ///
  /* Constructor
   * input param part_method - type of Isorropia partitioning to use
   */
  Isorropia_CrsGraph( Teuchos::ParameterList& partitionList )
  : partitionList_(partitionList)
  {}

  ///
  /* Generates the Isorropia partitioned Epetra_CrsGraph from the input object.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EPETRAEXT_ISORROPIA_CRSGRAPH_H
