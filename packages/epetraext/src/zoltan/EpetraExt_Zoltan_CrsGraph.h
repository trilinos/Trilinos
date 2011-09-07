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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_ZOLTAN_CRSGRAPH_H
#define EPETRAEXT_ZOLTAN_CRSGRAPH_H

#include <EpetraExt_Transform.h>

#include <string>

class Epetra_CrsGraph;
class Epetra_Map;

class Zoltan_LoadBalance;

namespace EpetraExt {

///
/** Generates an Epetra_CrsGraph based on the repartitioning algorithms of Zoltan
 */
class EPETRAEXT_DEPRECATED Zoltan_CrsGraph : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  const std::string partitionMethod_;

  Epetra_Map * NewRowMap_;

 public:

  ///
  /* Destructor
   */
  ~Zoltan_CrsGraph();

  ///
  /* Constructor
   * input param part_method - type of Zoltan partitioning to use
   */
  Zoltan_CrsGraph( const std::string & part_method = std::string("PartKway") )
  : partitionMethod_(part_method),
    NewRowMap_(0)
  {}

  ///
  /* Generates the Zoltan partitioned Epetra_CrsGraph from the input object.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EPETRAEXT_ZOLTAN_CRSGRAPH_H
