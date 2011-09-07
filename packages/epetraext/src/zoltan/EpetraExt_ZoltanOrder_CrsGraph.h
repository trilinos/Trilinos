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

#ifndef EPETRAEXT_ZOLTANORDER_CRSGRAPH_H
#define EPETRAEXT_ZOLTANORDER_CRSGRAPH_H

#ifdef ZOLTAN_ORDER

#include <EpetraExt_Transform.h>

namespace Zoltan {

class LoadBalance;

}

class Epetra_Map;
class Epetra_CrsGraph;

namespace EpetraExt {

///
/** Generates a local reordered Epetra_CrsGraph based on the Metis local ordering
 * from Zoltan
 */

class EPETRAEXT_DEPRECATED ZoltanOrder_CrsGraph : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  Epetra_Map * NewRowMap_;

 public:

  ///
  /** Destructor
   */
  ~ZoltanOrder_CrsGraph();

  ///
  /** Constructor
   */
  ZoltanOrder_CrsGraph()
  : NewRowMap_(0)
  {}

  ///
  /** Generates the reordered Epetra_CrsGraph from the input object.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //ZOLTAN_ORDER

#endif //EPETRAEXT_ZOLTANORDER_CRSGRAPH_H
