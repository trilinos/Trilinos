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
                                                                                                    
#ifndef EpetraExt_CRSGRAPH_OVERLAP_H
#define EpetraExt_CRSGRAPH_OVERLAP_H

#include <EpetraExt_Transform.h>

class Epetra_BlockMap;
class Epetra_CrsGraph;

namespace EpetraExt {

///
/** Given an input Epetra_CrsGraph, a "overlapped" Epetra_CrsGraph is generated
 * including rows associated with off processor contributions.
 */
class CrsGraph_Overlap : public StructuralSameTypeTransform<Epetra_CrsGraph> {

  const int levelOverlap_;

  const bool squareLocalBlock_;

  Epetra_BlockMap * OverlapMap_;

 public:

  ///
  /** Destructor
   */
  ~CrsGraph_Overlap();

  ///
  /** Constructor
   */
  CrsGraph_Overlap( int overlap, bool squareLocalBlock = false )
  : levelOverlap_(overlap),
    squareLocalBlock_(squareLocalBlock),
    OverlapMap_(0)
  {}

  ///
  /** Constructs "overlapped" Epetra_CrsGraph from original
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_OVERLAP_H
