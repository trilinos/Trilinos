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
//                                                                                                    
#ifndef EpetraExt_CRSGRAPH_BTF_H
#define EpetraExt_CRSGRAPH_BTF_H

#include <EpetraExt_Transform.h>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of Epetra_CrsGraph
 *
 * Uses Alex Pothien's BTF algorithm to find a block upper triangular
 * ordering for a Epetra_CrsGraph.
 */

class CrsGraph_BTF : public StructuralSameTypeTransform<Epetra_CrsGraph> {

 public:

  ///
  /** Destructor
   */
  ~CrsGraph_BTF();

  ///
  /** Default Constructor
   */
  CrsGraph_BTF()
  : NewRowMap_(0),
    NewDomainMap_(0),
    NewGraph_(0)
  {}

  ///
  /** Construction of BTF ordered Epetra_CrsGraph from <tt>orig</tt> object.
   *
   * Preconditions:<ul>
   * </ul>
   *
   * Invariants:<ul>
   * </ul>
   *
   * Postconditions:<ul>
   * </ul>
   *
   */
  NewTypeRef operator()( OriginalTypeRef orig );

 private:

  Epetra_Map * NewRowMap_;
  Epetra_Map * NewDomainMap_;
  
  Epetra_CrsGraph * NewGraph_;

};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_BTF_H
