//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#ifndef ET_CRSGRAPH_BTF_H
#define ET_CRSGRAPH_BTF_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;
class Epetra_Map;

namespace EpetraExt {

///
/** Block Triangular Factorization (Reordering) of Epetra_CrsGraph
 *
 * Uses Alex Pothien's BTF algorithm to find a block upper triangular
 * ordering for a Epetra_CrsGraph.
 */

class CrsGraph_BTF : public StructuralSameTypeTransform<Epetra_CrsGraph>
{

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

#endif //ET_CRSGRAPH_BTF_H
