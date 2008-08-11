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

#ifndef EpetraExt_LINEARPROBLEM_SOLVERMAP_H
#define EpetraExt_LINEARPROBLEM_SOLVERMAP_H

#include <EpetraExt_Transform.h>

#include <EpetraExt_SolverMap_CrsMatrix.h>

class Epetra_LinearProblem;

namespace EpetraExt {

///
/** Constructs a LinearProblem with a "fixed" Column Map for the CrsMatrix.
 *  Almost entirely a view except for the "fixed" Epetra_CrsGraph.
 */
class LinearProblem_SolverMap : public StructuralSameTypeTransform<Epetra_LinearProblem>
{
  CrsMatrix_SolverMap SolverMapTrans_;

 public:

  ///
  /** Destructor
   */
  ~LinearProblem_SolverMap();

  ///
  /** Constructor
   */
  LinearProblem_SolverMap() {}

  ///
  /** Constructs "fixed" Epetra_LinearProblem
   */
  NewTypeRef operator()( OriginalTypeRef orig );
};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_SOLVERMAP_H

