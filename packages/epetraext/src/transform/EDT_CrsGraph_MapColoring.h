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

#ifndef EDT_CRSGRAPH_MAPCOLORING_H
#define EDT_CRSGRAPH_MAPCOLORING_H

#include <Epetra_Transform.h>

class Epetra_CrsGraph;
class Epetra_MapColoring;

namespace EpetraExt {

///
/** Map Coloring of independent columns in a Graph
 *
 * Generates a Epetra_MapColoring object for which all column indices form
 * independent sets.
 */

class CrsGraph_MapColoring : public StructuralTransform<Epetra_CrsGraph,Epetra_MapColoring>
{

 public:

  ///
  /** Destructor
   */
  ~CrsGraph_MapColoring() {}

  ///
  /** Constructor
   */
  CrsGraph_MapColoring( bool verbose = false )
  : verbose_(verbose)
  {}

  ///
  /** Generates the Epetra_MapColoring object from an input Epetra_CrsGraph
   */
  NewTypeRef operator()( OriginalTypeRef orig );

 private:

  bool verbose_;

};

} //namespace EpetraExt

#endif //EDT_CRSGRAPH_MAPCOLORING_H
