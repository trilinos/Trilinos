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

#ifndef EpetraExt_CRSGRAPH_MAPCOLORINGINDEX_H
#define EpetraExt_CRSGRAPH_MAPCOLORINGINDEX_H

#include <EpetraExt_Transform.h>

#include <vector>

class Epetra_CrsGraph;
class Epetra_MapColoring;
class Epetra_IntVector;

namespace EpetraExt {

///
/** Generates a std::vector of Epetra_IntVector's to be used to map perturbation
 * contributions to a CrsGraph/CrsMatrix from a perturbed vector.
 */

class CrsGraph_MapColoringIndex
: public StructuralTransform< Epetra_CrsGraph,std::vector<Epetra_IntVector> > {

  const Epetra_MapColoring & ColorMap_;

 public:

  ///
  /** Destructor
   */
  ~CrsGraph_MapColoringIndex() {}

  ///
  /** Constructor
   * input param ColorMap defines the perturbation coloring
   */
  CrsGraph_MapColoringIndex( const Epetra_MapColoring & ColorMap )
  : ColorMap_( ColorMap )
  {}

  ///
  /** Generates a std::vector<Epetra_IntVector> from the input Epetra_CrsGraph
   */
  NewTypeRef operator()( OriginalTypeRef orig );
};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_MAPCOLORINGINDEX_H
