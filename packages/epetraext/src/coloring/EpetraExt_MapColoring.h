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

#ifndef EpetraExt_CRSGRAPH_MAPCOLORING_H
#define EpetraExt_CRSGRAPH_MAPCOLORING_H

#include <EpetraExt_Transform.h>

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

  enum ColoringAlgorithm{ ALGO_GREEDY, ALGO_LUBI };

  ///
  /** Destructor
   */
  ~CrsGraph_MapColoring() {}

  ///
  /** Constructor
   */
  CrsGraph_MapColoring( ColoringAlgorithm algo = ALGO_GREEDY,
                        bool verbose = false,
                        int reordering = 0,
		        bool colorParallel = false,
                        bool serialBoundaryColoring = true,
                        int verbosityLevel = 0 )
  : algo_(algo),
    verbose_(verbose),
    reordering_(reordering),
    colorParallel_(colorParallel),
    SerialBoundaryColoring_(serialBoundaryColoring),
    verbosityLevel_(verbosityLevel)
  {}

  ///
  /** Generates the Epetra_MapColoring object from an input Epetra_CrsGraph
   */
  NewTypeRef operator()( OriginalTypeRef orig );

 private:

  const bool verbose_;

  const ColoringAlgorithm algo_;

  const int reordering_;
  const bool colorParallel_;

  const bool SerialBoundaryColoring_;

  const int verbosityLevel_;

};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_MAPCOLORING_H
