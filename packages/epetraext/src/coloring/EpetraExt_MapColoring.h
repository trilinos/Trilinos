//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EpetraExt_CRSGRAPH_MAPCOLORING_H
#define EpetraExt_CRSGRAPH_MAPCOLORING_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

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

  enum ColoringAlgorithm{ GREEDY, LUBY, JONES_PLASSMAN, PSEUDO_PARALLEL  };

  ///
  /** Destructor
   */
  ~CrsGraph_MapColoring() {}

  ///
  /** Constructor
   */
  CrsGraph_MapColoring( ColoringAlgorithm algo = GREEDY,
                        int reordering = 0,
                        bool distance1 = false,
                        int verbosity = 0 )
  : algo_(algo),
    reordering_(reordering),
    distance1_(distance1),
    verbosity_(verbosity)
  {}

  ///
  /** Generates the Epetra_MapColoring object from an input Epetra_CrsGraph
   */
  CrsGraph_MapColoring::NewTypeRef operator()( CrsGraph_MapColoring::OriginalTypeRef orig );

 private:


  const ColoringAlgorithm algo_;

  const int reordering_;
  const bool distance1_;

  const int verbosity_;

  template<typename int_type>
  CrsGraph_MapColoring::NewTypeRef Toperator( CrsGraph_MapColoring::OriginalTypeRef orig );
};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_MAPCOLORING_H
