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

#ifndef EpetraExt_TCRSGRAPH_MAPCOLORINGINDEX_H
#define EpetraExt_TCRSGRAPH_MAPCOLORINGINDEX_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <EpetraExt_Transform.h>
#include <Epetra_GIDTypeVector.h>
#include <Epetra_MapColoring.h>

#include <vector>
#include <map>

class Epetra_CrsGraph;
class Epetra_IntVector;

namespace EpetraExt {

///
/** Generates a std::vector of Epetra_IntVector's to be used to map perturbation
 * contributions to a CrsGraph/CrsMatrix from a perturbed vector.
 */

template<typename int_type>
class TCrsGraph_MapColoringIndex
: public StructuralTransform< Epetra_CrsGraph,std::vector<typename Epetra_GIDTypeVector<int_type>::impl> > {

  const Epetra_MapColoring & ColorMap_;

 protected:

  ///
  /** Destructor
   */
  ~TCrsGraph_MapColoringIndex() {}

  ///
  /** Constructor
   * input param ColorMap defines the perturbation coloring
   */
  TCrsGraph_MapColoringIndex( const Epetra_MapColoring & ColorMap )
  : ColorMap_( ColorMap )
  {}
public:
  typedef StructuralTransform< Epetra_CrsGraph,std::vector<typename Epetra_GIDTypeVector<int_type>::impl> > Base;
  ///
  /** Generates a std::vector<Epetra_IntVector> from the input Epetra_CrsGraph
   */
  typedef typename Base::NewTypeRef NewTypeRef;
  typedef typename Base::OriginalTypeRef OriginalTypeRef;
  NewTypeRef operator()( OriginalTypeRef orig );
};

//////// Implementations ////////

template<typename int_type>
typename TCrsGraph_MapColoringIndex<int_type>::NewTypeRef
TCrsGraph_MapColoringIndex<int_type>::
operator()( OriginalTypeRef orig )
{
  if(!orig.RowMap(). template GlobalIndicesIsType<int_type>())
    throw "EpetraExt::TCrsGraph_MapColoringIndex::operator(): Global indices mismatch.";

  Base::origObj_ = &orig;

  const Epetra_BlockMap & RowMap = orig.RowMap();
  int nRows = RowMap.NumMyElements();

  int NumColors = ColorMap_.NumColors();
  int * ListOfColors = ColorMap_.ListOfColors();

  std::map<int,int> MapOfColors;
  for( int i = 0; i < NumColors; ++i ) MapOfColors[ ListOfColors[i] ] = i;

  //initial setup of stl vector of IntVectors for indexing
  std::vector<int_type> dummy( nRows, -1 );
  typename Base::NewTypePtr IndexVec = new typename Base::NewType( NumColors, typename Epetra_GIDTypeVector<int_type>::impl( Copy, RowMap, &dummy[0] ) );

  int MaxNumIndices = orig.MaxNumIndices();
  int NumIndices;
  std::vector<int_type> Indices( MaxNumIndices );

  for( int i = 0; i < nRows; ++i )
  {
    orig.ExtractGlobalRowCopy( (int_type) orig.GRID64(i), MaxNumIndices, NumIndices, &Indices[0] );

    for( int j = 0; j < NumIndices; ++j )
     (*IndexVec)[ MapOfColors[ColorMap_(Indices[j])] ][i] = Indices[j];
  }

  Base::newObj_ = IndexVec;

  return *IndexVec;
}

} //namespace EpetraExt

#endif // EpetraExt_TCRSGRAPH_MAPCOLORINGINDEX_H
