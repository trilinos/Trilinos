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

#ifndef EpetraExt_MAPCOLORINGINDEX_H
#define EpetraExt_MAPCOLORINGINDEX_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

#include <EpetraExt_TCrsGraph_MapColoringIndex.h>

namespace EpetraExt {

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

class CrsGraph_MapColoringIndex
: public TCrsGraph_MapColoringIndex<int> {

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
  : TCrsGraph_MapColoringIndex<int>( ColorMap )
  {}

  std::vector<Epetra_IntVector>& operator()( Epetra_CrsGraph& orig );
};
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

class CrsGraph_MapColoringIndex64
: public TCrsGraph_MapColoringIndex<long long> {

 public:

  ///
  /** Destructor
   */
  ~CrsGraph_MapColoringIndex64() {}

  ///
  /** Constructor
   * input param ColorMap defines the perturbation coloring
   */
  CrsGraph_MapColoringIndex64( const Epetra_MapColoring & ColorMap )
  : TCrsGraph_MapColoringIndex<long long>( ColorMap )
  {}

  std::vector<Epetra_LongLongVector>& operator()( Epetra_CrsGraph& orig );
};
#endif

// #ifndef EPETRA_NO_32BIT_GLOBAL_INDICES

// std::vector<Epetra_IntVector>&
// CrsGraph_MapColoringIndex::operator()( Epetra_CrsGraph& orig )
// {
//   return TCrsGraph_MapColoringIndex<int>::operator()(orig);
// }

// #endif

// #ifndef EPETRA_NO_64BIT_GLOBAL_INDICES

// std::vector<Epetra_LongLongVector>&
// CrsGraph_MapColoringIndex64::operator()( Epetra_CrsGraph& orig )
// {
//   return TCrsGraph_MapColoringIndex<long long>::operator()(orig);
// }

// #endif

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_MAPCOLORINGINDEX_H
