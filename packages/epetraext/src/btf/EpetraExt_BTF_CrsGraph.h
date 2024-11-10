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
//                                                                                                    
#ifndef EpetraExt_CRSGRAPH_BTF_H
#define EpetraExt_CRSGRAPH_BTF_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif

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
