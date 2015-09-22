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
                                                                                                    
#ifndef EpetraExt_LINEARPROBLEM_REINDEX_H
#define EpetraExt_LINEARPROBLEM_REINDEX_H

#include <EpetraExt_Transform.h>

class Epetra_Map;
class Epetra_LinearProblem;

namespace EpetraExt {

class CrsMatrix_Reindex;
class MultiVector_Reindex;

///
/** Given and input Epetra_LinearProblem, a "reindexed" version will be returned
 *  using the given NewRowMap.  If a null map is given, a lexigraphically indexed
 *  LP will be returned.  The data in the new E_LP is a "reindexed" view of the 
 *  original.
 */
class LinearProblem_Reindex : public ViewTransform<Epetra_LinearProblem>
{
  CrsMatrix_Reindex * MatTrans_;
  MultiVector_Reindex * LHSTrans_;
  MultiVector_Reindex * RHSTrans_;

  Epetra_Map * NewRowMap_;

  bool NewRowMapOwned_;

 public:

  ///
  /** Destructor
   */
  ~LinearProblem_Reindex();

  ///
  /** Constructor
   */
  LinearProblem_Reindex( Epetra_Map * NewRowMap )
  : MatTrans_(0),
    LHSTrans_(0),
    RHSTrans_(0),
    NewRowMap_(NewRowMap),
    NewRowMapOwned_(false)
  {}

  ///
  /** Constructs a new view the original LP, "reindexed" using the given NewRowMap.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_REINDEX_H

