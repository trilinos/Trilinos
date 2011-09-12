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

#ifndef EpetraExt_CRSGRAPH_SUBCOPY_H
#define EpetraExt_CRSGRAPH_SUBCOPY_H

#include <EpetraExt_Transform.h>
#include <Epetra_Map.h>

class Epetra_CrsMatrix;
class Epetra_Import;

namespace EpetraExt {

//! Generates a sub-block view of a Epetra_CrsMatrix
class CrsMatrix_SubCopy : public SameTypeTransform<Epetra_CrsMatrix> {

  Epetra_Map  newRowMap_;
  Epetra_Map  newColMap_;
  Epetra_Map  newDomainMap_;
  Epetra_Map  newRangeMap_;
  Epetra_Import * importer_;


 public:

  //! Destructor
  ~CrsMatrix_SubCopy();

  //! Constructor
  CrsMatrix_SubCopy( const Epetra_Map  & newMap )
  : newRowMap_(newMap),
    newColMap_(newMap),
    newDomainMap_(newMap),
    newRangeMap_(newMap)
  {}

  //! Constructor
  CrsMatrix_SubCopy( const Epetra_Map  & newRangeAndRowMap,
		    const Epetra_Map  & newDomainMap )
  : newRowMap_(newRangeAndRowMap),
    newColMap_(newDomainMap),
    newDomainMap_(newDomainMap),
    newRangeMap_(newRangeAndRowMap)
  {}

  //! Transformation Operator
  NewTypeRef operator()( OriginalTypeRef orig );

  ///
  bool fwd();

  ///
  bool rvs();
};

} //namespace EpetraExt

#endif //EpetraExt_CRSGRAPH_SUBCOPY_H
