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
                                                                                                    
#ifndef EpetraExt_VECTOR_DIRICHLET_H
#define EpetraExt_VECTOR_DIRICHLET_H

#include <EpetraExt_Transform.h>

#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>

namespace EpetraExt {

///
/** Given an input Epetra_Vector, apply given dirichlet conditions
 */
class Vector_Dirichlet : public InPlaceTransform<Epetra_Vector>
{
 public:

  ///
  /** Destructor
   */
  ~Vector_Dirichlet() {}

  ///
  /** Constructor
   @param Locations Integer Vector containing 1's for Dirichlet rows and 0's if not
   @param Values Vector containing values of the Dirichlet BC's
   */
  Vector_Dirichlet( const Epetra_IntVector & Locations,
                    const Epetra_Vector & Values )
  : locations_( Locations ),
    values_( Values )
  {}

  ///
  /** Applies Dirichlet BC's
   */
  bool fwd();

  ///
  /** NoOp
   */
  bool rvs();

 private:

  const Epetra_IntVector locations_; 
  const Epetra_Vector values_;
};

} //namespace EpetraExt

#endif //EpetraExt_VECTOR_DIRICHLET_H

