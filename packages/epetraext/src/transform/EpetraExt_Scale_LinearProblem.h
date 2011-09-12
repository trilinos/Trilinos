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
                                                                                                    
#ifndef EpetraExt_LINEARPROBLEM_SCALE_H
#define EpetraExt_LINEARPROBLEM_SCALE_H

#include <vector>

#include <EpetraExt_Transform.h>

class Epetra_LinearProblem;
class Epetra_Vector;

namespace EpetraExt {

///
/** Given an input Epetra_LinearProblem, recursive, left and right scaling are performed.
 */
class LinearProblem_Scale : public InPlaceTransform<Epetra_LinearProblem>
{
 public:

  enum ScaleType{ Sum, Max, Diag, None };

  ///
  /** Destructor
   */
  ~LinearProblem_Scale();

  ///
  /** Constructor
   */
  LinearProblem_Scale( ScaleType left = Sum,
                       ScaleType right = Sum,
                       double exp_fac = 1.0,
                       int iterations = 1 )
  : lScale_(left),
    rScale_(right),
    expFac_(exp_fac),
    iters_(iterations),
    scaled_(false)
  {}

  ///
  /** Applies forward scaling
   */
  bool fwd();

  ///
  /** Reverses scaling
   */
  bool rvs();

 private:

  const ScaleType lScale_;
  const ScaleType rScale_;

  const double expFac_;

  const int iters_;

  bool scaled_;

  std::vector<Epetra_Vector*> lScaleVecs_;
  std::vector<Epetra_Vector*> rScaleVecs_;
};

} //namespace EpetraExt

#endif //EpetraExt_LINEARPROBLEM_SCALE_H

