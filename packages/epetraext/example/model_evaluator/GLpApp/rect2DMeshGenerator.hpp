/*
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
*/


#include <iosfwd>

class Epetra_IntSerialDenseVector;
class Epetra_SerialDenseMatrix;
class Epetra_IntSerialDenseMatrix;

namespace GLpApp {

/** \brief Generate a simple rectangular 2D triangular mesh that is only
 * partitioned between processors in the <tt>y</tt> direction.
 *
 * ToDo: Finish documentation!
 */
void rect2DMeshGenerator(
  const int                      numProc
  ,const int                     procRank
  ,const double                  len_x
  ,const double                  len_y
  ,const int                     local_nx
  ,const int                     local_ny
  ,const int                     bndy_marker
  ,Epetra_IntSerialDenseVector   *ipindx_out
  ,Epetra_SerialDenseMatrix      *ipcoords_out
  ,Epetra_IntSerialDenseVector   *pindx_out
  ,Epetra_SerialDenseMatrix      *pcoords_out
  ,Epetra_IntSerialDenseMatrix   *t_out
  ,Epetra_IntSerialDenseMatrix   *e_out
  ,std::ostream                  *out
  ,const bool                    dumpAll
  );

} // namespace GLpApp
