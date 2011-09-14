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

#include "GLpApp_AdvDiffReactOptModel.hpp"
#include "Teuchos_CommandLineProcessor.hpp"

namespace GLpApp {

/** \brief A utility class for creating an <tt>AdvDiffReactOptModelCreator/tt>
 * object by reading from the command-line.
 */
class AdvDiffReactOptModelCreator {
public:

  /** \brief . */
  AdvDiffReactOptModelCreator();

  /** \brief . */
  void setupCLP( Teuchos::CommandLineProcessor *clp );

  /** \brief . */
  Teuchos::RefCountPtr<AdvDiffReactOptModel>
  createModel(
    const Teuchos::RefCountPtr<const Epetra_Comm>     &comm
    ,std::ostream                                     *out  = NULL
    ) const;

private:

  double              len_x_;
  double              len_y_;
  int                 local_nx_;
  int                 local_ny_;
  std::string         geomFileBase_;
  int                 np_;
  bool                normalizeBasis_;
  double              reactionRate_;
  double              beta_;
  double              x0_;
  double              p0_;
  bool                supportDerivatives_;
  
};

} // namespace GLpApp
