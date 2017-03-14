// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// @HEADER

#ifndef RBGEN_ISVD_MULTISDAUDV_H
#define RBGEN_ISVD_MULTISDAUDV_H

#include "RBGen_ISVDMultiSDA.h"
#include "RBGen_ISVDUDV.h"

namespace RBGen {

  //! IncSVD method implementing UDV storage with multiple steepest descent (variant A) passes.
  class ISVD_MultiSDAUDV : public virtual ISVDUDV, public virtual ISVDMultiSDA {
    public:
      //! @name Constructor/Destructor.
      //@{

      //! Default constructor.
      ISVD_MultiSDAUDV();

      //! Destructor.
      virtual ~ISVD_MultiSDAUDV() {};
      //@}

      //! @name Set Methods
      //@{

      //! Initialize the method with the given parameter list and snapshot set.
      void Initialize( const Teuchos::RCP< Teuchos::ParameterList >& params,
          const Teuchos::RCP< const Epetra_MultiVector >& init,
          const Teuchos::RCP< RBGen::FileIOHandler< Epetra_Operator > >& fileio = Teuchos::null );

      //@}
  };

} // end of RBGen namespace

#endif // RBGEN_ISVD_MULTISDAUDV_H
