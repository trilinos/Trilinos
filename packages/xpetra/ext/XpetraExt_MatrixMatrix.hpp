// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * XpetraExt_MatrixMatrix_decl.hpp
 *
 *  Created on: Oct 10, 2012
 *      Author: wiesner
 */

#ifndef XPETRAEXT_MATRIXMATRIX_DECL_HPP_
#define XPETRAEXT_MATRIXMATRIX_DECL_HPP_

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"
#include <string>

//#include "XpetraExt_MatrixMatrix.hpp"

#include "Xpetra_Matrix.hpp"

#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix_fwd.hpp>
#endif

#ifdef HAVE_XPETRA_EPETRAEXT
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_RowMatrixTransposer.h>
#endif // HAVE_XPETRA_EPETRAEXT

#ifdef HAVE_XPETRA_TPETRA
#include <TpetraExt_MatrixMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraVector.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif // HAVE_XPETRA_TPETRA

/*! \file Xpetra_MatrixMatrix.hpp

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Xpetra {

namespace Jacobi {



} // end namespace MatrixMatrix

} // end namespace Xpetra


#endif /* XPETRAEXT_MATRIXMATRIX_DECL_HPP_ */
