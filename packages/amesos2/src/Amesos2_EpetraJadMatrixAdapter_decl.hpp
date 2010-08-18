// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package 
//                  Copyright 2010 Sandia Corporation
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
//
// @HEADER

/**
  \file   Amesos2_EpetraJadMatrixAdapter_decl.hpp
  \author Eric Bavier <etbavier@sandia.gov>
  \date   Thu Jul 22 10:27:57 CDT 2010
  
  \brief  Amesos2 matrix adapter for the Epetra_JadMatrix class.
*/

#ifndef AMESOS2_EPETRA_JADMATRIX_ADAPTER_DECL_HPP
#define AMESOS2_EPETRA_JADMATRIX_ADAPTER_DECL_HPP

#include <Epetra_JadMatrix.h>

#include "Amesos2_MatrixAdapter_decl.hpp"
#include "Amesos2_EpetraRowMatrixAdapter.hpp"

namespace Amesos {


/**
 * \brief Amesos2 matrix adapter for the Epetra_JadMatrix class.
 *
 * This class inherits all methods and functionality from the
 * MatrixAdapter<Epetra_RowMatrix> class, but overrides the \c matrix_type
 * typedef and the \c name attribute.
 */
template<>
class MatrixAdapter< Epetra_JadMatrix > : public MatrixAdapter< Epetra_RowMatrix >
{
public:
  // override matrix typedef
  typedef Epetra_JadMatrix matrix_type;

  /// The name of this adapter class.
  static const char* name;


  /// Default constructor
  MatrixAdapter();


  /// Copy constructor
  MatrixAdapter(const MatrixAdapter<matrix_type>& adapter);


  /**
   * \brief Initialize an adapter from a matrix RCP
   *
   * \param m An RCP pointing to the matrix which is to be wrapped.
   */
  MatrixAdapter(const Teuchos::RCP<matrix_type>& m);

};


} // end namespace Amesos

#endif  // AMESOS2_EPETRA_JADMATRIX_ADAPTER_DECL_HPP
