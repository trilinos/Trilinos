/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef IFPACK2_DETAILS_ROWMATRIX_HPP
#define IFPACK2_DETAILS_ROWMATRIX_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Ifpack2_Details_throwBecauseDeprecated.hpp"

namespace Ifpack2 {
namespace Details {
    
/// \class RowMatrix
/// \brief All Ifpack2 implementations of Tpetra::RowMatrix must
///   inherit from this class.
/// \tparam MatrixType Tpetra::RowMatrix specialization.
///
/// \warning This class is an implementation detail of Ifpack2.  Users
///   should not rely on its interface.
///
/// This class exists to facilitate Tpetra interface changes.  See
/// e.g., GitHub Issue #2630.
template<class MatrixType>
class RowMatrix :
    public Tpetra::RowMatrix<typename MatrixType::scalar_type,
			     typename MatrixType::local_ordinal_type,
			     typename MatrixType::global_ordinal_type,
			     typename MatrixType::node_type> {
public:
  //! \name Typedefs
  //@{
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;

  //@}
  //! \name Destructor
  //@{

  //! Destructor (virtual for memory safety of derived classes)
  virtual ~RowMatrix () {}

  //@}
  //! @name Work-around implementations of deprecated virtual methods  
  //@{

  /// \brief The global number of diagonal entries.
  ///
  /// \warning This method is DEPRECATED and will be removed soon!
  Tpetra::global_size_t IFPACK2_DEPRECATED
  getGlobalNumDiags () const final
  {
    throwBecauseDeprecated ("getGlobalNumDiags");
    return Tpetra::global_size_t (0);    
  }

  /// \brief The local number of diagonal entries.
  ///
  /// \warning This method is DEPRECATED and will be removed soon!
  std::size_t IFPACK2_DEPRECATED
  getNodeNumDiags () const final
  {
    throwBecauseDeprecated ("getNodeNumDiags");
    return std::size_t (0);    
  }

  /// \brief Whether this graph is locally lower triangular.
  ///
  /// \warning This method is DEPRECATED and will be removed soon!
  bool IFPACK2_DEPRECATED
  isLowerTriangular () const final
  {
    throwBecauseDeprecated ("isLowerTriangular");
    return false;
  }

  /// \brief Whether this graph is locally upper triangular.
  ///
  /// \warning This method is DEPRECATED and will be removed soon!
  bool IFPACK2_DEPRECATED isUpperTriangular() const final
  {
    throwBecauseDeprecated ("isUpperTriangular");
    return false;
  }
};

} // namespace Details
} // namespace Ifpack2

#endif /* IFPACK2_DETAILS_ROWMATRIX_HPP */
