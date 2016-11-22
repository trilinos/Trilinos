/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_DETAILS_CONTAINERFACTORY_H
#define IFPACK2_DETAILS_CONTAINERFACTORY_H

#include "Ifpack2_Container.hpp"
#include "Ifpack2_TriDiContainer.hpp"
#include "Ifpack2_DenseContainer.hpp"
#include "Ifpack2_SparseContainer.hpp"
#include "Ifpack2_BandedContainer.hpp"
#include "Ifpack2_Partitioner.hpp"
#include "Ifpack2_ILUT_decl.hpp"
#ifdef HAVE_IFPACK2_AMESOS2
#  include "Ifpack2_Details_Amesos2Wrapper.hpp"
#endif
#include "Tpetra_RowMatrix_decl.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Ptr.hpp"
#include "Teuchos_ArrayViewDecl.hpp"
#include <string>


namespace Ifpack2 {
namespace Details {

//Based on the ETIs at the bottom of BlockRelaxation_def, these are all the supported container types:
//pass one of these as the container name

//Dense
//SparseILUT
//SparseAmesos
//TriDi
//Banded

template <typename MatrixType>
Teuchos::RCP< ::Ifpack2::Container< ::Tpetra::RowMatrix<typename MatrixType::scalar_type,
                                                        typename MatrixType::local_ordinal_type,
                                                        typename MatrixType::global_ordinal_type,
                                                        typename MatrixType::node_type> > >
createContainer (const std::string& containerName,
                 const Teuchos::RCP< const MatrixType>& A,
                 const Teuchos::Array< Teuchos::Array< typename MatrixType::local_ordinal_type> >& localRows,
                 const Teuchos::RCP<const Tpetra::Import<typename MatrixType::local_ordinal_type,
                                                         typename MatrixType::global_ordinal_type,
                                                         typename MatrixType::node_type> > importer,
                 int OverlapLevel,
                 typename MatrixType::scalar_type DampingFactor)
{
  using Teuchos::rcp;
  typedef Tpetra::RowMatrix<typename MatrixType::scalar_type,
    typename MatrixType::local_ordinal_type,
    typename MatrixType::global_ordinal_type,
    typename MatrixType::node_type> row_matrix_type;
  // Use std::decay to remove const-ness and references.
  static_assert (std::is_same<typename std::decay<MatrixType>::type, row_matrix_type>::value,
                 "MatrixType must be a Tpetra::RowMatrix specialization.");

  if (containerName == "TriDi") {
    return rcp (new TriDiContainer<MatrixType, typename MatrixType::scalar_type> (A, localRows, importer, OverlapLevel, DampingFactor));
  }
  else if (containerName == "Dense") {
    return rcp (new DenseContainer<MatrixType, typename MatrixType::scalar_type> (A, localRows, importer, OverlapLevel, DampingFactor));
  }
  else if (containerName == "SparseILUT") {
    return rcp (new SparseContainer<MatrixType, ILUT<MatrixType> > (A, localRows, importer, OverlapLevel, DampingFactor));
  }
#ifdef HAVE_IFPACK2_AMESOS2
  else if (containerName == "SparseAmesos2" || containerName == "SparseAmesos") {
    return rcp (new SparseContainer<MatrixType, Amesos2Wrapper<MatrixType> > (A, localRows, importer, OverlapLevel, DampingFactor));
  }
#endif
  else if (containerName == "Banded") {
    return rcp (new BandedContainer<MatrixType, typename MatrixType::scalar_type> (A, localRows, importer, OverlapLevel, DampingFactor));
  }

  TEUCHOS_TEST_FOR_EXCEPTION
    (true, std::invalid_argument, "Ifpack2::Details::createContainer: Input "
     "argument containerName=\"" << containerName << "\" is invalid.  Valid "
     "values include: \"TriDi\", \"Dense\", \"SparseILUT\", \"SparseAmesos2\", "
     "and \"Banded\".");
}

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CONTAINERFACTORY_H
