/*
//@HEADER
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

#ifndef _read_matrix_hpp_
#define _read_matrix_hpp_

#include <algorithm>
#include <iostream>
#include <fstream>

#include "Teuchos_Time.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "MatrixMarket_Tpetra.hpp"

/// \fn read_matrix_mm
/// \brief Read a sparse matrix from a Harwell-Boeing file
///
/// \param hb_file [in] Path of a Harwell-Boeing file.  To be opened
///   and read only by Process 0 of the given communicator.
/// \param comm [in] Communicator object, over which to distribute the
///   sparse matrix to return.
/// \param node [in] Node instance to be used by the returned sparse
///   matrix.
///
/// \return The sparse matrix, distributed over the given communicator.
///
/// \note This defers to Tpetra::Utils::readHBMatrix for reading the
///   sparse matrix from the Matrix Market file.
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_matrix_hb (const std::string& hb_file,
                const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                Teuchos::RCP<Node> node)
{
  using Teuchos::RCP;
  using std::cout;
  using std::endl;

  RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewCounter ("read_matrix");
  RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A;
  {
    Teuchos::TimeMonitor timeMon (*timer);
    Tpetra::Utils::readHBMatrix (hb_file, comm, node, A);
  }
  if (comm->getRank () == 0) {
    cout << "Proc 0: Time in seconds to read the Harwell-Boeing - format "
         << "sparse matrix and finish fillComplete(): "
         << timer->totalElapsedTime () << endl;
  }
  return A;
}

#endif // _read_matrix_hpp_

