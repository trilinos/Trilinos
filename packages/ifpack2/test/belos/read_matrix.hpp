// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
///   matrix.  This is an optional argument; if not provided or null,
///   we make one for you.
///
/// \return The sparse matrix, distributed over the given communicator.
///
/// \note This defers to Tpetra::Utils::readHBMatrix for reading the
///   sparse matrix from the Matrix Market file.
template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
read_matrix_hb (const std::string& hb_file,
                const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                const Teuchos::RCP<Node>& node = Teuchos::null)
{
  using Teuchos::RCP;
  using std::cout;
  using std::endl;

  RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::getNewCounter ("read_matrix");
  RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A;
  {
    Teuchos::TimeMonitor timeMon (*timer);
    Tpetra::Utils::readHBMatrix (hb_file, comm, A);
  }
  if (comm->getRank () == 0) {
    cout << "Proc 0: Time in seconds to read the Harwell-Boeing - format "
         << "sparse matrix and finish fillComplete(): "
         << timer->totalElapsedTime () << endl;
  }
  return A;
}

#endif // _read_matrix_hpp_

