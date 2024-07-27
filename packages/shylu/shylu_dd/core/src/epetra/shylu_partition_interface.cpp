// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_partition_interface.cpp
    
    \brief Epetra/Tpetra templated interface for calls to Zoltan(Isorropia)/Zoltans

    \author Joshua Dennis Booth
*/
#include "shylu_partition_interface_decl.hpp"
#include "shylu_partition_interface_def.hpp"

#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include "Tpetra_ConfigDefs.hpp"
#endif

namespace ShyLU {

  template class PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector>;

#if defined(HAVE_SHYLU_DDCORE_TPETRA)

#if defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)
  template class PartitionInterface<
    Tpetra::CrsMatrix<float, int, int>,
    Tpetra::MultiVector<float, int, int> >;
#endif // defined(HAVE_TPETRA_INST_FLOAT) && defined(HAVE_TPETRA_INST_INT_INT)

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
  template class PartitionInterface<
    Tpetra::CrsMatrix<double, int, int> ,
    Tpetra::MultiVector<double, int, int> >;
#endif // defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)

#endif // defined(HAVE_SHYLU_DDCORE_TPETRA)

} // namespace ShyLU
