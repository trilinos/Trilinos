// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_TPETRA_TYPEDEFS_HPP
#define NOX_TPETRA_TYPEDEFS_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Map_fwd.hpp"
#include "Tpetra_MultiVector.hpp" // fwd doesn't have default enabled template types
#include "Tpetra_Vector_fwd.hpp"
#include "Tpetra_Export_fwd.hpp"
#include "Tpetra_Import_fwd.hpp"
#include "Tpetra_RowGraph_fwd.hpp"
#include "Tpetra_RowMatrix_fwd.hpp"
#include "Tpetra_CrsGraph_fwd.hpp"
#include "Tpetra_CrsMatrix_fwd.hpp"
#include "Tpetra_Operator_fwd.hpp"

namespace NOX {

  using Scalar = Tpetra::MultiVector<>::scalar_type;
  using LocalOrdinal = Tpetra::Map<>::local_ordinal_type;
  using GlobalOrdinal = Tpetra::Map<>::global_ordinal_type;
  using GlobalSizeType = Tpetra::global_size_t;
  using DeviceSpace = Kokkos::DefaultExecutionSpace;
  using NodeType = Tpetra::KokkosCompat::KokkosDeviceWrapperNode<DeviceSpace>;
  
  using TMap = Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType>;
  using ConstTMap = const Tpetra::Map<LocalOrdinal, GlobalOrdinal, NodeType>;
  using TImport = Tpetra::Import<LocalOrdinal,GlobalOrdinal,NodeType>;
  using TExport = Tpetra::Export<LocalOrdinal,GlobalOrdinal,NodeType>;
  using TVector = Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,NodeType>;
  using TMultiVector = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,NodeType>;
  using TRowGraph = Tpetra::RowGraph< LocalOrdinal, GlobalOrdinal, NodeType>;
  using TRowMatrix = Tpetra::RowMatrix< Scalar, LocalOrdinal, GlobalOrdinal, NodeType>;
  using TCrsGraph = Tpetra::CrsGraph< LocalOrdinal, GlobalOrdinal, NodeType>;
  using TCrsMatrix = Tpetra::CrsMatrix< Scalar, LocalOrdinal, GlobalOrdinal, NodeType>;
  using TOperator =  Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, NodeType>;

}

#endif
