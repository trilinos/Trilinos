// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TEST_HELPERS_HO_H
#define MUELU_TEST_HELPERS_HO_H
#include "MueLu_ConfigDefs.hpp"

// Intrepid
#ifdef HAVE_MUELU_INTREPID2
#include "Kokkos_DynRankView.hpp"

#include "MueLu_TestHelpers.hpp"

#ifdef HAVE_MUELU_EPETRA
#include "Epetra_FECrsMatrix.h"
#endif

#include "MueLu_Utilities_def.hpp"

namespace MueLuTests {
namespace TestHelpers {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AllocateEpetraFECrsMatrix(RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& pn_rowmap, RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > pn_colmap, Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& B) {
  throw MueLu::Exceptions::RuntimeError("MueLuTests::TestHelpers::AllocateEpetraFECrsMatrix only works for Tpetra::KokkosCompat::KokkosSerialWrapperNode");
}

#if defined(HAVE_MUELU_EPETRA) &&                   \
    (!defined(HAVE_MUELU_EXPLICIT_INSTANTIATION) || \
     (defined(HAVE_MUELU_EXPLICIT_INSTANTIATION) && defined(HAVE_TPETRA_INST_SERIAL)))
template <>
void AllocateEpetraFECrsMatrix<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(RCP<const Xpetra::Map<int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& pn_rowmap, RCP<const Xpetra::Map<int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> > pn_colmap, Teuchos::RCP<Xpetra::Matrix<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode> >& B)

{
  // Epetra is hard
  const Epetra_Map& pn_rowmap_epetra = Xpetra::toEpetra(*pn_rowmap);
  const Epetra_Map& pn_colmap_epetra = Xpetra::toEpetra(*pn_colmap);
  RCP<Epetra_CrsMatrix> B_epetra     = rcp(new Epetra_FECrsMatrix(Copy, pn_rowmap_epetra, pn_colmap_epetra, 0));
  B                                  = MueLu::Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double, int, int, Tpetra::KokkosCompat::KokkosSerialWrapperNode>(B_epetra);
}
#endif

// Here nx is the number of nodes on the underlying (p=1) mesh.
// This mesh is then promoted up to degree
// Teuchos::RCP<Matrix>
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
Build1DPseudoPoissonHigherOrder(GlobalOrdinal nx, int degree,
                                Kokkos::DynRankView<LocalOrdinal, typename Node::device_type>& elem_to_node,
                                Xpetra::UnderlyingLib lib) {
#include "MueLu_UseShortNames.hpp"
  using Teuchos::arcp;
  using Teuchos::arcp_reinterpret_cast;
  using Teuchos::arcpFromArrayView;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcp_implicit_cast;
  using Teuchos::rcpFromRef;
  GO go_invalid = Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid();

  Teuchos::ParameterList matrixList;
  matrixList.set("nx", nx);
  matrixList.set("matrixType", "Laplace1D");
  // Build a lower order matrix
  RCP<Matrix> A                       = MueLuTests::TestHelpers::TestFactory<SC, LO, GO, NO>::BuildMatrix(matrixList, lib);
  RCP<const Teuchos::Comm<int> > comm = A->getRowMap()->getComm();
  int MyPID                           = comm->getRank();
  int Nproc                           = comm->getSize();

  // Get maps
  RCP<CrsMatrix> Acrs      = rcp_dynamic_cast<CrsMatrixWrap>(A)->getCrsMatrix();
  RCP<const Map> p1_colmap = Acrs->getColMap();
  RCP<const Map> p1_rowmap = Acrs->getRowMap();

  // Count edges.   For shared edges, lower PID gets the owning nodes
  GO global_num_nodes       = p1_rowmap->getGlobalNumElements();
  size_t local_num_nodes    = p1_rowmap->getLocalNumElements();
  GO global_num_elements    = global_num_nodes - 1;
  size_t local_num_elements = local_num_nodes;
  if (p1_rowmap->getGlobalElement(local_num_elements - 1) == global_num_nodes - 1) local_num_elements--;

  printf("[%d] P1 Problem Size: nodes=%d/%d elements=%d/%d\n", MyPID, (int)local_num_nodes, (int)global_num_nodes, (int)local_num_elements, (int)global_num_elements);

  int num_edge_dofs            = (degree - 1) * local_num_elements;
  size_t p1_num_ghost_col_dofs = p1_colmap->getLocalNumElements() - local_num_nodes;

  // Scansum owned edge counts
  int edge_start = 0;
  Teuchos::scan(*comm, Teuchos::REDUCE_SUM, 1, &num_edge_dofs, &edge_start);
  edge_start -= num_edge_dofs;
  GO go_edge_start = global_num_nodes + edge_start;

  // Build owned pn map
  Teuchos::Array<GO> pn_owned_dofs(local_num_nodes + num_edge_dofs);
  for (size_t i = 0; i < local_num_nodes; i++)
    pn_owned_dofs[i] = p1_rowmap->getGlobalElement(i);
  for (size_t i = 0; i < (size_t)num_edge_dofs; i++)
    pn_owned_dofs[local_num_nodes + i] = go_edge_start + i;
  RCP<const Map> pn_rowmap = MapFactory::Build(lib, go_invalid, pn_owned_dofs(), p1_rowmap->getIndexBase(), comm);

  // Build owned column map in [E|T]petra ordering - offproc nodes last
  size_t pn_num_col_dofs = pn_owned_dofs.size() + p1_num_ghost_col_dofs;
  if (MyPID != 0) pn_num_col_dofs += (degree - 1);  // pn ghosts; left side only
  Teuchos::Array<GO> pn_col_dofs(pn_num_col_dofs);
  for (size_t i = 0; i < (size_t)pn_owned_dofs.size(); i++)
    pn_col_dofs[i] = pn_owned_dofs[i];  // onproc

  //      printf("[%d/%d] DEBUG: degree = %d local_num_elements = %d local_num_nodes = %d num_edge_dofs = %d p1_num_ghost_col_dofs =%d  pn_owned_dofs.size() = %d,pn_num_col_dofs = %d mycount = %d\n",MyPID,Nproc,degree,(int)local_num_elements,(int)local_num_nodes,(int)num_edge_dofs,(int)p1_num_ghost_col_dofs, (int)pn_owned_dofs.size(),(int)pn_num_col_dofs,(int)pn_owned_dofs.size()+(MyPID!=0)*(degree) + (MyPID!=Nproc-1) );

  // We have to copy the ghosts from the p1_rowmap as well as the new edge dofs.
  // This needs to follow [E|T]petra ordering
  size_t idx = pn_owned_dofs.size();
  if (MyPID != 0) {
    // Left side nodal
    pn_col_dofs[idx] = p1_colmap->getGlobalElement(p1_rowmap->getLocalNumElements());
    idx++;
    // Left side, edge
    for (size_t i = 0; i < (size_t)(degree - 1); i++) {
      pn_col_dofs[idx] = go_edge_start - (degree - 1) + i;
      idx++;
    }
  }
  if (MyPID != Nproc - 1) {
    // Right side nodal
    pn_col_dofs[idx] = p1_colmap->getGlobalElement(p1_colmap->getLocalNumElements() - 1);
    idx++;
  }

  RCP<const Map> pn_colmap = MapFactory::Build(lib, go_invalid, pn_col_dofs(), p1_rowmap->getIndexBase(), comm);

#if 0
      {
        printf("[%d] TH P1 RowMap = ",MyPID);
        for(size_t i=0; i<p1_rowmap->getLocalNumElements(); i++)
          printf("%d ",(int)p1_rowmap->getGlobalElement(i));
        printf("\n");
        printf("[%d] TH P1 ColMap = ",MyPID);
        for(size_t i=0; i<p1_colmap->getLocalNumElements(); i++)
          printf("%d ",(int) p1_colmap->getGlobalElement(i));
        printf("\n");
        printf("[%d] TH Pn RowMap = ",MyPID);
        for(size_t i=0; i<pn_rowmap->getLocalNumElements(); i++)
          printf("%d ",(int) pn_rowmap->getGlobalElement(i));
        printf("\n");
        printf("[%d] TH Pn ColMap = ",MyPID);
        for(size_t i=0; i<pn_colmap->getLocalNumElements(); i++)
          printf("%d ",(int) pn_colmap->getGlobalElement(i));
        printf("\n");
        fflush(stdout);
      }
#endif

  // Fill elem_to_node using Kirby-style ordering
  // Ownership rule: I own the element if I own the left node in said element
  Kokkos::resize(elem_to_node, local_num_elements, degree + 1);
  auto elem_to_node_host = Kokkos::create_mirror_view(elem_to_node);
  for (size_t i = 0; i < local_num_elements; i++) {
    // End Nodes
    // NTS: This only works for lines
    GO row_gid                   = pn_colmap->getGlobalElement(i);
    GO col_gid                   = row_gid + 1;
    elem_to_node_host(i, 0)      = i;
    elem_to_node_host(i, degree) = pn_colmap->getLocalElement(col_gid);

    // Middle nodes (in local ids)
    for (size_t j = 0; j < (size_t)(degree - 1); j++)
      elem_to_node_host(i, 1 + j) = pn_colmap->getLocalElement(go_edge_start + i * (degree - 1) + j);
  }

  // Since we're inserting off-proc, we really need to use the Epetra_FECrsMatrix here if we're in Epetra mode
  RCP<Matrix> B;
  if (lib == Xpetra::UseEpetra) {
    AllocateEpetraFECrsMatrix(pn_rowmap, pn_colmap, B);
  } else {
    // Tpetra is easy
    B = rcp(new CrsMatrixWrap(pn_rowmap, pn_colmap, global_num_elements));
  }

  // Assemble pseudo-poisson matrix
  for (size_t i = 0; i < local_num_elements; i++) {
    // Fill in a fake stiffness matrix
    for (int j = 0; j < degree + 1; j++) {
      // Dirichlet check
      if ((j == 0 && pn_colmap->getGlobalElement(elem_to_node_host(i, j)) == 0) ||
          (j == degree && pn_colmap->getGlobalElement(elem_to_node_host(i, j)) == global_num_nodes - 1)) {
        // Stick a 1 on the diagonal
        GO row_gid = pn_colmap->getGlobalElement(elem_to_node_host(i, j));
        Teuchos::Array<GO> index(1);
        index[0] = row_gid;
        Teuchos::Array<SC> value(1);
        value[0] = 1.0;
        B->insertGlobalValues(row_gid, index(), value());
        continue;
      }

      GO rowj = pn_colmap->getGlobalElement(elem_to_node_host(i, j));
      for (int k = 0; k < degree + 1; k++) {
        GO rowk = pn_colmap->getGlobalElement(elem_to_node_host(i, k));
        Teuchos::Array<GO> index(1);
        index[0] = rowk;
        Teuchos::Array<SC> value(1);
        if (j == 0 && k == 0)
          value[0] = 1.0;
        else if (j == 0 && k == degree)
          value[0] = -1.0;
        else if (j == degree && k == 0)
          value[0] = -1.0;
        else if (j == degree && k == degree)
          value[0] = 1.0;
        else if (j == k)
          value[0] = (degree + 1) / 100.0;
        else
          value[0] = -1.0 / 100;
        B->insertGlobalValues(rowj, index(), value());
      }
    }
  }

  B->fillComplete(pn_rowmap, pn_rowmap);
  Kokkos::deep_copy(elem_to_node, elem_to_node_host);

#if 0
      std::cout<<"*** Pseudo Poisson ***"<<std::endl;
      Teuchos::FancyOStream ofs(rcp(&std::cout,false));
      B->describe(ofs,Teuchos::VERB_EXTREME);
      std::cout<<"**********************"<<std::endl;

      printf("\n[%d] Pn elem_to_node = \n***\n",MyPID);
      for(size_t i=0; i<(size_t)elem_to_node.dimension(0); i++) {
        for(size_t j=0; j<(size_t)elem_to_node.dimension(1); j++)
          printf("%d[%d] ",(int)elem_to_node(i,j),(int)pn_colmap->getGlobalElement(elem_to_node(i,j)));
        printf("\n");
        }
      printf("***\n");
#endif

  return B;
}  // Build1DPseudoPoissonHigherOrder()

}  // namespace TestHelpers

}  // namespace MueLuTests

#endif  // ifdef HAVE_MUELU_INTREPID2

#endif  // ifndef MUELU_TEST_HELPERS_HO_H
