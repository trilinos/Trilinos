// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//#include <Tpetra_CrsMatrix_decl.hpp>
//#include <Teuchos_RCPDecl.hpp>

#include "neighbours.hpp"

typedef Tpetra::CrsMatrix<double, int, int, KokkosClassic::DefaultNode::DefaultNodeType> tpetra_matrix_type;
typedef Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, int, int, KokkosClassic::DefaultNode::DefaultNodeType> tpetra_multivector_type;
typedef typename Teuchos::ArrayView<const int>::const_iterator iterator_type;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType node_type2;

void BAP1D(Teuchos::RCP<tpetra_matrix_type> BAP, Teuchos::RCP<tpetra_matrix_type> tpetra_prolong, Teuchos::RCP<tpetra_multivector_type> BAP_shrunk, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  // INPUT: BAP = matrix where the uncompressed version of B_DD * A_h * Ptentative must be stored
  // INPUT: tpetra_prolong = Ptentative
  // INPUT: BAP_shrunk = Tpetra:MultiVector contatining the shrunk version of B_DD * A_h * Ptentative resulting from domain decomposition with coloring
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)

  Teuchos::ArrayView<const int> myLocalElements = BAP->getRowMap()->getLocalElementList();
  int mypid                                     = comm->getRank();

  for (int color = 0; color < 3; ++color) {
    Teuchos::ArrayRCP<const double> localBAP = BAP_shrunk->getData(color);

    for (iterator_type it = myLocalElements.begin(); it != myLocalElements.end(); ++it) {
      const int i_local = *it;
      const int aux     = BAP->getRowMap()->getLocalElement(i_local);

      std::vector<int> BAP_inds;
      std::vector<double> BAP_vals;

      int aux2;

      if ((mypid - 1) % 3 == color && (mypid - 1) >= 0 && (mypid - 1) < tpetra_prolong->getGlobalNumCols())
        aux2 = BAP->getColMap()->getLocalElement(mypid - 1);
      else if ((mypid - 2) % 3 == color && (mypid - 2) >= 0 && (mypid - 2) < tpetra_prolong->getGlobalNumCols())
        aux2 = BAP->getColMap()->getLocalElement(mypid - 2);
      else if ((mypid) % 3 == color && (mypid) >= 0 && (mypid) < tpetra_prolong->getGlobalNumCols())
        aux2 = BAP->getColMap()->getLocalElement(mypid);

      if (aux2 >= 0) {
        BAP_inds.emplace_back(aux2);
        BAP_vals.emplace_back(localBAP[aux]);
        BAP->insertLocalValues(aux, BAP_inds, BAP_vals);
      }
    }
  }
}

void BAP2D(Teuchos::RCP<tpetra_matrix_type> BAP, Teuchos::RCP<tpetra_matrix_type> tpetra_prolong, Teuchos::RCP<tpetra_multivector_type> BAP_shrunk, Teuchos::RCP<const Teuchos::Comm<int> > comm, int ndx) {
  // INPUT: BAP = matrix where the uncompressed version of B_DD * A_h * Ptentative must be stored
  // INPUT: tpetra_prolong = Ptentative
  // INPUT: BAP_shrunk = Tpetra:MultiVector contatining the shrunk version of B_DD * A_h * Ptentative resulting from domain decomposition with coloring
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)
  // INPUT: ndx = number of domains along x-direction

  Teuchos::ArrayView<const int> myLocalElements = BAP->getRowMap()->getLocalElementList();
  int mypid                                     = comm->getRank();
  int brick_id                                  = mypid;
  int shifted_id                                = brick_id - 1;

  if (mypid > 0) {
    for (int color = 0; color < 9; ++color) {
      int neighbour                            = -1;
      Teuchos::ArrayRCP<const double> localBAP = BAP_shrunk->getData(color);

      // The following if statements control the neighbours of a subdomain in a 2D brick partitioned mesh
      if (coloring2D(brick_id, ndx) == color && (shifted_id) >= 0 && (shifted_id) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id);

      else if (coloring2D(brick_id - 1, ndx) == color && (shifted_id - 1) >= 0 && (shifted_id - 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - 1);

      else if (coloring2D(brick_id + 1, ndx) == color && (shifted_id + 1) >= 0 && (shifted_id + 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + 1);

      else if (coloring2D(brick_id - ndx, ndx) == color && (shifted_id - ndx) >= 0 && (shifted_id - ndx) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx);

      else if (coloring2D(brick_id - ndx - 1, ndx) == color && (shifted_id - ndx - 1) >= 0 && (shifted_id - ndx - 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx - 1);

      else if (coloring2D(brick_id - ndx + 1, ndx) == color && (shifted_id - ndx + 1) >= 0 && (shifted_id - ndx + 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx + 1);

      else if (coloring2D(brick_id + ndx, ndx) == color && (shifted_id + ndx) >= 0 && (shifted_id + ndx) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx);

      else if (coloring2D(brick_id + ndx - 1, ndx) == color && (shifted_id + ndx - 1) >= 0 && (shifted_id + ndx - 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx - 1);

      else if (coloring2D(brick_id + ndx + 1, ndx) == color && (shifted_id + ndx + 1) >= 0 && (shifted_id + ndx + 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx + 1);

      // in case neighbour>=0, it means that the current MPI processor (=subdomain) has a neighbour with the color analyzed at the current for-loop iteration
      // otherwise the current MPI processor sits on the boundary of the mesh and it has less than 8 neighbours.
      if (neighbour >= 0) {
        for (iterator_type it = myLocalElements.begin(); it != myLocalElements.end(); ++it) {
          const int i_local = *it;
          const int aux     = BAP->getRowMap()->getLocalElement(i_local);

          std::vector<int> BAP_inds;
          std::vector<double> BAP_vals;

          int aux2;

          aux2 = BAP->getColMap()->getLocalElement(neighbour);

          if (aux2 >= 0) {
            BAP_inds.emplace_back(aux2);
            BAP_vals.emplace_back(localBAP[aux]);
            BAP->insertLocalValues(aux, BAP_inds, BAP_vals);
          }
        }
      }
    }
  }
}

void BAP3D(Teuchos::RCP<tpetra_matrix_type> BAP, Teuchos::RCP<tpetra_matrix_type> tpetra_prolong, Teuchos::RCP<tpetra_multivector_type> BAP_shrunk, Teuchos::RCP<const Teuchos::Comm<int> > comm, int ndx, int ndy) {
  // INPUT: BAP = matrix where the uncompressed version of B_DD * A_h * Ptentative must be stored
  // INPUT: tpetra_prolong = Ptentative
  // INPUT: BAP_shrunk = Tpetra:MultiVector contatining the shrunk version of B_DD * A_h * Ptentative resulting from domain decomposition with coloring
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)
  // INPUT: ndx = number of domains along x-direction
  // INPUT: ndy = number of domains along y-direction

  Teuchos::ArrayView<const int> myLocalElements = BAP->getRowMap()->getLocalElementList();
  int mypid                                     = comm->getRank();
  int brick_id                                  = mypid;
  int shifted_id                                = brick_id - 1;

  if (mypid > 0) {
    for (int color = 0; color < 27; ++color) {
      int neighbour                            = -1;
      Teuchos::ArrayRCP<const double> localBAP = BAP_shrunk->getData(color);

      // The following if statements control the neighbours of a subdomain in a 3D brick partitioned mesh
      // Each subdomains is incorporated in a 3x3x3 cube which is sliced into 3 squares living on three different planes
      // The neighbours of a subdomain are checked plane by plane: in total there are three planes to span
      //
      // Identification of neighbours that are on the current plane
      if (coloring3D(brick_id, ndx, ndy) == color && (shifted_id) >= 0 && (shifted_id) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id);
      else if (coloring3D(brick_id - 1, ndx, ndy) == color && (shifted_id - 1) >= 0 && (shifted_id - 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - 1);
      else if (coloring3D(brick_id + 1, ndx, ndy) == color && (shifted_id + 1) >= 0 && (shifted_id + 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + 1);
      else if (coloring3D(brick_id - ndx, ndx, ndy) == color && (shifted_id - ndx) >= 0 && (shifted_id - ndx) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx);
      else if (coloring3D(brick_id - ndx - 1, ndx, ndy) == color && (shifted_id - ndx - 1) >= 0 && (shifted_id - ndx - 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx - 1);
      else if (coloring3D(brick_id - ndx + 1, ndx, ndy) == color && (shifted_id - ndx + 1) >= 0 && (shifted_id - ndx + 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx + 1);
      else if (coloring3D(brick_id + ndx, ndx, ndy) == color && (shifted_id + ndx) >= 0 && (shifted_id + ndx) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx);
      else if (coloring3D(brick_id + ndx - 1, ndx, ndy) == color && (shifted_id + ndx - 1) >= 0 && (shifted_id + ndx - 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx - 1);
      else if (coloring3D(brick_id + ndx + 1, ndx, ndy) == color && (shifted_id + ndx + 1) >= 0 && (shifted_id + ndx + 1) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx + 1);

      // Identification of the neighbours that are on the plane below
      else if (coloring3D(brick_id - ndx * ndy, ndx, ndy) == color && (shifted_id - ndx * ndy) >= 0 && (shifted_id - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx * ndy);
      else if (coloring3D(brick_id - 1 - ndx * ndy, ndx, ndy) == color && (shifted_id - 1 - ndx * ndy) >= 0 && (shifted_id - 1 - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - 1 - ndx * ndy);
      else if (coloring3D(brick_id + 1 - ndx * ndy, ndx, ndy) == color && (shifted_id + 1 - ndx * ndy) >= 0 && (shifted_id + 1 - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + 1 - ndx * ndy);
      else if (coloring3D(brick_id - ndx - ndx * ndy, ndx, ndy) == color && (shifted_id - ndx - ndx * ndy) >= 0 && (shifted_id - ndx - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx - ndx * ndy);
      else if (coloring3D(brick_id - ndx - 1 - ndx * ndy, ndx, ndy) == color && (shifted_id - ndx - 1 - ndx * ndy) >= 0 && (shifted_id - ndx - 1 - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx - 1 - ndx * ndy);
      else if (coloring3D(brick_id - ndx + 1 - ndx * ndy, ndx, ndy) == color && (shifted_id - ndx + 1 - ndx * ndy) >= 0 && (shifted_id - ndx + 1 - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx + 1 - ndx * ndy);
      else if (coloring3D(brick_id + ndx - ndx * ndy, ndx, ndy) == color && (shifted_id + ndx - ndx * ndy) >= 0 && (shifted_id + ndx - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx - ndx * ndy);
      else if (coloring3D(brick_id + ndx - 1 - ndx * ndy, ndx, ndy) == color && (shifted_id + ndx - 1 - ndx * ndy) >= 0 && (shifted_id + ndx - 1 - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx - 1 - ndx * ndy);
      else if (coloring3D(brick_id + ndx + 1 - ndx * ndy, ndx, ndy) == color && (shifted_id + ndx + 1 - ndx * ndy) >= 0 && (shifted_id + ndx + 1 - ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx + 1 - ndx * ndy);

      // Identification of the neighbours that are on the plane above
      else if (coloring3D(brick_id + ndx * ndy, ndx, ndy) == color && (shifted_id + ndx * ndy) >= 0 && (shifted_id + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx * ndy);
      else if (coloring3D(brick_id - 1 + ndx * ndy, ndx, ndy) == color && (shifted_id - 1 + ndx * ndy) >= 0 && (shifted_id - 1 + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - 1 + ndx * ndy);
      else if (coloring3D(brick_id + 1 + ndx * ndy, ndx, ndy) == color && (shifted_id + 1 + ndx * ndy) >= 0 && (shifted_id + 1 + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + 1 + ndx * ndy);
      else if (coloring3D(brick_id - ndx + ndx * ndy, ndx, ndy) == color && (shifted_id - ndx + ndx * ndy) >= 0 && (shifted_id - ndx + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx + ndx * ndy);
      else if (coloring3D(brick_id - ndx - 1 + ndx * ndy, ndx, ndy) == color && (shifted_id - ndx - 1 + ndx * ndy) >= 0 && (shifted_id - ndx - 1 + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx - 1 + ndx * ndy);
      else if (coloring3D(brick_id - ndx + 1 + ndx * ndy, ndx, ndy) == color && (shifted_id - ndx + 1 + ndx * ndy) >= 0 && (shifted_id - ndx + 1 + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id - ndx + 1 + ndx * ndy);
      else if (coloring3D(brick_id + ndx + ndx * ndy, ndx, ndy) == color && (shifted_id + ndx + ndx * ndy) >= 0 && (shifted_id + ndx + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx + ndx * ndy);
      else if (coloring3D(brick_id + ndx - 1 + ndx * ndy, ndx, ndy) == color && (shifted_id + ndx - 1 + ndx * ndy) >= 0 && (shifted_id + ndx - 1 + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx - 1 + ndx * ndy);
      else if (coloring3D(brick_id + ndx + 1 + ndx * ndy, ndx, ndy) == color && (shifted_id + ndx + 1 + ndx * ndy) >= 0 && (shifted_id + ndx + 1 + ndx * ndy) < tpetra_prolong->getGlobalNumCols())
        neighbour = (shifted_id + ndx + 1 + ndx * ndy);

      // in case neighbour>=0, it means that the current MPI processor (=subdomain) has a neighbour with the color analyzed at the current for-loop iteration
      // otherwise the current MPI processor sits on the boundary of the mesh and it has less than 26 neighbours.
      if (neighbour >= 0) {
        for (iterator_type it = myLocalElements.begin(); it != myLocalElements.end(); ++it) {
          const int i_local = *it;
          const int aux     = BAP->getRowMap()->getLocalElement(i_local);

          std::vector<int> BAP_inds;
          std::vector<double> BAP_vals;

          int aux2;

          aux2 = BAP->getColMap()->getLocalElement(neighbour);

          if (aux2 >= 0) {
            BAP_inds.emplace_back(aux2);
            BAP_vals.emplace_back(localBAP[aux]);
            BAP->insertLocalValues(aux, BAP_inds, BAP_vals);
          } else
            std::cout << "ID: " << mypid << " does not reach " << neighbour + 1 << std::endl;
        }
      }
    }
  }
}
