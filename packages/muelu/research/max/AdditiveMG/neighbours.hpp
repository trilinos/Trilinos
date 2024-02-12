//#include <Tpetra_CrsMatrix_decl.hpp>
//#include <Teuchos_RCPDecl.hpp>

#include "coloring.hpp"

typedef Tpetra::CrsMatrix<double, int, int, KokkosClassic::DefaultNode::DefaultNodeType> tpetra_matrix_type;
typedef Tpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, int, int, KokkosClassic::DefaultNode::DefaultNodeType> tpetra_multivector_type;
typedef typename Teuchos::ArrayView<const int>::const_iterator iterator_type;
typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType node_type2;

void neighbours1D(Teuchos::RCP<tpetra_matrix_type> tpetra_prolong, std::vector<int>& neighbours, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  // INPUT: tpetra_prolong is the tentative prolongator, here it is just used to get the Total number of aggregates given by the nubmer of columns
  // INPUT: neighoubrs is a vector where the neoghbours of the current MPI processor (=subdomain) are stored
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)

  neighbours.clear();

  int mypid = comm->getRank();

  for (int color = 0; color < 3; ++color) {
    if ((mypid - 1) % 3 == color && (mypid - 1) >= 0 && (mypid - 1) < tpetra_prolong->getGlobalNumCols())
      neighbours.emplace_back(mypid - 1);
    else if ((mypid - 2) % 3 == color && (mypid - 2) >= 0 && (mypid - 2) < tpetra_prolong->getGlobalNumCols())
      neighbours.emplace_back(mypid - 2);
    else if ((mypid) % 3 == color && (mypid) >= 0 && (mypid) < tpetra_prolong->getGlobalNumCols())
      neighbours.emplace_back(mypid);
  }

  neighbours.shrink_to_fit();
}

void neighbours2D(Teuchos::RCP<tpetra_matrix_type> tpetra_prolong, std::vector<int>& neighbours, Teuchos::RCP<const Teuchos::Comm<int> > comm, int ndx) {
  // INPUT: tpetra_prolong is the tentative prolongator, here it is just used to get the Total number of aggregates given by the nubmer of columns
  // INPUT: neighoubrs is a vector where the neoghbours of the current MPI processor (=subdomain) are stored
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)
  // INPUT: ndx is the number of domains along the x-direction

  neighbours.clear();
  int mypid      = comm->getRank();
  int brick_id   = mypid;
  int shifted_id = brick_id - 1;

  if (mypid > 0) {
    for (int color = 0; color < 9; ++color) {
      int neighbour = -1;

      // The following if statements control the neighbours of a subdomain in a 2D brick partitioned mesh
      // Each subdomains is incorporated in a 3x3 square which is sliced into 3 stripes
      // The neighbours of a subdomain are checked plane by plane: in total there are three planes to span
      //  In case the subdomains sit on a boundary, there are missing neighbours for a specific color
      //  (this is what the last two conditions of each if statement take care of)
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

      if (neighbour >= 0)
        neighbours.emplace_back(neighbour);
    }
  }

  neighbours.shrink_to_fit();
}

void neighbours3D(Teuchos::RCP<tpetra_matrix_type> tpetra_prolong, std::vector<int>& neighbours, Teuchos::RCP<const Teuchos::Comm<int> > comm, int ndx, int ndy) {
  // INPUT: tpetra_prolong is the tentative prolongator, here it is just used to get the Total number of aggregates given by the nubmer of columns
  // INPUT: neighoubrs is a vector where the neoghbours of the current MPI processor (=subdomain) are stored
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)
  // INPUT: ndx is the number of domains along the x-direction
  // INPUT: ndy is the nubmer of domains along the y-direction

  neighbours.clear();
  int mypid      = comm->getRank();
  int brick_id   = mypid;
  int shifted_id = brick_id - 1;

  if (mypid > 0) {
    for (int color = 0; color < 27; ++color) {
      int neighbour = -1;

      // The following if statements control the neighbours of a subdomain in a 3D brick partitioned mesh
      // Each subdomains is incorporated in a 3x3x3 cube which is sliced into 3 squares living on three different planes
      // The neighbours of a subdomain are checked plane by plane: in total there are three planes to span
      //  In case the subdomains sit on a boundary, there are missing neighbours for a specific color
      //  (this is what the last two conditions of each if statement take care of)
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

      if (neighbour >= 0)
        neighbours.emplace_back(neighbour);
    }
  }

  neighbours.shrink_to_fit();
}
