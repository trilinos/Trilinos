// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Core.hpp>
// This is the only header file you need to include for the "core"
// part of Kokkos.  That includes Kokkos::View, Kokkos::parallel_*,
// and atomic updates.
#include <Kokkos_Core.hpp>

int main (int argc, char* argv[]) {
  using std::cout;
  using std::endl;
  // Once we fill into Tpetra data structures, we will need to respect
  // Tpetra's choices of local (LO) and global (GO) indices.  For now,
  // we define these to int.
  typedef int LO;
  typedef int GO;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    LO numLclElements = 10000;
    GO numGblElements = numProcs * numLclElements;

    // We happened to choose a discretization with the same number of
    // nodes and degrees of freedom as elements.  That need not always
    // be the case.
    LO numLclNodes = numLclElements;
    LO numLclRows = numLclNodes;

    // Describe the physical problem (heat equation with
    // nonhomogeneous Dirichlet boundary conditions) and its
    // discretization.
    Kokkos::View<double*> temperature ("temperature", numLclNodes);
    const double diffusionCoeff = 1.0;
    const double x_left = 0.0; // position of the left boundary
    const double T_left = 0.0; // temperature at the left boundary
    const double x_right = 1.0; // position of the right boundary
    const double T_right = 1.0; // temperature at the right boundary
    const double dx = (x_right - x_left) / numGblElements;
    Kokkos::View<double*> forcingTerm ("forcing term", numLclNodes);

    // Set the forcing term.  We picked it so that we can know the
    // exact solution of the heat equation, namely
    //
    // u(x) = -4.0 * (x - 0.5) * (x - 0.5) + 1.0 + (T_left - T_right)x.
    Kokkos::parallel_for ("Set forcing term", numLclNodes,
      KOKKOS_LAMBDA (const LO node) {
	//const double x = x_left + node * dx;

	forcingTerm(node) = -8.0;
	// We multiply dx*dx into the forcing term, so the matrix's
	// entries don't need to know it.
	forcingTerm(node) *= dx * dx;
      });

    // Do a reduction over local elements to count the total number of
    // (local) entries in the graph.  While doing so, count the number
    // of (local) entries in each row, using Kokkos' atomic updates.
    Kokkos::View<size_t*> rowCounts ("row counts", numLclRows);
    size_t numLclEntries = 0;
    Kokkos::parallel_reduce ("Count graph", numLclElements,
      KOKKOS_LAMBDA (const LO elt, size_t& curNumLclEntries) {
        const LO lclRows = elt;

	// Always add a diagonal matrix entry.
	Kokkos::atomic_fetch_add (&rowCounts(lclRows), size_t(1));
	curNumLclEntries++;

	// Each neighboring MPI process contributes an entry to the
	// current row.  In a more realistic code, we might handle
	// this either through a global assembly process (requiring
	// MPI communication), or through ghosting a layer of elements
	// (no MPI communication).

	// MPI process to the left sends us an entry
	if (myRank > 0 && lclRows == 0) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows), size_t(1));
	  curNumLclEntries++;
	}
	// MPI process to the right sends us an entry
	if (myRank + 1 < numProcs && lclRows + 1 == numLclRows) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows), size_t(1));
	  curNumLclEntries++;
	}

	// Contribute a matrix entry to the previous row.
	if (lclRows > 0) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows-1), size_t(1));
	  curNumLclEntries++;
	}
	// Contribute a matrix entry to the next row.
	if (lclRows + 1 < numLclRows) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows+1), size_t(1));
	  curNumLclEntries++;
	}
      }, numLclEntries /* reduction result */);

    // Use a parallel scan (prefix sum) over the array of row counts,
    // to compute the array of row offsets for the sparse graph.
    Kokkos::View<size_t*> rowOffsets ("row offsets", numLclRows+1);
    Kokkos::parallel_scan ("Row offsets", numLclRows+1,
      KOKKOS_LAMBDA (const LO lclRows, size_t& update, const bool final) {
	if (final) {
	  // Kokkos uses a multipass algorithm to implement scan.
	  // Only update the array on the final pass.  Updating the
	  // array before changing 'update' means that we do an
	  // exclusive scan.  Update the array after for an inclusive
	  // scan.
	  rowOffsets[lclRows] = update;
	}
	if (lclRows < numLclRows) {
	  update += rowCounts(lclRows);
	}
      });

    // Use the array of row counts to keep track of where to put each
    // new column index, when filling the graph.  Updating the entries
    // of rowCounts atomically lets us parallelize over elements
    // (which may touch multiple rows at a time -- esp. in 2-D or 3-D,
    // or with higher-order discretizations), rather than rows.
    //
    // We leave as an exercise to the reader how to use this array
    // without resetting its entries.
    Kokkos::deep_copy (rowCounts, static_cast<size_t> (0));

    Kokkos::View<LO*> colIndices ("column indices", numLclEntries);
    Kokkos::View<double*> matrixValues ("matrix values", numLclEntries);

    // Iterate over elements in parallel to fill the graph, matrix,
    // and right-hand side (forcing term).  The latter gets the
    // boundary conditions (a trick for nonzero Dirichlet boundary
    // conditions).
    Kokkos::parallel_for ("Assemble", numLclElements,
      KOKKOS_LAMBDA (const LO elt) {
	// We multiply dx*dx into the forcing term, so the matrix's
	// entries don't need to know it.
	const double offCoeff = -diffusionCoeff / 2.0;
	const double midCoeff =  diffusionCoeff;
	// In this discretization, every element corresponds to a
	// degree of freedom, and to a row of the matrix.  (Boundary
	// conditions are Dirichlet, so they don't count as degrees of
	// freedom.)
	const int lclRows = elt;

	// Always add a diagonal matrix entry.
	{
	  const size_t count = Kokkos::atomic_fetch_add (&rowCounts(lclRows), size_t(1));
	  colIndices(rowOffsets(lclRows) + count) = lclRows;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows) + count), midCoeff);
	}

	// Each neighboring MPI process contributes an entry to the
	// current row.  In a more realistic code, we might handle
	// this either through a global assembly process (requiring
	// MPI communication), or through ghosting a layer of elements
	// (no MPI communication).

	// MPI process to the left sends us an entry
	if (myRank > 0 && lclRows == 0) {
	  const size_t count = Kokkos::atomic_fetch_add (&rowCounts(lclRows), size_t(1));
	  colIndices(rowOffsets(lclRows) + count) = numLclRows;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows) + count), offCoeff);
	}
	// MPI process to the right sends us an entry
	if (myRank + 1 < numProcs && lclRows + 1 == numLclRows) {
	  const size_t count = Kokkos::atomic_fetch_add (&rowCounts(lclRows), size_t(1));

	  // Give this entry the right local column index, depending
	  // on whether the MPI process to the left has already sent
	  // us an entry.
	  const int colInd = (myRank > 0) ? numLclRows + 1 : numLclRows;
	  colIndices(rowOffsets(lclRows) + count) = colInd;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows) + count), offCoeff);
	}

	// Contribute a matrix entry to the previous row.
	if (lclRows > 0) {
	  const size_t count = Kokkos::atomic_fetch_add (&rowCounts(lclRows-1), size_t(1));
	  colIndices(rowOffsets(lclRows-1) + count) = lclRows;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows-1) + count), offCoeff);
	}
	// Contribute a matrix entry to the next row.
	if (lclRows + 1 < numLclRows) {
	  const size_t count = Kokkos::atomic_fetch_add (&rowCounts(lclRows+1), size_t(1));
	  colIndices(rowOffsets(lclRows+1) + count) = lclRows;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows+1) + count), offCoeff);
	}
      });

    //
    // TODO put the discretization into Tpetra data structures, solve
    // the linear system, and correct the solution for the
    // nonhomogenous Dirichlet boundary conditions.
    //
  }
}
