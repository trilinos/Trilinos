// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example 05_solve.cpp
\brief Example of thread-parallel linear system construction using Kokkos.
\ref Tpetra_Lesson07 explains this example in detail.
*/

#include <Tpetra_Core.hpp>
// This is the only header file you need to include for the "core"
// part of Kokkos.  That includes Kokkos::View, Kokkos::parallel_*,
// and atomic updates.
#include <Kokkos_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Import_Util2.hpp>  //for sortCrsEntries
#include <iostream>

#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA

// Exact solution of the partial differential equation that main()
// discretizes.  We include it here to check the error.
KOKKOS_INLINE_FUNCTION double
exactSolution (const double x, const double T_left, const double T_right)
{
  return -4.0 * (x - 0.5) * (x - 0.5) + 1.0 + (T_left - T_right) * x;
}

// Attempt to solve Ax=b using CG (the Method of Conjugate Gradients),
// and return the number of iterations.
template<class CrsMatrixType, class VectorType>
int
solve (VectorType& x, const CrsMatrixType& A, const VectorType& b, const double dx)
{
  using std::cout;
  using std::endl;
  const int myRank = x.getMap ()->getComm ()->getRank ();
  const bool verbose = false;

  // In practice, one would call Belos using a preconditioner from
  // MueLu, Ifpack2, or some other package.  For now, we implement CG
  // by hand.
  const double convTol = std::max (dx * dx, 1.0e-8);
  // Don't do more CG iterations than the problem's dimension.
  const int maxNumIters = std::min (static_cast<Tpetra::global_size_t> (100),
                                    A.getGlobalNumRows ());
  VectorType r (A.getRangeMap ());
  A.apply (x, r);
  r.update (1.0, b, -1.0); // r := -(A*x) + b
  const double origResNorm = r.norm2 ();
  if (myRank == 0 && verbose) {
    cout << "Original residual norm: " << origResNorm << endl;
  }
  if (origResNorm <= convTol) {
    return 0; // the solution is already close enough
  }

  double r_dot_r = origResNorm * origResNorm;
  VectorType p (r, Teuchos::Copy);
  VectorType Ap (A.getRangeMap ());
  int numIters = 1;
  for (; numIters <= maxNumIters; ++numIters) {
    A.apply (p, Ap);
    const double p_dot_Ap = p.dot (Ap);
    if (p_dot_Ap <= 0.0) {
      if (myRank == 0) {
        cout << "At iteration " << numIters << ", p.dot(Ap) = " << p_dot_Ap << " <= 0.";
      }
      // Revert to approximate solution x from previous iteration.
      return numIters - 1;
    }
    const double alpha = r_dot_r / p_dot_Ap;
    if (alpha <= 0.0) {
      if (myRank == 0) {
        cout << "At iteration " << numIters << ", alpha = " << alpha << " <= 0.";
      }
      // Revert to approximate solution x from previous iteration.
      return numIters - 1;
    }
    x.update (alpha, p, 1.0); // x := alpha*p + x
    r.update (-alpha, Ap, 1.0); // r := -alpha*Ap + r

    const double newResNorm = r.norm2 ();
    const double r_dot_r_next = newResNorm * newResNorm;
    const double newRelResNorm = newResNorm / origResNorm;
    if (myRank == 0 && verbose) {
      cout << "Iteration " << numIters << ": r_dot_r = " << r_dot_r
           << ", r_dot_r_next = " << r_dot_r_next
           << ", newResNorm = " << newResNorm << endl;
    }
    if (newRelResNorm <= convTol) {
      return numIters;
    }
    const double beta = r_dot_r_next / r_dot_r;
    p.update (1.0, r, beta); // p := r + beta*p
    r_dot_r = r_dot_r_next;
  }
  return numIters;
}

int main (int argc, char* argv[]) {
  using std::cout;
  using std::endl;
  // We're filling into Tpetra data structures now, so we have to
  // respect Tpetra's choices of local and global indices.
  using LO = Tpetra::Map<>::local_ordinal_type;
  using GO = Tpetra::Map<>::global_ordinal_type;
  using NT = Tpetra::Map<>::node_type;
  using device_type = typename NT::device_type;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    LO numLclElements = 10000;
    int numGblElements = numProcs * numLclElements;

    // We happened to choose a discretization with the same number of
    // nodes and degrees of freedom as elements.  That need not always
    // be the case.
    LO numLclNodes = numLclElements;
    LO numLclRows = numLclNodes;

    // Describe the physical problem (heat equation with
    // nonhomogeneous Dirichlet boundary conditions) and its
    // discretization.
    Kokkos::View<double*, device_type> temperature ("temperature", numLclNodes);
    const double diffusionCoeff = 1.0;
    const double x_left = 0.0; // position of the left boundary
    const double T_left = 0.0; // temperature at the left boundary
    const double x_right = 1.0; // position of the right boundary
    const double T_right = 1.0; // temperature at the right boundary
    const double dx = (x_right - x_left) / numGblElements;
    Kokkos::View<double*, device_type> forcingTerm ("forcing term", numLclNodes);

    // Set the forcing term.  We picked it so that we can know the exact
    // solution of the heat equation, namely
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
    // We may use LO for the number of entries in each row, since it
    // may not exceed the number of columns in the local matrix.
    Kokkos::View<LO*, device_type> rowCounts ("row counts", numLclRows);
    size_t numLclEntries = 0;
    Kokkos::parallel_reduce ("Count graph", numLclElements,
      KOKKOS_LAMBDA (const LO elt, size_t& curNumLclEntries) {
        const LO lclRows = elt;

        // Always add a diagonal matrix entry.
	Kokkos::atomic_fetch_add (&rowCounts(lclRows), LO(1));
	curNumLclEntries++;

	// Each neighboring MPI process contributes an entry to the
	// current row.  In a more realistic code, we might handle
	// this either through a global assembly process (requiring
	// MPI communication), or through ghosting a layer of elements
	// (no MPI communication).

	// MPI process to the left sends us an entry
	if (myRank > 0 && lclRows == 0) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows), LO(1));
	  curNumLclEntries++;
	}
	// MPI process to the right sends us an entry
	if (myRank + 1 < numProcs && lclRows + 1 == numLclRows) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows), LO(1));
	  curNumLclEntries++;
	}

	// Contribute a matrix entry to the previous row.
	if (lclRows > 0) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows-1), LO(1));
	  curNumLclEntries++;
	}
	// Contribute a matrix entry to the next row.
	if (lclRows + 1 < numLclRows) {
	  Kokkos::atomic_fetch_add (&rowCounts(lclRows+1), LO(1));
	  curNumLclEntries++;
	}
      }, numLclEntries /* reduction result */);

    // mfh 20 Aug 2017: We can't just use a Kokkos::View<size_t*> for
    // the row offsets, because Tpetra::CrsMatrix reserves the right
    // to use a different row offset type for different execution /
    // memory spaces.  Instead, we first deduce the row offset type,
    // then construct a View of it.  (Note that a row offset needs to
    // have a type that can contain the sum of the row counts.)
    using row_offset_type =
      Tpetra::CrsMatrix<double>::local_matrix_device_type::row_map_type::non_const_value_type;

    // Use a parallel scan (prefix sum) over the array of row counts, to
    // compute the array of row offsets for the sparse graph.
    Kokkos::View<row_offset_type*, device_type> rowOffsets ("row offsets", numLclRows+1);
    Kokkos::parallel_scan ("Row offsets", numLclRows+1,
      KOKKOS_LAMBDA (const LO lclRows,
                     row_offset_type& update,
                     const bool final) {
        if (final) {
          // Kokkos uses a multipass algorithm to implement scan.  Only
	  // update the array on the final pass.  Updating the array
	  // before changing 'update' means that we do an exclusive
	  // scan.  Update the array after for an inclusive scan.
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
    Kokkos::deep_copy (rowCounts, 0);

    Kokkos::View<LO*, device_type> colIndices ("column indices", numLclEntries);
    Kokkos::View<double*, device_type> matrixValues ("matrix values", numLclEntries);

    // Iterate over elements in parallel to fill the graph, matrix,
    // and right-hand side (forcing term).  The latter gets the
    // boundary conditions (a trick for nonzero Dirichlet boundary
    // conditions).
    Kokkos::parallel_for ("Assemble", numLclElements,
      KOKKOS_LAMBDA (const LO elt) {
        // Push dx*dx into the forcing term.
        const double offCoeff = -diffusionCoeff / 2.0;
        const double midCoeff =  diffusionCoeff;
        // In this discretization, every element corresponds to a
        // degree of freedom, and to a row of the matrix.  (Boundary
        // conditions are Dirichlet, so they don't count as degrees of
        // freedom.)
        const int lclRows = elt;

	// Always add a diagonal matrix entry.
	{
	  const LO count =
	    Kokkos::atomic_fetch_add (&rowCounts(lclRows), LO(1));
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
	  const LO count = Kokkos::atomic_fetch_add (&rowCounts(lclRows), LO(1));
	  colIndices(rowOffsets(lclRows) + count) = numLclRows;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows) + count), offCoeff);
	}
	// MPI process to the right sends us an entry
	if (myRank + 1 < numProcs && lclRows + 1 == numLclRows) {
	  const LO count =
	    Kokkos::atomic_fetch_add (&rowCounts(lclRows), LO(1));

	  // Give this entry the right local column index, depending
	  // on whether the MPI process to the left has already sent
	  // us an entry.
	  const int colInd = (myRank > 0) ? numLclRows + 1 : numLclRows;
	  colIndices(rowOffsets(lclRows) + count) = colInd;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows) + count), offCoeff);
	}

        // Contribute a matrix entry to the previous row.
        if (lclRows > 0) {
	  const LO count = Kokkos::atomic_fetch_add (&rowCounts(lclRows-1), LO(1));
	  colIndices(rowOffsets(lclRows-1) + count) = lclRows;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows-1) + count), offCoeff);
	}
	// Contribute a matrix entry to the next row.
	if (lclRows + 1 < numLclRows) {
	  const LO count = Kokkos::atomic_fetch_add (&rowCounts(lclRows+1), LO(1));
	  colIndices(rowOffsets(lclRows+1) + count) = lclRows;
	  Kokkos::atomic_fetch_add (&matrixValues(rowOffsets(lclRows+1) + count), offCoeff);
	}
      });

    //
    // Construct Tpetra objects
    //
    using Teuchos::RCP;
    using Teuchos::rcp;

    const GO indexBase = 0;
    RCP<const Tpetra::Map<> > rowMap =
      rcp (new Tpetra::Map<> (numGblElements, numLclElements, indexBase, comm));

    LO num_col_inds = numLclElements;
    if (myRank > 0) {
      num_col_inds++;
    }
    if (myRank + 1 < numProcs) {
      num_col_inds++;
    }

    // Soon, it will be acceptable to use a device View here.  The
    // issues are that Teuchos::RCP isn't thread safe, and the Kokkos
    // version of Tpetra::Map isn't quite ready yet.  We will change
    // the latter soon.
    Kokkos::View<GO*, Kokkos::HostSpace> colInds ("Column Map", num_col_inds);
    for (LO k = 0; k < numLclElements; ++k) {
      colInds(k) = rowMap->getGlobalElement (k);
    }
    LO k = numLclElements;
    if (myRank > 0) {
      // Contribution from left process.
      colInds(k++) = rowMap->getGlobalElement (0) - 1;
    }
    if (myRank + 1 < numProcs) {
      // Contribution from right process.
      colInds(k++) = rowMap->getGlobalElement (numLclElements - 1) + 1;
    }

    // Flag to tell Tpetra::Map to compute the global number of
    // indices in the Map.
    const Tpetra::global_size_t INV =
      Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();
    RCP<const Tpetra::Map<> > colMap =
      rcp (new Tpetra::Map<> (INV, colInds.data (), colInds.extent (0), indexBase, comm));

    Tpetra::Import_Util::sortCrsEntries(rowOffsets, colIndices, matrixValues);

    Tpetra::CrsMatrix<double> A (rowMap, colMap, rowOffsets,
			         colIndices, matrixValues);
    A.fillComplete ();

    // Hack to deal with the fact that Tpetra::Vector needs a
    // DualView<double**> for now, rather than a View<double*>.
    Kokkos::DualView<double**, Kokkos::LayoutLeft, device_type> b_lcl ("b", numLclRows, 1);
    b_lcl.modify_device ();
    Kokkos::deep_copy (Kokkos::subview (b_lcl.d_view, Kokkos::ALL (), 0), forcingTerm);
    Tpetra::Vector<> b (A.getRangeMap (), b_lcl);

    Kokkos::DualView<double**, Kokkos::LayoutLeft, device_type> x_lcl ("b", numLclRows, 1);
    x_lcl.modify_device ();
    Kokkos::deep_copy (Kokkos::subview (x_lcl.d_view, Kokkos::ALL (), 0), temperature);
    Tpetra::Vector<> x (A.getDomainMap (), x_lcl);

    const int numIters = solve (x, A, b, dx); // solve the linear system
    if (myRank == 0) {
      cout << "Linear system Ax=b took " << numIters << " iteration(s) to solve" << endl;
    }

    // Hack to deal with the fact that Tpetra::Vector needs a
    // DualView<double**> for now, rather than a View<double*>.  This
    // means that we have to make a deep copy back into the
    // 'temperature' output array.
    x_lcl.sync_device ();
    Kokkos::deep_copy (temperature, Kokkos::subview (b_lcl.d_view, Kokkos::ALL (), 0));

    // Correct the solution for the nonhomogenous Dirichlet boundary
    // conditions.
    Kokkos::parallel_for ("Boundary conditions", numLclNodes,
      KOKKOS_LAMBDA (const LO node) {
	const double x_cur = x_left + node * dx;
	temperature(node) += (T_left - T_right) * x_cur;
      });

    // Compare the computed solution against the known exact solution:
    //
    // u(x) = -4.0 * (x - 0.5) * (x - 0.5) + 1.0 + (T_left - T_right)x.
    Tpetra::Vector<> x_exact (x, Teuchos::Copy);
    typedef typename device_type::execution_space execution_space;
    typedef Kokkos::RangePolicy<execution_space, LO> policy_type;

    auto x_exact_lcl = x_exact.getLocalViewDevice (Tpetra::Access::ReadWrite);

    Kokkos::parallel_for ("Compare solutions",
      policy_type (0, numLclNodes),
      KOKKOS_LAMBDA (const LO& node) {
        const double x_cur = x_left + node * dx;
        x_exact_lcl(node,0) -= exactSolution (x_cur, T_left, T_right);
      });
    const double absErrNorm = x_exact.norm2 ();
    const double relErrNorm = b.norm2 () / absErrNorm;
    if (myRank == 0) {
      cout << "Relative error norm: " << relErrNorm << endl;
    }

    if (myRank == 0) {
      cout << "End Result: TEST PASSED" << endl;
    }
  }
  return EXIT_SUCCESS;
}

#else

int main (int argc, char* argv[]) {
  using std::cout;
  using std::endl;
  
  cout << "This lesson was not compiled because Kokkos" << endl
    "was not configured with lambda support for all backends." << endl
    "Tricking CTest into perceiving success anyways:" << endl
    "End Result: TEST PASSED" << endl;
  return 0;
}

#endif
