// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TrilinosCouplings_TpetraIntrepidPoissonExample_hpp
#define __TrilinosCouplings_TpetraIntrepidPoissonExample_hpp

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "TrilinosCouplings_config.h"

#ifdef HAVE_TRILINOSCOUPLINGS_MUELU
#  include "MueLu_TpetraOperator.hpp"
#endif // HAVE_TRILINOSCOUPLINGS_MUELU

#include "TrilinosCouplings_IntrepidPoissonExample_SolveWithBelos.hpp"

namespace TrilinosCouplings {
/// \namespace TpetraIntrepidPoissonExample
/// \brief Tpetra version of the Intrepid Poisson test problem example.
///
/// The Intrepid Poisson test problem uses Pamgen to construct a 3-D
/// mesh (a simple rectangular prism with hex elements) in parallel,
/// Sacado automatic differentiation to construct a right-hand side of
/// the PDE corresponding to a given exact solution, and Intrepid to
/// build a discretization.
///
/// We provide two variants of the Intrepid Poisson test: one that
/// fills Epetra objects, and one that fills Tpetra objects.  The two
/// variants do exactly the same things otherwise, so you can use them
/// to compare the performance of Epetra and Tpetra fill.
///
/// This namespace contains the Tpetra variant.  It defines typedefs
/// which you can use when writing a main() driver to run the test.
/// The makeMatrixAndRightHandSide() function does all the work.  You
/// can use the exactResidualNorm() function to test correctness of
/// the discretization.  In particular, if the continuous exact
/// solution was chosen from the space of finite element polynomials,
/// the exact solution of the discrete linear system AX=B should match
/// the continuous exact solution exactly, modulo rounding error when
/// assembling the discretization.
///
/// The solveWithBelos() function solves the given linear system using
/// a Belos iterative solver.  You can provide a left and/or right
/// preconditioner if you want.
namespace TpetraIntrepidPoissonExample {

//
// mfh 19 Apr 2012: Leave these typedefs up here for use by main() and
// the other functions.  This example probably only works for ST =
// double and LO,GO = int, but it probably works for other Node types
// besides the default.
//
typedef double ST;
typedef int    LO;
typedef Tpetra::Map<>::global_ordinal_type GO;
typedef Tpetra::Map<>::node_type           Node;

//
// mfh 19 Apr 2012: If you want to change the template parameters of
// these typedefs, modify the typedefs (ST, LO, GO, Node) above.
//
typedef Tpetra::CrsMatrix<ST, LO, GO, Node>    sparse_matrix_type;
typedef Tpetra::Operator<ST, LO, GO, Node>     operator_type;
typedef Tpetra::MultiVector<ST, LO, GO, Node>  multivector_type;
typedef Tpetra::Vector<ST, LO, GO, Node>       vector_type;

/// \brief Create the mesh and build the linear system to solve.
///
/// \param A [out] The sparse matrix.
/// \param B [out] The right-hand side(s).
/// \param X_exact [out] The exact solution of the PDE, projected onto
///   the discrete mesh.  This may not necessarily be the same as the
///   exact solution of the discrete linear system.
/// \param X [out] The approximate solution(s).
/// \param err [out] Output stream for errors.
/// \param out [out] Output stream for verbose output.
/// \param comm [in] Communicator.
/// \param node [in/out] Kokkos Node instance.
/// \param meshInput [in] Pamgen mesh specification string.
///
/// Pamgen is a parallel mesh generation library (that nevertheless
/// performs no communication, since it limits itself to generate
/// simple meshes).  Here is its technical report:
///
/// @techreport{hensinger2008pamgen,
/// author = "David M. Hensinger and Richard R. Drake and James G. Foucar and Thomas A. Gardiner",
/// title = "Pamgen, a Library for Parallel Generation of Simple Finite Element Meshes",
/// institution = "Sandia National Laboratories",
/// number = "SAND2008-1933",
/// month = "April",
/// year = "2008"
/// }
///
/// Quote from its abstract: "\textsc{Pamgen} is a parallel mesh
/// generation library that allows on-the-fly scalable generation of
/// hexahedral and quadrilateral finite element meshes for several
/// simple geometries.  It has been used to generate more than 1.1
/// billion elements on 17,576 processors."
///
/// Pamgen takes a string of commands as input.  Read the
/// "Poisson.xml" file in this directory for an example.  The example
/// uses a "rectilinear" mesh, which has the following fields:
///
/// nx: Number of cells in X direction
/// ny: Number of cells in Y direction
/// nz: Number of cells in Z direction
/// bx: Number of blocks in X direction
/// by: Number of blocks in Y direction
/// bz: Number of blocks in Z direction
/// gmin: Minimum domain coordinates x y z
/// gmax: Maximum domain coordinates x y z
///
/// Poisson.xml specifies a cube that lives in the unit positive
/// octant (that is, (x,y,z) with 0 <= x,y,z <= 1), with 20 cells
/// along each dimension (thus, 20*20*20 = 8000 cells total).  If you
/// want to make the problem bigger, you can scale up nx, ny, and/or
/// nz.  It's probably better for communication balance to scale them
/// equally.
///
/// The "set assign ... end" statement specifies boundaries
/// ("nodesets" and "sidesets") of the domain.  The Poisson.xml
/// example names the exterior faces of the cube with IDs.
///
/// The Poisson.xml example does not include a (parallel)
/// "decomposition strategy ... end" statement.  The default strategy
/// is "bisection," which attempts to assign an equal number of cells
/// to each parallel process.
void
makeMatrixAndRightHandSide (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<vector_type>& B,
                            Teuchos::RCP<vector_type>& X_exact,
                            Teuchos::RCP<vector_type>& X,
                            Teuchos::RCP<multivector_type> & coords,
                            Teuchos::RCP<vector_type>& node_sigma,
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const std::string& meshInput,
                            Teuchos::ParameterList & inputList,
                            Teuchos::ParameterList & problemStatistics,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose = false,
                            const bool debug = false);

//! Just like above, but with multivector_type output arguments.
void
makeMatrixAndRightHandSide (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<multivector_type>& B,
                            Teuchos::RCP<multivector_type>& X_exact,
                            Teuchos::RCP<multivector_type>& X,
                            Teuchos::RCP<multivector_type> & coords,
                            Teuchos::RCP<vector_type>& node_sigma,
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const std::string& meshInput,
                            Teuchos::ParameterList & inputList,
                            Teuchos::ParameterList & problemStatistics,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose = false,
                            const bool debug = false);

//! Return \f$\|B - A X_{\text{exact}}\|_2\f$, \f$\|B\|\f$, and \f$\|A\|_F\f$.
std::vector<Teuchos::ScalarTraits<ST>::magnitudeType>
exactResidualNorm (const Teuchos::RCP<const sparse_matrix_type>& A,
                   const Teuchos::RCP<const vector_type>& B,
                   const Teuchos::RCP<const vector_type>& X_exact);

/// \brief Solve the linear system(s) AX=B with Belos.
///
/// X and B are both multivectors, meaning that you may ask Belos to
/// solve more than one linear system at a time.  X and B must have
/// the same number of columns.
///
/// This interface will change in the future to accept the name of the
/// Belos solver to use.  For now, the solver is hard-coded to
/// Pseudoblock CG (implemented by Belos::PseudoBlockCGSolMgr).
///
/// \param converged [out] Whether Belos reported that the iterative
///   method converged (solved the linear system to the desired
///   tolerance).
///
/// \param numItersPerformed [out] Number of iterations that the Belos
///   solver performed.
///
/// \param solverName [in] Which iterative linear solver to use.
///   Any name that Belos::SolverFactory knows will work here.
///
/// \param tol [in] Convergence tolerance for the iterative method.
///   The meaning of this depends on the particular iterative method.
///
/// \param maxNumIters [in] Maximum number of iterations that the
///   iterative method should perform, regardless of whether it
///   converged.
///
/// \param num_steps [in] Number of "time steps", i.e., the number of
//    times the solver is called in a fake time-step loop.
///
/// \param X [in/out] On input: the initial guess(es) for the iterative
///   method.  On output: the computed approximate solution.
///
/// \param A [in] The matrix in the linear system(s) AX=B to solve.
///
/// \param B [in] The right-hand side(s) in the linear system AX=B to solve.
///
/// \param M_left [in] If nonnull, a left preconditioner that the
///   iterative method may use.  If null, the iterative method will
///   not use a left preconditioner.
///
/// \param M_right [in] If nonnull, a right preconditioner that the
///   iterative method may use.  If null, the iterative method will
///   not use a right preconditioner.
void
solveWithBelos (bool& converged,
                int& numItersPerformed,
                const std::string& solverName,
                const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
                const int maxNumIters,
                const int num_steps,
                const Teuchos::RCP<multivector_type>& X,
                const Teuchos::RCP<const sparse_matrix_type>& A,
                const Teuchos::RCP<const multivector_type>& B,
                const Teuchos::RCP<const operator_type>& M_left=Teuchos::null,
                const Teuchos::RCP<const operator_type>& M_right=Teuchos::null);

/// \brief Solve the linear system(s) AX=B with Belos on the GPU.
///
/// In addition to the parameters taken by solveWithBelos(), this function
/// additionally takes the following parameters to determine which GPU
/// device should be attached to which MPI rank:
///
/// \param ranks_per_node [in]  Number of MPI ranks per node.  For OpenMPI and
///   MVAPICH MPI libraries, this is determined automatically and thus needn't
///   be specified.  For other MPI libraries it must be supplied by the user.
///
/// \param gpu_ranks_per_node [in] Number of MPI ranks per node associated
///   with GPUs.  This must be less than or equal to the number of devices
///   per node.
///
/// \param device_offset [in] Offset for computing the CUDA device ID for
///   each MPI rank associated with a GPU.  If there are n GPU's per node,
///   the first n MPI ranks for each node will be associated with a GPU starting
///   at device number device_offset.
///
/// \param prec_type [in] The preconditioner type (e.g., "MueLu").
void
solveWithBelosGPU (
  bool& converged,
  int& numItersPerformed,
  const Teuchos::ScalarTraits<ST>::magnitudeType& tol,
  const int maxNumIters,
  const int num_steps,
  const int ranks_per_node,
  const int gpu_ranks_per_node,
  const int device_offset,
  const std::string& prec_type,
  const Teuchos::RCP<multivector_type>& X,
  const Teuchos::RCP<const sparse_matrix_type>& A,
  const Teuchos::RCP<const multivector_type>& B,
  const Teuchos::RCP<const operator_type>& M_left=Teuchos::null,
  const Teuchos::RCP<const operator_type>& M_right=Teuchos::null);

} // namespace TpetraIntrepidPoissonExample
} // namespace TrilinosCouplings

#endif // __TrilinosCouplings_TpetraIntrepidPoissonExample_hpp
