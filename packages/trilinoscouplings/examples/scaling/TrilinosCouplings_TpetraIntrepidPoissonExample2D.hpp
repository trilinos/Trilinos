// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

#ifndef __TrilinosCouplings_TpetraIntrepidPoissonExample2D_hpp
#define __TrilinosCouplings_TpetraIntrepidPoissonExample2D_hpp

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
typedef int    GO;
typedef Tpetra::Map<>::node_type  Node;

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
/// "Poisson2D.xml" file in this directory for an example.  The example
/// uses a "rectilinear" mesh, which has the following fields:
///
/// nx: Number of cells in X direction
/// ny: Number of cells in Y direction
/// bx: Number of blocks in X direction
/// by: Number of blocks in Y direction
/// gmin: Minimum domain coordinates x y
/// gmax: Maximum domain coordinates x y
///
/// Poisson.xml specifies a cube that lives in the unit positive
/// octant (that is, (x,y,z) with 0 <= x,y <= 1), with 20 cells
/// along each dimension (thus, 20*20 = 400 cells total).  If you
/// want to make the problem bigger, you can scale up nx and/or ny.
/// It's probably better for communication balance to scale them
/// equally.
///
/// The "set assign ... end" statement specifies boundaries
/// ("nodesets" and "sidesets") of the domain.  The Poisson.xml
/// example names the exterior faces of the cube with IDs.
///
/// The Poisson2D.xml example does not include a (parallel)
/// "decomposition strategy ... end" statement.  The default strategy
/// is "bisection," which attempts to assign an equal number of cells
/// to each parallel process.
void
makeMatrixAndRightHandSide2D (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<vector_type>& B,
                            Teuchos::RCP<vector_type>& X_exact,
                            Teuchos::RCP<vector_type>& X,
                            Teuchos::RCP<multivector_type> & coords,
                            Teuchos::RCP<vector_type>& node_sigma,
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Node>& node,
                            const std::string& meshInput,
                            Teuchos::ParameterList & inputList,
                            Teuchos::ParameterList & problemStatistics,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose = false,
                            const bool debug = false);

//! Just like above, but with multivector_type output arguments.
void
makeMatrixAndRightHandSide2D (Teuchos::RCP<sparse_matrix_type>& A,
                            Teuchos::RCP<multivector_type>& B,
                            Teuchos::RCP<multivector_type>& X_exact,
                            Teuchos::RCP<multivector_type>& X,
                            Teuchos::RCP<multivector_type> & coords,
                            Teuchos::RCP<vector_type>& node_sigma,
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Node>& node,
                            const std::string& meshInput,
                            Teuchos::ParameterList & inputList,
                            Teuchos::ParameterList & problemStatistics,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose = false,
                            const bool debug = false);

} // namespace TpetraIntrepidPoissonExample2D
} // namespace TrilinosCouplings

#endif // __TrilinosCouplings_TpetraIntrepidPoissonExample2D_hpp
