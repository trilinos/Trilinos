#ifndef __TrilinosCouplings_TpetraIntrepidPoissonExample_hpp
#define __TrilinosCouplings_TpetraIntrepidPoissonExample_hpp

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ScalarTraits.hpp"

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
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  Node;

//
// mfh 19 Apr 2012: If you want to change the template parameters of
// these typedefs, modify the typedefs (ST, LO, GO, Node) above.
//
typedef Tpetra::CrsMatrix<ST, LO, GO, Node>    sparse_matrix_type;
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
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Node>& node,
                            const std::string& meshInput,
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
                            const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                            const Teuchos::RCP<Node>& node,
                            const std::string& meshInput,
                            const Teuchos::RCP<Teuchos::FancyOStream>& out,
                            const Teuchos::RCP<Teuchos::FancyOStream>& err,
                            const bool verbose = false,
                            const bool debug = false);

//! Return ||B - A*X_exact||_2, ||B||.
std::pair<Teuchos::ScalarTraits<ST>::magnitudeType,
          Teuchos::ScalarTraits<ST>::magnitudeType>
exactResidualNorm (const Teuchos::RCP<const sparse_matrix_type>& A,
                   const Teuchos::RCP<const vector_type>& B,
                   const Teuchos::RCP<const vector_type>& X_exact);

} // namespace TpetraIntrepidPoissonExample

#endif // __TrilinosCouplings_TpetraIntrepidPoissonExample_hpp
