#ifndef ML_WORKSPACE_H
#define ML_WORKSPACE_H

/*!
\file ml_Workspace.h
*/

/*!
 * \namespace ML_FiniteElements
 *
 * \brief Default namespace for ML finite element example.
 *
 * The ML_FiniteElements namespace contains all the classes and the
 * examples that show a possible structure for a finite element code
 * using Epetra and ML.
 *
 * Scalar second-order, symmetric and non-symmetric PDEs of type
 *
 *  \f[ - \eps \Delta u + b * \grad u + \sigma u = f  on \Omega \f]
 *  \f[                                        u = g  on \partial \Omega \f]
 *
 * can be discretized using this example. Neumann boundary conditions
 * require minimal changes to the code. \f$\Omega\f$ must be a 2D or a 3D 
 * domain, discretized using triangles, quadrilaterals, tetrahedra or
 * hexahedra. The code can be quite easily extended to tackle vector
 * problems, using the same finite element space for all unknowns.
 *
 * Two discretizations are available:
 * - GalerkinVariational is a pure Galerkin approach. Since this approach is
 *   unstable for large Peclet numbers, the class defines b = 0.
 * - SUPGVariational is a SUPG discretization, only for triangles and
 *   tetrahedra, and coth formula for tau. It should be used for
 *   advective problems.
 *
 * The code is based on the following pure virtual classes:
 * - AbstractGrid defines the query methods that each grid class must
 *   implement. AbstractGrid is based on the getrow() concept: instead
 *   of prescribing a given format, the user must define a set of methods,
 *   that will be used by the code to query for the vertices in a local
 *   elements, the coordinates of each vertex, the boundary faces, and so on.
 * - AbstractQuadrature defines the basis and test functions and
 *   the quadrature strategy.
 * - AbstractVariational defines a set of abstract interfaces for the
 *   variational form. An implementation of this class constructs the
 *   elemental finite element matrix and right-hand side.
 * - AbstractProblem takes an AbstractGrid, an AbstractVariational,
 *   and an AbstractQuadrature in input,
 *   then builds the global finite element matrix, and imposes the
 *   boundary conditions as specified by the user. It can also compute
 *   the norm of the numerical solution, exact solution, and error, using
 *   H1, semi-H1 or L2 norms.
 * 
 * The solution can be visualized using MEDIT (see web page
 * http://www.ann.jussieu.fr/~frey/logiciels/medit.html for details and
 * download). 
 * 
 * A nice feature of this small example is that some grids can be created
 * on-the-fly, with no need of input files. The following domains can be
 * triangulated:
 * - In 2D, a rectangle, using classes TriangleRectangleGrid and
 *   QuadRectangleGrid.
 * - In 3D, a cube, using classes TetCubeGrid and HexCubeGrid.
 * 
 * These classes are valid for both serial and parallel runs.
 *
 * An interface to triangle (see http://www-2.cs.cmu.edu/~quake/triangle.html
 * for more details and download) is under development; the interested user
 * can modify the file ml_TRIANGLEGrid.h. This interface is currently working
 * for serial runs only.
 *
 * The interested user should look to the following examples, all located
 * in ml/examples/FiniteElements:
 * - \ref Laplacian3D describes how to solve the Laplace equation in
 *   a unitary cube, discretized using hexahedra. After solution, the norm
 *   of the numerical solution, exact solution and error is computed, and the
 *   numerical solution is visualized using MEDIT.
 * - \ref AdvDiff2D describes how to solve an advection-diffusion problem
 *   in a 2D domain, using SUPG.
 * - \ref MLAPI shows how to interface with the MLAPI classes.
 * - TestAll is for testing purposed only.
 *
 * \author Marzio Sala, SNL 9214.
 *
 * \date Last updated on Apr-05.
 */

namespace ML_FiniteElements
{
  const int ML_INTERNAL = 0;
  const int ML_BOTTOM   = -1;
  const int ML_RIGHT    = -2;
  const int ML_TOP      = -3;
  const int ML_LEFT     = -4;
  const int ML_FRONT    = -5;
  const int ML_REAR     = -6;

  const int ML_DIRICHLET = 1000;
  const int ML_NEUMANN   = 1001;
}

/*!
 * \page Laplacian3D Laplacian3D.cpp
 * \include Laplacian3D.cpp
 */

/*!
 * \page AdvDiff2D AdvDiff2D.cpp
 * \include AdvDiff2D.cpp
 */

/*!
 * \page ml_MLAPI MLAPI.cpp
 * \include MLAPI.cpp
 */

#endif
