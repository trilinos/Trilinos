// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_FINITEDIFFERENCECOLORING_H
#define NOX_EPETRA_FINITEDIFFERENCECOLORING_H

#include "Epetra_RowMatrix.h"            // base class
#include "NOX_Epetra_FiniteDifference.H"    // base class
#include "NOX_Common.H"
#include "Teuchos_RCP.hpp"
#include <vector>
#include <list>
#include <map>

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_MapColoring;
class Epetra_Import;
class Epetra_Vector;
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
class Epetra_IntVector;
#else
class Epetra_LongLongVector;
#endif
class Epetra_CrsGraph;
class Epetra_CrsMatrix;

namespace NOX {

namespace Epetra {

/*! \brief Concrete implementation for creating an Epetra_RowMatrix Jacobian via finite differencing of the residual using coloring.

The Jacobian entries are calculated via 1st or 2nd order finite
differencing.  This requires \f$ N  + 1 \f$ or \f$ 2N + 1 \f$ calls to
computeF(), respectively, where \f$ N \f$ is the number of colors.

  \f[ J_{ij} = \frac{\partial F_i}{\partial x_j} = \frac{F_i(x+\delta\mathbf{e}_j) - F_i(x)}{\delta}  \f]

where \f$ J\f$ is the Jacobian, \f$ F\f$ is the function evaluation, \f$ x\f$
is the solution vector, and \f$\delta\f$ is a small perturbation to the
\f$ x_j\f$ entry.

Instead of perturbing each \f$ N_{dof} \f$ problem degrees of freedom
sequentially and then evaluating all \f$ N_{dof} \f$ functions for each
perturbation, coloring allows several degrees of freedom (all belonging
to the same color) to be perturbed at the same time.  This reduces the
total number of function evaluations needed to compute \f$\mathbf{J}\f$
from \f$ N_{dof}^2 \f$ as is required using FiniteDifference to \f$
N\cdot N_{dof} \f$, often representing substantial computational savings.

Coloring is based on a user-supplied color map generated using an
appropriate algorithm, eg greedy-algorithm - Y. Saad, "Iterative
Methods for Sparse %Linear Systems, 2<sup>nd</sup> ed.," chp. 3,
SIAM, 2003..  Use can be made of the coloring algorithm provided by
the EpetraExt package in Trilinos.  The 1Dfem_nonlinearColoring and
Brusselator example problems located in the nox/epetra-examples
subdirectory demonstrate use of the EpetraExt package, and the
1Dfem_nonlinearColoring directory also contains a stand-alone coloring
algorithm very similar to that in EpetraExt.

The perturbation, \f$ \delta \f$, is calculated using the following equation:

\f[ \delta = \alpha * | x_j | + \beta \f]

where \f$ \alpha \f$ is a scalar value (defaults to 1.0e-4) and \f$
\beta \f$ is another scalar (defaults to 1.0e-6).

Since both %FiniteDifferenceColoring and FiniteDifference inherit from
the Epetra_RowMatrix class, they can be used as preconditioning matrices
for AztecOO preconditioners.

As for FiniteDifference, 1st order accurate Forward and Backward differences
as well as 2nd order accurate Centered difference can be specified using
setDifferenceMethod with the appropriate enumerated type passed as
the argument.

<B>Using %FiniteDifferenceColoring in Parallel</B>

Two ways of using this class in a distributed parallel environment are
currently supported.  From an application standpoint, the two approaches
differ only in the status of the solution iterate used in the residual
fill.  If an object of this class is contructed with %parallelColoring = true
the solution iterate will be passe back in a non-ghosted form.  On the
contrary, setting this parameter to false in the constructor will
cause the solution iterate to be in a ghosted form when calling back for
a residual fill.  When using the second approach, the user should be
aware that the perturbed vector used to compute residuals has already
been scattered to a form consistent with the column space of the
Epetra_CrsGraph.  In practice, this means that the perturbed vector used
by computeF() has already been scattered to a ghosted or overlapped state.
The application should then not perform this step but rather simply use
the vector provided with the possible exception of requiring a local
index reordering to bring the column-space based vector in sync with a
potentially different ghosted index ordering.  See the %Brusselator and
%1Dfem_nonlinearColoring example problems for details.

<B>Special Case for Approximate Jacobian Construction</B>

Provision is made for a simplified and cheaper use of coloring that
currently provides only for the diagonal of the Jacobian to be computed.
This is based on using a first-neighbors coloring of the original
Jacobian graph using the Epetra_Ext MapColoring class with the distance1
argument set to true.  This same argument should also be set to true in
the constructor to this class.  The result will be a diagonal Jacobian
filled in a much more efficient manner.

*/

class FiniteDifferenceColoring : public FiniteDifference {

 protected:

  enum ColoringType { NOX_SERIAL, NOX_PARALLEL };

 public:

  //! Constructor with output control
  FiniteDifferenceColoring(
      Teuchos::ParameterList& printingParams,
      const Teuchos::RCP<Interface::Required>& i,
      const NOX::Epetra::Vector& initialGuess,
      const Teuchos::RCP<Epetra_MapColoring>& colorMap,
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
      const Teuchos::RCP<std::vector<Epetra_IntVector> >& columns,
#else
      const Teuchos::RCP<std::vector<Epetra_LongLongVector> >& columns,
#endif
      bool parallelColoring = false,
      bool distance1 = false,
      double beta = 1.0e-6, double alpha = 1.0e-4);

  //! Constructor with output control
  FiniteDifferenceColoring(
      Teuchos::ParameterList& printingParams,
      const Teuchos::RCP<Interface::Required>& i,
      const NOX::Epetra::Vector& initialGuess,
      const Teuchos::RCP<Epetra_CrsGraph>& rawGraph,
      const Teuchos::RCP<Epetra_MapColoring>& colorMap,
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
     const Teuchos::RCP<std::vector<Epetra_IntVector> >& columns,
#else
     const Teuchos::RCP<std::vector<Epetra_LongLongVector> >& columns,
#endif
      bool parallelColoring = false,
      bool distance1 = false,
      double beta = 1.0e-6, double alpha = 1.0e-4);

  //! Pure virtual destructor
  virtual ~FiniteDifferenceColoring();

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x);

  //! Output the coloring map, index map and underlying matrix
  //virtual void Print(ostream&);

  //! Create containers for using color and index maps in parallel coloring
  virtual void createColorContainers();

protected:

  //! Enum flag for type of coloring being used
  ColoringType coloringType;

  //! bool flag for specifying special case of distance-1 coloring
  bool distance1;

  //! Color Map created by external algorithm
  Teuchos::RCP<const Epetra_MapColoring> colorMap;

  //! vector of Epetra_IntVectors containing columns corresponding to a given row and color
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  Teuchos::RCP<std::vector<Epetra_IntVector> > columns;
#else
  Teuchos::RCP<std::vector<Epetra_LongLongVector> > columns;
#endif

  //! Number of colors in Color Map
  int numColors;

  //! Max Number of colors on all procs in Color Map
  int maxNumColors;

  //! List of colors in Color Map
  int* colorList;

  //! List of colors in Color Map
  std::list<int> listOfAllColors;

  //! Inverse mapping from color id to its position in colorList
  std::map<int, int> colorToNumMap;

  //! Coloring Map created by external algorithm
  Epetra_Map* cMap;

  //! Importer needed for mapping Color Map to unColored Map
  Epetra_Import* Importer;

  //! Color vector based on Color Map containing perturbations
  Epetra_Vector* colorVect;

  //! Color vector based on Color Map containing beta value(s)
  Epetra_Vector* betaColorVect;

  //! Color vector based on unColorred Map containing perturbations
  Epetra_Vector* mappedColorVect;

  //! Perturbed solution vector based on column map
  Epetra_Vector* xCol_perturb;

  //! Overlap Map (Column Map of Matrix Graph) needed for parallel
  const Epetra_BlockMap* columnMap;

  //! An Import object needed in parallel to map from row-space to column-space
  Epetra_Import* rowColImporter;

};
}  // namespace Epetra
}  // namespace NOX

#endif /* _EPETRA_FINITEDIFFERENCECOLORING_H_ */
