
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_FINITEDIFFERENCECOLORING_H_
#define _EPETRA_FINITEDIFFERENCECOLORING_H_

#include "Epetra_RowMatrix.h"    	    // base class
#include "NOX_Epetra_FiniteDifference.H"    // base class

#include "NOX_Common.H"          // for <string>

#include "vector"                // for column indices from coloring

// The following using statement facilitates proper wrapping by SWIG,
// the simple wrapper interface generator, when creating a python
// interface to Trilinos/NOX.  The statement can be wrapped with
// #ifdefs if it creates a porting problem, but must be in this file.
using namespace std;

// Forward Declarations
class Epetra_Comm;
class Epetra_Map;
class Epetra_MapColoring;
class Epetra_Import;
class Epetra_Vector;
class Epetra_IntVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;

namespace NOX {
  class Utils;
  namespace Epetra {
    class Interface;
  }
}

namespace NOX {

namespace Epetra {

/*! \brief Concrete implementation for creating an Epetra_RowMatrix Jacobian via finite differencing of the residual using coloring.

The Jacobian entries are calculated via 1st or 2nd order finite
differencing.  This requires \f$ N  + 1 \f$ or \f$ 2N + 1 \f$ calls to
computeF(), respectively, where \f$ N \f$ is the number of colors.

  \f[ J_{ij} = \frac{\partial F_i}{\partial x_j} = \frac{F_i(x+\delta\mathbf{e}_j) - F_i(x)}{\delta}  \f]

where \f$J\f$ is the Jacobian, \f$F\f$ is the function evaluation, \f$x\f$
is the solution vector, and \f$\delta\f$ is a small perturbation to the
\f$x_j\f$ entry.

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

This class can be used to numerically approximate \f$ J \f$ for
distributed parallel problems provided a very minimal of conditions are
met.  Namely, an Epetra_CrsGraph reflecting a row-based distribution must
be used in the constructor.  Next, the application must be aware that the
perturbed vector used to compute residuals has already been scattered
to a form consistent with the column space of the Epetra_CrsGraph.
In practice, this means that the perturbed vector used by computeF()
has already been scattered to a ghosted or overlapped state.  The
application should then not perform this step but rather simply use
the vector provided with the possible exception of requiring a local
index reordering to bring the column-space based vector in sync with a
potentially different ghosted index ordering.  See the %Brusselator and
%1Dfem_nonlinearColoring example problems for details.


*/

class FiniteDifferenceColoring : public FiniteDifference {
      
 public:

  //! Constructor

  FiniteDifferenceColoring(Interface& i, const Epetra_Vector& initialGuess, 
                           Epetra_MapColoring& colorMap,
                           vector<Epetra_IntVector>& columns,
                           double beta = 1.0e-6, double alpha = 1.0e-4);

  //! Constructor

  FiniteDifferenceColoring(Interface& i, const Epetra_Vector& initialGuess, 
                           Epetra_CrsGraph& rawGraph,
                           Epetra_MapColoring& colorMap,
                           vector<Epetra_IntVector>& columns,
                           double beta = 1.0e-6, double alpha = 1.0e-4);

  //! Pure virtual destructor
  virtual ~FiniteDifferenceColoring();

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  //! Compute Jacobian given the specified input vector, x. Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x);

protected:

  //! Color Map created by external algorithm 
  const Epetra_MapColoring* colorMap;

  //! vector of Epetra_IntVectors containing columns corresponding to a given row and color
  vector<Epetra_IntVector>* columns;

  //! Number of colors in Color Map 
  int numColors;

  //! List of colors in Color Map 
  int* colorList;

  //! Coloring Map created by external algorithm
  Epetra_Map* cMap;

  //! Importer needed for mapping Color Map to unColored Map
  Epetra_Import* Importer;

  //! Color vector based on Color Map containing perturbations
  Epetra_Vector* colorVect;

  //! Color vector based on Color Map containing beta value(s)
  Epetra_Vector* betaColorVect;

  //! Color vector based on unColorred Map containing perturbations
  Epetra_Vector mappedColorVect;

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
