// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_SCALING_H
#define NOX_EPETRA_SCALING_H

#include "NOX_Common.H"           // for <string>
#include "NOX_Abstract_Vector.H"  //
#include "Teuchos_RCP.hpp"
#include <vector>
#include <ostream>

// Forward Declarations
class Epetra_Vector;
class Epetra_RowMatrix;
class Epetra_LinearProblem;
namespace NOX {
  class Utils;
}

namespace NOX {

namespace Epetra {

  /*! \brief Object to control scaling of vectors and linear systems.

  Currently this assumes a diagonal scaling only!  Once epetra can
  do matrix-matrix multiplies we will expand this class.

  */
class Scaling {

public:

  //! Describes where the scaling vector comes from.
  enum SourceType {None, RowSum, ColSum, UserDefined};

  //! Describes the type of scaling to apply.
  enum ScaleType {Left, Right};

  //! Constructor.
  Scaling();

  //! Virtual destructor
  virtual ~Scaling();

  //! Add a user supplied diagonal scale vector to the scaling object.
  virtual void addUserScaling(ScaleType type, const Teuchos::RCP<Epetra_Vector>& D);

  //! Add "Row Sum" scaling to the scaling object.  The supplied vector is used to store the current row sum vector.
  virtual void addRowSumScaling(ScaleType type, const Teuchos::RCP<Epetra_Vector>& D);

  //! Add "Col Sum" scaling to the scaling object.  The supplied vector is used to store the current column sum vector.
  virtual void addColSumScaling(ScaleType type, const Teuchos::RCP<Epetra_Vector>& D);

  //! Computes Row Sum scaling diagonal vectors.  Only needs to be called if a row or column sum scaling has been requested.
  virtual void computeScaling(const Epetra_LinearProblem& problem);

  //! Scales the linear system.
  virtual void scaleLinearSystem(Epetra_LinearProblem& problem);

  //! Remove the scaling from the linear system.
  virtual void unscaleLinearSystem(Epetra_LinearProblem& problem);

  //! Applies any RIGHT scaling vectors to an input vector.
  virtual void applyRightScaling(const Epetra_Vector& input,
                 Epetra_Vector& result);

  //! Applies any LEFT scaling vectors to an input vector.
  virtual void applyLeftScaling(const Epetra_Vector& input,
                Epetra_Vector& result);

  //! Printing
  virtual void print(std::ostream& os);

private:

  //! Scaling type
  std::vector<ScaleType> scaleType;

  //! Source type
  std::vector<SourceType> sourceType;

  //! Scaling vector pointer
  std::vector< Teuchos::RCP<Epetra_Vector> > scaleVector;

  //! Temporary vector
  Teuchos::RCP<Epetra_Vector> tmpVectorPtr;

};

  std::ostream& operator<<(std::ostream& os, NOX::Epetra::Scaling& scalingObject);

}  // namespace Epetra
}  // namespace NOX

#endif /* NOX_EPETRA_SCALING_H */
