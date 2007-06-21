//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#ifndef NOX_EPETRA_SCALING_H
#define NOX_EPETRA_SCALING_H

#include "NOX_Common.H"           // for <string>
#include "NOX_Abstract_Vector.H"  // 
#include "Teuchos_RCP.hpp"

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
  virtual void print(ostream& os);

private:
  
  //! Scaling type
  vector<ScaleType> scaleType;

  //! Source type
  vector<SourceType> sourceType;

  //! Scaling vector pointer
  vector< Teuchos::RCP<Epetra_Vector> > scaleVector;

  //! Temporary vector
  Teuchos::RCP<Epetra_Vector> tmpVectorPtr;

};
}  // namespace Epetra
}  // namespace NOX

ostream& operator<<(ostream& os, NOX::Epetra::Scaling& scalingObject);

#endif /* NOX_EPETRA_SCALING_H */
