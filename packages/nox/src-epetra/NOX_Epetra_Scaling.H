
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

#ifndef NOX_EPETRA_SCALING_H
#define NOX_EPETRA_SCALING_H

#include "NOX_Common.H"           // for <string>
#include "NOX_Abstract_Vector.H"  // 

// Forward Declarations
class Epetra_Vector;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

namespace NOX {
  class Utils;
  namespace Epetra {
    class Interface;
  }
}

namespace NOX {

namespace Epetra {

  /*! \brief Object to control scaling of a linear system.  

  Scaling of the linear system can affect the Nonlinear solves if using 
  an inexact Newton method.  We must account for computed norms.

  Currently this assumes a diagonal scaling only!  Once epetra can 
  do matrix-matrix multiplies we will expand this class.
  
  */
class Scaling {
      
public:
  
  //! Describes where the scaling vector comes from.
  enum SourceType {None, RowSum, UserDefined};
  
  //! Describes the type of scaling to apply.
  enum ScaleType {Left, Right};
  
  //! Constructor. 
  Scaling();
  
  //! Virtual destructor
  virtual ~Scaling();
  
  //! Add a user supplied diagonal scale vector to the scaling object.
  virtual void addUserScaling(ScaleType type, Epetra_Vector& D);

  //! Add "Row Sum" scaling to the scaling object.  The supplied vector is used to store the current row sum vector.
  virtual void addRowSumScaling(ScaleType type, Epetra_Vector& D);

  //! Computes Row Sum scaling diagonal vectors.  Only needs to be called if a row sum scaling has been rewuested.
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
  
  vector<ScaleType> scaleType;
  vector<SourceType> sourceType;
  vector<Epetra_Vector*> scaleVector;

  Epetra_Vector* tmpVectorPtr;

};
}  // namespace Epetra
}  // namespace NOX

ostream& operator<<(ostream& os, NOX::Epetra::Scaling& scalingObject);

#endif /* NOX_EPETRA_SCALING_H */
