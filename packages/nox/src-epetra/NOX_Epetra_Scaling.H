//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
