/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_CRSILUT_H
#define IFPACK_CRSILUT_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_ScalingType.h"
#include "Ifpack_OverlapGraph.h"
#include "Ifpack_OverlapFactorObject.h"
#include "Ifpack_OverlapSolveObject.h"
#include "Epetra_CompObject.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Object.h"
class Epetra_Comm;
class Epetra_Map;
class Epetra_RowMatrix;
class Epetra_Vector;
class Epetra_MultiVector;

//! Ifpack_CrsIlut: ILUT preconditioner of a given Epetra_RowMatrix.



class Ifpack_CrsIlut: public Epetra_Object, public Epetra_CompObject, public Ifpack_OverlapFactorObject, public Ifpack_OverlapSolveObject {
  
 public:
  //@{ \name Constructors/Destructor

  //! Constructor using Ifpack_OverlapGraph.
  /*! Creates an object from the overlap graph. 
    \param OverlapGraph (In) - Graph describing the graph that should be used for the factors.
    \param DropTol (In/Default) - Drop tolerance used by ILUT algorithm.
    \param FillTol (In/Default) - Fill tolerance used by ILUT algorithm.

  */
  Ifpack_CrsIlut(const Ifpack_OverlapGraph * OverlapGraph, double DropTol = 1.0E-4, 
		 double FillTol = 1.0);

  //! Constructor using Epetra_RowMatrix.
  /*! Creates an Ifpack_Graph object from the user graph implicitly defined by the
	 Epetra_RowMatrix interface. 
    \param RowMatrix (In) - An object that has implemented the Epetra_RowMatrix interface.
    \param DropTol (In/Default) - Drop tolerance used by ILUT algorithm.
    \param FillTol (In/Default) - Fill tolerance used by ILUT algorithm.

  */
  Ifpack_CrsIlut(const Epetra_RowMatrix * UserMatrix, double DropTol = 1.0E-4, 
		 double FillTol = 1.0);
  
  //! Copy constructor.
  Ifpack_CrsIlut(const Ifpack_CrsIlut & Source);
  
  //! Ifpack_CrsIlut Destructor
  virtual ~Ifpack_CrsIlut();
  //@}

  //@{ \name Initialization methods.

  //! Set Drop tolerance value as defined by the ILUT algorithm.
  int SetDropTol(double DropTol) {DropTol_ = DropTol; return(0);};

  //! Set fill tolerance value as defined by the ILUT algorithm.
  int SetFillTol(double FillTol) {FillTol_ = FillTol; return(0);};

#ifdef HAVE_IFPACK_TEUCHOS
  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes two parameter names: drop_tolerance and
     fill_tolerance. These names are case insensitive. For both, the
     ParameterEntry must have type double.
  */
  int SetParameters(const Teuchos::ParameterList& parameterlist,
                    bool cerr_warning_if_unused=false);
#endif
  //@}
  
  //@{ \name Attribute access methods.

  //! Set Drop tolerance value as defined by the ILUT algorithm.
  double DropTol() const {return(DropTol_);};

  //! Set fill tolerance value as defined by the ILUT algorithm.
  double FillTol() const {return(FillTol_);};
  //@}
  
 protected:
  //@{ \name Methods needed to implement OverlapFactorObject.

  //! Processes the overlapped user matrix for computing the ILUT preconditioner.
  int ProcessOverlapMatrix(const Epetra_RowMatrix &A);
  //! Compute ILUT factors L and U: WARNING: THIS ROUTINE IS NOT USER CALLABLE, CALL Factor().
  int DerivedFactor();
  //@}

 private:

 double DropTol_;
 double FillTol_;

};

#endif /* IFPACK_CRSILUT_H */
