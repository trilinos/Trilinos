/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * October 20, 2002, the United States Government is granted for itself and
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
