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

#ifndef IFPACK_PERTURBEDMATRIX_H
#define IFPACK_PERTURBEDMATRIX_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Object.h"
class Epetra_Comm;
class Epetra_Map;
//! Ifpack_PerturbedMatrix: Supports the use of diagonal perturbations in Ifpack.



class Ifpack_PerturbedMatrix: public Epetra_Object {

 public:
  //@{ \name Constructors/Destructor

  //! Constructor using Epetra_RowMatrix.
  /*! Creates an Ifpack_PerturbedMatrix object from the
	 Epetra_RowMatrix interface. 
    \param In
            RowMatrix - An object that has implemented the Epetra_RowMatrix interface.
  */
  Ifpack_PerturbedMatrix(const Epetra_RowMatrix & UserMatrix);
  
  //! Copy constructor.
  Ifpack_PerturbedMatrix(const Ifpack_PerturbedMatrix & Source);

  //! Ifpack_PerturbedMatrix Destructor
  virtual ~Ifpack_PerturbedMatrix();
  //@}

  //@{ \name Initialization methods.

  //! Set absolute threshold value
  void SetAbsoluteThreshold( double Athresh) {Athresh_ = Athresh; return;}

  //! Set relative threshold value
  void SetRelativeThreshold( double Rthresh) {Rthresh_ = Rthresh; return;}


  //@}
#endif // IFPACK_PERTURBEDMATRIX_H
