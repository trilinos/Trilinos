//@HEADER
/*!
 * \file Epetra_Multi_CrsMatrix.h
 *
 * \class Epetra_Multi_CrsMatrix
 *
 * \brief Class that encapsulates multiple CRS matrices and supports
 * MatrixMatrix Multiplies (as a derived class of Epetra_Operator_With_MatMat).
 *
 * \date Last update to Doxygen: 25-Jan-07
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef ML_EPETRA_MULTI_CRS_MATRIX_H
#define ML_EPETRA_MULTI_CRS_MATRIX_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
#if defined(HAVE_ML_EPETRA)
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "ml_Preconditioner.h"
#include "Epetra_Operator_With_MatMat.h"
#include "ml_include.h"

namespace ML_Epetra{

/*! The Epetra_MultiCrs_Matrix is an Epetra_Operator_With_MatMat that
  encapsulates and operator equal to the product of a number of
  Epetra_CrsMatrices.  However, this operator is never explicitly formed.
  Instead, the matrices are stored individually and the matvec/matmat operations
  are applied one at a time.

  This class is useful in applications where matvecs/matmats with an operator
  representing the product of matrices is used a number of times, but not enough
  to warrant forming an entire matrix.  Maintaining sparse representations of
  complicated operators is another potential use for this class.
*/

class Epetra_Multi_CrsMatrix: public Epetra_Operator_With_MatMat{
public:
  //! @name Constructor
  //@{
  //! Constructor - CRS matrices are applied from right to left.
  /*! WARNING: These are shallow links, so be sure these matrices do not get
    deallocated before you're done using the Epetra_Multi_CrsMatrix
  */
  Epetra_Multi_CrsMatrix(int NumMatrices,Epetra_CrsMatrix ** CrsMatrices);
  //@}

  //! @name Destructor
  //@{
  //! Destructor
  virtual ~Epetra_Multi_CrsMatrix();
  //@}

  //! @name Attribute set methods
  //@{

  //! Sets use transpose (not implemented).
  virtual int SetUseTranspose(bool /* useTranspose */ ){return(-1);}


  //@}

  //! @name Mathematical functions
  //@{

  //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*!
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a Epetra_Operator inverse applied to an
  //Epetra_MultiVector X in Y. NOT IMPLEMEMENTED!
  virtual int ApplyInverse(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const {return -1;}


  //! Computes C= <me> * A
  virtual int MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, Epetra_CrsMatrix **C) const;

  //! Computes C= <me> * A
  virtual int MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, ML_Comm *comm, ML_Operator **C) const;


  //@}

  //! @name Attribute access functions
  //@{

  //! Returns the infinity norm (not implemented).
  virtual double NormInf() const {return(0.0);};

  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const {return(false);};

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const{return(false);};

  //! Prints label associated to this object.
  virtual const char* Label() const{return(Label_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm& Comm() const{return(*Comm_);};

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map& OperatorDomainMap() const {return(*DomainMap_);};

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map& OperatorRangeMap() const {return(*RangeMap_);};
  //@}

private:
  //! Private Data
  //@{

  //! Number of Crs matrices to store
  int NumMatrices_;

  //! Pointers to the CRS matrices.  They are applied from right to left
  Epetra_CrsMatrix ** CrsMatrices_;

  //! Label for this object
  char* Label_;

  //! Domain Map
  const Epetra_Map* DomainMap_;
  //! Range Map
  const Epetra_Map* RangeMap_;
  //! Epetra communicator object
  const Epetra_Comm* Comm_;


  //@}

};//end Epetra_Multi_CrsMatrix

}/*end namespace*/

#endif
#endif /*ML_EPETRA_MULTI_CRS_MATRIX_H*/
