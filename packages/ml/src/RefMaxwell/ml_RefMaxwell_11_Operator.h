/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
//@HEADER
/*!
 * \file ml_RefMaxwell_11_Operator.h
 *
 * \class ML_RefMaxwell_11_Operator
 *
 * \brief Epetra_Operator that encapsulates the (1,1) block of the reformulated
 * operator for Maxwell.
 *
 *
 * \date Last update to Doxygen: 25-Jan-07
 *
 */

#ifndef ML_REFMAXWELL_11_OPERATOR_H
#define ML_REFMAXWELL_11_OPERATOR_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
#if defined(HAVE_ML_EPETRA) && defined (HAVE_ML_EPETRAEXT)
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "ml_Preconditioner.h"
#include "Epetra_Operator_With_MatMat.h"
#include "Epetra_Multi_CrsMatrix.h"
#include "EpetraExt_Reindex_CrsMatrix.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "EpetraExt_SolverMap_CrsMatrix.h"
namespace ML_Epetra{

/*! ML_RefMaxwell_11_Operator encapsulates the reformulated (1,1) block
  operator of the system described in Bochev, Hu, Siefert and Tuminaro, 2007.
  It inherits from Epetra_Operator_With_MatMat, and provides encapsulation for
  the operator:
     S + M + M1 D0 M0^{-1} DO^T M1
*/

class ML_RefMaxwell_11_Operator: public Epetra_Operator_With_MatMat{
public:
  //! @name Constructor
  //@{
  //! Constructor - All the matrices needed for Maxwell.  OptimizeStorage *must*
  // be called for all of these matrices before you try to use the matmat functions.
  // WARNING:  All these matrices will be shallow pointed to.  Please be sure
  // they stick around until after ML_RefMaxwell_11_Operator is done.
  ML_RefMaxwell_11_Operator(const Epetra_CrsMatrix& SM_Matrix,    //S+M
                            const Epetra_CrsMatrix& D0_Matrix,    //T or D0
                            const Epetra_CrsMatrix& M0inv_Matrix, //M0^{-1}
                            const Epetra_CrsMatrix& M1_Matrix);   //M1(1)
  //@}

  //! @name Destructor
  //@{
  //! Destructor
  virtual ~ML_RefMaxwell_11_Operator();
  //@}

  //! @name Attribute set methods
  //@{

  //! Sets use transpose (not implemented).
  virtual int SetUseTranspose(bool /* UseTranspose */){return(-1);}


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


  //! Computes C= <me> * A.  OptimizeStorage *must* be called for both A and the
  //matrices in *this, before this routine can work.
  virtual int MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, Epetra_CrsMatrix **C) const;

  //! Computes C= <me> * A.  OptimizeStorage *must* be called for both A and the
  //matrices in *this, before this routine can work.
  virtual int MatrixMatrix_Multiply(const Epetra_CrsMatrix & A, ML_Comm *comm, ML_Operator **C) const;

  //! Computes C= A^T * <me> * A.  OptimizeStorage *must* be called for both A and the
  //matrices in *this, before this routine can work.
  virtual int PtAP(const Epetra_CrsMatrix & A, ML_Comm *comm, ML_Operator **C) const;


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

  //! EXPERIMENTAL: Return SM Matrix
  virtual const Epetra_CrsMatrix & SM_Matrix(){return *SM_Matrix_;}
  //@}

private:
  //! Private Data
  //@{
  //! Matrix: S+M1(sigma)
  const Epetra_CrsMatrix * SM_Matrix_;
  //! Matrix: M1 D0^T M0inv D0 M1
  Epetra_CrsMatrix ** Addon_Matrix_;

  //! Matrix: D0^T
  Epetra_CrsMatrix * D0T_Matrix_;
  EpetraExt::RowMatrix_Transpose * D0_Matrix_Transposer_;
  EpetraExt::CrsMatrix_SolverMap D0T_Matrix_Trans_;

  //! Multi_Crs_Matrix
  Epetra_Multi_CrsMatrix *Addon_;

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

}//end namespace
#endif
#endif /*ML_REFMAXWELL_11_OPERATOR_H*/
