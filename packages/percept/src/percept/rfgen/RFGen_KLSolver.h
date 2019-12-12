/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef RFGen_KLSolver_h
#define RFGen_KLSolver_h

#include "Epetra_Operator.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"

#include "RFGen_API_KLSolver.h"
#include "RFGen_CovarianceFunction.h"

class Epetra_CrsMatrix;

namespace RFGen
{

class KLSolver 
  : public Epetra_Operator
{
 public:
  explicit KLSolver(const unsigned spatialDim);

  ~KLSolver() {}

  void setAPI(const Teuchos::RCP<API_KLSolver> &appInterface)
  {
    appInterface_ = appInterface;
    mpiComm_      = appInterface_->getParallelComm();
    epetraComm_   = Epetra_MpiComm(mpiComm_);
  }

  void setCovarianceFunction(const Teuchos::RCP<const CovarianceFunction> &covarFunc)
    {covarFunc_ = covarFunc;}

  void setPL(const Teuchos::ParameterList &pl) 
    {pl_=pl;}

  void setMatrixFree(const bool useMatrixFree) 
    {useMatrixFree_=useMatrixFree;}

  void solve(const int maxNev);

  // virtuals from Epetra_Operator
  int SetUseTranspose(bool UseTranspose);

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual double NormInf() const;
    //! Returns a character string describing the operator
  virtual const char * Label() const;
    //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const;

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const;
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const;
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const;
  
  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const;

  // END virtuals from Epetra_Operator

 private:
  //void assembleMatrixFreeCovariance();
  void assembleFullCovariance();

  void assembleMassMat();

  void partialAssemble(
    const shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> &nonlocalIntgPtCoords,	       
    const shards::Array<double,shards::NaturalOrder,Cell,Point> nonlocalVolumeWeights);

  void partialMatVec(
    const shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> &nonlocalIntgPtCoords,
    const shards::Array<double,shards::NaturalOrder,Cell,Point> &nonlocalVolumeWeights,
    const shards::Array<double,shards::NaturalOrder,Eigen,Cell> &X_vec,
    shards::Array<double,shards::NaturalOrder,Eigen,Cell> &Y_vec) const;

  // covariance operator
  Teuchos::RCP<const Epetra_Operator> covarianceOp_;

  // vector space
  Teuchos::RCP<const Epetra_Map> epetra_Map_;

  // application class to: 
  //   1) read in intg pt coords and volume weights
  //   2) provide parallel communicator
  //   3) receive KL solution
  Teuchos::RCP<API_KLSolver> appInterface_;

  MPI_Comm mpiComm_;
  Epetra_MpiComm epetraComm_;

  Teuchos::ParameterList pl_;

  int localNumElem_, globalNumElem_;
  int localMaxIntgPts_, maxIntgPts_;

  const unsigned spatialDim_;

  bool useMatrixFree_;

  std::vector<int> localElemIDs_, nonlocalElemIDs_;

  Teuchos::RCP<const CovarianceFunction> covarFunc_;

  Teuchos::RCP<Epetra_CrsMatrix> epetra_covarianceMat_, epetra_massMat_;

  std::vector<double> lambda_mem_;

  shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> localIntgPtCoords_;
  shards::Array<double,shards::NaturalOrder,Cell,Point> localVolumeWeights_;

  bool log_debug_;
};

Teuchos::RCP<KLSolver> 
  buildKLSolver(
    const Teuchos::RCP<API_KLSolver> &appInterface,
    const Teuchos::RCP<const CovarianceFunction> &covarFunc,
    const Teuchos::ParameterList &pl,
    const bool useMatrixFree = false
    );


}

#endif // RFGen_KLSolver_h
