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

#include "Tpetra_Operator.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_ParameterList.hpp"

#include "RFGen_API_KLSolver.h"
#include "RFGen_CovarianceFunction.h"

namespace RFGen
{

class KLSolver 
  : public Tpetra::Operator<double>
{
private:
  using Scalar = double;
  using Operator = Tpetra::Operator<Scalar>;
  using LocalOrdinal = Operator::local_ordinal_type;
  using GlobalOrdinal = Operator::global_ordinal_type;
  using Node = Operator::node_type;
  using CrsMatrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using MultiVector = Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using Vector = Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

public:
  explicit KLSolver(const unsigned spatialDim);

  ~KLSolver() {}

  void setAPI(const Teuchos::RCP<API_KLSolver> &appInterface)
  {
    appInterface_ = appInterface;
    mpiComm_      = appInterface_->getParallelComm();
  }

  void setCovarianceFunction(const Teuchos::RCP<const CovarianceFunction> &covarFunc)
    {covarFunc_ = covarFunc;}

  void setPL(const Teuchos::ParameterList &pl) 
    {pl_=pl;}

  void setMatrixFree(const bool useMatrixFree) 
    {useMatrixFree_=useMatrixFree;}

  void solve(const int maxNev);

  // virtuals from Epetra_Operator
  bool hasTransposeApply() const override { return false; }

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const override { return map_; }
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const override { return map_; }

  void apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
          Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
          Teuchos::ETransp mode,
          Scalar alpha,
          Scalar beta) const override;

 private:
  //void assembleMatrixFreeCovariance();
  void assembleFullCovariance();

  void assembleMassMat();

  void partialAssemble(
    const shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> &nonlocalIntgPtCoords,	       
    const shards::Array<double,shards::NaturalOrder,Cell,Point> nonlocalVolumeWeights);

  void partialMatVec(
    const double alpha,
    const shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> &nonlocalIntgPtCoords,
    const shards::Array<double,shards::NaturalOrder,Cell,Point> &nonlocalVolumeWeights,
    const shards::Array<double,shards::NaturalOrder,Eigen,Cell> &X_vec,
    shards::Array<double,shards::NaturalOrder,Eigen,Cell> &Y_vec) const;

  // covariance operator
  Teuchos::RCP<const Operator> covarianceOp_;

  // vector space
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> map_;

  // application class to: 
  //   1) read in intg pt coords and volume weights
  //   2) provide parallel communicator
  //   3) receive KL solution
  Teuchos::RCP<API_KLSolver> appInterface_;

  MPI_Comm mpiComm_;

  Teuchos::ParameterList pl_;

  int localNumElem_, globalNumElem_;
  int localMaxIntgPts_, maxIntgPts_;

  const unsigned spatialDim_;

  bool useMatrixFree_;

  std::vector<GlobalOrdinal> localElemIDs_, nonlocalElemIDs_;

  Teuchos::RCP<const CovarianceFunction> covarFunc_;

  Teuchos::RCP<CrsMatrix> tpetra_covarianceMat_, tpetra_massMat_;

  std::vector<double> lambda_mem_;

  shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> localIntgPtCoords_;
  shards::Array<double,shards::NaturalOrder,Cell,Point> localVolumeWeights_;
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
