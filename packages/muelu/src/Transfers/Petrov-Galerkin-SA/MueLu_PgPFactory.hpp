/*
 * MueLu_PgPFactory.hpp
 *
 *  Created on: 23.09.2011
 *      Author: tobias
 */

#ifndef MUELU_PGPFACTORY_HPP_
#define MUELU_PGPFACTORY_HPP_

#include <Teuchos_Assert.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraVector.hpp>
#endif
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
    @class PgPFactory class.
    @brief Factory for building Petrov-Galerkin Smoothed Aggregation prolongators.
 */

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class PgPFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, PgPFactory<AA,BB,CC,DD, EE> &factory);

public:

  /* Options defining how to pick-up the next root node in the local aggregation procedure */
  enum MinimizationNorm {
    ANORM = 0, /* A norm   */
    L2NORM = 1, /* L2 norm */
    DINVANORM  = 2, /* Dinv A norm */
    ATDINVTPLUSDINVANORM = 3 /* most expensive variant */
  };

  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.
      User can supply a factory for generating the tentative prolongator.
   */
  PgPFactory(RCP<PFactory> InitialPFact = Teuchos::null, RCP<SingleLevelFactoryBase> AFact = Teuchos::null)
  : initialPFact_(InitialPFact), AFact_(AFact),
    diagonalView_("current") {
    if(initialPFact_ == Teuchos::null)
    {
      // use tentative P factory as default
      initialPFact_ = rcp(new TentativePFactory());
    }

    min_norm_ = DINVANORM;
  }

  //! Destructor.
  virtual ~PgPFactory() {}

  //@}

  //! @name Set methods.
  //@{

  //! Change view of diagonal.
  void SetDiagonalView(std::string const& diagView) {
    diagonalView_ = diagView;
  }

  //! Set minimization mode (L2NORM for cheapest, ANORM more expensive, DINVANORM = default, ATDINVTPLUSDINVANORM most expensive method)
  inline void SetMinimizationMode(MinimizationNorm minnorm) { min_norm_ = minnorm; }
  //@}

  //! @name Get methods.
  //@{

  //! Returns current view of diagonal.
  std::string GetDiagonalView() {
    return diagonalView_;
  }

  //@}

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("A",AFact_.get());
    coarseLevel.DeclareInput("P",initialPFact_.get());
  };
  //@}

  //! @name Build methods.
  //@{

  /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
   */
  void Build(Level& fineLevel, Level &coarseLevel) const {
    Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(std::cout));
    std::ostringstream buf; buf << coarseLevel.GetLevelID();
    RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("PgPFactory::BuildP_"+buf.str()));
    timer->start(true);

    // Level Get
    RCP<Operator> Ptent = coarseLevel.Get< RCP<Operator> >("P", initialPFact_.get());
    RCP<Operator> A     = fineLevel.  Get< RCP<Operator> >("A", AFact_.get());

    /////////////////// switch from A to A^T in restriction mode (necessary as long as implicit transpose not working for Epetra)
    if(restrictionMode_)
      A = Utils2<Scalar,LocalOrdinal,GlobalOrdinal>::Transpose(A,true); // build transpose of A explicitely

    Monitor m(*this, "prolongator smoothing (PG-AMG)");

    /////////////////// calculate D^{-1} A Ptent (needed for smoothing)
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Operator> DinvAP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

    doFillComplete=true;
    optimizeStorage=false;
    Teuchos::ArrayRCP<Scalar> diag = Utils::GetMatrixDiagonal(A);
    Utils::MyOldScaleMatrix(DinvAP0,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    /////////////////// calculate local damping factors omega

    // create array for row based omegas
    Teuchos::ArrayRCP<Scalar> RowBasedOmega( DinvAP0->getNodeNumRows(), Teuchos::ScalarTraits<Scalar>::zero() );

    // calculate row based omegas
    ComputeRowBasedOmegas(A,Ptent,DinvAP0,diag,RowBasedOmega);


    /////////////////// prolongator smoothing using local damping parameters omega
    RCP<Operator> P_smoothed = Teuchos::null;
    Utils::MyOldScaleMatrix(DinvAP0,RowBasedOmega,false,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    Utils::TwoMatrixAdd(Ptent, false, Teuchos::ScalarTraits<Scalar>::one(),
                                         DinvAP0, false, -Teuchos::ScalarTraits<Scalar>::one(),
                                         P_smoothed);
    P_smoothed->fillComplete(Ptent->getDomainMap(), Ptent->getRangeMap());

    //////////////////// store results in Level

    // Level Set
    if(!restrictionMode_)
    {
        // prolongation factory is in prolongation mode
        coarseLevel.Set("P", P_smoothed, this);
    }
    else
    {
        // prolongation factory is in restriction mode
        RCP<Operator> R = Utils2<SC,LO,GO>::Transpose(P_smoothed,true); // use Utils2 -> specialization for double
        coarseLevel.Set("R", R, this);
    }

    timer->stop();
    MemUtils::ReportTimeAndMemory(*timer, *(P_smoothed->getRowMap()->getComm()));
  }

  void BuildP(Level &fineLevel, Level &coarseLevel) const {
    std::cout << "TODO: remove me" << std::endl;
  }

  //@}

private:

  //! @brief Compute row based omegas
  //!
  //! This routine computes the row based omegas for the current processor
  //! The input parameters are the matrix A, its diagonal, the tentative prolongator
  //! and DinvAPtent. These are already available from the calling function, so we
  //! can use them without recalculating.
  //! First the column based omegas are calculated and then transformed to row-based omegas.
  //!
  //! @param[in] const RCP<Operator>& A: operator A
  //! @param[in] const RCP<Operator>& Ptent: tentative prolongation operator
  //! @param[in] const RCP<Operator>& DinvAPtent: scaled product of A and Ptent
  //! @param[in] const ArrayRCP<Scalar>& diagA: diagonal of matrix A
  //! @param[out] ArrayRCP<Scalar>& RowBasedOmegas: vector of row based omegas
  void ComputeRowBasedOmegas(const RCP<Operator>& A, const RCP<Operator>& Ptent, const RCP<Operator>& DinvAPtent,const Teuchos::ArrayRCP<Scalar>& diagA, Teuchos::ArrayRCP<Scalar>& RowBasedOmegas) const
  {

    std::map<GlobalOrdinal, GlobalOrdinal> GID2localgid;
    BuildLocalReplicatedColMap(DinvAPtent,GID2localgid);

    RCP<Teuchos::Array<Scalar> > Numerator   = Teuchos::null;
    RCP<Teuchos::Array<Scalar> > Denominator = Teuchos::null;

    switch (min_norm_)
    {
    case ANORM: {
      // MUEMAT mode (=paper)
      // Minimize with respect to the (A)' A norm.
      // Need to be smart here to avoid the construction of A' A
      //
      //                   diag( P0' (A' A) D^{-1} A P0)
      //   omega =   ------------------------------------------
      //             diag( P0' A' D^{-1}' ( A'  A) D^{-1} A P0)
      //
      // expensive, since we have to recalculate AP0 due to the lack of an explicit scaling routine for DinvAP0

      // calculate A * Ptent
      bool doFillComplete=true;
      bool optimizeStorage=false;
      RCP<Operator> AP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

      // compute A * D^{-1} * A * P0
      RCP<Operator> ADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);

      Numerator = MultiplyAll(AP0, ADinvAP0, GID2localgid);
      Denominator = MultiplySelfAll(ADinvAP0, GID2localgid);
      }
      break;
    case L2NORM: {
      // ML mode 1 (cheapest)
      // Minimize with respect to L2 norm
      //                  diag( P0' D^{-1} A P0)
      //   omega =   -----------------------------
      //             diag( P0' A' D^{-1}' D^{-1} A P0)
      //
      Numerator = MultiplyAll(Ptent, DinvAPtent, GID2localgid);
      Denominator = MultiplySelfAll(DinvAPtent, GID2localgid);
      }
      break;
    case DINVANORM: {
      // ML mode 2
      // Minimize with respect to the (D^{-1} A)' D^{-1} A norm.
      // Need to be smart here to avoid the construction of A' A
      //
      //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
      //   omega =   --------------------------------------------------------
      //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
      //

      // compute D^{-1} * A * D^{-1} * A * P0
      bool doFillComplete=true;
      bool optimizeStorage=false;
      RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);
      Utils::MyOldScaleMatrix(DinvADinvAP0,diagA,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

      Numerator = MultiplyAll(DinvAPtent, DinvADinvAP0, GID2localgid);
      Denominator = MultiplySelfAll(DinvADinvAP0, GID2localgid);
      }
      break;
    case ATDINVTPLUSDINVANORM: {
      // ML mode 3 (most expensive)
      //             diag( P0' ( A'D' + DA) D A P0)
      //   omega =   -----------------------------
      //             diag( P0'A'D' ( A'D' + DA) D A P0)
      //
      //             diag( DinvAP0'DinvAP0 + P0'DinvADinvAP0)
      //         =   -----------------------------
      //                2*diag( DinvADinvAP0'DinvAP0)
      //
      //

      // compute D^{-1} * A * D^{-1} * A * P0
      bool doFillComplete=true;
      bool optimizeStorage=false;
      RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);
      Utils::MyOldScaleMatrix(DinvADinvAP0,diagA,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

      Numerator = MultiplyAll(Ptent, DinvADinvAP0, GID2localgid);
      RCP<Teuchos::Array<Scalar> > Numerator2= MultiplySelfAll(DinvAPtent, GID2localgid);
      TEUCHOS_TEST_FOR_EXCEPTION(Numerator->size() != Numerator2->size(), Exceptions::RuntimeError, "PgPFactory::ComputeRowBasedOmegas: size of Numerator and Numerator2 different. Error");
      for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
        (*Numerator)[i] += (*Numerator2)[i];
      Denominator = MultiplyAll(DinvAPtent,DinvADinvAP0, GID2localgid);
      for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
        (*Denominator)[i] *= 2.;

      }
      break;
    }

    /////////////////// DEBUG: check for zeros in denominator
    size_t zeros_in_denominator = 0;
    for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
    {
      if((*Denominator)[i] == Teuchos::ScalarTraits<Scalar>::zero()) zeros_in_denominator ++;
    }
    if(zeros_in_denominator>Teuchos::ScalarTraits<Scalar>::zero())
      GetOStream(Warnings0, 0) << "Found " << zeros_in_denominator<< " zeros in Denominator. very suspicious!" << std::endl;

    /////////////////// build ColBasedOmegas
    RCP<Teuchos::ArrayRCP<Scalar> > ColBasedOmegas = Teuchos::rcp(new Teuchos::ArrayRCP<Scalar>(Numerator->size(),Teuchos::ScalarTraits<Scalar>::zero()));
    for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
    {
      (*ColBasedOmegas)[i] = (*Numerator)[i]/(*Denominator)[i];
      if((*ColBasedOmegas)[i] < Teuchos::ScalarTraits<Scalar>::zero())
        (*ColBasedOmegas)[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }

    /////////////////// transform ColBasedOmegas to row based omegas (local ids)
    TransformCol2RowBasedOmegas(ColBasedOmegas, GID2localgid, DinvAPtent, RowBasedOmegas);
  }


  //! @brief Transform column based omegas to row based omegas
  //!
  //! This routine transforms the local replicated column based omegas to local row based omegas.
  //! @note all processors have the same ColBasedOmegas (locally replicated). Each processor then uses
  //! the information from Op to calculate its (local) row based omegas.
  //!
  //! @param[in] const RCP<ArrayRCP<Scalar> >& ColBasedOmegas: array with column based omegas
  //! @param[in] const map<GlobalOrdinal,GlobalOrdinal>& GID2localgid: maps global IDs to local (=global) ids of ColBasedOmegas
  //! @param[in] const RCP<const Operator>& Op: operator with information how to transform col based data to row based data
  //! @param[out] ArrayRCP<Scalar>& RowBasedOmegas: vector of row based omegas
  void TransformCol2RowBasedOmegas(const RCP<Teuchos::ArrayRCP<Scalar> >& ColBasedOmegas, const std::map<GlobalOrdinal,GlobalOrdinal>& GID2localgid, const RCP<const Operator>& Op, Teuchos::ArrayRCP<Scalar>& RowBasedOmegas) const
  {
    // build RowBasedOmegas

    // mark all local row based omegas as "not assigned"
    RowBasedOmegas.assign(Op->getNodeNumRows(), -666*Teuchos::ScalarTraits<Scalar>::one());
    //TEUCHOS_TEST_FOR_EXCEPTION(RowBasedOmegas.size() != Teuchos::as<size_t>(Op->getNodeNumRows()), Exceptions::RuntimeError, "MueLu::PgPFactory::TransformCol2RowBasedOmegas: initial assignement of RowBasedOmegas failed. Error.");

    // loop over local rows of Op and extract row information
    for(size_t row=0; row<Op->getNodeNumRows(); row++)
    {
      size_t nnz = Op->getNumEntriesInLocalRow(row);

      Teuchos::ArrayView<const LocalOrdinal> indices;
      Teuchos::ArrayView<const Scalar> vals;
      Op->getLocalRowView(row, indices, vals);

      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<size_t>(indices.size()) != nnz, Exceptions::RuntimeError, "MueLu::PgPFactory::TransformCol2RowBasedOmegas: number of nonzeros not equal to number of indices? Error.");

      if(nnz == 0)
      {
        RowBasedOmegas[row] = Teuchos::ScalarTraits<Scalar>::zero();
      }
      else if(nnz == 1)
      {
        // dirichlet dofs
        RowBasedOmegas[row] = Teuchos::ScalarTraits<Scalar>::zero();
      }
      else
      {
        for(size_t j=0; j<nnz; j++)
        {
          GlobalOrdinal col_gid = Op->getColMap()->getGlobalElement(indices[j]);
          GlobalOrdinal localgid = GID2localgid.at(col_gid);
          Scalar omega = (*ColBasedOmegas)[localgid];
          if(RowBasedOmegas[row] == -666.0)
            RowBasedOmegas[row] = omega;
          else if (omega < RowBasedOmegas[row])
            RowBasedOmegas[row] = omega;
        }
      }
    }
    // check for negative entries in RowBasedOmegas
    for(size_t i=0; i<Teuchos::as<size_t>(RowBasedOmegas.size()); i++)
    {
      if(RowBasedOmegas[i] < Teuchos::ScalarTraits<Scalar>::zero()) RowBasedOmegas[i] = Teuchos::ScalarTraits<Scalar>::zero();
    }
  }

  RCP<Teuchos::Array<Scalar> > MultiplySelfAll(const RCP<Operator>& Op, const std::map<GlobalOrdinal,GlobalOrdinal>& GID2localgid) const
  {
    Teuchos::Array<Scalar> InnerProd_local(GID2localgid.size(),Teuchos::ScalarTraits<Scalar>::zero());

    for(size_t n=0; n<Op->getNodeNumRows(); n++)
    {
      Teuchos::ArrayView<const LocalOrdinal> lindices;
      Teuchos::ArrayView<const Scalar> lvals;
      Op->getLocalRowView(n, lindices, lvals);

      for(size_t i=0; i<Teuchos::as<size_t>(lindices.size()); i++)
      {
        GlobalOrdinal gid = Op->getColMap()->getGlobalElement(lindices[i]);
        GlobalOrdinal localgid = GID2localgid.at(gid);
        InnerProd_local[localgid] += lvals[i]*lvals[i];
      }
    }
    RCP<Teuchos::Array<Scalar> > InnerProd = rcp(new Teuchos::Array<Scalar>(GID2localgid.size(),Teuchos::ScalarTraits<Scalar>::zero()));

    /////////////// sum up all values to global
    Teuchos::reduceAll(*(Op->getRowMap()->getComm()),Teuchos::REDUCE_SUM, Teuchos::as<const GlobalOrdinal>(InnerProd->size()) ,&InnerProd_local[0], &(*InnerProd)[0]);

    return InnerProd;
  }


  RCP<Teuchos::Array<Scalar> > MultiplyAll(const RCP<Operator>& left, const RCP<Operator>& right, const std::map<GlobalOrdinal,GlobalOrdinal>& GID2localgid) const
  {


    TEUCHOS_TEST_FOR_EXCEPTION(!left->getDomainMap()->isSameAs(*right->getDomainMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: domain maps of left and right do not match. Error.");
    TEUCHOS_TEST_FOR_EXCEPTION(!left->getRowMap()->isSameAs(*right->getRowMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: row maps of left and right do not match. Error.");

    // test
    /*Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::RCP<Xpetra::Operator<double> > leftT = MueLu::Utils2<double>::Transpose(left,false);
    Teuchos::RCP<Xpetra::Operator<double> > leftTright = MueLu::Utils<double>::TwoMatrixMultiply(leftT,false,right,false);
    leftTright->describe(*fos,Teuchos::VERB_EXTREME);*/

    Teuchos::Array<Scalar> InnerProd_local(GID2localgid.size(),Teuchos::ScalarTraits<Scalar>::zero());

    for(size_t n=0; n<left->getNodeNumRows(); n++)
    {
      Teuchos::ArrayView<const LocalOrdinal> lindices_left;
      Teuchos::ArrayView<const Scalar> lvals_left;
      left->getLocalRowView(n, lindices_left, lvals_left);

      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;
      right->getLocalRowView(n, lindices_right, lvals_right);

      for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
      {
        for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
        {
          GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
          GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
          if(left_gid == right_gid)
          {
            GlobalOrdinal left_localgid = GID2localgid.at(left_gid);
            InnerProd_local[left_localgid] += lvals_left[i]*lvals_right[j];
            break; // skip remaining gids of right operator
          }
        }
      }
    }

    RCP<Teuchos::Array<Scalar> > InnerProd = rcp(new Teuchos::Array<Scalar>(GID2localgid.size(),Teuchos::ScalarTraits<Scalar>::zero()));

    /////////////// sum up all values to global
    Teuchos::reduceAll(*(left->getRowMap()->getComm()),Teuchos::REDUCE_SUM, Teuchos::as<const GlobalOrdinal>(InnerProd->size()) ,&InnerProd_local[0], &(*InnerProd)[0]);

    return InnerProd;
  }


  /*Teuchos::RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >*/
  void BuildLocalReplicatedColMap(const RCP<Operator>& DinvAP0,std::map<GlobalOrdinal, GlobalOrdinal>& GID2localgid) const
  {
    //Teuchos::RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::RCP< const Teuchos::Comm< int > > comm = DinvAP0->getRangeMap()->getComm();

    // generate allreduced column-based array vector for local damping factors
    Teuchos::Array<GlobalOrdinal> localreplicatedcolgids;
    reduceAllXpetraMap(localreplicatedcolgids,*(DinvAP0->getDomainMap()));

    // map GID -> localgid
    GlobalOrdinal curLocalGid = 0;

    for (size_t i = 0; i<Teuchos::as<size_t>(localreplicatedcolgids.size()); i++)
    {
      GID2localgid[localreplicatedcolgids[i]] = curLocalGid;
      //std::cout << "GID " << localreplicatedcolgids[i] << " localgid " << GID2localgid[localreplicatedcolgids[i]] << std::endl;
      curLocalGid++;
    }

  }

  //! @name helper function for allreducing a Xpetra::Map
  //@{

  GlobalOrdinal FindMyPos(size_t nummyelements, const Teuchos::Comm<int>& comm) const
  {
    const int myrank = comm.getRank();
    const int numproc = comm.getSize();

    std::vector<size_t> snum (numproc,0); // local vector of no of elements
    std::vector<size_t> rnum (numproc);   // gobal vector of no of elements on proc
    snum[myrank] = nummyelements;

    Teuchos::reduceAll(comm, Teuchos::REDUCE_SUM, numproc, &snum[0], &rnum[0]);

    return std::accumulate(&rnum[0], &rnum[myrank], 0);
  }

  void reduceAllXpetraMap(Teuchos::Array<GlobalOrdinal>& rredundant, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const
  {
    const GlobalOrdinal mynodepos = FindMyPos(map.getNodeNumElements(), *(map.getComm()));

    std::vector<GlobalOrdinal> sredundant(map.getGlobalNumElements(),0);

    Teuchos::ArrayView< const GlobalOrdinal > gids = map.getNodeElementList();
    std::copy(gids.getRawPtr(), gids.getRawPtr() + map.getNodeNumElements(), &sredundant[mynodepos]); // check me

    rredundant.resize(map.getGlobalNumElements());
    GlobalOrdinal numEle = map.getGlobalNumElements();
    Teuchos::reduceAll(*(map.getComm()),Teuchos::REDUCE_SUM, numEle ,&sredundant[0], &rredundant[0]);
  }

  //@}


private:

  //! Input factories
  RCP<PFactory> initialPFact_;        //! Ptentative Factory
  RCP<SingleLevelFactoryBase> AFact_; //! A Factory

  //! Factory parameters
  std::string diagonalView_;

  //! minimization norm
  MinimizationNorm min_norm_;
};

//! Friend print function.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, PgPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
  os << "Printing an PgPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_PGPFACTORY_SHORT

#endif /* MUELU_PGPFACTORY_HPP_ */
