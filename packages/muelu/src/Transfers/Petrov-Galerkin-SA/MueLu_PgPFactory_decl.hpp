#ifndef MUELU_PGPFACTORY_DECL_HPP
#define MUELU_PGPFACTORY_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

/*
 * MueLu_PgPFactory.hpp
 *
 *  Created on: 23.09.2011
 *      Author: tobias
 */

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

public:

  /* Options defining how to pick-up the next root node in the local aggregation procedure */
  enum MinimizationNorm ;;

  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.
      User can supply a factory for generating the tentative prolongator.
   */
  PgPFactory(RCP<PFactory> InitialPFact = Teuchos::null, RCP<SingleLevelFactoryBase> AFact = Teuchos::null)
  : initialPFact_(InitialPFact), AFact_(AFact),
    diagonalView_("current") ;

  //! Destructor.
  virtual ~PgPFactory() ;

  //@}

  //! @name Set methods.
  //@{

  //! Change view of diagonal.
  void SetDiagonalView(std::string const& diagView) ;

  //! Set minimization mode (L2NORM for cheapest, ANORM more expensive, DINVANORM = default, ATDINVTPLUSDINVANORM most expensive method)
  inline void SetMinimizationMode(MinimizationNorm minnorm) ;
  //@}

  //! @name Get methods.
  //@{

  //! Returns current view of diagonal.
  std::string GetDiagonalView() ;

  //@}

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const ;
  //@}

  //! @name Build methods.
  //@{

  /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
   */
  void Build(Level& fineLevel, Level &coarseLevel) const ;

  void BuildP(Level &fineLevel, Level &coarseLevel) const ;

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
  ;


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
  ;

  RCP<Teuchos::Array<Scalar> > MultiplySelfAll(const RCP<Operator>& Op, const std::map<GlobalOrdinal,GlobalOrdinal>& GID2localgid) const
  ;


  RCP<Teuchos::Array<Scalar> > MultiplyAll(const RCP<Operator>& left, const RCP<Operator>& right, const std::map<GlobalOrdinal,GlobalOrdinal>& GID2localgid) const
  ;


  /*Teuchos::RCP<const Xpetra::Map< LocalOrdinal, GlobalOrdinal, Node > >*/
  void BuildLocalReplicatedColMap(const RCP<Operator>& DinvAP0,std::map<GlobalOrdinal, GlobalOrdinal>& GID2localgid) const
  ;

  //! @name helper function for allreducing a Xpetra::Map
  //@{

  GlobalOrdinal FindMyPos(size_t nummyelements, const Teuchos::Comm<int>& comm) const
  ;

  void reduceAllXpetraMap(Teuchos::Array<GlobalOrdinal>& rredundant, const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& map) const
  ;

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

} //namespace MueLu

#define MUELU_PGPFACTORY_SHORT
#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_PGPFACTORY_DECL_HPP
