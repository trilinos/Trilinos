// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BRICKAGGREGATIONFACTORY_DECL_HPP_
#define MUELU_BRICKAGGREGATIONFACTORY_DECL_HPP_

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_BrickAggregationFactory_fwd.hpp"

#include "MueLu_LWGraph_fwd.hpp"

#include "MueLu_LWGraph_fwd.hpp"
#include "MueLu_Level_fwd.hpp"
#include "MueLu_Aggregates_fwd.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities_fwd.hpp"

/*!
  @class BrickAggregationFactory
  @brief Aggregation method for generating "brick" aggregates.  It also does "hotdogs" and "pancakes."

  This factory can generate aggregates of size 1, 2 or 3 in each dimension, in any combination.

*/

namespace MueLu {

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class BrickAggregationFactory : public SingleLevelFactoryBase {
#undef MUELU_BRICKAGGREGATIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"
 private:
  typedef Teuchos::ScalarTraits<Scalar> STS;

  // Comparator for doubles
  // Generally, the coordinates for coarser levels would come out of averaging of fine level coordinates
  // It is possible that the result of the averaging differs slightly between clusters, as we might have
  // 3x2 and 2x2 cluster which would result in averaging 6 and 4 y-coordinates respectively, leading to
  // slightly different results.
  // Therefore, we hardcode a constant so that close points are considered the same.
  class compare {
   public:
    bool operator()(const Scalar& x, const Scalar& y) const {
      if (STS::magnitude(x - y) < 1e-14)
        return false;
      return STS::real(x) < STS::real(y);
    }
  };
  typedef std::map<Scalar, GlobalOrdinal, compare> container;

 public:
  //! @name Constructors/Destructors.
  //@{

  //! Constructor.
  BrickAggregationFactory()
    : nDim_(-1)
    , nx_(-1)
    , ny_(-1)
    , nz_(-1)
    , bx_(-1)
    , by_(-1)
    , bz_(-1){};

  //! Destructor.
  virtual ~BrickAggregationFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  // Options shared by all aggregation algorithms

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*! @brief Build aggregates. */
  void Build(Level& currentLevel) const;

  //@}

 private:
  void Setup(const RCP<const Teuchos::Comm<int> >& comm, const RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::magnitudeType, LO, GO, NO> >& coords, const RCP<const Map>& map) const;
  RCP<container> Construct1DMap(const RCP<const Teuchos::Comm<int> >& comm, const ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType>& x) const;

  void BuildGraph(Level& currentLevel, const RCP<Matrix>& A) const;

  bool isDirichlet(LocalOrdinal LID) const;
  bool isRoot(LocalOrdinal LID) const;
  GlobalOrdinal getRoot(LocalOrdinal LID) const;
  GlobalOrdinal getAggGID(LocalOrdinal LID) const;

  void getIJK(LocalOrdinal LID, int& i, int& j, int& k) const;
  void getAggIJK(LocalOrdinal LID, int& i, int& j, int& k) const;

  mutable int nDim_;
  mutable RCP<container> xMap_, yMap_, zMap_;
  mutable ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType> x_, y_, z_;
  mutable int nx_, ny_, nz_;
  mutable int bx_, by_, bz_;
  mutable bool dirichletX_, dirichletY_, dirichletZ_;
  mutable int naggx_, naggy_, naggz_;

  mutable std::map<GlobalOrdinal, GlobalOrdinal> revMap_;
};  // class BrickAggregationFactory

}  // namespace MueLu

#define MUELU_BRICKAGGREGATIONFACTORY_SHORT
#endif /* MUELU_BRICKAGGREGATIONFACTORY_DECL_HPP_ */
