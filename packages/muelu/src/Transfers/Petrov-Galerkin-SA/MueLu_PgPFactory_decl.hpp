// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_PGPFACTORY_DECL_HPP_
#define MUELU_PGPFACTORY_DECL_HPP_

#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_Import_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>
#include <Xpetra_Export_fwd.hpp>
#include <Xpetra_ExportFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_PerfUtils_fwd.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities_fwd.hpp"

namespace MueLu {

/* Options defining how to pick-up the next root node in the local aggregation procedure */
enum MinimizationNorm {
  ANORM     = 0, /* A norm   */
  L2NORM    = 1, /* L2 norm */
  DINVANORM = 2  /* Dinv A norm */
};

/*!
  @class PgPFactory class.
  @brief Factory for building Petrov-Galerkin Smoothed Aggregation prolongators.
  @ingroup MueLuTransferClasses
*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class PgPFactory : public PFactory {
#undef MUELU_PGPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.
    User can supply a factory for generating the tentative prolongator.
  */
  PgPFactory() {}

  //! Destructor.
  virtual ~PgPFactory() {}

  //@}

  //! @name Set methods.
  //@{

  RCP<const ParameterList> GetValidParameterList() const;

  //! Set minimization mode (L2NORM for cheapest, ANORM more expensive, DINVANORM = default)
  void SetMinimizationMode(MinimizationNorm minnorm);

  //! return minimization mode
  MinimizationNorm GetMinimizationMode();
  //@}

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;
  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
  */
  void Build(Level& fineLevel, Level& coarseLevel) const;

  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

  void ReUseDampingParameters(bool bReuse);

 private:
  void MultiplySelfAll(const RCP<Matrix>& Op, Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& InnerProdVec) const;

  void MultiplyAll(const RCP<Matrix>& left, const RCP<Matrix>& right, Teuchos::RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& InnerProdVec) const;

  void ComputeRowBasedOmega(Level& fineLevel, Level& coarseLevel, const RCP<Matrix>& A, const RCP<Matrix>& P0, const RCP<Matrix>& DinvAP0, RCP<Xpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& RowBasedOmega) const;

 private:
  //! Factory parameters
  std::string diagonalView_;  // TODO do we need this?
};

}  // namespace MueLu

#define MUELU_PGPFACTORY_SHORT
#endif /* MUELU_PGPFACTORY_DECL_HPP_ */
