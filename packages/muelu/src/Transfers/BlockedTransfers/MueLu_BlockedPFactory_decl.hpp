// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BLOCKEDPFACTORY_DECL_HPP_
#define MUELU_BLOCKEDPFACTORY_DECL_HPP_

#include <set>
#include <Xpetra_MapExtractorFactory_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_TentativePFactory_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"

namespace MueLu {

/*!
  @class BlockedPFactory class.
  @brief Factory for building blocked, segregated prolongation operators.

  Factory for building blocked, segregated prolongation operators of the form
  \f$ P=\diag(P_{11},P_{22},\ldots)\f$, where \f$ P_{ii}\f$ are prolongation operators
  for the corresponding subblocks \f$A_{ii}\f$ in the blocked operator \f$ A \f$.

  @param RCP<FactoryBase> AFact = Teuchos::null: factory for generating blocked operator \f$ A\f$.

  The blocked operator \f$ A \f$ is needed for accessing the underlaying blocked maps (MapExtractors).
  Use the AddFactoryManager function to define block rows in the blocked operator. For each block row you have to define
  a special factory manager that handles the prolongation and restriction operator for the corresponding blocks. The
  SubBlockAFactory class provides a factory interface for accessing specific blocks in the blocked operator \f$ A \f$.

  \code{.cpp}
  // define SubBlockAFactory objects for accessing the diagonal blocks in the blocked operator A
  RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
  RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

  // define transfer operators for first block row
  RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());
  RCP<TransPFactory> R11Fact = rcp(new TransPFactory(P11Fact));

  // define transfer operators for second block row
  RCP<TentativePFactory> P22TentFact = rcp(new TentativePFactory());
  RCP<PgPFactory> P22Fact = rcp(new PgPFactory(P22TentFact));
  RCP<GenericRFactory> R22Fact = rcp(new GenericRFactory(P22Fact));

  // setup facatory managers for first and second block row
  RCP<FactoryManager> M11 = rcp(new FactoryManager());
  M11->SetFactory("A", A11Fact);
  M11->SetFactory("P", P11Fact);
  M11->SetFactory("R", R11Fact);

  RCP<FactoryManager> M22 = rcp(new FactoryManager());
  M22->SetFactory("A", A22Fact);
  M22->SetFactory("P", P22Fact);
  M22->SetFactory("R", R22Fact);

  // setup blocked transfer operator
  // use Teuchos::null as AFactory -> standard "A" variable (in this case the blocked operator A)
  RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory());
  PFact->SetFactory("A",Teuchos::null); // not necessary (use blocked operator A)
  PFact->AddFactoryManager(M11); // add factory manager for first block row
  PFact->AddFactoryManager(M22); // add factory manager for second block row

  // define restriction operator
  // Note: always use GenericRFactory here.
  // The simple TransPFactory class cannot handle blocked operators, yet.
  RCP<GenericRFactory> RFact = rcp(new GenericRFactory(PFact));

  // RAP factory
  RCP<RAPFactory> AcFact = rcp(new RAPFactory(PFact, RFact));
  \endcode

*/

template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class BlockedPFactory : public PFactory {
#undef MUELU_BLOCKEDPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"
 public:
  //! @name Constructors/Destructors.
  //@{

  /*! @brief Constructor.
    User can supply a factory for generating the tentative prolongator.
  */
  BlockedPFactory(/*RCP<FactoryBase> AFact = Teuchos::null*/)
    : diagonalView_("current") {}

  //! Destructor.
  virtual ~BlockedPFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! @name Set methods.
  //@{

  //! Change view of diagonal.
  void SetDiagonalView(std::string const &diagView) { diagonalView_ = diagView; }

  //! Add a factory manager
  void AddFactoryManager(RCP<const FactoryManagerBase> FactManager) { FactManager_.push_back(FactManager); }

  //@}

  //! @name Get methods.
  //@{

  //! Returns current view of diagonal.
  std::string GetDiagonalView() { return diagonalView_; }

  //@}

  //! Input
  //@{

  void DeclareInput(Level &fineLevel, Level &coarseLevel) const;
  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Build segregated blocked prolongator by assembly from prolongators of each subblock in A
    and return it in \c coarseLevel .
  */
  void Build(Level &fineLevel, Level &coarseLevel) const;

  void BuildP(Level &fineLevel, Level &coarseLevel) const;

  //@}

 private:
  bool areGidsUnique(const std::vector<GO> &X) const {
    std::set<GO> Y(X.begin(), X.end());
    return X.size() == Y.size();
  }

  //! Input factories
  std::vector<Teuchos::RCP<const FactoryManagerBase> > FactManager_;

  //! Factory parameters
  std::string diagonalView_;
};

}  // namespace MueLu

#define MUELU_BLOCKEDPFACTORY_SHORT
#endif /* MUELU_BLOCKEDPFACTORY_DECL_HPP_ */
