// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GlobalEvaluationData.hpp"

// Thyra
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_VectorBase.hpp"

namespace panzer
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  Default Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  BlockedVector_ReadOnly_GlobalEvaluationData::
  BlockedVector_ReadOnly_GlobalEvaluationData()
    : isInitialized_(false)
  {
  } // end of Default Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Copy Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  BlockedVector_ReadOnly_GlobalEvaluationData::
  BlockedVector_ReadOnly_GlobalEvaluationData(
    const BlockedVector_ReadOnly_GlobalEvaluationData& src)
    :
    isInitialized_(false)
  {
    initialize(src.ghostedSpace_, Teuchos::null, src.gedBlocks_);
  } // end of Copy Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Initializing Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  BlockedVector_ReadOnly_GlobalEvaluationData::
  BlockedVector_ReadOnly_GlobalEvaluationData(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ghostedSpace,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ownedSpace,
    const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>>&
      gedBlocks)
    :
    isInitialized_(false)
  {
    initialize(ghostedSpace, ownedSpace, gedBlocks);
  } // end of Initializing Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  initialize()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  BlockedVector_ReadOnly_GlobalEvaluationData::
  initialize(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>>& ghostedSpace,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>>& /* ownedSpace */,
    const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>>&
      gedBlocks)
  {
    using std::logic_error;
    using std::size_t;
    using Teuchos::rcp_dynamic_cast;
    using Thyra::DefaultProductVectorSpace;

    // Assert that all the gedBlocks are initialized.
    for (size_t i(0); i < gedBlocks.size(); ++i)
      TEUCHOS_TEST_FOR_EXCEPTION(not gedBlocks[i]->isInitialized(),
        logic_error, "BlockedVector_ReadOnly_GlobalEvaluationData::"          \
        "initialize:  GED block " << i << " is not initialized.")
    gedBlocks_    = gedBlocks;
    ghostedSpace_ =
      rcp_dynamic_cast<const DefaultProductVectorSpace<double>>(ghostedSpace);
    TEUCHOS_TEST_FOR_EXCEPTION(ghostedSpace_.is_null(), logic_error,
      "BlockedVector_ReadOnly_GlobalEvaluationData::initialize():  Ghosted "  \
      "space must be a Thyra::DefaultProductVectorSpace");
    isInitialized_ = true;
  } // end of initialize()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  globalToGhost()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  BlockedVector_ReadOnly_GlobalEvaluationData::
  globalToGhost(
    int mem)
  {
    using std::logic_error;
    using std::size_t;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "BlockedVector_ReadOnly_GlobalEvaluationData has not been "             \
      "initialized; cannot call \"globalToGhost()\"!");
    for (size_t i(0); i < gedBlocks_.size(); ++i)
      gedBlocks_[i]->globalToGhost(mem);
  } // end of globalToGhost()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  initializeData()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  BlockedVector_ReadOnly_GlobalEvaluationData::
  initializeData()
  {
    using std::logic_error;
    using std::size_t;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "BlockedVector_ReadOnly_GlobalEvaluationData has not been "             \
      "initialized; cannot call \"initializeData()\"!");
    for (size_t i(0); i < gedBlocks_.size(); ++i)
      gedBlocks_[i]->initializeData();
  } // end of initializeData()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  setOwnedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  BlockedVector_ReadOnly_GlobalEvaluationData::
  setOwnedVector(
    const Teuchos::RCP<const Thyra::VectorBase<double>>& ownedVector)
  {
    using std::logic_error;
    using std::size_t;
    using Teuchos::as;
    using Teuchos::RCP;
    using Thyra::castOrCreateProductVectorBase;
    using Thyra::ProductVectorBase;
    ownedVector_ = ownedVector;
    RCP<const ProductVectorBase<double>> blocks =
      castOrCreateProductVectorBase(ownedVector_);
    TEUCHOS_TEST_FOR_EXCEPTION(blocks->productSpace()->numBlocks() !=
      as<int>(gedBlocks_.size()), logic_error,
      "BlockedVector_ReadOnly_GlobalEvaluationData owned vector has the "     \
      "wrong number of blocks!");
    for (size_t i(0); i < gedBlocks_.size(); ++i)
      gedBlocks_[i]->setOwnedVector(blocks->getVectorBlock(i));
  } // end of setOwnedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getOwnedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<const Thyra::VectorBase<double>>
  BlockedVector_ReadOnly_GlobalEvaluationData::
  getOwnedVector() const
  {
    return ownedVector_;
  } // end of getOwnedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getGhostedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Thyra::VectorBase<double>>
  BlockedVector_ReadOnly_GlobalEvaluationData::
  getGhostedVector() const
  {
    using std::logic_error;
    using std::size_t;
    using std::vector;
    using Teuchos::arrayViewFromVector;
    using Teuchos::RCP;
    using Thyra::defaultProductVector;
    using Thyra::VectorBase;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "BlockedVector_ReadOnly_GlobalEvaluationData has not been "             \
      "initialized; cannot call \"getGhostedVector()\"!");
    vector<RCP<VectorBase<double>>> blocks;
    for (size_t i(0); i < gedBlocks_.size(); ++i)
      blocks.push_back(gedBlocks_[i]->getGhostedVector());
    const vector<RCP<VectorBase<double>>>& constBlocks = blocks;
    return defaultProductVector(ghostedSpace_,
      arrayViewFromVector(constBlocks));
  } // end of getGhostedVector()

} // end of namespace panzer

// end of Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.cpp