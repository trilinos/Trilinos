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
#include "Panzer_BlockedVector_Write_GlobalEvaluationData.hpp"
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
  BlockedVector_Write_GlobalEvaluationData::
  BlockedVector_Write_GlobalEvaluationData()
    : isInitialized_(false)
  {
  } // end of Default Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Copy Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  BlockedVector_Write_GlobalEvaluationData::
  BlockedVector_Write_GlobalEvaluationData(
    const BlockedVector_Write_GlobalEvaluationData& src)
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
  BlockedVector_Write_GlobalEvaluationData::
  BlockedVector_Write_GlobalEvaluationData(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ghostedSpace,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ownedSpace,
    const std::vector<Teuchos::RCP<WriteVector_GlobalEvaluationData>>&
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
  BlockedVector_Write_GlobalEvaluationData::
  initialize(
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>>& ghostedSpace,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double>>& /* ownedSpace */,
    const std::vector<Teuchos::RCP<WriteVector_GlobalEvaluationData>>&
      gedBlocks)
  {
    using std::logic_error;
    using std::size_t;
    using Teuchos::rcp_dynamic_cast;
    using Thyra::DefaultProductVectorSpace;

    // Assert that all the gedBlocks are initialized.
    for (size_t i(0); i < gedBlocks.size(); ++i)
      TEUCHOS_TEST_FOR_EXCEPTION(not gedBlocks[i]->isInitialized(),
        logic_error, "BlockedVector_Write_GlobalEvaluationData::"             \
        "initialize:  GED block " << i << " is not initialized.")
    gedBlocks_    = gedBlocks;
    ghostedSpace_ =
      rcp_dynamic_cast<const DefaultProductVectorSpace<double>>(ghostedSpace);
    TEUCHOS_TEST_FOR_EXCEPTION(ghostedSpace_.is_null(), logic_error,
      "BlockedVector_Write_GlobalEvaluationData::initialize():  Ghosted "     \
      "space must be a Thyra::DefaultProductVectorSpace");
    isInitialized_ = true;
  } // end of initialize()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  ghostToGlobal()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  BlockedVector_Write_GlobalEvaluationData::
  ghostToGlobal(
    int mem)
  {
    using std::logic_error;
    using std::size_t;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "BlockedVector_Write_GlobalEvaluationData has not been "                \
      "initialized; cannot call \"ghostToGlobal()\"!");
    for (size_t i(0); i < gedBlocks_.size(); ++i)
      gedBlocks_[i]->ghostToGlobal(mem);
  } // end of ghostToGlobal()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  initializeData()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  BlockedVector_Write_GlobalEvaluationData::
  initializeData()
  {
    using std::logic_error;
    using std::size_t;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "BlockedVector_Write_GlobalEvaluationData has not been "                \
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
  BlockedVector_Write_GlobalEvaluationData::
  setOwnedVector(
    const Teuchos::RCP<Thyra::VectorBase<double>>& ownedVector)
  {
    using std::logic_error;
    using std::size_t;
    using Teuchos::as;
    using Teuchos::RCP;
    using Thyra::castOrCreateNonconstProductVectorBase;
    using Thyra::ProductVectorBase;
    ownedVector_ = ownedVector;
    RCP<ProductVectorBase<double>> blocks =
      castOrCreateNonconstProductVectorBase(ownedVector_);
    TEUCHOS_TEST_FOR_EXCEPTION(blocks->productSpace()->numBlocks() !=
      as<int>(gedBlocks_.size()), logic_error,
      "BlockedVector_Write_GlobalEvaluationData owned vector has the "        \
      "wrong number of blocks!");
    for (size_t i(0); i < gedBlocks_.size(); ++i)
      gedBlocks_[i]->setOwnedVector(blocks->getNonconstVectorBlock(i));
  } // end of setOwnedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getOwnedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Thyra::VectorBase<double>>
  BlockedVector_Write_GlobalEvaluationData::
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
  BlockedVector_Write_GlobalEvaluationData::
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
      "BlockedVector_Write_GlobalEvaluationData has not been "                \
      "initialized; cannot call \"getGhostedVector()\"!");
    vector<RCP<VectorBase<double>>> blocks;
    for (size_t i(0); i < gedBlocks_.size(); ++i)
      blocks.push_back(gedBlocks_[i]->getGhostedVector());
    const vector<RCP<VectorBase<double>>>& constBlocks = blocks;
    return defaultProductVector(ghostedSpace_,
      arrayViewFromVector(constBlocks));
  } // end of getGhostedVector()

} // end of namespace panzer

// end of Panzer_BlockedVector_Write_GlobalEvaluationData.cpp