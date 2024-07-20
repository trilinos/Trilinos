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

// Epetra
#include "Epetra_Import.h"

// Panzer
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"

// Thyra
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

namespace panzer
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  useConstantValues()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  useConstantValues(
    const std::vector<int>& indices,
    double                  value)
  {
    using std::logic_error;
    TEUCHOS_TEST_FOR_EXCEPTION(isInitialized_, logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData has been initialized; "     \
      "cannot call \"useConstantValues()\"!");

    // Add this specification to the filtered pairs vector.
    FilteredPair fp;
    fp.first  = indices;
    fp.second = value;
    filteredPairs_.push_back(fp);
  } // end of useConstantValues()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  initialize()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  initialize(
    const Teuchos::RCP<const Epetra_Import>& importer,
    const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
    const Teuchos::RCP<const Epetra_Map>&    ownedMap)
  {
    using panzer::kokkos_utils::getView;
    using std::size_t;
    using std::vector;
    using Teuchos::rcp;
    using Thyra::create_Vector;
    using Thyra::create_VectorSpace;

    // Save the input.
    importer_   = importer;
    ghostedMap_ = ghostedMap;
    ownedMap_   = ownedMap;

    // Build up the Thyra conversion data structures.
    ghostedSpace_ = create_VectorSpace(ghostedMap_);
    ownedSpace_   = create_VectorSpace(ownedMap_  );

    // Allocate the vectors.
    ghostedVector_ = rcp(new Epetra_Vector(*ghostedMap_));
    auto ownedVector = rcp(new Epetra_Vector(*ownedMap_));
    ownedVector_ = create_Vector(ownedVector, ownedSpace_);

    // Translate filtered pair GIDs to LIDs and initialize some ghosted values
    // to the user-specified values.
    for (size_t i(0); i < filteredPairs_.size(); ++i)
    {
      vector<int> lids;
      const vector<int>& gids = filteredPairs_[i].first;
      for (size_t j(0); j < gids.size(); ++j)
      {
        // Add legitimate LIDs to the list.
        int lid = ghostedMap->LID(gids[j]);
        if (lid >= 0)
          lids.push_back(lid);
      } // end loop over gids

      // Overwrite the original GID vector with the new LID vector.
      filteredPairs_[i].first = lids;
    } // end loop over filteredPairs_
    isInitialized_ = true;

    // Get the PHX::Views corresponding to the owned and ghosted vectors.
    ownedView_   = getView<const Epetra_Vector>(*ownedVector_);
    ghostedView_ = getView<Epetra_Vector>(*getGhostedVector());
  } // end of initialize()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  globalToGhost()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  globalToGhost(
    int /* mem */)
  {
    using std::logic_error;
    using Teuchos::RCP;
    using Thyra::get_Epetra_Vector;
    TEUCHOS_TEST_FOR_EXCEPTION(ownedVector_.is_null(), logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::globalToGhost():  Owned "  \
      "vector has not been set; can't perform the halo exchange!")

    // Initialize the ghosted data, zeroing out things, and filling in
    // specified constants.
    initializeData();
    RCP<const Epetra_Vector> ownedVector_ep =
      get_Epetra_Vector(*ownedMap_, ownedVector_);

    // Do the global distribution.
    ghostedVector_->Import(*ownedVector_ep, *importer_, Insert);
  } // end of globalToGhost()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  initializeData()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  initializeData()
  {
    using std::logic_error;
    using std::size_t;
    using std::vector;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData has not been initialized, " \
      "cannot call \"initializeData()\"!");
    ghostedVector_->PutScalar(0);

    // Initialize some ghosted values to the user-specified values.
    for (size_t i(0); i < filteredPairs_.size(); ++i)
    {
      const vector<int>& lids = filteredPairs_[i].first;
      for (size_t j(0); j < lids.size(); ++j)
        (*ghostedVector_)[lids[j]] = filteredPairs_[i].second;
    } // end loop over filteredPairs_
  } // end of initializeData()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  ghostToGlobal()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  ghostToGlobal(
    int /* mem = 0 */)
  {
  } // end of ghostToGlobal()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  setOwnedVector_Epetra()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  setOwnedVector_Epetra(
    const Teuchos::RCP<const Epetra_Vector>& ownedVector)
  {
    using panzer::kokkos_utils::getView;
    using std::logic_error;
    using Thyra::create_Vector;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::"                          \
      "setOwnedVector_Epetra():  This object hasn't yet been initialized.")
    ownedVector_ = create_Vector(ownedVector, ownedSpace_);
    ownedView_   = getView<const Epetra_Vector>(*ownedVector_);
  } // end of setOwnedVector_Epetra()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getGhostedVector_Epetra()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Epetra_Vector>
  EpetraVector_ReadOnly_GlobalEvaluationData::
  getGhostedVector_Epetra() const
  {
    using std::logic_error;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::"                          \
      "getGhostedVector_Epetra():  This object hasn't yet been initialized.")
    TEUCHOS_TEST_FOR_EXCEPTION(ghostedVector_.is_null(), logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::"                          \
      "getGhostedVector_Epetra():  The ghosted vector is just a null RCP.")
    return ghostedVector_;
  } // end of getGhostedVector_Epetra()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  setOwnedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  setOwnedVector(
    const Teuchos::RCP<const Thyra::VectorBase<double>>& ownedVector)
  {
    using panzer::kokkos_utils::getView;
    using std::logic_error;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::setOwnedVector():  This "  \
      "object hasn't yet been initialized.")
    ownedVector_ = ownedVector;
    ownedView_   = getView<const Epetra_Vector>(*ownedVector_);
  } // end of setOwnedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getOwnedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<const Thyra::VectorBase<double>>
  EpetraVector_ReadOnly_GlobalEvaluationData::
  getOwnedVector() const
  {
    using std::logic_error;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::getOwnedVector():  This "  \
      "object hasn't yet been initialized.")
    return ownedVector_;
  } // end of getOwnedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getGhostedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Thyra::VectorBase<double>>
  EpetraVector_ReadOnly_GlobalEvaluationData::
  getGhostedVector() const
  {
    using std::logic_error;
    using Thyra::create_Vector;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::getGhostedVector():  "     \
      "This object hasn't yet been initialized.")
    TEUCHOS_TEST_FOR_EXCEPTION(ghostedVector_.is_null(), logic_error,
      "EpetraVector_ReadOnly_GlobalEvaluationData::getGhostedVector():  The " \
      "ghosted vector is just a null RCP.")
    return create_Vector(ghostedVector_, ghostedSpace_);
  } // end of getGhostedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  print()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_ReadOnly_GlobalEvaluationData::
  print(
    std::ostream& os) const
  {
    using std::endl;
    using std::string;
    const string tab("    ");
    os << endl
       << tab << "EpetraVector_ReadOnly_GlobalEvaluationData" << endl
       << tab << "  init    = " << isInitialized_             << endl
       << tab << "  owned   = " << ownedVector_               << endl
       << tab << "  ghosted = " << ghostedVector_             << endl;
  } // end of print()

} // end of namespace panzer

// end of Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.cpp