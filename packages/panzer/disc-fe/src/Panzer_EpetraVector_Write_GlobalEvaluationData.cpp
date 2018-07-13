// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Epetra
#include "Epetra_Export.h"

// Panzer
#include "Panzer_EpetraVector_Write_GlobalEvaluationData.hpp"

// Thyra
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

namespace panzer
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  initialize()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_Write_GlobalEvaluationData::
  initialize(
    const Teuchos::RCP<const Epetra_Export>& exporter,
    const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
    const Teuchos::RCP<const Epetra_Map>&    ownedMap)
  {
    using panzer::kokkos_utils::getView;
    using Teuchos::rcp;
    using Thyra::create_Vector;
    using Thyra::create_VectorSpace;

    // Save the input.
    exporter_   = exporter;
    ghostedMap_ = ghostedMap;
    ownedMap_   = ownedMap;
  
    // Build up the Thyra conversion data structures.
    ghostedSpace_ = create_VectorSpace(ghostedMap_);
    ownedSpace_   = create_VectorSpace(ownedMap_);
  
    // Allocate the vectors.
    ghostedVector_ = rcp(new Epetra_Vector(*ghostedMap_));
    auto ownedVector = rcp(new Epetra_Vector(*ownedMap_));
    ownedVector_ = create_Vector(ownedVector, ownedSpace_);
    isInitialized_ = true;

    // Get the Kokkos::View corresponding to the ghosted vector.
    ownedView_   = getView<Epetra_Vector>(*ownedVector_);
    ghostedView_ = getView<Epetra_Vector>(*getGhostedVector());
  } // end of initialize()
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //  ghostToGlobal()
  //
  /////////////////////////////////////////////////////////////////////////////
  void 
  EpetraVector_Write_GlobalEvaluationData::
  ghostToGlobal(
    int /* mem */)
  {
    using std::invalid_argument;
    using std::logic_error;
    using Teuchos::RCP;
    using Thyra::get_Epetra_Vector;
    TEUCHOS_TEST_FOR_EXCEPTION(ownedVector_.is_null(), logic_error,
      "EpetraVector_Write_GlobalEvaluationData::ghostToGlobal():  Owned "     \
      "vector has not been set; can't perform the halo exchange!")
  
    // Set different combine modes.
    Epetra_CombineMode cm = Add;
    switch (getCombineMode())
    {
      case CM_Sum:
        cm = Add;
        break;
      case CM_Min:
        cm = Epetra_Min;
        break;
      case CM_Max:
        cm = Epetra_Max;
        break;
      case CM_Insert:
        cm = Insert;
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPTION(true, invalid_argument,
          "EpetraVector_Write_GlobalEvaluationData::ghostToGlobal():  "       \
          "Invalid CombineMode.  Valid modes are CM_Sum, CM_Max, CM_Min, "    \
          "and CM_Insert.")
    }; // end switch (getCombineMode())
    
    // Do the global distribution.
    RCP<Epetra_Vector> ownedVector_ep = get_Epetra_Vector(*ownedMap_,
      ownedVector_);
    ownedVector_ep->Export(*ghostedVector_, *exporter_, cm);
  } // end of ghostToGlobal()
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //  initializeData()
  //
  /////////////////////////////////////////////////////////////////////////////
  void 
  EpetraVector_Write_GlobalEvaluationData::
  initializeData()
  {
    using std::logic_error;
    using Thyra::put_scalar;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_Write_GlobalEvaluationData has not been initialized; "    \
      "cannot call \"initializeData()\"!")
    put_scalar(0.0, ownedVector_.ptr());
  } // end of initializeData()
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //  setOwnedVector_Epetra()
  //
  /////////////////////////////////////////////////////////////////////////////
  void 
  EpetraVector_Write_GlobalEvaluationData::
  setOwnedVector_Epetra(
    const Teuchos::RCP<Epetra_Vector>& ownedVector)
  {
    using panzer::kokkos_utils::getView;
    using std::logic_error;
    using Thyra::create_Vector;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_Write_GlobalEvaluationData::setOwnedVector_Epetra():  "   \
      "This object hasn't yet been initialized.")
    ownedVector_ = create_Vector(ownedVector, ownedSpace_);
    ownedView_   = getView<Epetra_Vector>(*ownedVector_);
  } // end of setOwnedVector_Epetra()
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //  getGhostedVector_Epetra()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Epetra_Vector> 
  EpetraVector_Write_GlobalEvaluationData::
  getGhostedVector_Epetra() const
  {
    using std::logic_error;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_Write_GlobalEvaluationData::setGhostedVector_Epetra():  " \
      "This object hasn't yet been initialized.")
    TEUCHOS_TEST_FOR_EXCEPTION(ghostedVector_.is_null(), logic_error,
      "EpetraVector_Write_GlobalEvaluationData::setGhostedVector_Epetra():  " \
      "The ghosted vector is just a null RCP.")
    return ghostedVector_;
  } // end of getGhostedVector_Epetra()
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //  setOwnedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  void 
  EpetraVector_Write_GlobalEvaluationData::
  setOwnedVector(
    const Teuchos::RCP<Thyra::VectorBase<double>>& ownedVector)
  {
    using panzer::kokkos_utils::getView;
    using std::logic_error;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_Write_GlobalEvaluationData::setOwnedVector():  This "     \
      "object hasn't yet been initialized.")
    ownedVector_ = ownedVector;
    ownedView_   = getView<Epetra_Vector>(*ownedVector_);
  } // end of setOwnedVector()
  
  /////////////////////////////////////////////////////////////////////////////
  //
  //  getOwnedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Thyra::VectorBase<double>> 
  EpetraVector_Write_GlobalEvaluationData::
  getOwnedVector() const
  {
    using std::logic_error;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_Write_GlobalEvaluationData::getOwnedVector():  This "     \
      "object hasn't yet been initialized.")
    return ownedVector_;
  } // end of getOwnedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getGhostedVector()
  //
  /////////////////////////////////////////////////////////////////////////////
  Teuchos::RCP<Thyra::VectorBase<double>> 
  EpetraVector_Write_GlobalEvaluationData::
  getGhostedVector() const
  {
    using std::logic_error;
    using Thyra::create_Vector;
    TEUCHOS_TEST_FOR_EXCEPTION(not isInitialized_, logic_error,
      "EpetraVector_Write_GlobalEvaluationData::getGhostedVector():  This "   \
      "object hasn't yet been initialized.")
    TEUCHOS_TEST_FOR_EXCEPTION(ghostedVector_.is_null(), logic_error,
      "EpetraVector_Write_GlobalEvaluationData::getGhostedVector():  The "    \
      "ghosted vector is just a null RCP.")
    return create_Vector(ghostedVector_, ghostedSpace_);
  } // end of getGhostedVector()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  print()
  //
  /////////////////////////////////////////////////////////////////////////////
  void
  EpetraVector_Write_GlobalEvaluationData::
  print(
    std::ostream& os) const
  {
    using std::string;
    const string tab("    ");
    os << "\n";
    os << tab << "EpetraVector_Write_GlobalEvaluationData\n"
       << tab << "  init    = " << isInitialized_ << "\n"
       << tab << "  owned   = " << ownedVector_   << "\n"
       << tab << "  ghosted = " << ghostedVector_ << "\n";
  } // end of print()

} // end of namespace panzer

// end of Panzer_EpetraVector_Write_GlobalEvaluationData.cpp
