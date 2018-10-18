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

#ifndef   __Panzer_Integerator_BasisTimesVector_impl_hpp__
#define   __Panzer_Integerator_BasisTimesVector_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

namespace panzer
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  Main Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_BasisTimesVector<EvalT, Traits>::
  Integrator_BasisTimesVector(
    const panzer::EvaluatorStyle&   evalStyle,
    const std::string&              resName,
    const std::string&              valName,
    const panzer::BasisIRLayout&    basis,
    const panzer::IntegrationRule&  ir,
    const double&                   multiplier, /* = 1 */
    const std::vector<std::string>& fmNames     /* =
      std::vector<std::string>() */)
    :
    evalStyle_(evalStyle),
    useDescriptors_(false),
    multiplier_(multiplier),
    basisName_(basis.name())
  {
    using Kokkos::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using panzer::PureBasis;
    using PHX::MDField;
    using PHX::typeAsString;
    using std::invalid_argument;
    using std::logic_error;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(resName == "", invalid_argument, "Error:  "    \
      "Integrator_BasisTimesVector called with an empty residual name.")
    TEUCHOS_TEST_FOR_EXCEPTION(valName == "", invalid_argument, "Error:  "    \
      "Integrator_BasisTimesVector called with an empty value name.")
    RCP<const PureBasis> tmpBasis = basis.getBasis();
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->isVectorBasis(), logic_error,
      "Error:  Integrator_BasisTimesVector:  Basis of type \""
      << tmpBasis->name() << "\" is not a vector basis.");
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->requiresOrientations(),
      logic_error, "Integrator_BasisTimesVector:  Basis of type \""
      << tmpBasis->name() << "\" does not require orientations.  This seems " \
      "very strange, so I'm failing.");

    // Create the field for the vector-valued quantity we're integrating.
    vector_ = MDField<const ScalarT, Cell, IP, Dim>(valName, ir.dl_vector);
    this->addDependentField(vector_);

    // Create the field that we're either contributing to or evaluating
    // (storing).
    field_ = MDField<ScalarT, Cell, BASIS>(resName, basis.functional);
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      this->addContributedField(field_);
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      this->addEvaluatedField(field_);

    // Add the dependent field multipliers, if there are any.
    int i(0);
    fieldMults_.resize(fmNames.size());
    kokkosFieldMults_ =
      View<View<const ScalarT**,typename PHX::DevLayout<ScalarT>::type,PHX::Device>*>("BasisTimesVector::KokkosFieldMultipliers",
      fmNames.size());
    for (const auto& name : fmNames)
    {
      fieldMults_[i++] = MDField<const ScalarT, Cell, IP>(name, ir.dl_scalar);
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_BasisTimesVector (");
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      n += "Cont";
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      n += "Eval";
    n += ", " + typeAsString<EvalT>() + "):  " + field_.fieldTag().name();
    this->setName(n);
  } // end of Main Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Descriptor Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_BasisTimesVector<EvalT, Traits>::
  Integrator_BasisTimesVector(
    const panzer::EvaluatorStyle&                         evalStyle,
    const PHX::FieldTag&                                  resTag,
    const PHX::FieldTag&                                  valTag,
    const BasisDescriptor&                                bd,
    const IntegrationDescriptor&                          id,
    const double&                                         multiplier, /* = 1 */
    const std::vector<PHX::FieldTag>&                     multipliers /* =
      std::vector<PHX::FieldTag>() */)
    :
    evalStyle_(evalStyle),
    useDescriptors_(true),
    bd_(bd),
    id_(id),
    multiplier_(multiplier)
  {
    using Kokkos::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using panzer::PureBasis;
    using PHX::MDField;
    using PHX::typeAsString;
    using std::invalid_argument;
    using std::logic_error;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(not ((bd_.getType() == "HCurl") or
      (bd_.getType() == "HDiv")), logic_error, "Error:  "                     \
      "Integrator_BasisTimesVector:  Basis of type \"" << bd_.getType()
      << "\" is not a vector basis.")

    // Create the field for the vector-valued quantity we're integrating.
    vector_ = valTag;
    this->addDependentField(vector_);

    // Create the field that we're either contributing to or evaluating
    // (storing).
    field_ = resTag;
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      this->addContributedField(field_);
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      this->addEvaluatedField(field_);

    // Add the dependent field multipliers, if there are any.
    int i(0);
    fieldMults_.resize(multipliers.size());
    kokkosFieldMults_ =
      View<View<const ScalarT**,typename PHX::DevLayout<ScalarT>::type,PHX::Device>*>("BasisTimesVector::KokkosFieldMultipliers",
      multipliers.size());
    for (const auto& fm : multipliers)
    {
      fieldMults_[i++] = fm;
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_BasisTimesVector (");
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      n += "Cont";
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      n += "Eval";
    n += ", " + typeAsString<EvalT>() + "):  " + field_.fieldTag().name();
    this->setName(n);
  } // end of Main Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  ParameterList Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_BasisTimesVector<EvalT, Traits>::
  Integrator_BasisTimesVector(
    const Teuchos::ParameterList& p)
    :
    Integrator_BasisTimesVector(
      panzer::EvaluatorStyle::EVALUATES,
      p.get<std::string>("Residual Name"),
      p.get<std::string>("Value Name"),
      (*p.get<Teuchos::RCP<panzer::BasisIRLayout>>("Basis")),
      (*p.get<Teuchos::RCP<panzer::IntegrationRule>>("IR")),
      p.get<double>("Multiplier"),
      p.isType<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers") ?
        (*p.get<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers")) : std::vector<std::string>())
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // Ensure that the input ParameterList didn't contain any bogus entries.
    RCP<ParameterList> validParams = this->getValidParameters();
    p.validateParameters(*validParams);
  } // end of ParameterList Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  postRegistrationSetup()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_BasisTimesVector<EvalT, Traits>::
  postRegistrationSetup(
    typename Traits::SetupData sd,
    PHX::FieldManager<Traits>& /* fm */)
  {
    using panzer::getBasisIndex;
    using std::size_t;

    // Get the Kokkos::Views of the field multipliers.
    for (size_t i(0); i < fieldMults_.size(); ++i)
      kokkosFieldMults_(i) = fieldMults_[i].get_static_view();

    // Determine the number of quadrature points and the dimensionality of the
    // vector that we're integrating.
    numQP_  = vector_.extent(1);
    numDim_ = vector_.extent(2);

    // Determine the index in the Workset bases for our particular basis name.
    if (not useDescriptors_)
      basisIndex_ = getBasisIndex(basisName_, (*sd.worksets_)[0], this->wda);
  } // end of postRegistrationSetup()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  operator()()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  template<int NUM_FIELD_MULT>
  KOKKOS_INLINE_FUNCTION
  void
  Integrator_BasisTimesVector<EvalT, Traits>::
  operator()(
    const FieldMultTag<NUM_FIELD_MULT>& /* tag */,
    const size_t&                       cell) const
  {
    using panzer::EvaluatorStyle;

    // Initialize the evaluated field.
    const int numBases(basis_.extent(1));
    if (evalStyle_ == EvaluatorStyle::EVALUATES)
      for (int basis(0); basis < numBases; ++basis)
        field_(cell, basis) = 0.0;

    // The following if-block is for the sake of optimization depending on the
    // number of field multipliers.
    ScalarT tmp;
    if (NUM_FIELD_MULT == 0)
    {
      // Loop over the quadrature points and dimensions of our vector fields,
      // scale the integrand by the multiplier, and then perform the
      // integration, looping over the bases.
      for (int qp(0); qp < numQP_; ++qp)
      {
        for (int dim(0); dim < numDim_; ++dim)
        {
          tmp = multiplier_ * vector_(cell, qp, dim);
          for (int basis(0); basis < numBases; ++basis)
            field_(cell, basis) += basis_(cell, basis, qp, dim) * tmp;
        } // end loop over the dimensions of the vector field
      } // end loop over the quadrature points
    }
    else if (NUM_FIELD_MULT == 1)
    {
      // Loop over the quadrature points and dimensions of our vector fields,
      // scale the integrand by the multiplier and the single field multiplier,
      // and then perform the actual integration, looping over the bases.
      for (int qp(0); qp < numQP_; ++qp)
      {
        for (int dim(0); dim < numDim_; ++dim)
        {
          tmp = multiplier_ * vector_(cell, qp, dim) *
            kokkosFieldMults_(0)(cell, qp);
          for (int basis(0); basis < numBases; ++basis)
            field_(cell, basis) += basis_(cell, basis, qp, dim) * tmp;
        } // end loop over the dimensions of the vector field
      } // end loop over the quadrature points
    }
    else
    {
      // Loop over the quadrature points and pre-multiply all the field
      // multipliers together.  Then loop over the dimensions of our vector
      // fields, scale the integrand by the multiplier and the combination of
      // the field multipliers, and then perform the actual integration,
      // looping over the bases.
      const int numFieldMults(kokkosFieldMults_.extent(0));
      for (int qp(0); qp < numQP_; ++qp)
      {
        ScalarT fieldMultsTotal(1);
        for (int fm(0); fm < numFieldMults; ++fm)
          fieldMultsTotal *= kokkosFieldMults_(fm)(cell, qp);
        for (int dim(0); dim < numDim_; ++dim)
        {
          tmp = multiplier_ * vector_(cell, qp, dim) * fieldMultsTotal;
          for (int basis(0); basis < numBases; ++basis)
            field_(cell, basis) += basis_(cell, basis, qp, dim) * tmp;
        } // end loop over the dimensions of the vector field
      } // end loop over the quadrature points
    } // end if (NUM_FIELD_MULT == something)
  } // end of operator()()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evaluateFields()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_BasisTimesVector<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
    using Kokkos::parallel_for;
    using Kokkos::RangePolicy;

    // Grab the basis information.
    const panzer::BasisValues2<double>& bv = useDescriptors_ ?
      this->wda(workset).getBasisValues(bd_,id_) :
      *this->wda(workset).bases[basisIndex_];
    basis_ = bv.weighted_basis_vector;

    // The following if-block is for the sake of optimization depending on the
    // number of field multipliers.  The parallel_fors will loop over the cells
    // in the Workset and execute operator()() above.
    if (fieldMults_.size() == 0)
      parallel_for(RangePolicy<FieldMultTag<0>>(0, workset.num_cells), *this);
    else if (fieldMults_.size() == 1)
      parallel_for(RangePolicy<FieldMultTag<1>>(0, workset.num_cells), *this);
    else
      parallel_for(RangePolicy<FieldMultTag<-1>>(0, workset.num_cells), *this);
  } // end of evaluateFields()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getValidParameters()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Teuchos::RCP<Teuchos::ParameterList>
  Integrator_BasisTimesVector<EvalT, Traits>::
  getValidParameters() const
  {
    using panzer::BasisIRLayout;
    using panzer::IntegrationRule;
    using std::string;
    using std::vector;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Create a ParameterList with all the valid keys we support.
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<string>("Residual Name", "?");
    p->set<string>("Value Name", "?");
    RCP<BasisIRLayout> basis;
    p->set("Basis", basis);
    RCP<IntegrationRule> ir;
    p->set("IR", ir);
    p->set<double>("Multiplier", 1.0);
    RCP<const vector<string> > fms;
    p->set("Field Multipliers", fms);
    return p;
  } // end of getValidParameters()

} // end of namespace panzer

#endif // __Panzer_Integerator_BasisTimesVector_impl_hpp__
