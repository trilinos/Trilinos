// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   PANZER_INTEGRATOR_BASISTIMESSCALAR_IMPL_HPP
#define   PANZER_INTEGRATOR_BASISTIMESSCALAR_IMPL_HPP

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
  Integrator_BasisTimesScalar<EvalT, Traits>::
  Integrator_BasisTimesScalar(
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
    multiplier_(multiplier),
    basisName_(basis.name())
  {
    using PHX::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using panzer::PureBasis;
    using PHX::MDField;
    using std::invalid_argument;
    using std::logic_error;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(resName == "", invalid_argument, "Error:  "    \
      "Integrator_BasisTimesScalar called with an empty residual name.")
    TEUCHOS_TEST_FOR_EXCEPTION(valName == "", invalid_argument, "Error:  "    \
      "Integrator_BasisTimesScalar called with an empty value name.")
    RCP<const PureBasis> tmpBasis = basis.getBasis();
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->isScalarBasis(), logic_error,
      "Error:  Integrator_BasisTimesScalar:  Basis of type \""
      << tmpBasis->name() << "\" is not a scalar basis.")

    // Create the field for the scalar quantity we're integrating.
    scalar_ = MDField<const ScalarT, Cell, IP>(valName, ir.dl_scalar);
    this->addDependentField(scalar_);

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
      View<Kokkos::View<const ScalarT**, typename PHX::DevLayout<ScalarT>::type, Kokkos::MemoryUnmanaged>*>("BasisTimesScalar::KokkosFieldMultipliers",
      fmNames.size());
    for (const auto& name : fmNames)
    {
      fieldMults_[i++] = MDField<const ScalarT, Cell, IP>(name, ir.dl_scalar);
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_BasisTimesScalar (");
    if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
      n += "CONTRIBUTES";
    else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
      n += "EVALUATES";
    n += "):  " + field_.fieldTag().name();
    this->setName(n);
  } // end of Main Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  ParameterList Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_BasisTimesScalar<EvalT, Traits>::
  Integrator_BasisTimesScalar(
    const Teuchos::ParameterList& p)
    :
    Integrator_BasisTimesScalar(
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
  Integrator_BasisTimesScalar<EvalT, Traits>::
  postRegistrationSetup(
    typename Traits::SetupData sd,
    PHX::FieldManager<Traits>& /* fm */)
  {
    using panzer::getBasisIndex;
    using std::size_t;

    auto kokkosFieldMults_h = Kokkos::create_mirror_view(kokkosFieldMults_);

    // Get the PHX::Views of the field multipliers.
    for (size_t i(0); i < fieldMults_.size(); ++i)
      kokkosFieldMults_h(i) = fieldMults_[i].get_static_view();

    Kokkos::deep_copy(kokkosFieldMults_, kokkosFieldMults_h);

    // Determine the index in the Workset bases for our particular basis name.
    basisIndex_ = getBasisIndex(basisName_, (*sd.worksets_)[0], this->wda);

    // Allocate temporary
    if (Sacado::IsADType<ScalarT>::value) {
      const auto fadSize = Kokkos::dimension_scalar(field_.get_view());
      tmp_ = PHX::View<ScalarT*>("IntegratorBasisTimesScalar::tmp_",field_.extent(0),fadSize);
      if (fieldMults_.size() > 1)
	tmp2_ = PHX::View<ScalarT*>("IntegratorBasisTimesScalar::tmp_",field_.extent(0),fadSize);
    } else {
      tmp_ = PHX::View<ScalarT*>("IntegratorBasisTimesScalar::tmp_",field_.extent(0));
      if (fieldMults_.size() > 1)
	tmp2_ = PHX::View<ScalarT*>("IntegratorBasisTimesScalar::tmp_",field_.extent(0));
    }
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
  Integrator_BasisTimesScalar<EvalT, Traits>::
  operator()(
    const FieldMultTag<NUM_FIELD_MULT>& /* tag */,
    const size_t&                       cell) const
  {
    using panzer::EvaluatorStyle;

    // Initialize the evaluated field.
    const int numQP(scalar_.extent(1)), numBases(basis_.extent(1));
    if (evalStyle_ == EvaluatorStyle::EVALUATES)
      for (int basis(0); basis < numBases; ++basis)
        field_(cell, basis) = 0.0;

    // The following if-block is for the sake of optimization depending on the
    // number of field multipliers.
    if (NUM_FIELD_MULT == 0)
    {
      // Loop over the quadrature points, scale the integrand by the
      // multiplier, and then perform the actual integration, looping over the
      // bases.
      for (int qp(0); qp < numQP; ++qp)
      {
        tmp_(cell) = multiplier_ * scalar_(cell, qp);
        for (int basis(0); basis < numBases; ++basis)
          field_(cell, basis) += basis_(cell, basis, qp) * tmp_(cell);
      } // end loop over the quadrature points
    }
    else if (NUM_FIELD_MULT == 1)
    {
      // Loop over the quadrature points, scale the integrand by the multiplier
      // and the single field multiplier, and then perform the actual
      // integration, looping over the bases.
      for (int qp(0); qp < numQP; ++qp)
      {
        tmp_(cell) = multiplier_ * scalar_(cell, qp) * kokkosFieldMults_(0)(cell, qp);
        for (int basis(0); basis < numBases; ++basis)
          field_(cell, basis) += basis_(cell, basis, qp) * tmp_(cell);
      } // end loop over the quadrature points
    }
    else
    {
      // Loop over the quadrature points and pre-multiply all the field
      // multipliers together, scale the integrand by the multiplier and
      // the combination of the field multipliers, and then perform the actual
      // integration, looping over the bases.
      const int numFieldMults(kokkosFieldMults_.extent(0));
      for (int qp(0); qp < numQP; ++qp)
      {
	tmp2_(cell) = 1.0;
        for (int fm(0); fm < numFieldMults; ++fm)
          tmp2_(cell) *= kokkosFieldMults_(fm)(cell, qp);
        tmp_(cell) = multiplier_ * scalar_(cell, qp) * tmp2_(cell);
        for (int basis(0); basis < numBases; ++basis)
          field_(cell, basis) += basis_(cell, basis, qp) * tmp_(cell);
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
  Integrator_BasisTimesScalar<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
    using Kokkos::parallel_for;
    using Kokkos::RangePolicy;

    // Grab the basis information.
    basis_ = this->wda(workset).bases[basisIndex_]->weighted_basis_scalar;

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
  template<typename EvalT, typename TRAITS>
  Teuchos::RCP<Teuchos::ParameterList>
  Integrator_BasisTimesScalar<EvalT, TRAITS>::
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
    RCP<const vector<string>> fms;
    p->set("Field Multipliers", fms);
    return p;
  } // end of getValidParameters()

} // end of namespace panzer

#endif // PANZER_INTEGRATOR_BASISTIMESSCALAR_IMPL_HPP
