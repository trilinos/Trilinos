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

#ifndef   __Panzer_Integrator_GradBasisDotVector_impl_hpp__
#define   __Panzer_Integrator_GradBasisDotVector_impl_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Kokkos
#include "Kokkos_ViewFactory.hpp"

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
  Integrator_GradBasisDotVector<EvalT, Traits>::
  Integrator_GradBasisDotVector(
    const panzer::EvaluatorStyle&        evalStyle,
    const std::string&                   resName,
    const std::string&                   fluxName,
    const panzer::BasisIRLayout&         basis,
    const panzer::IntegrationRule&       ir,
    const double&                        multiplier, /* = 1 */
    const std::vector<std::string>&      fmNames,    /* =
      std::vector<std::string>() */
    const Teuchos::RCP<PHX::DataLayout>& vecDL       /* = Teuchos::null */)
    :
    evalStyle_(evalStyle),
    multiplier_(multiplier),
    numDim_(static_cast<int>(ir.dl_vector->extent(2))),
    basisName_(basis.name())
  {
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using PHX::DataLayout;
    using PHX::MDField;
    using PHX::typeAsString;
    using std::invalid_argument;
    using std::logic_error;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(resName == "", invalid_argument, "Error:  "   \
      "Integrator_GradBasisDotVector called with an empty residual name.")
    TEUCHOS_TEST_FOR_EXCEPTION(fluxName == "", invalid_argument, "Error:  "   \
      "Integrator_GradBasisDotVector called with an empty flux name.")
    RCP<const PureBasis> tmpBasis = basis.getBasis();
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->supportsGrad(), logic_error,
      "Integrator_GradBasisDotVector:  Basis of type \"" << tmpBasis->name() <<
      "\" does not support the gradient operator.");
    RCP<DataLayout> tmpVecDL = ir.dl_vector;
    if (not vecDL.is_null())
    {
      tmpVecDL = vecDL;
      TEUCHOS_TEST_FOR_EXCEPTION(
        static_cast<int>(tmpVecDL->extent(2)) < numDim_, logic_error,
        "Integrator_GradBasisDotVector:  Dimension of space exceeds "         \
        "dimension of Vector Data Layout.");
    } // end if (not vecDL.is_null())

    // Create the field for the vector-valued function we're integrating.
    flux_ = MDField<const ScalarT, Cell, IP, Dim>(fluxName, tmpVecDL);
    this->addDependentField(flux_);

    // Create the field that we're either contributing to or evaluating
    // (storing).
    field_ = MDField<ScalarT, Cell, BASIS>(resName, basis.functional);
    if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
      this->addContributedField(field_);
    else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
      this->addEvaluatedField(field_);
 
    // Add the dependent field multipliers, if there are any.
    for (const auto& name : fmNames)
      fieldMults_.push_back(MDField<const ScalarT, Cell, IP>(name,
        ir.dl_scalar));
    for (const auto& mult : fieldMults_)
      this->addDependentField(mult);

    // Set the name of this object.
    string n("Integrator_GradBasisDotVector (");
    if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
      n += "Cont";
    else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
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
  Integrator_GradBasisDotVector<EvalT, Traits>::
  Integrator_GradBasisDotVector(
    const Teuchos::ParameterList& p)
    :
    Integrator_GradBasisDotVector(
      panzer::EvaluatorStyle::EVALUATES,
      p.get<std::string>("Residual Name"),
      p.get<std::string>("Flux Name"),
      (*p.get<Teuchos::RCP<panzer::BasisIRLayout>>("Basis")),
      (*p.get<Teuchos::RCP<panzer::IntegrationRule>>("IR")),
      p.get<double>("Multiplier"),
      p.isType<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers") ?
        (*p.get<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers")) : std::vector<std::string>(),
      p.isType<Teuchos::RCP<PHX::DataLayout>>("Vector Data Layout") ?
        p.get<Teuchos::RCP<PHX::DataLayout>>("Vector Data Layout") :
        Teuchos::null)
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
  Integrator_GradBasisDotVector<EvalT, Traits>::
  postRegistrationSetup(
    typename Traits::SetupData sd,
    PHX::FieldManager<Traits>& fm)
  {
    using Kokkos::createDynRankView;
    using panzer::getBasisIndex;

    // Determine the number of nodes and quadrature points.
    numNodes_ = static_cast<int>(field_.extent(1));
    numQP_    = static_cast<int>(flux_.extent(1));

    // Determine the index in the Workset bases for our particular basis name.
    basisIndex_ = getBasisIndex(basisName_, (*sd.worksets_)[0], this->wda);

    // Create a temporary View that we'll use in computing the integral.
    tmp = createDynRankView(field_.get_static_view(), "tmp",
      flux_.dimension(0), numQP_, numDim_);
  } // end of postRegistrationSetup()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evaluateFields()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_GradBasisDotVector<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
    using panzer::index_t;
    using PHX::MDField;

    // Initialize the evaluated field.
    if (evalStyle_ == EvaluatorStyle::EVALUATES)
      field_.deep_copy(ScalarT(0));

    // Scale the integrand by the multiplier, and any field multipliers, out in
    // front of the integral.
    for (index_t cell(0); cell < workset.num_cells; ++cell)
    {
      for (int qp(0); qp < numQP_; ++qp)
      {
        for (int dim(0); dim < numDim_; ++dim)
        {
          tmp(cell, qp, dim) = multiplier_ * flux_(cell, qp, dim);
          for (const auto& mult : fieldMults_)
            tmp(cell, qp, dim) *= mult(cell, qp);
        } // end loop over the dimensions
      } // end loop over the quadrature points
    } // end loop over the cells in the workset

    // Perform integration and vector dot product via looping over the cells in
    // the worksheet, the bases, the quadrature points, and the number of
    // dimensions of the flux vector.
    const BasisValues2<double> & bv = *this->wda(workset).bases[basisIndex_];
    ScalarT quantity(0);
    for (index_t cell(0); cell < workset.num_cells; ++cell)
      for (int basis(0); basis < numNodes_; ++basis)
        for (int qp(0); qp < numQP_; ++qp)
          for (int dim(0); dim < numDim_; ++dim)
            field_(cell, basis) += tmp(cell, qp, dim) *
              bv.weighted_grad_basis(cell, basis, qp, dim);
  } // end of evaluateFields()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getValidParameters()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename TRAITS>
  Teuchos::RCP<Teuchos::ParameterList>
  Integrator_GradBasisDotVector<EvalT, TRAITS>::
  getValidParameters() const
  {
    using panzer::BasisIRLayout;
    using panzer::IntegrationRule;
    using PHX::DataLayout;
    using std::string;
    using std::vector;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Create a ParameterList with all the valid keys we support.
    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<string>("Residual Name", "?");
    p->set<string>("Flux Name", "?");
    RCP<BasisIRLayout> basis;
    p->set("Basis", basis);
    RCP<IntegrationRule> ir;
    p->set("IR", ir);
    p->set<double>("Multiplier", 1.0);
    RCP<const vector<string>> fms;
    p->set("Field Multipliers", fms);
    RCP<DataLayout> vecDL;
    p->set("Vector Data Layout", vecDL);
    return p;
  } // end of getValidParameters()

} // end of namespace panzer

#endif // __Panzer_Integrator_GradBasisDotVector_impl_hpp__
