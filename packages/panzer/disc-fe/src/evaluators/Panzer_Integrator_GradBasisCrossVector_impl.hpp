// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   PANZER_EVALUATOR_GRADBASISCROSSVECTOR_IMPL_HPP
#define   PANZER_EVALUATOR_GRADBASISCROSSVECTOR_IMPL_HPP

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Intrepid2
#include "Intrepid2_FunctionSpaceTools.hpp"

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
  //  Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_GradBasisCrossVector<EvalT, Traits>::
  Integrator_GradBasisCrossVector(
    const panzer::EvaluatorStyle&        evalStyle,
    const std::vector<std::string>&      resNames,
    const std::string&                   vecName,
    const panzer::BasisIRLayout&         basis,
    const panzer::IntegrationRule&       ir,
    const double&                        multiplier, /* = 1 */
    const std::vector<std::string>&      fmNames,    /* =
      std::vector<std::string>() */
    const Teuchos::RCP<PHX::DataLayout>& vecDL       /* = Teuchos::null */)
    :
    evalStyle_(evalStyle),
    multiplier_(multiplier),
    numDims_(resNames.size()),
    numGradDims_(ir.dl_vector->extent(2)),
    basisName_(basis.name())
  {
    using PHX::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using PHX::DataLayout;
    using PHX::Device;
    using PHX::DevLayout;
    using PHX::MDField;
    using std::invalid_argument;
    using std::logic_error;
    using std::size_t;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(numDims_ != 3, invalid_argument, "Error:  "    \
      "Integrator_GradBasisCrossVector called with the number of residual "   \
      "names not equal to three.")
    for (const auto& name : resNames)
      TEUCHOS_TEST_FOR_EXCEPTION(name == "", invalid_argument, "Error:  "     \
        "Integrator_GradBasisCrossVector called with an empty residual name.")
    TEUCHOS_TEST_FOR_EXCEPTION(vecName == "", invalid_argument, "Error:  "    \
      "Integrator_GradBasisCrossVector called with an empty vector name.")
    RCP<const PureBasis> tmpBasis = basis.getBasis();
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->supportsGrad(), logic_error,
      "Error:  Integrator_GradBasisCrossVector:  Basis of type \""
      << tmpBasis->name() << "\" does not support the gradient operator.")
    RCP<DataLayout> tmpVecDL = ir.dl_vector;
    if (not vecDL.is_null())
    {
      tmpVecDL = vecDL;
      TEUCHOS_TEST_FOR_EXCEPTION(
        tmpVecDL->extent(2) < ir.dl_vector->extent(2), logic_error,
        "Error:  Integrator_GradBasisCrossVector:  Dimension of space "       \
        "exceeds dimension of Vector Data Layout.")
      TEUCHOS_TEST_FOR_EXCEPTION(numDims_ !=
        static_cast<int>(vecDL->extent(2)), logic_error, "Error:  "           \
        "Integrator_GradBasisCrossVector:  The vector must be the same "      \
        "length as the number of residuals.")
    } // end if (not vecDL.is_null())
    TEUCHOS_TEST_FOR_EXCEPTION(numGradDims_ > numDims_, logic_error,
      "Error:  Integrator_GradBasisCrossVector:  The vector must have at "    \
      "least as many components as there are dimensions in the mesh.")

    // Create the field for the vector-valued function we're integrating.
    vector_ = MDField<const ScalarT, Cell, IP, Dim>(vecName, tmpVecDL);
    this->addDependentField(vector_);

    // Create the fields that we're either contributing to or evaluating
    // (storing).
    fields_host_.resize(resNames.size());
    fields_ = OuterView("Integrator_GradBasisCrossVector::fields_", resNames.size());
    {
      int i=0;
      for (const auto& name : resNames)
        fields_host_[i++] = MDField<ScalarT, Cell, BASIS>(name, basis.functional);
    } // end loop over resNames

    for (std::size_t i=0; i< fields_.extent(0); ++i) {
      const auto& field = fields_host_[i];
      if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
        this->addContributedField(field);
      else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
        this->addEvaluatedField(field);
    }

    // Add the dependent field multipliers, if there are any.
    int i = 0;
    fieldMults_.resize(fmNames.size());
    kokkosFieldMults_ = PHX::View<PHX::UnmanagedView<const ScalarT**>*>("GradBasisCrossVector::KokkosFieldMultipliers", fmNames.size());
    for (const auto& name : fmNames)
    {
      fieldMults_[i++] = MDField<const ScalarT, Cell, IP>(name, ir.dl_scalar);
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_GradBasisCrossVector (");
    if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
      n += "CONTRIBUTES";
    else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
      n += "EVALUATES";
    n += "):  {";
    for (size_t j=0; j < fields_host_.size() - 1; ++j)
      n += resNames[j] + ", ";
    n += resNames[resNames.size()-1] + "}";
    this->setName(n);
  } // end of Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  ParameterList Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_GradBasisCrossVector<EvalT, Traits>::
  Integrator_GradBasisCrossVector(
    const Teuchos::ParameterList& p)
    :
    Integrator_GradBasisCrossVector(
      panzer::EvaluatorStyle::EVALUATES,
      p.get<const std::vector<std::string>>("Residual Names"),
      p.get<std::string>("Vector Name"),
      (*p.get<Teuchos::RCP<panzer::BasisIRLayout>>("Basis")),
      (*p.get<Teuchos::RCP<panzer::IntegrationRule>>("IR")),
      p.get<double>("Multiplier"),
      p.isType<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers") ?
        (*p.get<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers")) : std::vector<std::string>(),
      p.isType<Teuchos::RCP<PHX::DataLayout>>("Data Layout Vector") ?
        p.get<Teuchos::RCP<PHX::DataLayout>>("Data Layout Vector") :
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
  Integrator_GradBasisCrossVector<EvalT, Traits>::
  postRegistrationSetup(
    typename Traits::SetupData sd,
    PHX::FieldManager<Traits>& /* fm */)
  {
    using Kokkos::createDynRankView;
    using panzer::getBasisIndex;

    // Get the PHX::Views of the fields.
    auto fields_host_mirror_ = Kokkos::create_mirror_view(fields_);
    for (size_t i=0; i < fields_host_.size(); ++i) {
      fields_host_mirror_(i) = fields_host_[i].get_static_view();
    }
    Kokkos::deep_copy(fields_,fields_host_mirror_);

    // Get the PHX::Views of the field multipliers.
    auto field_mults_host_mirror_ = Kokkos::create_mirror_view(kokkosFieldMults_);
    for (size_t i=0; i < fieldMults_.size(); ++i)
      field_mults_host_mirror_(i) = fieldMults_[i].get_static_view();
    Kokkos::deep_copy(kokkosFieldMults_,field_mults_host_mirror_);

    // Determine the index in the Workset bases for our particular basis name.
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
  Integrator_GradBasisCrossVector<EvalT, Traits>::
  operator()(
    const FieldMultTag<NUM_FIELD_MULT>& /* tag */,
    const size_t&                       cell) const
  {
    using panzer::EvaluatorStyle;

    // Initialize the evaluated fields.
    const int numBases(fields_[0].extent(1)), numQP(vector_.extent(1));
    if (evalStyle_ == EvaluatorStyle::EVALUATES)
      for (int dim(0); dim < numDims_; ++dim)
        for (int basis(0); basis < numBases; ++basis)
          fields_[dim](cell, basis) = 0.0;

    // The following if-block is for the sake of optimization depending on the
    // number of field multipliers.
    ScalarT tmp[3];
    const int X(0), Y(1), Z(2);
    if (NUM_FIELD_MULT == 0)
    {
      if (numGradDims_ == 1)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[Y] = multiplier_ * vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * vector_(cell, qp, Z);
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[Y](cell, basis) +=  tmp[Z] * basis_(cell, basis, qp, X);
            fields_[Z](cell, basis) += -tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      }
      else if (numGradDims_ == 2)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[X] = multiplier_ * vector_(cell, qp, X);
          tmp[Y] = multiplier_ * vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * vector_(cell, qp, Z);
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[X](cell, basis) += -tmp[Z] * basis_(cell, basis, qp, Y);
            fields_[Y](cell, basis) +=  tmp[Z] * basis_(cell, basis, qp, X);
            fields_[Z](cell, basis) +=  tmp[X] * basis_(cell, basis, qp, Y) -
                                        tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      }
      else if (numGradDims_ == 3)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[X] = multiplier_ * vector_(cell, qp, X);
          tmp[Y] = multiplier_ * vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * vector_(cell, qp, Z);
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[X](cell, basis) += tmp[Y] * basis_(cell, basis, qp, Z) -
                                       tmp[Z] * basis_(cell, basis, qp, Y);
            fields_[Y](cell, basis) += tmp[Z] * basis_(cell, basis, qp, X) -
                                       tmp[X] * basis_(cell, basis, qp, Z);
            fields_[Z](cell, basis) += tmp[X] * basis_(cell, basis, qp, Y) -
                                       tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      } // end if (numGradDims_ == something)
    }
    else if (NUM_FIELD_MULT == 1)
    {
      if (numGradDims_ == 1)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[Y] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, Z);
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[Y](cell, basis) +=  tmp[Z] * basis_(cell, basis, qp, X);
            fields_[Z](cell, basis) += -tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      }
      else if (numGradDims_ == 2)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[X] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, X);
          tmp[Y] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, Z);
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[X](cell, basis) += -tmp[Z] * basis_(cell, basis, qp, Y);
            fields_[Y](cell, basis) +=  tmp[Z] * basis_(cell, basis, qp, X);
            fields_[Z](cell, basis) +=  tmp[X] * basis_(cell, basis, qp, Y) -
                                        tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      }
      else if (numGradDims_ == 3)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[X] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, X);
          tmp[Y] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * kokkosFieldMults_(0)(cell, qp) *
            vector_(cell, qp, Z);
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[X](cell, basis) += tmp[Y] * basis_(cell, basis, qp, Z) -
                                       tmp[Z] * basis_(cell, basis, qp, Y);
            fields_[Y](cell, basis) += tmp[Z] * basis_(cell, basis, qp, X) -
                                       tmp[X] * basis_(cell, basis, qp, Z);
            fields_[Z](cell, basis) += tmp[X] * basis_(cell, basis, qp, Y) -
                                       tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      } // end if (numGradDims_ == something)
    }
    else
    {
      const int numFieldMults(kokkosFieldMults_.extent(0));
      if (numGradDims_ == 1)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[Y] = multiplier_ * vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * vector_(cell, qp, Z);
          for (int fm(0); fm < numFieldMults; ++fm)
          {
            tmp[Y] *= kokkosFieldMults_(fm)(cell, qp);
            tmp[Z] *= kokkosFieldMults_(fm)(cell, qp);
          } // end loop over the field multipliers
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[Y](cell, basis) +=  tmp[Z] * basis_(cell, basis, qp, X);
            fields_[Z](cell, basis) += -tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      }
      else if (numGradDims_ == 2)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[X] = multiplier_ * vector_(cell, qp, X);
          tmp[Y] = multiplier_ * vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * vector_(cell, qp, Z);
          for (int fm(0); fm < numFieldMults; ++fm)
          {
            tmp[X] *= kokkosFieldMults_(fm)(cell, qp);
            tmp[Y] *= kokkosFieldMults_(fm)(cell, qp);
            tmp[Z] *= kokkosFieldMults_(fm)(cell, qp);
          } // end loop over the field multipliers
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[X](cell, basis) += -tmp[Z] * basis_(cell, basis, qp, Y);
            fields_[Y](cell, basis) +=  tmp[Z] * basis_(cell, basis, qp, X);
            fields_[Z](cell, basis) +=  tmp[X] * basis_(cell, basis, qp, Y) -
                                        tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      }
      else if (numGradDims_ == 3)
      {
        for (int qp(0); qp < numQP; ++qp)
        {
          tmp[X] = multiplier_ * vector_(cell, qp, X);
          tmp[Y] = multiplier_ * vector_(cell, qp, Y);
          tmp[Z] = multiplier_ * vector_(cell, qp, Z);
          for (int fm(0); fm < numFieldMults; ++fm)
          {
            tmp[X] *= kokkosFieldMults_(fm)(cell, qp);
            tmp[Y] *= kokkosFieldMults_(fm)(cell, qp);
            tmp[Z] *= kokkosFieldMults_(fm)(cell, qp);
          } // end loop over the field multipliers
          for (int basis(0); basis < numBases; ++basis)
          {
            fields_[X](cell, basis) += tmp[Y] * basis_(cell, basis, qp, Z) -
                                       tmp[Z] * basis_(cell, basis, qp, Y);
            fields_[Y](cell, basis) += tmp[Z] * basis_(cell, basis, qp, X) -
                                       tmp[X] * basis_(cell, basis, qp, Z);
            fields_[Z](cell, basis) += tmp[X] * basis_(cell, basis, qp, Y) -
                                       tmp[Y] * basis_(cell, basis, qp, X);
          } // end loop over the bases
        } // end loop over the quadrature points
      } // end if (numGradDims_ == something)
    } // end if (NUM_FIELD_MULT == something)
  } // end of operator()()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evaluateFields()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_GradBasisCrossVector<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
    using Kokkos::parallel_for;
    using Kokkos::RangePolicy;

    // Grab the basis information.
    basis_ = this->wda(workset).bases[basisIndex_]->weighted_grad_basis;

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
  Integrator_GradBasisCrossVector<EvalT, TRAITS>::
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

    RCP<const vector<string>> resNames;
    p->set("Residual Names", resNames);
    p->set<string>("Vector Name", "?");
    RCP<BasisIRLayout> basis;
    p->set("Basis", basis);
    RCP<IntegrationRule> ir;
    p->set("IR", ir);
    p->set<double>("Multiplier", 1.0);
    RCP<const vector<string>> fms;
    p->set("Field Multipliers", fms);
    RCP<DataLayout> vecDL;
    p->set("Data Layout Vector", vecDL);

    return p;
  } // end of getValidParameters()

} // end of namespace panzer

#endif // PANZER_EVALUATOR_GRADBASISCROSSVECTOR_IMPL_HPP
