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
    using Kokkos::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using PHX::DataLayout;
    using PHX::MDField;
    using std::invalid_argument;
    using std::logic_error;
    using std::size_t;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(resNames.size() != 3, invalid_argument,
      "Error:  Integrator_GradBasisCrossVector called with the number of "    \
      "residual names not equal to three.")
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
        TEUCHOS_TEST_FOR_EXCEPTION(numDims_ != static_cast<int>(vecDL->extent(2)), logic_error,
        "Error:  Integrator_GradBasisCrossVector:  The vector must be the "   \
        "same length as the number of residuals.")
    } // end if (not vecDL.is_null())
    TEUCHOS_TEST_FOR_EXCEPTION(numGradDims_ > numDims_, logic_error,
      "Error:  Integrator_GradBasisCrossVector:  The vector must have at "    \
      "least as many components as there are dimensions in the mesh.")

    // Create the field for the vector-valued function we're integrating.
    vector_ = MDField<const ScalarT, Cell, IP, Dim>(vecName, tmpVecDL);
    this->addDependentField(vector_);

    // Create the fields that we're either contributing to or evaluating
    // (storing).
    fields_.clear();
    for (const auto& name : resNames)
    {
      MDField<ScalarT, Cell, BASIS> res(name, basis.functional);
      fields_.push_back(res);
    } // end loop over resNames
    for (const auto& field : fields_)
      if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
        this->addContributedField(field);
      else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
        this->addEvaluatedField(field);

    // Add the dependent field multipliers, if there are any.
    int i = 0;
    fieldMults_.resize(fmNames.size());
    kokkosFieldMults_ = View<View<const ScalarT**,typename PHX::DevLayout<ScalarT>::type,PHX::Device>*>(
      "GradBasisCrossVector::KokkosFieldMultipliers", fmNames.size());
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
    for (size_t j=0; j < fields_.size() - 1; ++j)
      n += fields_[j].fieldTag().name() + ", ";
    n += fields_.back().fieldTag().name() + "}";
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
    using PHX::Device;

    // Get the Kokkos::Views of the field multipliers.
    for (size_t i(0); i < fieldMults_.size(); ++i)
      kokkosFieldMults_(i) = fieldMults_[i].get_static_view();
    Device::fence();

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
