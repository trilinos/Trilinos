// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Panzer_HierarchicParallelism.hpp"

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
    using PHX::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using panzer::PureBasis;
    using PHX::MDField;
    using PHX::print;
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
    TEUCHOS_TEST_FOR_EXCEPTION(fmNames.size() > 3,std::runtime_error,
                               "ERROR: Integrator_BasisTimesVector only supports up to three multipliers!");

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
    for (const auto& name : fmNames) {
      fieldMults_[i++] = MDField<const ScalarT, Cell, IP>(name, ir.dl_scalar);
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_BasisTimesVector<");
    n += std::to_string(fmNames.size()) + ">(";
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      n += "Cont";
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      n += "Eval";
    n += ", " + print<EvalT>() + "):  " + field_.fieldTag().name();
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
    using PHX::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using panzer::PureBasis;
    using PHX::MDField;
    using PHX::print;
    using std::invalid_argument;
    using std::logic_error;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(not ((bd_.getType() == "HCurl") or
      (bd_.getType() == "HDiv")), logic_error, "Error:  "                     \
      "Integrator_BasisTimesVector:  Basis of type \"" << bd_.getType()
      << "\" is not a vector basis.")

    TEUCHOS_TEST_FOR_EXCEPTION(multipliers.size() > 3,std::runtime_error,
                               "ERROR: Integrator_BasisTimesVector only supports up to three multipliers!");

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
    for (const auto& fm : multipliers)
    {
      fieldMults_[i++] = fm;
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_BasisTimesVector<");
    n += std::to_string(multipliers.size()) + ">(";
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      n += "Cont";
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      n += "Eval";
    n += ", " + print<EvalT>() + "):  " + field_.fieldTag().name();
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

    // Get the PHX::Views of the field multipliers.
    for (size_t i=0; i < fieldMults_.size(); ++i)
      kokkosFieldMults_[i] = fieldMults_[i].get_static_view();

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
    const typename Kokkos::TeamPolicy<FieldMultTag<NUM_FIELD_MULT>,PHX::exec_space>::member_type& team) const
  {
    using panzer::EvaluatorStyle;
    const int cell = team.league_rank();
    const int numBases = basis_.extent(1);

    // Initialize the evaluated field.
    if (evalStyle_ == EvaluatorStyle::EVALUATES) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
        field_(cell, basis) = 0.0;
      });
    }

    for (int qp(0); qp < numQP_; ++qp) {
      for (int dim(0); dim < numDim_; ++dim) {
        if constexpr (NUM_FIELD_MULT == 0) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
            field_(cell, basis) += basis_(cell, basis, qp, dim) * multiplier_ * vector_(cell, qp, dim);
          });
        }
        else if constexpr (NUM_FIELD_MULT == 1) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
            field_(cell, basis) += basis_(cell, basis, qp, dim) * multiplier_ * vector_(cell, qp, dim)
              * kokkosFieldMults_[0](cell, qp);
          });
        }
        else if constexpr (NUM_FIELD_MULT == 2) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
            field_(cell, basis) += basis_(cell, basis, qp, dim) * multiplier_ * vector_(cell, qp, dim) 
              * kokkosFieldMults_[0](cell, qp) * kokkosFieldMults_[1](cell, qp);
          });
        }
        else if constexpr (NUM_FIELD_MULT == 3) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
            field_(cell, basis) += basis_(cell, basis, qp, dim) * multiplier_ * vector_(cell, qp, dim)
              * kokkosFieldMults_[0](cell, qp) * kokkosFieldMults_[1](cell, qp) * kokkosFieldMults_[2](cell, qp);
          });
        }
        else {
          Kokkos::abort("Panzer_Integrator_BasisTimesVector: NUM_FIELD_MULT out of range!");
        }
      } // end loop over the dimensions of the vector field
    } // end loop over the quadrature points
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
    using Array=typename BasisValues2<double>::ConstArray_CellBasisIPDim;
    basis_ = useDescriptors_ ? bv.getVectorBasisValues(true) : Array(bv.weighted_basis_vector);

    // The following if-block is for the sake of optimization depending on the
    // number of field multipliers.  The parallel_fors will loop over the cells
    // in the Workset and execute operator()() above.
    if (fieldMults_.size() == 0) {
      auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag<0>,PHX::exec_space>(workset.num_cells);
      parallel_for("Panzer_Integrator_BasisTimesVector<0>", policy, *this);
    }
    else if (fieldMults_.size() == 1) {
      auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag<1>,PHX::exec_space>(workset.num_cells);
      parallel_for("Panzer_Integrator_BasisTimesVector<1>", policy, *this);
    }
    else  if (fieldMults_.size() == 2) {
      auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag<2>,PHX::exec_space>(workset.num_cells);
      parallel_for("Panzer_Integrator_BasisTimesVector<2>", policy, *this);
    }
    else  if (fieldMults_.size() == 3) {
      auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag<3>,PHX::exec_space>(workset.num_cells);
      parallel_for("Panzer_Integrator_BasisTimesVector<3>", policy, *this);
    }
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
