// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_Integrator_DivBasisTimesScalar_impl_hpp__
#define   __Panzer_Integrator_DivBasisTimesScalar_impl_hpp__

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
#include "Panzer_HierarchicParallelism.hpp"

namespace panzer
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  Main Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_DivBasisTimesScalar<EvalT, Traits>::
  Integrator_DivBasisTimesScalar(
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
      "Integrator_DivBasisTimesScalar called with an empty residual name.")
    TEUCHOS_TEST_FOR_EXCEPTION(valName == "", invalid_argument, "Error:  "    \
      "Integrator_DivBasisTimesScalar called with an empty value name.")
    RCP<const PureBasis> tmpBasis = basis.getBasis();
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->supportsDiv(), logic_error,
      "Error:  Integrator_DivBasisTimesScalar:  Basis of type \""
      << tmpBasis->name() << "\" does not support DIV.")
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->requiresOrientations(), logic_error,
      "Error:  Integration_DivBasisTimesScalar:  Basis of type \""
      << tmpBasis->name() << "\" should require orientations.")

    // Create the field for the scalar quantity we're integrating.
    scalar_ = MDField<const ScalarT, Cell, IP>(valName, ir.dl_scalar);
    this->addDependentField(scalar_);

    // Create the field that we're either contributing to or evaluating
    // (storing).
    field_ = MDField<ScalarT, Cell, BASIS>(resName, basis.functional);
    if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
      this->addContributedField(field_);
    else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
      this->addEvaluatedField(field_);

    // Add the dependent field multipliers, if there are any.
    int i(0);
    fieldMults_.resize(fmNames.size());
    kokkosFieldMults_ =
      View<PHX::UnmanagedView<const ScalarT**>*>("DivBasisTimesScalar::KokkosFieldMultipliers",
      fmNames.size());
    for (const auto& name : fmNames)
    {
      fieldMults_[i++] = MDField<const ScalarT, Cell, IP>(name, ir.dl_scalar);
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_DivBasisTimesScalar (");
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
  Integrator_DivBasisTimesScalar<EvalT, Traits>::
  Integrator_DivBasisTimesScalar(
    const Teuchos::ParameterList& p)
    :
    Integrator_DivBasisTimesScalar(
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
  Integrator_DivBasisTimesScalar<EvalT, Traits>::
  postRegistrationSetup(
    typename Traits::SetupData sd,
    PHX::FieldManager<Traits>& /* fm */)
  {
    using Kokkos::createDynRankView;
    using panzer::getBasisIndex;

    auto kokkosFieldMults_h = Kokkos::create_mirror_view(kokkosFieldMults_);

    // Get the PHX::Views of the field multipliers.
    for (size_t i(0); i < fieldMults_.size(); ++i)
      kokkosFieldMults_h(i) = fieldMults_[i].get_static_view();

    Kokkos::deep_copy(kokkosFieldMults_, kokkosFieldMults_h);

    // Allocate temporary if not using shared memory
    bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();
    if (!use_shared_memory) {
      if (Sacado::IsADType<ScalarT>::value) {
	const auto fadSize = Kokkos::dimension_scalar(field_.get_view());
	tmp_ = PHX::View<ScalarT*>("panzer::Integrator::DivBasisTimesScalar::tmp_",field_.extent(0),fadSize);
      } else {
	tmp_ = PHX::View<ScalarT*>("panzer::Integrator::DivBasisTimesScalar::tmp_",field_.extent(0));
      }
    }

    // Determine the index in the Workset bases for our particular basis name.
    basisIndex_ = getBasisIndex(basisName_, (*sd.worksets_)[0], this->wda);
  } // end of postRegistrationSetup()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  No shared memory operator()()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  template<int NUM_FIELD_MULT>
  KOKKOS_INLINE_FUNCTION
  void
  Integrator_DivBasisTimesScalar<EvalT, Traits>::
  operator()(
    const FieldMultTag<NUM_FIELD_MULT>& /* tag */,
    const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    using panzer::EvaluatorStyle;
    const int cell = team.league_rank();

    // Initialize the evaluated field.
    const int numQP(scalar_.extent(1)), numBases(basis_.extent(1));
    if (evalStyle_ == EvaluatorStyle::EVALUATES) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
	field_(cell, basis) = 0.0;
      });
    }

    // The following if-block is for the sake of optimization depending on the
    // number of field multipliers.
    if (NUM_FIELD_MULT == 0)
    {
      // Loop over the quadrature points, scale the integrand by the
      // multiplier, and then perform the actual integration, looping over the
      // bases.
      for (int qp(0); qp < numQP; ++qp)
      {
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
          field_(cell, basis) += basis_(cell, basis, qp) * multiplier_ * scalar_(cell, qp);
	});
      } // end loop over the quadrature points
    }
    else if (NUM_FIELD_MULT == 1)
    {
      // Loop over the quadrature points, scale the integrand by the multiplier
      // and the single field multiplier, and then perform the actual
      // integration, looping over the bases.
      for (int qp(0); qp < numQP; ++qp)
      {
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
          field_(cell, basis) += basis_(cell, basis, qp) * multiplier_ * scalar_(cell, qp) * kokkosFieldMults_(0)(cell, qp);
	});
      } // end loop over the quadrature points
    }
    else
    {
      // Loop over the quadrature points and pre-multiply all the field
      // multipliers together, scale the integrand by the multiplier and the
      // combination of field multipliers, and then perform the actual
      // integration, looping over the bases.
      const int numFieldMults(kokkosFieldMults_.extent(0));
      for (int qp(0); qp < numQP; ++qp)
      {
	team.team_barrier();
	tmp_(cell) = 1.0;
        for (int fm(0); fm < numFieldMults; ++fm)
          tmp_(cell) *= kokkosFieldMults_(fm)(cell, qp);
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int& basis) {
	    field_(cell, basis) += basis_(cell, basis, qp) * multiplier_ * scalar_(cell, qp) * tmp_(cell);
	  });
      } // end loop over the quadrature points
    } // end if (NUM_FIELD_MULT == something)
  } // end of operator()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Shared memory operator()()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  template<int NUM_FIELD_MULT>
  KOKKOS_INLINE_FUNCTION
  void
  Integrator_DivBasisTimesScalar<EvalT, Traits>::
  operator()(
    const SharedFieldMultTag<NUM_FIELD_MULT>& /* tag */,
    const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    using panzer::EvaluatorStyle;
    const int cell = team.league_rank();
    const int numQP = scalar_.extent(1);
    const int numBases = basis_.extent(1);
    const int fadSize = Kokkos::dimension_scalar(field_.get_view());

    scratch_view tmp;
    scratch_view tmp_field;
    if (Sacado::IsADType<ScalarT>::value) {
      tmp = scratch_view(team.team_shmem(),1,fadSize);
      tmp_field = scratch_view(team.team_shmem(),numBases,fadSize);
    }
    else {
      tmp = scratch_view(team.team_shmem(),1);
      tmp_field = scratch_view(team.team_shmem(),numBases);
    }

    // Initialize the evaluated field.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int& basis) {
      tmp_field(basis) = 0.0;
    });

    // The following if-block is for the sake of optimization depending on the
    // number of field multipliers.
    if (NUM_FIELD_MULT == 0)
    {
      // Loop over the quadrature points, scale the integrand by the
      // multiplier, and then perform the actual integration, looping over the
      // bases.
      for (int qp(0); qp < numQP; ++qp)
      {
	team.team_barrier();
	tmp(0) = multiplier_ * scalar_(cell, qp);
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int& basis) {
	  tmp_field(basis) += basis_(cell, basis, qp) * tmp(0);
	});
      } // end loop over the quadrature points
    }
    else if (NUM_FIELD_MULT == 1)
    {
      // Loop over the quadrature points, scale the integrand by the multiplier
      // and the single field multiplier, and then perform the actual
      // integration, looping over the bases.
      for (int qp(0); qp < numQP; ++qp)
      {
	team.team_barrier();
        tmp(0) = multiplier_ * scalar_(cell, qp) * kokkosFieldMults_(0)(cell, qp);
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int& basis) {
	  tmp_field(basis) += basis_(cell, basis, qp) * tmp(0);
	});
      } // end loop over the quadrature points
    }
    else
    {
      // Loop over the quadrature points and pre-multiply all the field
      // multipliers together, scale the integrand by the multiplier and the
      // combination of field multipliers, and then perform the actual
      // integration, looping over the bases.
      const int numFieldMults = kokkosFieldMults_.extent(0);
      for (int qp(0); qp < numQP; ++qp)
      {
	team.team_barrier();
        ScalarT fieldMultsTotal(1); // need shared mem here
        for (int fm(0); fm < numFieldMults; ++fm)
          fieldMultsTotal *= kokkosFieldMults_(fm)(cell, qp);
        tmp(0) = multiplier_ * scalar_(cell, qp) * fieldMultsTotal;
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int& basis) {
	  tmp_field(basis) += basis_(cell, basis, qp) * tmp(0);
	});
      } // end loop over the quadrature points
    } // end if (NUM_FIELD_MULT == something)

    // Put final values into target field
    if (evalStyle_ == EvaluatorStyle::EVALUATES) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int& basis) {
	field_(cell,basis) = tmp_field(basis);
      });
    }
    else { // Contributed
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int& basis) {
	field_(cell,basis) += tmp_field(basis);
      });
    }

  } // end of operator()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evaluateFields()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_DivBasisTimesScalar<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
    using Kokkos::parallel_for;
    using Kokkos::TeamPolicy;

    // Grab the basis information.
    basis_ = this->wda(workset).bases[basisIndex_]->weighted_div_basis;

    bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();

    if (use_shared_memory) {
      int bytes;
      if (Sacado::IsADType<ScalarT>::value) {
	const int fadSize = Kokkos::dimension_scalar(field_.get_view());
	bytes = scratch_view::shmem_size(1,fadSize) + scratch_view::shmem_size(basis_.extent(1),fadSize);
      }
      else
	bytes = scratch_view::shmem_size(1) + scratch_view::shmem_size(basis_.extent(1));

      // The following if-block is for the sake of optimization depending on the
      // number of field multipliers.  The parallel_fors will loop over the cells
      // in the Workset and execute operator()() above.
      if (fieldMults_.size() == 0) {
	auto policy = panzer::HP::inst().teamPolicy<ScalarT,SharedFieldMultTag<0>,PHX::Device>(workset.num_cells).set_scratch_size(0,Kokkos::PerTeam(bytes));
	parallel_for(this->getName(), policy, *this);
      } else if (fieldMults_.size() == 1) {
	auto policy = panzer::HP::inst().teamPolicy<ScalarT,SharedFieldMultTag<1>,PHX::Device>(workset.num_cells).set_scratch_size(0,Kokkos::PerTeam(bytes));
	parallel_for(this->getName(), policy, *this);
      } else {
	auto policy = panzer::HP::inst().teamPolicy<ScalarT,SharedFieldMultTag<-1>,PHX::Device>(workset.num_cells).set_scratch_size(0,Kokkos::PerTeam(bytes));
	parallel_for(this->getName(), policy, *this);
      }
    }
    else {
      // The following if-block is for the sake of optimization depending on the
      // number of field multipliers.  The parallel_fors will loop over the cells
      // in the Workset and execute operator()() above.
      if (fieldMults_.size() == 0) {
	auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag<0>,PHX::Device>(workset.num_cells);
	parallel_for(this->getName(), policy, *this);
      } else if (fieldMults_.size() == 1) {
	auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag<1>,PHX::Device>(workset.num_cells);
	parallel_for(this->getName(), policy, *this);
      } else {
	auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag<-1>,PHX::Device>(workset.num_cells);
	parallel_for(this->getName(), policy, *this);
      }
    }
  } // end of evaluateFields()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getValidParameters()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename TRAITS>
  Teuchos::RCP<Teuchos::ParameterList>
  Integrator_DivBasisTimesScalar<EvalT, TRAITS>::getValidParameters() const
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
    p->set<string>("Test Field Name", "?");
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

#endif // __Panzer_Integrator_DivBasisTimesScalar_impl_hpp__
