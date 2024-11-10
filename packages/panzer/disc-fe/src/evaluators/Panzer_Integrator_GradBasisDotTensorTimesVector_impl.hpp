// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_Integrator_GradBasisDotTensorTimesVector_impl_hpp__
#define   __Panzer_Integrator_GradBasisDotTensorTimesVector_impl_hpp__

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
  Integrator_GradBasisDotTensorTimesVector<EvalT, Traits>::
  Integrator_GradBasisDotTensorTimesVector(
    const panzer::EvaluatorStyle&        evalStyle,
    const std::string&                   resName,
    const std::string&                   fluxName,
    const panzer::BasisIRLayout&         basis,
    const panzer::IntegrationRule&       ir,
    const std::string&                   tensorName,
    const Teuchos::RCP<PHX::DataLayout>& vecDL       /* = Teuchos::null */)
    :
    evalStyle_(evalStyle),
    basisName_(basis.name())
  {
    using PHX::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using PHX::DataLayout;
    using PHX::MDField;
    using std::invalid_argument;
    using std::logic_error;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(resName == "", invalid_argument, "Error:  "   \
      "Integrator_GradBasisDotTensorTimesVector called with an empty residual name.")
    TEUCHOS_TEST_FOR_EXCEPTION(fluxName == "", invalid_argument, "Error:  "   \
      "Integrator_GradBasisDotTensorTimesVector called with an empty flux name.")
    RCP<const PureBasis> tmpBasis = basis.getBasis();
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->supportsGrad(), logic_error,
      "Error:  Integrator_GradBasisDotTensorTimesVector:  Basis of type \""
      << tmpBasis->name() << "\" does not support the gradient operator.")
    RCP<DataLayout> tmpVecDL = ir.dl_vector;
    if (not vecDL.is_null())
    {
      tmpVecDL = vecDL;
      TEUCHOS_TEST_FOR_EXCEPTION(
        tmpVecDL->extent(2) < ir.dl_vector->extent(2), logic_error,
        "Integrator_GradBasisDotTensorTimesVector:  Dimension of space exceeds "         \
        "dimension of Vector Data Layout.");
    } // end if (not vecDL.is_null())

    // Create the field for the vector-valued function we're integrating.
    vector_ = MDField<const ScalarT, Cell, IP, Dim>(fluxName, tmpVecDL);
    this->addDependentField(vector_);

    // Create the field that we're either contributing to or evaluating
    // (storing).
    field_ = MDField<ScalarT, Cell, BASIS>(resName, basis.functional);
    if (evalStyle_ == EvaluatorStyle::CONTRIBUTES)
      this->addContributedField(field_);
    else // if (evalStyle_ == EvaluatorStyle::EVALUATES)
      this->addEvaluatedField(field_);

    // Add the tensor.
    tensor_ = MDField<const ScalarT, Cell, IP, Dim, Dim>(tensorName, ir.dl_tensor);
    this->addDependentField(tensor_);

    // Set the name of this object.
    string n("Integrator_GradBasisDotTensorTimesVector (");
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
  Integrator_GradBasisDotTensorTimesVector<EvalT, Traits>::
  Integrator_GradBasisDotTensorTimesVector(
    const Teuchos::ParameterList& p)
    :
    Integrator_GradBasisDotTensorTimesVector(
      panzer::EvaluatorStyle::EVALUATES,
      p.get<std::string>("Residual Name"),
      p.get<std::string>("Flux Name"),
      (*p.get<Teuchos::RCP<panzer::BasisIRLayout>>("Basis")),
      (*p.get<Teuchos::RCP<panzer::IntegrationRule>>("IR")),
      p.get<std::string>("Tensor Name"),
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
  Integrator_GradBasisDotTensorTimesVector<EvalT, Traits>::
  postRegistrationSetup(
    typename Traits::SetupData sd,
    PHX::FieldManager<Traits>& /* fm */)
  {
    using panzer::getBasisIndex;
    using std::size_t;

    // Get the PHX::View of the tensor.
    kokkosTensor_ = tensor_.get_static_view();

    // Determine the index in the Workset bases for our particular basis name.
    basisIndex_ = getBasisIndex(basisName_, (*sd.worksets_)[0], this->wda);

    // Allocate temporary if not using shared memory
    bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();
    if (!use_shared_memory) {
      if (Sacado::IsADType<ScalarT>::value) {
	const auto fadSize = Kokkos::dimension_scalar(field_.get_view());
	tmp_ = PHX::View<ScalarT*>("GradBasisDotTensorTimesVector::tmp_",field_.extent(0),fadSize);
      } else {
	tmp_ = PHX::View<ScalarT*>("GradBasisDotTensorTimesVector::tmp_",field_.extent(0));
      }
    }
  } // end of postRegistrationSetup()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  operator()() NO shared memory
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  KOKKOS_INLINE_FUNCTION
  void
  Integrator_GradBasisDotTensorTimesVector<EvalT, Traits>::
  operator()(
    const FieldMultTag& /* tag */,
    const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    using panzer::EvaluatorStyle;
    const int cell = team.league_rank();

    // Initialize the evaluated field.
    const int numQP(vector_.extent(1)), numDim(vector_.extent(2)),
              numBases(basis_.extent(1));
    if (evalStyle_ == EvaluatorStyle::EVALUATES)
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
        field_(cell, basis) = 0.0;
      });

    {
      // Loop over the quadrature points and dimensions of our vector fields,
      // scale the integrand by the multiplier, and then perform the
      // integration, looping over the bases.
      for (int qp(0); qp < numQP; ++qp)
      {
        for (int dim(0); dim < numDim; ++dim)
        {
          tmp_(cell) = 0.0;
          for (int dim2(0); dim2 < numDim; ++dim2)
            tmp_(cell) += kokkosTensor_(cell, qp, dim, dim2) * vector_(cell, qp, dim2);
	  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
	    field_(cell, basis) += basis_(cell, basis, qp, dim) * tmp_(cell);
	  });
        } // end loop over the dimensions of the vector field
      } // end loop over the quadrature points
    }
  } // end of operator()()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  operator()() With shared memory
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  KOKKOS_INLINE_FUNCTION
  void
  Integrator_GradBasisDotTensorTimesVector<EvalT, Traits>::
  operator()(
    const SharedFieldMultTag& /* tag */,
    const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    using panzer::EvaluatorStyle;
    const int cell = team.league_rank();
    const int numQP = vector_.extent(1);
    const int numDim = vector_.extent(2);
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
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
      tmp_field(basis) = 0.0;
    });

    {
      // Loop over the quadrature points and dimensions of our vector fields,
      // scale the integrand by the multiplier, and then perform the
      // integration, looping over the bases.
      for (int qp(0); qp < numQP; ++qp)
      {
        for (int dim(0); dim < numDim; ++dim)
        {
          tmp_(cell) = 0.0;
          for (int dim2(0); dim2 < numDim; ++dim2)
            tmp_(cell) += kokkosTensor_(cell, qp, dim, dim2) * vector_(cell, qp, dim2);
	  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases), [&] (const int basis) {
	    tmp_field(basis) += basis_(cell, basis, qp, dim) * tmp_(cell);
	  });
        } // end loop over the dimensions of the vector field
      } // end loop over the quadrature points
    }

    // Put final values into target field
    if (evalStyle_ == EvaluatorStyle::EVALUATES) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int basis) {
	field_(cell,basis) = tmp_field(basis);
      });
    }
    else { // Contributed
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,numBases),[&] (const int basis) {
	field_(cell,basis) += tmp_field(basis);
      });
    }

  } // end of operator()()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evaluateFields()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_GradBasisDotTensorTimesVector<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
    using Kokkos::parallel_for;
    using Kokkos::TeamPolicy;

    // Grab the basis information.
    basis_ = this->wda(workset).bases[basisIndex_]->weighted_grad_basis;

    bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();
    if (use_shared_memory) {
      int bytes;
      if (Sacado::IsADType<ScalarT>::value) {
	const int fadSize = Kokkos::dimension_scalar(field_.get_view());
	bytes = scratch_view::shmem_size(1,fadSize) + scratch_view::shmem_size(basis_.extent(1),fadSize);
      }
      else
	bytes = scratch_view::shmem_size(1) + scratch_view::shmem_size(basis_.extent(1));

      auto policy = panzer::HP::inst().teamPolicy<ScalarT,SharedFieldMultTag,PHX::Device>(workset.num_cells).set_scratch_size(0,Kokkos::PerTeam(bytes));
      parallel_for(this->getName(), policy, *this);
    }
    else {
      auto policy = panzer::HP::inst().teamPolicy<ScalarT,FieldMultTag,PHX::Device>(workset.num_cells);
      parallel_for(this->getName(), policy, *this);
    }
  } // end of evaluateFields()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getValidParameters()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename TRAITS>
  Teuchos::RCP<Teuchos::ParameterList>
  Integrator_GradBasisDotTensorTimesVector<EvalT, TRAITS>::
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
    p->set<std::string>("Tensor Name", "?");
    RCP<DataLayout> vecDL;
    p->set("Vector Data Layout", vecDL);
    return p;
  } // end of getValidParameters()

} // end of namespace panzer

#endif // __Panzer_Integrator_GradBasisDotTensorTimesVector_impl_hpp__
