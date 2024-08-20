// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DOF_GRADIENT_IMPL_HPP
#define PANZER_DOF_GRADIENT_IMPL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_HierarchicParallelism.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

namespace panzer {

namespace {

  template <typename ScalarT,typename ArrayT>
  struct evaluateGrad_withSens {
    using scratch_view = Kokkos::View<ScalarT*,typename PHX::DevLayout<ScalarT>::type,typename PHX::exec_space::scratch_memory_space,Kokkos::MemoryUnmanaged>;
    using team_policy = Kokkos::TeamPolicy<PHX::exec_space>::member_type;

    PHX::MDField<ScalarT>  dof_grad_;
    PHX::MDField<const ScalarT,Cell,Point>  dof_value_;
    const ArrayT grad_basis_;
    const int num_fields_;
    const int num_points_;
    const int space_dim_;
    const int fad_size_;
    const bool use_shared_memory_;

    evaluateGrad_withSens(PHX::MDField<ScalarT> dof_grad,
			  PHX::MDField<const ScalarT,Cell,Point> dof_value,
			  const ArrayT  grad_basis,
			  const bool use_shared_memory) :
      dof_grad_(dof_grad),
      dof_value_(dof_value),
      grad_basis_(grad_basis),
      num_fields_(grad_basis.extent(1)),
      num_points_(grad_basis.extent(2)),
      space_dim_(grad_basis.extent(3)),
      fad_size_(Kokkos::dimension_scalar(dof_grad.get_view())),
      use_shared_memory_(use_shared_memory) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const team_policy& team) const
    {
      const int cell = team.league_rank();

      if (not use_shared_memory_) {
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points_), [&] (const int& pt) {
	  for (int d=0; d<space_dim_; ++d) {
	    // first initialize to the right thing (prevents over writing with 0)
	    // then loop over one less basis function
	    dof_grad_(cell,pt,d) = dof_value_(cell, 0) * grad_basis_(cell, 0, pt, d);
	    for (int bf=1; bf<num_fields_; ++bf)
	      dof_grad_(cell,pt,d) += dof_value_(cell, bf) * grad_basis_(cell, bf, pt, d);
	  }
	});
      } else {
        scratch_view dof_values;
        scratch_view point_values;
        if (Sacado::IsADType<ScalarT>::value) {
          dof_values = scratch_view(team.team_shmem(),num_fields_,fad_size_);
          point_values = scratch_view(team.team_shmem(),num_points_,fad_size_);
        }
        else {
          dof_values = scratch_view(team.team_shmem(),num_fields_);
          point_values = scratch_view(team.team_shmem(),num_points_);
        }

	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_fields_), [&] (const int& dof) {
	  dof_values(dof) = dof_value_(cell,dof);
	});
	team.team_barrier();

	for (int dim=0; dim < space_dim_; ++dim) {
	  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points_), [&] (const int& pt) {
	    point_values(pt) = 0.0;
	  });
	  // Perform contraction
	  for (int dof=0; dof<num_fields_; ++dof) {
	    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points_), [&] (const int& pt) {
	      point_values(pt) += dof_values(dof) * grad_basis_(cell,dof,pt,dim);
	    });
	  }
	  // Copy to main memory
	  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points_), [&] (const int& pt) {
	    dof_grad_(cell,pt,dim) = point_values(pt);
	  });
	} // loop over dim
      }
    }

    size_t team_shmem_size(int /* team_size */ ) const
    {
      if (not use_shared_memory_)
	return 0;

      size_t bytes;
      if (Sacado::IsADType<ScalarT>::value)
	bytes = scratch_view::shmem_size(num_fields_,fad_size_) + scratch_view::shmem_size(num_points_,fad_size_);
      else
	bytes = scratch_view::shmem_size(num_fields_) + scratch_view::shmem_size(num_points_);
      return bytes;
    }
  };

} // anonymous namespace

//**********************************************************************
template<typename EvalT, typename Traits>
DOFGradient<EvalT, Traits>::
DOFGradient(
  const Teuchos::ParameterList& p) :
  use_descriptors_(false),
  dof_value( p.get<std::string>("Name"),
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  dof_gradient( p.get<std::string>("Gradient Name"),
		p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector ),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsGrad(),std::logic_error,
                             "DOFGradient: Basis of type \"" << basis->name() << "\" does not support GRAD");

  this->addEvaluatedField(dof_gradient);
  this->addDependentField(dof_value);

  std::string n = "DOFGradient: " + dof_gradient.fieldTag().name() + " ("+PHX::print<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
DOFGradient<EvalT, TRAITS>::
DOFGradient(const PHX::FieldTag & input,
            const PHX::FieldTag & output,
            const panzer::BasisDescriptor & bd,
            const panzer::IntegrationDescriptor & id)
  : use_descriptors_(true)
  , bd_(bd)
  , id_(id)
  , dof_value(input)
  , dof_gradient(output)
{
  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(bd_.getType()=="HGrad",std::logic_error,
                             "DOFGradient: Basis of type \"" << bd_.getType() << "\" does not support GRAD");

  this->addEvaluatedField(dof_gradient);
  this->addDependentField(dof_value);

  std::string n = "DOFGradient: " + dof_gradient.fieldTag().name() + " ("+PHX::print<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DOFGradient<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& /* fm */)
{
  if(not use_descriptors_)
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DOFGradient<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  if (workset.num_cells == 0)
    return;

  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  using Array=typename panzer::BasisValues2<double>::ConstArray_CellBasisIPDim;
  Array grad_basis = use_descriptors_ ? basisValues.getGradBasisValues(false) : Array(basisValues.grad_basis);

  bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();
  evaluateGrad_withSens<ScalarT, Array> eval(dof_gradient,dof_value,grad_basis,use_shared_memory);
  auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::Device>(workset.num_cells);
  Kokkos::parallel_for("panzer::DOFGradient::evaluateFields", policy, eval);
}

//**********************************************************************

}

#endif
