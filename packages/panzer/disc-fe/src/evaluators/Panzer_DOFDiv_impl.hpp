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

#ifndef PANZER_DOF_DIV_IMPL_HPP
#define PANZER_DOF_DIV_IMPL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_HierarchicParallelism.hpp"

namespace panzer {

//**********************************************************************
template<typename ScalarT,typename ArrayT>
class EvaluateDOFDiv_withSens {
  PHX::MDField<ScalarT,Cell,IP> dof_div;
  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  const ArrayT div_basis;
  const int num_fields;
  const int num_points;
  const int fad_size;
  const bool use_shared_memory;

public:
  using scratch_view = Kokkos::View<ScalarT* ,typename PHX::DevLayout<ScalarT>::type,typename PHX::exec_space::scratch_memory_space,Kokkos::MemoryUnmanaged>;

  EvaluateDOFDiv_withSens(PHX::MDField<ScalarT,Cell,IP> & in_dof_div,
                          PHX::MDField<const ScalarT,Cell,Point> & in_dof_value,
                          const ArrayT & in_div_basis,
			  const bool in_use_shared_memory = false)
    : dof_div(in_dof_div),
      dof_value(in_dof_value),
      div_basis(in_div_basis),
      num_fields(static_cast<int>(div_basis.extent(1))),
      num_points(static_cast<int>(div_basis.extent(2))),
      fad_size(static_cast<int>(Kokkos::dimension_scalar(in_dof_div.get_view()))),
      use_shared_memory(in_use_shared_memory)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
  {
    const int cell = team.league_rank();

    if (use_shared_memory) {
      // Copy reused data into fast scratch space
      scratch_view dof_values;
      scratch_view point_values;
      if (Sacado::IsADType<ScalarT>::value) {
        dof_values = scratch_view(team.team_shmem(),num_fields,fad_size);
        point_values = scratch_view(team.team_shmem(),num_points,fad_size);
      }
      else {
        dof_values = scratch_view(team.team_shmem(),num_fields);
        point_values = scratch_view(team.team_shmem(),num_points);
      }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_fields), [&] (const int& dof) {
	dof_values(dof) = dof_value(cell,dof);
      });

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int& pt) {
	point_values(pt) = 0.0;
      });

      team.team_barrier();

      // Perform contraction
      for (int dof=0; dof<num_fields; ++dof) {
	Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int& pt) {
	  point_values(pt) += dof_values(dof) * div_basis(cell,dof,pt);
        });
      }

      // Copy to main memory
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int& pt) {
	dof_div(cell,pt) = point_values(pt);
      });
    }
    else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,0,num_points), [&] (const int& pt) {
	// first initialize to the right thing (prevents over writing with 0)
	// then loop over one less basis function
	// ScalarT & div = dof_div(cell,pt);
	dof_div(cell,pt) = dof_value(cell, 0) * div_basis(cell, 0, pt);
	for (int bf=1; bf<num_fields; bf++)
	  dof_div(cell,pt) += dof_value(cell, bf) * div_basis(cell, bf, pt);
      });
    }
  }
  size_t team_shmem_size(int /* team_size */ ) const
  {
    if (not use_shared_memory)
      return 0;

    size_t bytes;
    if (Sacado::IsADType<ScalarT>::value)
      bytes = scratch_view::shmem_size(num_fields,fad_size) + scratch_view::shmem_size(num_points,fad_size);
    else
      bytes = scratch_view::shmem_size(num_fields) + scratch_view::shmem_size(num_points);
    return bytes;
  }

};

//**********************************************************************
// MOST EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename TRAITS>
DOFDiv<EvalT, TRAITS>::
DOFDiv(const Teuchos::ParameterList & p) :
  use_descriptors_(false),
  dof_value( p.get<std::string>("Name"),
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the div operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsDiv(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" does not support DIV");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" in DOF Div should require orientations. So we are throwing.");

  // build dof_div
  dof_div = PHX::MDField<ScalarT,Cell,IP>(p.get<std::string>("Div Name"),
      	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );

  // add to evaluation graph
  this->addEvaluatedField(dof_div);
  this->addDependentField(dof_value);

  std::string n = "DOFDiv: " + dof_div.fieldTag().name() + " ("+PHX::print<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
DOFDiv<EvalT, TRAITS>::
DOFDiv(const PHX::FieldTag & input,
       const PHX::FieldTag & output,
       const panzer::BasisDescriptor & bd,
       const panzer::IntegrationDescriptor & id)
  : use_descriptors_(true)
  , bd_(bd)
  , id_(id)
  , dof_value(input)
{
  TEUCHOS_ASSERT(bd.getType()=="HDiv");

  // build dof_div
  dof_div = output;

  // add to evaluation graph
  this->addEvaluatedField(dof_div);
  this->addDependentField(dof_value);

  std::string n = "DOFDiv: " + dof_div.fieldTag().name() + " ("+PHX::print<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
void DOFDiv<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  if(not use_descriptors_)
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
void DOFDiv<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{
  if (workset.num_cells == 0)
    return;

  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  const bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();
  auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::exec_space>(workset.num_cells);
  auto f = EvaluateDOFDiv_withSens<ScalarT,typename panzer::BasisValues2<double>::Array_CellBasisIP>(dof_div,dof_value,basisValues.div_basis,use_shared_memory);
  Kokkos::parallel_for(policy,f,this->getName());
}

//**********************************************************************

//**********************************************************************
// JACOBIAN EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename TRAITS>
DOFDiv<panzer::Traits::Jacobian, TRAITS>::
DOFDiv(const Teuchos::ParameterList & p) :
  use_descriptors_(false),
  dof_value( p.get<std::string>("Name"),
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // do you specialize because you know where the basis functions are and can
  // skip a large number of AD calculations?
  if(p.isType<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector")) {
    offsets = *p.get<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector");
    accelerate_jacobian = true;  // short cut for identity matrix
  }
  else
    accelerate_jacobian = false; // don't short cut for identity matrix
  accelerate_jacobian = false; // don't short cut for identity matrix

  // Verify that this basis supports the div operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsDiv(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" does not support DIV");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "DOFDiv: Basis of type \"" << basis->name() << "\" in DOF Div should require orientations. So we are throwing.");

  // build dof_div
  dof_div = PHX::MDField<ScalarT,Cell,IP>(p.get<std::string>("Div Name"),
      	                                  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );

  // add to evaluation graph
  this->addEvaluatedField(dof_div);
  this->addDependentField(dof_value);

  std::string n = "DOFDiv: " + dof_div.fieldTag().name() + " ("+PHX::print<panzer::Traits::Jacobian>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>
DOFDiv<panzer::Traits::Jacobian, TRAITS>::
DOFDiv(const PHX::FieldTag & input,
       const PHX::FieldTag & output,
       const panzer::BasisDescriptor & bd,
       const panzer::IntegrationDescriptor & id)
  : use_descriptors_(true)
  , bd_(bd)
  , id_(id)
  , dof_value(input)
{
  TEUCHOS_ASSERT(bd.getType()=="HDiv");

  // build dof_div
  dof_div = output;

  accelerate_jacobian = false; // don't short cut for identity matrix

  // add to evaluation graph
  this->addEvaluatedField(dof_div);
  this->addDependentField(dof_value);

  std::string n = "DOFDiv: " + dof_div.fieldTag().name() + " ("+PHX::print<panzer::Traits::Jacobian>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>
void DOFDiv<panzer::Traits::Jacobian, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_div,fm);

  if(not use_descriptors_)
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

template<typename TRAITS>
void DOFDiv<panzer::Traits::Jacobian,TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{
  if (workset.num_cells == 0)
    return;

  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  if(!accelerate_jacobian) {
    const bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();
    auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::exec_space>(workset.num_cells);
    auto f = EvaluateDOFDiv_withSens<ScalarT,typename panzer::BasisValues2<double>::Array_CellBasisIP>(dof_div,dof_value,basisValues.div_basis,use_shared_memory);
    Kokkos::parallel_for(policy,f,this->getName());
    return;
  }

  TEUCHOS_ASSERT(false);
}

}

#endif
