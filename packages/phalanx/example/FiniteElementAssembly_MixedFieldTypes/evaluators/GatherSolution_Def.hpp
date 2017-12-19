// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
GatherSolution<PHX::MyTraits::Residual, Traits>::
GatherSolution(const std::string& field_name,
               const Teuchos::RCP<PHX::DataLayout>& layout,
               const int& in_num_equations,
               const int& in_field_index,
               const Kokkos::View<double*,PHX::Device>& in_x) :
  num_equations(in_num_equations),
  field_index(in_field_index),
  x(in_x)
{
  PHX::Tag<ScalarT> field_tag(field_name,layout);
  this->addEvaluatedField(field_tag,field);
  this->setName("Gather Solution: "+field_tag.name());
}

// **********************************************************************
template<typename Traits>
void GatherSolution<PHX::MyTraits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  gids = workset.gids_;
  cell_global_offset_index = workset.first_cell_global_index_;
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

// **********************************************************************
template<typename Traits>
void GatherSolution<PHX::MyTraits::Residual,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();
  if (team.team_rank() == 0) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,field.extent(1)), [&] (const int& node) {
      field(cell,node) = x( gids(cell_global_offset_index+cell,node) * num_equations + field_index);
    });
  }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template<typename Traits>
GatherSolution<PHX::MyTraits::Jacobian, Traits>::
GatherSolution(const std::string& field_name,
               const Teuchos::RCP<PHX::DataLayout>& layout,
               const int& in_num_equations,
               const int& in_field_index,
               const Kokkos::View<double*,PHX::Device>& in_x) :
  num_equations(in_num_equations),
  field_index(in_field_index),
  x(in_x)
{
  PHX::Tag<ScalarT> field_tag(field_name,layout);
  this->addEvaluatedField(field_tag,field);
  this->setName("Gather Solution: "+field_tag.name());
}

// **********************************************************************
template<typename Traits>
void GatherSolution<PHX::MyTraits::Jacobian, Traits>::
evaluateFields(typename Traits::EvalData workset)
{
  gids = workset.gids_;
  cell_global_offset_index = workset.first_cell_global_index_;
  Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::exec_space>(workset.num_cells_,workset.team_size_,workset.vector_size_),*this);
}

// **********************************************************************
template<typename Traits>
void GatherSolution<PHX::MyTraits::Jacobian,Traits>::
operator()(const Kokkos::TeamPolicy<PHX::exec_space>::member_type& team) const
{
  const int cell = team.league_rank();
  if (team.team_rank() == 0) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,field.extent(1)), [&] (const int& node) {
      field(cell,node).val() = x(gids(cell_global_offset_index+cell,node) * num_equations + field_index);
      field(cell,node).fastAccessDx(num_equations * node + field_index) = 1.0;
    });
  }
}

/*
// **********************************************************************
// Specialization: Jv
// **********************************************************************

template<typename Traits>
GatherSolution<PHX::MyTraits::Jv, Traits>::
GatherSolution(const std::string& field_name,
               const Teuchos::RCP<PHX::DataLayout>& layout,
               const int& in_field_index,
               const Kokkos::View<double*,PHX::Device>& in_x)
{
  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Solution Names"));

  Teuchos::RCP<PHX::DataLayout> dl = 
    p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout");

  x = p.get< Teuchos::RCP<Epetra_Vector> >("Solution Vector");

  val.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<ScalarT,Cell,Node> f(names[eq],dl);
    val[eq] = f;
    this->addEvaluatedField(val[eq]);
  }

  this->setName("Gather Solution");
}

// **********************************************************************
template<typename Traits> 
void GatherSolution<PHX::MyTraits::Jv, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  for (std::size_t eq = 0; eq < val.size(); ++eq)
    this->utils.setFieldData(val[eq],fm);

  num_nodes = val[0].dimension(1);
}

// **********************************************************************
template<typename Traits>
void GatherSolution<PHX::MyTraits::Jv, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
  const std::size_t num_eq = val.size();
  std::vector<Element_Linear2D>::iterator element = workset.begin;

  std::size_t cell = 0;
  for (; element != workset.end; ++element,++cell) {
    
    for (std::size_t node = 0; node < num_nodes; node++) {
      int node_GID = element->globalNodeId(node);
      int firstDOF = x->Map().LID(static_cast<int>(node_GID * num_eq));

      for (std::size_t eq = 0; eq < val.size(); eq++) {
	(val[eq])(cell,node) = 
	  Sacado::Fad::SFad<double,1>((*x)[firstDOF + eq]);
	(val[eq])(cell,node).fastAccessDx(0) = (*(workset.v))[firstDOF + eq];
	//(val[eq])(cell,node).fastAccessDx(num_eq * node + eq) = 1.0;

	// std::cout << "val[" << eq << "](" << cell << "," << node << ") = "
// 		  << (val[eq])(cell,node) << std::endl;
      }
    }

  }

}

// **********************************************************************
*/
