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


#ifndef PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP
#define PHX_EXAMPLE_VP_NONLINEAR_SOURCE_DEF_HPP

//**********************************************************************
template<typename EvalT, typename Traits> NonlinearSource<EvalT, Traits>::
NonlinearSource(const Teuchos::ParameterList& p) :
  source("Nonlinear Source", 
	 p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  density("Density", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout")),
  temp("Temperature", p.get< Teuchos::RCP<PHX::DataLayout> >("Data Layout"))
{ 
  this->addEvaluatedField(source);
  this->addDependentField(density);
  this->addDependentField(temp);

  this->setName("NonlinearSource");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& vm)
{
  // NOTE: We no longer need manually call setFieldData(). It happens
  // automatically. We set them here to make sure backwards
  // compatiblity still works. We demonstrate both valid calling
  // pathways.
  this->utils.setFieldData(source,vm);
  vm.template getFieldData<ScalarT,EvalT>(density);
  vm.template getFieldData<ScalarT,EvalT>(temp);

  cell_data_size = source.size() / source.dimension(0);
}
//*********************************************************************
template<typename EvalT, typename Traits>
KOKKOS_INLINE_FUNCTION
void NonlinearSource<EvalT, Traits>::operator () (const int i) const
{
  for (PHX::index_size_type ip = 0; ip < static_cast<PHX::index_size_type>(density.dimension_1()); ++ip)
    source(i,ip) =  density(i,ip) * temp(i,ip) * temp(i,ip);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
evaluateFields(typename Traits::EvalData d)
{ 
 Kokkos::parallel_for(d.num_cells, *this);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
preEvaluate(typename Traits::PreEvalData d)
{ 
  using namespace std;
  cout << "In Source Pre Op" << endl;
}

//**********************************************************************
template< typename EvalT, typename Traits>
void NonlinearSource<EvalT, Traits>::
postEvaluate(typename Traits::PostEvalData d)
{ 
  using namespace std;
  cout << "In Source Post Op" << endl;
}

//**********************************************************************

#endif
