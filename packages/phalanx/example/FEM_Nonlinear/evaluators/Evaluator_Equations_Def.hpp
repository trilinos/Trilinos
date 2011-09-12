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


#include "Teuchos_TestForException.hpp"
#include "Phalanx_DataLayout.hpp"

//**********************************************************************
PHX_EVALUATOR_CTOR(Equations,p) :
  temp("Temperature", 
       p.get< Teuchos::RCP<PHX::DataLayout> >("QP Data Layout")), 
  vel("Velocity X", 
      p.get< Teuchos::RCP<PHX::DataLayout> >("QP Data Layout")),
  grad_temp("Temperature Gradient", 
	    p.get< Teuchos::RCP<PHX::DataLayout> >("Gradient QP Data Layout")), 
  grad_vel("Velocity X Gradient", 
	   p.get< Teuchos::RCP<PHX::DataLayout> >("Gradient QP Data Layout")),
  residual_temp("Residual Temperature", 
		p.get< Teuchos::RCP<PHX::DataLayout> >("Node Data Layout")), 
  residual_vel("Residual Velocity X", 
	       p.get< Teuchos::RCP<PHX::DataLayout> >("Node Data Layout"))
{ 
  this->addDependentField(temp);
  this->addDependentField(vel);
  this->addDependentField(grad_temp);
  this->addDependentField(grad_vel);

  this->addEvaluatedField(residual_temp);
  this->addEvaluatedField(residual_vel);
  
  this->setName("Equations");
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Equations,data,fm)
{
  this->utils.setFieldData(temp,fm);
  this->utils.setFieldData(vel,fm);
  this->utils.setFieldData(grad_temp,fm);
  this->utils.setFieldData(grad_vel,fm);
  this->utils.setFieldData(residual_temp,fm);
  this->utils.setFieldData(residual_vel,fm);

  num_qp = grad_temp.dimension(1);
  num_dim = grad_temp.dimension(2);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Equations,workset)
{ 
  for (int i=0; i < residual_temp.size(); ++i)
    residual_temp[i] = 0.0;

  for (int i=0; i < residual_vel.size(); ++i)
    residual_vel[i] = 0.0;

  std::vector<Element_Linear2D>::iterator element = workset.begin;
  
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell,++element) {
    
    const shards::Array<double,shards::NaturalOrder,QuadPoint,Node>& phi = 
      element->basisFunctions();

    const shards::Array<double,shards::NaturalOrder,QuadPoint,Node,Dim>& 
      grad_phi = element->basisFunctionGradientsRealSpace();

    const shards::Array<double,shards::NaturalOrder,QuadPoint>& det_jac = 
      element->detJacobian();

    const shards::Array<double,shards::NaturalOrder,QuadPoint>& weights = 
      element->quadratureWeights();

    for (int node = 0; node < element->numNodes(); ++node) {
      
      for (int qp = 0; qp < num_qp; ++qp) {

	residual_temp(cell,node) += det_jac(qp) * weights(qp) *
	  ( grad_phi(qp,node,0) * grad_temp(cell,qp,0) 
	    + grad_phi(qp,node,1) * grad_temp(cell,qp,1)
	    + 1000.0 * phi(qp,node) * temp(cell,qp) * vel(cell,qp) );
	
	residual_vel(cell,node) += det_jac(qp) * weights(qp) *
	  ( grad_phi(qp,node,0) * grad_vel(cell,qp,0) 
	    + grad_phi(qp,node,1) * grad_vel(cell,qp,1)
	    + 1000.0 * phi(qp,node) * temp(cell,qp) * vel(cell,qp) );
	
      }
      
    }

  }
   
//   std::cout << "Temp Residual" << std::endl;
//   residual_temp.print(std::cout,true);
//   std::cout << "Vx Residual" << std::endl;
//   residual_vel.print(std::cout,true);
 
}

//**********************************************************************
