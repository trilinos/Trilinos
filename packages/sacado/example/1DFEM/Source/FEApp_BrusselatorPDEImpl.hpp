// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

template <typename ScalarT>
FEApp::BrusselatorPDE<ScalarT>::
BrusselatorPDE(double alpha_, double beta_, double D1_, double D2_) : 
  num_qp(0),
  num_nodes(0),
  phi(),
  dphi(),
  jac(),
  T(),
  C(),
  dT(),
  dC(),
  alpha(alpha_),
  beta(beta_),
  D1(D1_),
  D2(D2_)
{
}

template <typename ScalarT>
FEApp::BrusselatorPDE<ScalarT>::
~BrusselatorPDE()
{
}

template <typename ScalarT>
unsigned int 
FEApp::BrusselatorPDE<ScalarT>::
numEquations() const
{
  return 2;
}

template <typename ScalarT>
void
FEApp::BrusselatorPDE<ScalarT>::
init(unsigned int numQuadPoints, unsigned int numNodes)
{
  num_qp = numQuadPoints;
  num_nodes = numNodes;

  phi.resize(num_qp);
  dphi.resize(num_qp);
  jac.resize(num_qp);
  T.resize(num_qp);
  C.resize(num_qp);
  dT.resize(num_qp);
  dC.resize(num_qp);

  for (unsigned int i=0; i<num_qp; i++) {
    phi[i].resize(num_nodes);
    dphi[i].resize(num_nodes);
  }
}

template <typename ScalarT>
void
FEApp::BrusselatorPDE<ScalarT>::
evaluateElementResidual(const FEApp::AbstractQuadrature& quadRule,
			const FEApp::AbstractElement& element,
			const std::vector<ScalarT>& solution,
			std::vector<ScalarT>& residual)
 {
  
  // Quadrature points
  const std::vector<double>& xi = quadRule.quadPoints();

  // Weights
  const std::vector<double>& w = quadRule.weights();

  // Evaluate shape functions
  element.evaluateShapes(xi, phi);

  // Evaluate shape function derivatives
  element.evaluateShapeDerivs(xi, dphi);

  // Evaluate Jacobian of transformation to standard element
  element.evaluateJacobian(xi, jac);

  // Compute discretizes solution
  for (unsigned int qp=0; qp<num_qp; qp++) {
    T[qp] = 0.0;
    C[qp] = 0.0;
    dT[qp] = 0.0;
    dC[qp] = 0.0;

    for (unsigned int node=0; node<num_nodes; node++) {
      T[qp] += solution[2*node] * phi[qp][node];
      C[qp] += solution[2*node+1] * phi[qp][node];
      dT[qp] += solution[2*node] * dphi[qp][node];
      dC[qp] += solution[2*node+1] * dphi[qp][node];
    }

  }

  // Evaluate residual
  for (unsigned int node=0; node<num_nodes; node++) {
    residual[2*node] = 0.0;
    residual[2*node+1] = 0.0;

    for (unsigned int qp=0; qp<num_qp; qp++) {
      residual[2*node] += 
	w[qp]*jac[qp]*((1.0/(jac[qp]*jac[qp]))*D1*dT[qp]*dphi[qp][node] + 
		       phi[qp][node]*(-alpha + (beta+1.0)*T[qp] - 
				      T[qp]*T[qp]*C[qp]));
      residual[2*node+1] += 
	w[qp]*jac[qp]*((1.0/(jac[qp]*jac[qp]))*D2*dC[qp]*dphi[qp][node] + 
		       phi[qp][node]*(-beta*T[qp] + T[qp]*T[qp]*C[qp]));
    }
  }

}
