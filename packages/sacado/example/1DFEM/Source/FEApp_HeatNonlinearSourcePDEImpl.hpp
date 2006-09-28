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
FEApp::HeatNonlinearSourcePDE<ScalarT>::
HeatNonlinearSourcePDE(const Teuchos::RefCountPtr< const FEApp::AbstractSourceFunction<ScalarT> >& src_func) : 
  source(src_func),
  num_qp(0),
  num_nodes(0),
  phi(),
  dphidxi(),
  jac(),
  u(),
  f()
{
}

template <typename ScalarT>
FEApp::HeatNonlinearSourcePDE<ScalarT>::
~HeatNonlinearSourcePDE()
{
}

template <typename ScalarT>
unsigned int 
FEApp::HeatNonlinearSourcePDE<ScalarT>::
numEquations() const
{
  return 1;
}

template <typename ScalarT>
void
FEApp::HeatNonlinearSourcePDE<ScalarT>::
init(unsigned int numQuadPoints, unsigned int numNodes)
{
  num_qp = numQuadPoints;
  num_nodes = numNodes;

  phi.resize(num_qp);
  dphidxi.resize(num_qp);
  jac.resize(num_qp);
  u.resize(num_qp);
  f.resize(num_qp);

  for (unsigned int i=0; i<num_qp; i++) {
    phi[i].resize(num_nodes);
    dphidxi[i].resize(num_nodes);
  }
}

template <typename ScalarT>
void
FEApp::HeatNonlinearSourcePDE<ScalarT>::
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
  element.evaluateShapeDerivs(xi, dphidxi);

  // Evaluate Jacobian of transformation to standard element
  element.evaluateJacobian(xi, jac);

  // Compute u
  for (unsigned int i=0; i<num_qp; i++) {
    u[i] = 0.0;
    for (unsigned int j=0; j<num_nodes; j++)
      u[i] += solution[j] * phi[i][j];
  }

  // Evaluate source function
  source->evaluate(u, f);

  // Evaluate residual
  ScalarT tmp;
  for (unsigned int i=0; i<num_nodes; i++) {
    residual[i] = 0.0;
    for (unsigned int k=0; k<num_qp; k++) {
      tmp = 0.0;
      for (unsigned int j=0; j<num_nodes; j++)
	tmp += solution[j] * dphidxi[k][j];
      residual[i] += 
	w[k] * (f[k]*phi[k][i]*jac[k] + tmp*dphidxi[k][i] / jac[k]);
    }
  }

}
