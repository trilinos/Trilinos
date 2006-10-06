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

#ifndef FEAPP_BRUSSELATORPDE_HPP
#define FEAPP_BRUSSELATORPDE_HPP

#include "Teuchos_RefCountPtr.hpp"

#include "FEApp_AbstractPDE.hpp"

namespace FEApp {

  template <typename ScalarT>
  class BrusselatorPDE : public FEApp::AbstractPDE<ScalarT> {
  public:
  
    //! Constructor
    BrusselatorPDE(double alpha, double beta, double D1, double D2);

    //! Destructor
    virtual ~BrusselatorPDE();

    //! Number of discretized equations
    virtual unsigned int numEquations() const;

    //! Initialize PDE
    virtual void init(unsigned int numQuadPoints, unsigned int numNodes);

    //! Evaluate discretized PDE element-level residual
    virtual void
    evaluateElementResidual(const FEApp::AbstractQuadrature& quadRule,
			    const FEApp::AbstractElement& element,
			    const std::vector<ScalarT>& solution,
			    std::vector<ScalarT>& residual);

  private:
    
    //! Private to prohibit copying
    BrusselatorPDE(const BrusselatorPDE&);
    
    //! Private to prohibit copying
    BrusselatorPDE& operator=(const BrusselatorPDE&);
    
  protected:
    
    //! Number of quad points
    unsigned int num_qp;
    
    //! Number of nodes
    unsigned int num_nodes;
    
    //! Shape function values
    std::vector< std::vector<double> > phi;

    //! Shape function derivatives
    std::vector< std::vector<double> > dphi;

    //! Element transformation Jacobian
    std::vector<double> jac;

    //! Discretized solution
    std::vector<ScalarT> T;

    //! Discretized solution
    std::vector<ScalarT> C;

    //! Discretized solution derivative
    std::vector<ScalarT> dT;

    //! Discretized solution derivative
    std::vector<ScalarT> dC;

    //! Model parameters
    double alpha, beta, D1, D2;

  };

  class BrusselatorPDE_TemplateBuilder {
  public:
    BrusselatorPDE_TemplateBuilder(double alpha_, double beta_, double D1_, 
				   double D2_) :
      alpha(alpha_), beta(beta_), D1(D1_), D2(D2_) {}
    template <typename T>
    Teuchos::RefCountPtr<FEApp::AbstractPDE_NTBase> build() const {
      return Teuchos::rcp( new BrusselatorPDE<T>(alpha, beta, D1, D2));
    }
  protected:
    double alpha, beta, D1, D2;
  };

}

// Include implementation
#include "FEApp_BrusselatorPDEImpl.hpp"

#endif // FEAPP_BRUSSELATORPDE_HPP
