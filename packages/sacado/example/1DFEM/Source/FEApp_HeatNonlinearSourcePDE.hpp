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

#ifndef FEAPP_HEATNONLINEARSOURCEPDE_HPP
#define FEAPP_HEATNONLINEARSOURCEPDE_HPP

#include "Teuchos_RefCountPtr.hpp"

#include "FEApp_AbstractPDE.hpp"
#include "FEApp_SourceFunctionFactory.hpp"

namespace FEApp {

  template <typename ScalarT>
  class HeatNonlinearSourcePDE : public FEApp::AbstractPDE<ScalarT> {
  public:
  
    //! Constructor
    HeatNonlinearSourcePDE(const Teuchos::RefCountPtr< const FEApp::AbstractSourceFunction<ScalarT> >& src_func);

    //! Destructor
    virtual ~HeatNonlinearSourcePDE();

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
    HeatNonlinearSourcePDE(const HeatNonlinearSourcePDE&);

    //! Private to prohibit copying
    HeatNonlinearSourcePDE& operator=(const HeatNonlinearSourcePDE&);

  protected:
    
    //! Pointer to source function
    Teuchos::RefCountPtr< const FEApp::AbstractSourceFunction<ScalarT> > source;

    //! Number of quad points
    unsigned int num_qp;
    
    //! Number of nodes
    unsigned int num_nodes;

    //! Shape function values
    std::vector< std::vector<double> > phi;

    //! Shape function derivatives
    std::vector< std::vector<double> > dphidxi;

    //! Element transformation Jacobian
    std::vector<double> jac;

    //! Discretized solution
    std::vector<ScalarT> u;

    //! Source function values
    std::vector<ScalarT> f;

  };

  class HeatNonlinearSourcePDE_TemplateBuilder {
  public:
    HeatNonlinearSourcePDE_TemplateBuilder(
		const Teuchos::RefCountPtr<Teuchos::ParameterList>& params_) :
      params(Teuchos::rcp(&(params_->sublist("Source Function")),false)) {}
    template <typename T>
    Teuchos::RefCountPtr<FEApp::AbstractPDE_NTBase> build() const {
      FEApp::SourceFunctionFactory<T> factory(params);
      Teuchos::RefCountPtr< FEApp::AbstractSourceFunction<T> > source =
	factory.create();
      return Teuchos::rcp( new FEApp::HeatNonlinearSourcePDE<T>(source));
    }
  protected:
    Teuchos::RefCountPtr<Teuchos::ParameterList> params;
  };

}

// Include implementation
#include "FEApp_HeatNonlinearSourcePDEImpl.hpp"

#endif // FEAPP_HEATNONLINERASOURCEPDE_HPP
