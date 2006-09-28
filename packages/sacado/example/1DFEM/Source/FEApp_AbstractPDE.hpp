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

#ifndef FEAPP_ABSTRACTPDE_HPP
#define FEAPP_ABSTRACTPDE_HPP

#include <vector>

#include "FEApp_AbstractElement.hpp"
#include "FEApp_AbstractQuadrature.hpp"

namespace FEApp {

  /*!
   * \brief Abstract interface for representing a discretized 1-D PDE.
   */
  template <typename ScalarT>
  class AbstractPDE {
  public:
  
    //! Default constructor
    AbstractPDE() {};

    //! Destructor
    virtual ~AbstractPDE() {};

    //! Number of discretized equations
    virtual unsigned int numEquations() const = 0;

    //! Initialize PDE
    virtual void init(unsigned int numQuadPoints, unsigned int numNodes) = 0;

    //! Evaluate discretized PDE element-level residual
    virtual void
    evaluateElementResidual(const FEApp::AbstractQuadrature& quadRule,
			    const FEApp::AbstractElement& element,
			    const std::vector<ScalarT>& solution,
			    std::vector<ScalarT>& residual) = 0;

  private:
    
    //! Private to prohibit copying
    AbstractPDE(const AbstractPDE&);

    //! Private to prohibit copying
    AbstractPDE& operator=(const AbstractPDE&);

  };

}

#endif // FEAPP_ABSTRACTPDE_HPP
