// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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

#include "FEApp_AbstractPDE_NTBase.hpp"
#include "FEApp_AbstractElement.hpp"
#include "FEApp_AbstractQuadrature.hpp"

#include "Sacado_TemplateManager.hpp"

namespace FEApp {

  /*!
   * \brief Abstract interface for representing a discretized 1-D PDE.
   */
  template <typename ScalarT>
  class AbstractPDE : public FEApp::AbstractPDE_NTBase {
  public:
  
    //! Default constructor
    AbstractPDE() {};

    //! Destructor
    virtual ~AbstractPDE() {};

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

  template <typename TypeSeq>
  class AbstractPDE_TemplateManager : 
    public Sacado::TemplateManager<TypeSeq, FEApp::AbstractPDE_NTBase,
				   AbstractPDE> {
  public:
    AbstractPDE_TemplateManager() {}
    ~AbstractPDE_TemplateManager() {}
  };

}

#endif // FEAPP_ABSTRACTPDE_HPP
