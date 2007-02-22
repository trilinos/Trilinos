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

#ifndef FEAPP_ABSTRACTSOURCEFUNCTION_HPP
#define FEAPP_ABSTRACTSOURCEFUNCTION_HPP

#include <vector>

namespace FEApp {

  /*!
   * \brief Abstract interface for representing a PDE source function
   */
  template <typename ScalarT>
  class AbstractSourceFunction {
  public:
  
    //! Default constructor
    AbstractSourceFunction() {};

    //! Destructor
    virtual ~AbstractSourceFunction() {};

    //! Evaluate source function
    virtual void
    evaluate(const std::vector<ScalarT>& solution,
	     std::vector<ScalarT>& value) const = 0;

  private:
    
    //! Private to prohibit copying
    AbstractSourceFunction(const AbstractSourceFunction&);

    //! Private to prohibit copying
    AbstractSourceFunction& operator=(const AbstractSourceFunction&);

  };

}

#endif // FEAPP_ABSTRACTSOURCEFUNCTION_HPP
