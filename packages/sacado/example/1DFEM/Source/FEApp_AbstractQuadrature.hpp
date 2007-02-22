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

#ifndef FEAPP_ABSTRACTQUADRATURE_HPP
#define FEAPP_ABSTRACTQUADRATURE_HPP

#include <vector>

namespace FEApp {

  /*!
   * \brief A simple abstract quadrature class for integrating functions 
   over the interval \f$(-1,1)\f$.
   */
  class AbstractQuadrature {
  public:

    //! Default constructor
    AbstractQuadrature() {};

    //! Destructor
    virtual ~AbstractQuadrature() {}

    //! Return the number of quadrature points
    virtual unsigned int numPoints() const = 0;

    //! Return the quadrature points
    virtual const std::vector<double>& quadPoints() const = 0;

    //! Return the weights
    virtual const std::vector<double>& weights() const = 0;

  private:

    //! Private to prohibit copying
    AbstractQuadrature(const AbstractQuadrature&);

    //! Private to prohibit copying
    AbstractQuadrature& operator=(const AbstractQuadrature&);

  };

}

#endif // FEAPP_ABSTRACTQUADRATURE_HPP
