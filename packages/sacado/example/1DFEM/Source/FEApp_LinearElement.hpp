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

#ifndef FEAPP_LINEARELEMENT_HPP
#define FEAPP_LINEARELEMENT_HPP

#include "FEApp_AbstractElement.hpp"

namespace FEApp {

  /*!
   * \brief A 1-D linear finite element.
   */
  class LinearElement : public FEApp::AbstractElement {
  public:
    
    //! Default constructor
    LinearElement();

    //! Destructor
    virtual ~LinearElement();

    //! Get the number of nodes the element requires
    virtual unsigned int numNodes() const;

    //! Create the nodes for this element
    virtual void createNodes(double x_left, double x_right,
			     unsigned int first_node_gid);

    //! Return GID of ith node
    virtual unsigned int nodeGID(unsigned int i) const;
  
    //! Evaluate all shape functions at a set of points in (-1,1)
    virtual void 
    evaluateShapes(const std::vector<double>& xi,
		   std::vector< std::vector<double> >& phi) const;

    //! Evaluate all shape function derivatives at a set of points in (-1,1)
    virtual void
    evaluateShapeDerivs(const std::vector<double>& xi,
			std::vector< std::vector<double> >& dphidxi) const;

    /*
     * \brief Evaluate Jacobian of element transformation at a set of 
     * points in (-1,1)
     */
    virtual void 
    evaluateJacobian(const std::vector<double>& xi,
		     std::vector<double>& jac) const;

  private:

    //! Private to prohibit copying
    LinearElement(const LinearElement&);

    //! Private to prohibit copying
    LinearElement& operator=(const LinearElement&);

  protected:

    //! Coordinate of left node
    double xl;

    //! Coordinate of right node
    double xr;

    //! GID of left node
    unsigned int left_GID;

    //! GID of right node
    unsigned int right_GID;

  };

}

#endif // FEAPP_LINEARELEMENT_HPP
