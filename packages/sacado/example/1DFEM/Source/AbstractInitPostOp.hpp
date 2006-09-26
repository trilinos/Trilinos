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

#ifndef ABSTRACTINITPOSTOP_HPP
#define ABSTRACTINITPOSTOP_HPP

#include <vector>

#include "AbstractElement.hpp"

template <typename ScalarT>
class AbstractInitPostOp {
public:

  //! Fill type
  typedef ScalarT fill_type;

  //! Constructor
  AbstractInitPostOp() {};

  //! Destructor
  virtual ~AbstractInitPostOp() {};

  //! Evaulate init operator
  virtual void evalInit(const AbstractElement& e,
			unsigned int neqn,
			std::vector<ScalarT>& elem_x) = 0;

  //! Evaluate post operator
  virtual void evalPost(const AbstractElement& e,
			unsigned int neqn,
			std::vector<ScalarT>& elem_f) = 0;

private:

  //! Private to prohibit copying
  AbstractInitPostOp(const AbstractInitPostOp&);

  //! Private to prohibit copying
  AbstractInitPostOp& operator=(const AbstractInitPostOp&);

};

#endif // ABSTRACTINITPOSTOP_HPP
