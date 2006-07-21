// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOP_SERVER_HPP
#define RTOP_SERVER_HPP

//#include <typeinfo>
//#include <ostream>
//#include <iomanip>

#include "RTOpPack_RTOpServerDecl.hpp"

namespace RTOpPack {

template<class Scalar>
void RTOpServer<Scalar>::add_op_factory(
  const Teuchos::RefCountPtr<Teuchos::AbstractFactory<RTOpPack::RTOpT<Scalar> > > &op_factory
  )
{
  // ToDo: RAB: 20030620: Validate op_factory properly!
  op_factories_[op_factory->create()->op_name()] = op_factory;
}

template<class Scalar>
Teuchos::RefCountPtr<Teuchos::AbstractFactory<RTOpPack::RTOpT<Scalar> > >
RTOpServer<Scalar>::get_op_factory( const char op_name[] ) const
{
  typename op_factories_t::const_iterator itr = op_factories_.find(op_name);
  TEST_FOR_EXCEPTION(
    itr == op_factories_.end(), std::logic_error
    ,"RTOpServer<Scalar>::get_op_factory(...): Error, an operator factory with the "
    "operator name \'" << op_name << "\' does not exist!"
    );
  return itr->second;
}

template<class Scalar>
void RTOpServer<Scalar>::print_op_factories(std::ostream& o) const
{
  using std::setw;
  const int w = 30;
  o << "\nRTOpServer<Scalar>::print_op_factories(...): RTOp operator factories currently registered\n\n" << std::left;
  o << setw(w) << "Operator name" << "Operator type" << std::endl;
  o << setw(w) << "-------------" << "-------------" << std::endl;
  for( typename op_factories_t::const_iterator itr = op_factories_.begin(); itr != op_factories_.end(); ++itr ) {
    o << setw(w) << itr->first << typeid(*itr->second->create()).name() << std::endl;
  }
  o << std::endl;
}

} // namespace RTOpPack

#endif // RTOP_SERVER_HPP
