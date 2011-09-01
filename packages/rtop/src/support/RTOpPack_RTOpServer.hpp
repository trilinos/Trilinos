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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
  const Teuchos::RCP<Teuchos::AbstractFactory<RTOpPack::RTOpT<Scalar> > > &op_factory
  )
{
  // ToDo: RAB: 20030620: Validate op_factory properly!
  op_factories_[op_factory->create()->op_name()] = op_factory;
}

template<class Scalar>
Teuchos::RCP<Teuchos::AbstractFactory<RTOpPack::RTOpT<Scalar> > >
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
    o << setw(w) << itr->first << typeName(*itr->second->create()) << std::endl;
  }
  o << std::endl;
}

} // namespace RTOpPack

#endif // RTOP_SERVER_HPP
