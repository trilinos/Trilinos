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

#ifndef RTOP_SERVER_DECL_HPP
#define RTOP_SERVER_DECL_HPP

#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace RTOpPack {

/** \brief Server for creating <tt>RTOpT</tt> objects given just the
 * operators name.
 *
 * \ingroup RTOpPack_Client_Server_support_grp
 */
template<class Scalar>
class RTOpServer {
public:

  /** \brief Add a new abstract factory for an <tt>RTOpT</tt> operator.
   *
   * @param  op_factory  [in] Smart pointer to factory object that will create
   *                     <tt>RTOpT</tt> objects of a given type.  The name of
   *                     the operator type will be extracted from
   *                     <tt>op_factory->create()->op_name()</tt>.
   *
   * Preconditions:<ul>
   * \item <tt>op_factory.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:<ul>
   * \item <tt>this->get_op_factory(op_factory->create()->op_name()).get() != op_factory.get()</tt>
   * </ul>
   */
  void add_op_factory( const Teuchos::RCP<Teuchos::AbstractFactory<RTOpPack::RTOpT<Scalar> > > &op_factory );

  /** \brief Get an operator factory given the name of the operator.
   *
   * @param  op_name  [in] Null-terminated C-style string for the name of an
   *                  operator type.
   *
   * Preconditions:<ul>
   * \item An operator factory with this name must have been added in a call to <tt>this->add_op_factory(...)</tt>
   * </ul>
   */
  Teuchos::RCP<Teuchos::AbstractFactory<RTOpPack::RTOpT<Scalar> > > get_op_factory( const char op_name[] ) const;

  /** \brief Print out all of the operator factories that have been added.
   */
  void print_op_factories(std::ostream& o) const;

private:
  typedef std::map< std::string, Teuchos::RCP<Teuchos::AbstractFactory<RTOpT<Scalar> > > >  op_factories_t;
  op_factories_t  op_factories_;

}; // class RTOpServer

} // namespace RTOpPack

#endif // RTOP_SERVER_DECL_HPP
