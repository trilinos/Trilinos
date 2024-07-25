// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
