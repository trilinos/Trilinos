// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef EVIL_BASE_DECL_HPP
#define EVIL_BASE_DECL_HPP


#include "Teuchos_RCP.hpp"


namespace EvilPack {


/** \brief . */
using Teuchos::RCP;


/** \brief Evil base class that people often write with a factory function
 * to all of the subclasses in the interface.
 */
template<class T>
class EvilBase {
public:

  /** \brief. Required virtual destructor. */
  virtual ~EvilBase();

  /** \brief The virtual function. */
  virtual void soundOff(const T& obj) const = 0;

  /** \brief The factory in the interface. */
  static RCP<EvilBase<T> >
  createEvil(const std::string& concreteEvilName);

};


} // namespace EvilPack


#endif // EVIL_BASE_DECL_HPP
