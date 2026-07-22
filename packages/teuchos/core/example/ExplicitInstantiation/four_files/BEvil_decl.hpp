// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef B_EVIL_DECL_HPP
#define B_EVIL_DECL_HPP

// Only include the declaration, not any implementations in case of cicular
// dependencies!
#include "EvilBase_decl.hpp"


namespace EvilPack {


// Need a forward for B to declare function callAEvil(...)
template<class T> class AEvil;


/** \brief A subclass of EvilBase that calls AEvil.
 */
template<class T>
class BEvil : public EvilBase<T> {
public:
  /** \brief . */
  void callAEvil(const AEvil<T> &aEvil, const T& obj) const;
  /** \brief . */
  void soundOff(const T& obj) const;
};


/** \brief Nonmember constructor.
 *
 * \relates BEvil
 **/
template<class T>
inline
RCP<BEvil<T> > bEvil()
{
  return Teuchos::rcp(new BEvil<T>);
}


} // namespace EvilPack


#endif // B_EVIL_DECL_HPP
