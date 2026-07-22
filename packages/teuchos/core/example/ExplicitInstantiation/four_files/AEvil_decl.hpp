// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef A_EVIL_DECL_HPP
#define A_EVIL_DECL_HPP

// Only include the declaration, not any implementations in case of cicular
// dependencies!
#include "EvilBase_decl.hpp"


namespace EvilPack {


// Need a forward for B to declare function callBEvil(...)
template<class T> class BEvil;


/** \brief A subclass of EvilBase that calls BEvil.
 */
template<class T>
class AEvil : public EvilBase<T> {
public:
  /** \brief . */
  void callBEvil(const BEvil<T> &bEvil, const T& obj) const;
  /** \brief . */
  void soundOff(const T& obj) const;
};


/** \brief Nonmember constructor.
 *
 * \relates AEvil
 */
template<class T>
inline
RCP<AEvil<T> > aEvil()
{
  return Teuchos::rcp(new AEvil<T>);
}


} // namespace EvilPack


#endif // A_EVIL_DECL_HPP
