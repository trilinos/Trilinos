// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_PROVIDERBASE_H
#define PIRO_PROVIDERBASE_H

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace Piro {

/*! \brief Generic abstract base class for an auxiliary object factory
 *
 *  Implementing the ProviderBase interface is the most flexible approach to create concrete Provider objects.
 */
template <typename T>
class ProviderBase {
public:
  /*! \brief Returns an owning pointer to an object of a subclass of T
   *
   *  The returned value may refer to a new instance, to an already existing -- therefore possibly shared -- one,
   *  or simply be the null pointer.
   */
  virtual Teuchos::RCP<T> getInstance(const Teuchos::RCP<Teuchos::ParameterList> &params) = 0;

  //! \name Constructor and destructor
  //@{
  //! \brief Empty default constructor
  ProviderBase() {}

  //! \brief Virtual empty destructor
  virtual ~ProviderBase() {}
  //@}

private:
  // Disallow copy & assignment
  ProviderBase(const ProviderBase &);
  ProviderBase &operator=(const ProviderBase &);
};

} // namespace Piro

#endif /*PIRO_PROVIDERBASE_H*/
