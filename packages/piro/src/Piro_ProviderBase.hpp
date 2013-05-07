// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
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
