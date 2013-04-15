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

#ifndef PIRO_OBSERVERBASE_HPP
#define PIRO_OBSERVERBASE_HPP

#include "Thyra_VectorBase.hpp"

namespace Piro {

template <typename Scalar>
class ObserverBase {
public:
  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution);

  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Scalar &stamp);

  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Thyra::VectorBase<Scalar> &solution_dot,
      const Scalar &stamp);

  virtual ~ObserverBase() {}
};

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Scalar &/*stamp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Thyra::VectorBase<Scalar> &/*solution_dot*/,
    const Scalar &/*stamp*/)
{
  // Nothing to do by default
}

} // namespace Piro

#endif /* PIRO_OBSERVERBASE_HPP */
