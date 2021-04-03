// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER


#ifndef ROL_SIMCONTROLLER_H
#define ROL_SIMCONTROLLER_H

#include "ROL_Vector.hpp"
#include "ROL_UpdateType.hpp"

namespace ROL {

template <class Real, class Key=std::vector<Real>>
class VectorController {
private:
  // Storage
  std::map<Key,int>               indices_, indices_trial_, indices_temp_;
  std::vector<bool>               flags_, flags_trial_, flags_temp_;
  std::vector<Ptr<Vector<Real>>>  vectors_, vectors_trial_, vectors_temp_;
  int maxIndex_, maxIndex_trial_, maxIndex_temp_;

  // Update flags
  bool trial_, temp_;
  bool objUpdated_, conUpdated_;

public:
  /** \brief Constructor.
  */
  VectorController(void);

  void reset(bool flag = true);

  /** \brief Objective function update for VectorController storage.
  */
  void objectiveUpdate(bool flag = true);

  /** \brief Equality constraint update for VectorController storage.
  */
  void constraintUpdate(bool flag = true);

  /** \brief Objective function update for VectorController storage.
  */
  void objectiveUpdate(UpdateType type);

  /** \brief Constraint update for VectorController storage.
  */
  void constraintUpdate(UpdateType type);

  /** \brief Check if vector associated with provided key is allocated.
  */
  bool isNull(const Key &param) const;

  /** \brief Check if vector has been computed.
  */
  bool isComputed(const Key &param) const;

  /** \brief Allocate the vector associated with provided key.
  */
  void allocate(const Vector<Real> &x, const Key &param);

  /** \brief Set the vector associated with provided key.  This assumes
      the vector data will be changed.
  */
  const Ptr<Vector<Real>> set(const Key &param);

  /** \brief Return the vector associated with provided key.
  */
  const Ptr<const Vector<Real>> get(const Key &param) const;

  /** \brief Return vector corresponding to input parameter.
  */
  bool get(Vector<Real> &x, const Key &param);

  /** \brief Set vector corresponding to input parameter.
  */
  void set(const Vector<Real> &x, const Key &param);

  /** \brief Push the contents of *this into another VectorController.
  */
  void push(VectorController<Real,Key> &to) const;

private:

  void resetTrial(void);

  void resetTemp(void);

  bool isNull(const Key &param, const std::map<Key,int> &indices) const;

  bool isComputed(const Key &param, const std::map<Key,int> &indices,
           const std::vector<bool> &flags) const;

  void allocate(const Vector<Real> &x, const Key &param,
           std::map<Key,int> &indices, std::vector<bool> &flags,
           std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const;

  const Ptr<const Vector<Real>> get(const Key &param,
           const std::map<Key,int> &indices, const std::vector<bool> &flags,
           const std::vector<Ptr<Vector<Real>>> &vectors, const int &maxIndex) const;

  const Ptr<Vector<Real>> set(const Key &param,
           std::map<Key,int> &indices, std::vector<bool> &flags,
           std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const;

  bool get(Vector<Real> &x, const Key &param,
           std::map<Key,int> &indices, std::vector<bool> &flags,
           std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const;

  void set(const Vector<Real> &x, const Key &param,
           std::map<Key,int> &indices, std::vector<bool> &flags,
           std::vector<Ptr<Vector<Real>>> &vectors, int &maxIndex) const;

  void accept(void);
}; // class VectorController

} // namespace ROL

#include "ROL_VectorController_Def.hpp"

#endif
