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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(ROL_MiniTensor_Function_hpp)
#define ROL_MiniTensor_Function_hpp

#include "MiniTensor_Solvers.h"
#include "ROL_Objective.hpp"
#include "ROL_MiniTensor_Vector.hpp"

namespace ROL {

using Index = minitensor::Index;

///
/// Function base class that defines the interface to Mini Solvers.
///
template<typename MSFN, typename S, Index M>
class MiniTensor_Objective : public Objective<S>
{
public:

  MiniTensor_Objective(MSFN & msfn);

  MiniTensor_Objective() = delete;

  virtual
  ~MiniTensor_Objective();

  // ROL interface
  virtual
  S
  value(Vector<S> const & x, S & tol) final;

  virtual
  void
  gradient(Vector<S> & g, Vector<S> const & x, S & tol) final;

  virtual
  void
  hessVec(Vector<S> & hv, Vector<S> const & v,
      Vector<S> const & x, S & tol) final;

  virtual
  void
  invHessVec(Vector<S> & hv, Vector<S> const & v,
      Vector<S> const & x, S & tol) final;

private:
  MSFN
  minisolver_fn_;
};

} // namespace ROL

#include "ROL_MiniTensor_Function_Def.hpp"

#endif // ROL_MiniTensor_Function_hpp
