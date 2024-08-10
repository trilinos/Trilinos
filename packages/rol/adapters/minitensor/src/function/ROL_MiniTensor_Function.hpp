// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
