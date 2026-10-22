// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_ROL_AKRI_ROLDISTRIBUTEDVECTOR_HPP_
#define KRINO_KRINO_ROL_AKRI_ROLDISTRIBUTEDVECTOR_HPP_

#include <Akri_DistributedVector.hpp>
#include <ROL_Vector.hpp>
#include <cmath>

namespace krino {

// ROL::Vector implementation that wraps a krino::DistributedVector so that the
// ROL optimizer computes parallel-consistent inner products and norms. Unlike
// ROL::StdVector (which sums over all local entries with no MPI reduction, thus
// double-counting shared/ghosted DOFs), every reduction here routes through the
// krino::Dot free function, which sums owned entries only and then does an
// all_reduce across ranks.
class ROLDistributedVector : public ROL::Vector<double>
{
public:
  explicit ROLDistributedVector(const DistributedVector & v) : myVec(v) {}
  explicit ROLDistributedVector(DistributedVector && v) : myVec(std::move(v)) {}

  DistributedVector & get() { return myVec; }
  const DistributedVector & get() const { return myVec; }

  void plus(const ROL::Vector<double> & x) override { myVec = xpby(myVec, 1.0, vec(x)); }

  void scale(const double alpha) override { myVec = scalar_times_vector(alpha, myVec); }

  void axpy(const double alpha, const ROL::Vector<double> & x) override { myVec = xpby(myVec, alpha, vec(x)); }

  // Parallel-consistent inner product: owned entries only, then all_reduce.
  double dot(const ROL::Vector<double> & x) const override { return krino::Dot(myVec, vec(x)); }

  double norm() const override { return std::sqrt(this->dot(*this)); }

  void set(const ROL::Vector<double> & x) override { myVec = vec(x); }

  void zero() override { for (auto & entry : myVec) entry = 0.0; }

  ROL::Ptr<ROL::Vector<double>> clone() const override
  {
    DistributedVector z = myVec;
    z.assign(myVec.comm(), myVec.size(), myVec.local_size(), 0.0); // same shape and comm, zeroed
    return ROL::makePtr<ROLDistributedVector>(std::move(z));
  }

  // Not called on the configured L-BFGS/backtracking path; returns the local
  // size only so it never triggers a collective. Do not use for global sizing.
  int dimension() const override { return static_cast<int>(myVec.local_size()); }

private:
  // clone() only ever produces ROLDistributedVector, so this downcast is safe.
  static const DistributedVector & vec(const ROL::Vector<double> & x)
  {
    return static_cast<const ROLDistributedVector &>(x).myVec;
  }

  DistributedVector myVec;
};

}

#endif /* KRINO_KRINO_ROL_AKRI_ROLDISTRIBUTEDVECTOR_HPP_ */
