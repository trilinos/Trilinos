// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_TPETRA_CG_HPP
#define STOKHOS_TPETRA_CG_HPP

namespace Stokhos {

// A simple CG solver for Tpetra-like objects
template <typename Matrix,
          typename Vector,
          typename Ordinal>
bool CG_Solve(const Matrix& A, Vector& x, const Vector& b,
              typename Vector::mag_type tol, Ordinal max_its,
              std::ostream* out = 0)
{
  typedef typename Vector::mag_type mag_type;
  typedef typename Vector::dot_type dot_type;

  Vector r(b.getMap());
  Vector p(x.getMap());
  Vector w(x.getMap());
  A.apply(x, r);
  r.update(1.0, b, -1.0);

  dot_type rho = r.dot(r);
  dot_type rho_old = rho;
  mag_type nrm = std::sqrt(rho);
  Ordinal k=0;
  dot_type alpha, beta, pAp;
  while (k < max_its && nrm > tol) {
    if (k == 0) {
      p.update(1.0, r, 0.0);
    }
    else {
      beta = rho / rho_old;
      p.update(1.0, r, beta);
    }
    A.apply(p, w);
    pAp = p.dot(w);
    alpha = rho/pAp;
    x.update( alpha, p, 1.0);
    r.update(-alpha, w, 1.0);

    rho_old = rho;
    rho = r.dot(r);
    nrm = std::sqrt(rho);
    ++k;

    if (out)
      *out << k << ":  " << nrm << std::endl;
  }

  if (nrm <= tol)
    return true;
  return false;
}

// A simple preconditioned CG solver based for Tpetra-like objects
template <typename Matrix,
          typename Vector,
          typename Prec,
          typename Ordinal>
bool PCG_Solve(const Matrix& A, Vector& x, const Vector& b, const Prec& M,
               typename Vector::mag_type tol, Ordinal max_its,
               std::ostream* out = 0)
{
  typedef typename Vector::mag_type mag_type;
  typedef typename Vector::dot_type dot_type;

  Vector r(b.getMap());
  Vector p(x.getMap());
  Vector w(x.getMap());
  A.apply(x, r);
  r.update(1.0, b, -1.0);

  mag_type nrm = r.norm2();
  dot_type rho = 1.0;
  Ordinal k=0;
  dot_type rho_old, alpha, beta, pAp;
  while (k < max_its && nrm > tol) {
    M.apply(r, w);
    rho_old = rho;
    rho = r.dot(w);

    if (k == 0) {
      p.update(1.0, w, 0.0);
    }
    else {
      beta = rho / rho_old;
      p.update(1.0, w, beta);
    }
    A.apply(p, w);
    pAp = p.dot(w);
    alpha = rho/pAp;
    x.update( alpha, p, 1.0);
    r.update(-alpha, w, 1.0);

    nrm = r.norm2();
    ++k;

    if (out)
      *out << k << ":  " << nrm << std::endl;
  }

  if (nrm <= tol)
    return true;
  return false;
}

}

#endif // STOKHOS_TPETRA_CG
