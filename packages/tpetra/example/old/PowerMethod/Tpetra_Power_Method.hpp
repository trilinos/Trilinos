/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef TPETRA_POWER_METHOD_HPP
#define TPETRA_POWER_METHOD_HPP

#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_ScalarTraits.hpp>

namespace TpetraExamples {

  /** \brief Simple power iteration eigensolver for a Tpetra::Operator.
   */
  template <class Scalar, class Ordinal>
  Scalar powerMethod(const Teuchos::RCP<const Tpetra::Operator<Scalar,Ordinal> > &A, int niters, typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance, bool verbose)
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    typedef Tpetra::Vector<Scalar,Ordinal> Vector;

    if ( A->getRangeMap() != A->getDomainMap() ) {
      throw std::runtime_error("TpetraExamples::powerMethod(): operator must have domain and range maps that are equivalent.");
    }
    // create three vectors, fill z with random numbers
    Teuchos::RCP<Vector> z, q, r;
    q = Tpetra::createVector<Scalar>(A->getRangeMap());
    r = Tpetra::createVector<Scalar>(A->getRangeMap());
    z = Tpetra::createVector<Scalar>(A->getRangeMap());
    z->randomize();
    //
    Scalar lambda = 0.0;
    Magnitude normz, residual = 0.0;
    // power iteration
    for (int iter = 0; iter < niters; ++iter) {
      normz = z->norm2();                             // Compute 2-norm of z
      q->scale(1.0/normz, *z);                        // Set q = z / normz
      A->apply(*q, *z);                               // Compute z = A*q
      lambda = q->dot(*z);                            // Approximate maximum eigenvalue: lamba = dot(q,z)
      if ( iter % 100 == 0 || iter + 1 == niters ) {
        r->update(1.0, *z, -lambda, *q, 0.0);         // Compute A*q - lambda*q
        residual = Teuchos::ScalarTraits<Scalar>::magnitude(r->norm2() / lambda);
        if (verbose) {
          std::cout << "Iter = " << iter
                    << "  Lambda = " << lambda
                    << "  Residual of A*q - lambda*q = " << residual
                    << std::endl;
        }
      }
      if (residual < tolerance) {
        break;
      }
    }
    return lambda;
  }

} // end of namespace TpetraExamples

/** \example Tpetra_Power_Method_From_File.cpp
    \brief Power method example with reading a file.

    Read a Harwell-Boeing file into a Tpetra::CrsMatrix sparse matrix,
    and compute the matrix's leading eigenvalue using
    TpetraExamples::powerMethod().
  */

#endif
