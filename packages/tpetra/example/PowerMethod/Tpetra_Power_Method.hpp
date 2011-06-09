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
    Demonstrate reading a Harwell-Boeing file into a Tpetra::CrsMatrix and computing its leading eigenvalue using the TpetraExamples::powerMethod() function.
  */

/** \example Tpetra_Power_Method.cpp
    Demonstrate building a simple sparse matrix and computing its leading eigenvalue using the TpetraExamples::powerMethod() function.
  */

#endif
