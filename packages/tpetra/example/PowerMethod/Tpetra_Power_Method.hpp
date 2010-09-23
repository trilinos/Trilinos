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
    if ( A->getRangeMap() != A->getDomainMap() ) throw std::runtime_error("TpetraExamples::powerMethod(): operator must have domain and range maps that are equivalent.");
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    const bool NO_INITIALIZE_TO_ZERO = false;
    // create three vectors; do not bother initializing q to zero, as we will fill it with random below
    Tpetra::Vector<Scalar,Ordinal> z(A->getRangeMap(), NO_INITIALIZE_TO_ZERO),
                                   q(A->getRangeMap(), NO_INITIALIZE_TO_ZERO),
                                   r(A->getRangeMap(), NO_INITIALIZE_TO_ZERO);
    // Fill z with random numbers
    z.randomize();
    // Variables needed for iteration
    const Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();
    const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
    Scalar lambda = static_cast<Scalar>(0.0);
    Magnitude normz, residual = static_cast<Magnitude>(0.0);
    // power iteration
    for (int iter = 0; iter < niters; ++iter) {
      normz = z.norm2();                            // Compute 2-norm of z
      q.scale(ONE/normz, z);                        // Set q = z / normz
      A->apply(q, z);                               // Compute z = A*q
      lambda = q.dot(z);                            // Approximate maximum eigenvalue: lamba = dot(q,z)
      if ( iter % 100 == 0 || iter + 1 == niters ) {
        r.update(ONE, z, -lambda, q, ZERO);     // Compute A*q - lambda*q
        residual = Teuchos::ScalarTraits<Scalar>::magnitude(r.norm2() / lambda);
        if (verbose) {
          std::cout << "Iter = " << iter << "  Lambda = " << lambda 
                    << "  Residual of A*q - lambda*q = " 
                    << residual << std::endl;
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
