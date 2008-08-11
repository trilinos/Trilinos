#ifndef GENSQP_VECTOR_H
#define GENSQP_VECTOR_H

#include "Teuchos_RefCountPtr.hpp"

/** \class GenSQP::Vector
    \brief Provides the interface to generic abstract vector libraries.

    The interfaced functionality is very basic and includes routines for:\n
    \li linear combinations,
    \li inner products,
    \li scaling and copying operations,
    \li the cloning of a vector.
*/


namespace GenSQP {

class Vector {
public:

  virtual ~Vector() {}

  /** \brief Returns inner(*this,x).
  */
  virtual double innerProd( const Vector &x ) const = 0;

  /** \brief <tt>y = alpha*x + beta*y</tt> where <tt>y == *this</tt>.
  */
  virtual void linComb( const double &alpha, const Vector &x, const double &beta = 1.0 ) = 0;

  /** \brief <tt>y = alpha*y</tt> where <tt>y == *this</tt>.
  */
  virtual void Scale( const double &alpha ) = 0;

  /** \brief <tt>y = alpha</tt> where <tt>y == *this</tt>.
  */
  virtual void Set( const double &alpha ) = 0;

  /** \brief <tt>y = alpha*x</tt> where <tt>y == *this</tt>.
  */
  virtual void Set( const double &alpha, const Vector &x ) = 0;

  /** Clone to make a new (uninitialized) vector.
  */
  virtual Teuchos::RefCountPtr<Vector> createVector() const = 0;

}; // class Vector

} // namespace GenSQP

#endif
