//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_VECTOR_H
#define ROL_VECTOR_H

#include "Teuchos_RefCountPtr.hpp"

/** \class ROL::Vector
    \brief Provides the vector space interface.

    The basic interface to be supplied by the user includes:\n
    \li vector addition,
    \li scalar multiplication,
    \li dot (scalar) product of vectors,
    \li vector norm,
    \li cloning of vectors.

    The dot product can represent an inner product (in Hilbert space) or
    a duality pairing (in general Banach space).

    There are additional virtual member functions that the user
    may want to reimplement for added efficiency. 
*/


namespace ROL {

template <class Real>
class Vector {
public:

  virtual ~Vector() {}

  /** \brief Compute \f$y \leftarrow x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void plus( const Vector &x ) = 0;

  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void scale( const Real alpha ) = 0;

  /** \brief Returns \f$ \langle y,x \rangle \f$ where \f$y = \mbox{*this}\f$.
  */
  virtual Real dot( const Vector &x ) const = 0;

  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mbox{*this}\f$.
  */
  virtual Real norm() const = 0;

  /** \brief Clone to make a new (uninitialized) vector.
  */
  virtual Teuchos::RCP<Vector> clone() const = 0;


  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void axpy( const Real alpha, const Vector &x ) {
    Teuchos::RCP<Vector> ax = x.clone();
    ax->set(x);
    ax->scale(alpha);
    this->plus(*ax);
  }

  /**  \brief Set to zero vector.
  */
  virtual void zero() {
    this->scale( (Real)0 );
  }

  /** \brief Return i-th basis vector: define if finite-difference gradients
             and Hessians are used.
  */
  virtual Teuchos::RCP<Vector> basis( const int i ) const {return Teuchos::null;}

  /** \brief Return dimension of the vector space.
  */
  virtual int dimension() {return 0;}

  /**  \brief Set \f$y \leftarrow x\f$ where \f$y = \mbox{*this}\f$.
  */
  virtual void set( const Vector &x ) {
    this->zero();
    this->plus(x);
  }

}; // class Vector

} // namespace ROL

#endif
