
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_VECTOR_H
#define ROL_VECTOR_H

#define ROL_UNUSED(x) (void) x

#include <ostream>
#include <vector>
#include <algorithm>

#include "ROL_Elementwise_Function.hpp"

#include "ROL_Ptr.hpp"
#include "ROL_Stream.hpp"

/** @ingroup la_group
    \class ROL::Vector
    \brief Defines the linear algebra or vector space interface.

    The basic linear algebra interface, to be implemented by the user, includes:\n
    \li #plus -- vector addition;
    \li #scale -- scalar multiplication;
    \li #dot -- dot (scalar) product of vectors;
    \li #norm -- vector norm;
    \li #clone -- cloning of existing vectors.

    The dot product can represent an inner product (in Hilbert space) or
    a duality pairing (in general Banach space).

    There are additional virtual member functions that can be
    overloaded for computational efficiency. 
*/

namespace ROL {

template <class Real>
class Vector
#ifdef ENABLE_PYROL
 : public std::enable_shared_from_this<Vector<Real>> 
#endif 
{
public:

  virtual ~Vector() {}


  /** \brief Compute \f$y \leftarrow y + x\f$, where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector to be added to \f$\mathtt{*this}\f$.

             On return \f$\mathtt{*this} = \mathtt{*this} + x\f$.

             ---
  */
  virtual void plus( const Vector &x ) = 0;


  /** \brief Compute \f$y \leftarrow \alpha y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of \f$\mathtt{*this}\f$.

             On return \f$\mathtt{*this} = \alpha (\mathtt{*this}) \f$.

             ---
  */
  virtual void scale( const Real alpha ) = 0;


  /** \brief Compute \f$ \langle y,x \rangle \f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x  is the vector that forms the dot product with \f$\mathtt{*this}\f$.
             @return         The number equal to \f$\langle \mathtt{*this}, x \rangle\f$.

             ---
  */
  virtual Real dot( const Vector &x ) const = 0;


  /** \brief Returns \f$ \| y \| \f$ where \f$y = \mathtt{*this}\f$.

             @return         A nonnegative number equal to the norm of \f$\mathtt{*this}\f$.

             ---
  */
  virtual Real norm() const = 0;


  /** \brief Clone to make a new (uninitialized) vector.

             @return         A reference-counted pointer to the cloned vector.

             Provides the means of allocating temporary memory in ROL.

             ---             
  */
  virtual ROL::Ptr<Vector> clone() const = 0;


  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of @b x.
             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = \mathtt{*this} + \alpha x \f$.
             Uses #clone, #set, #scale and #plus for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void axpy( const Real alpha, const Vector &x ) {
    ROL::Ptr<Vector> ax = x.clone();
    ax->set(x);
    ax->scale(alpha);
    this->plus(*ax);
  }

  /** \brief Set to zero vector.

             Uses #scale by zero for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void zero() {
    this->scale( (Real)0 );
  }


  /** \brief Return i-th basis vector.

             @param[in] i is the index of the basis function.
             @return A reference-counted pointer to the basis vector with index @b i.

             Overloading the basis is only required if the default gradient implementation
             is used, which computes a finite-difference approximation.

             ---
  */
  virtual ROL::Ptr<Vector> basis( const int i ) const {
    ROL_UNUSED(i);
    return ROL::nullPtr;
  }


  /** \brief Return dimension of the vector space.

             @return The dimension of the vector space, i.e., the total number of basis vectors.

             Overload if the basis is overloaded.

             ---
  */
  virtual int dimension() const {return 0;}


  /** \brief Set \f$y \leftarrow x\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = x\f$.
             Uses #zero and #plus methods for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void set( const Vector &x ) {
    this->zero();
    this->plus(x);
  }


  /** \brief Return dual representation of \f$\mathtt{*this}\f$, for example,
             the result of applying a Riesz map, or change of basis, or
             change of memory layout.

             @return         A const reference to dual representation.

             By default, returns the current object.
             Please overload if you need a dual representation.

             ---
  */
  virtual const Vector & dual() const {
    return *this;
  }

  /** \brief Apply \f$\mathtt{*this}\f$ to a dual vector.  This is equivalent
             to the call \f$\mathtt{this->dot(x.dual())}\f$.

             @param[in]     x      is a vector
             @return         The number equal to \f$\langle \mathtt{*this}, x \rangle\f$.

             ---
  */
  virtual Real apply(const Vector<Real> &x) const {
    return this->dot(x.dual());
  }

  virtual void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    ROL_UNUSED(f);
    ROL_TEST_FOR_EXCEPTION( true, std::logic_error,
      "The method applyUnary was called, but not implemented" << std::endl);
  }

  virtual void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector &x ) {
    ROL_UNUSED(f);
    ROL_UNUSED(x);
    ROL_TEST_FOR_EXCEPTION( true, std::logic_error,
      "The method applyBinary was called, but not implemented" << std::endl);
  }

  virtual Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    ROL_UNUSED(r);
    ROL_TEST_FOR_EXCEPTION( true, std::logic_error,
      "The method reduce was called, but not implemented" << std::endl); 
  }

  virtual void print( std::ostream &outStream ) const {
    outStream << "The method print was called, but not implemented" << std::endl;
  }

  /** \brief Set \f$y \leftarrow C\f$ where \f$C\in\mathbb{R}\f$.

             @param[in]      C     is a scalar.

             On return \f$\mathtt{*this} = C\f$.
             Uses #applyUnary methods for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void setScalar( const Real C ) {
    this->applyUnary(Elementwise::Fill<Real>(C));
  }

  /** \brief Set vector to be uniform random between [l,u].

             @param[in]      l     is a the lower bound.
             @param[in]      u     is a the upper bound.

             On return the components of \f$\mathtt{*this}\f$ are uniform
             random numbers on the interval \f$[l,u]\f$.
       	     The default implementation uses #applyUnary methods for the
       	     computation. Please overload if a more efficient implementation is
             needed.

             ---
  */
  virtual void randomize( const Real l = 0.0, const Real u = 1.0 ) {
    Elementwise::UniformlyRandom<Real> ur(l,u);
    this->applyUnary(ur);
  }

  /** \brief Verify vector-space methods.

             @param[in]      x     is a vector.
             @param[in]      y     is a vector.

             Returns a vector of Reals, all of which should be close to zero.
             They represent consistency errors in the vector space properties,
             as follows:

             - Commutativity of addition: \f$\mathtt{*this} + x = x + \mathtt{*this}\f$.
             - Associativity of addition: \f$\mathtt{*this} + (x + y) = (\mathtt{*this} + x) + y\f$.
             - Identity element of addition: \f$0 + x = x\f$.
             - Inverse elements of addition: \f$\mathtt{*this} + (- \mathtt{*this}) = 0\f$.
             - Identity element of scalar multiplication: \f$ 1 \cdot \mathtt{*this} = \mathtt{*this} \f$.
             - Consistency of scalar multiplication with field multiplication: \f$a (b \cdot \mathtt{*this}) = (a b) \cdot \mathtt{*this}\f$.
             - Distributivity of scalar multiplication with respect to field addition: \f$(a+b) \cdot \mathtt{*this} = a \cdot \mathtt{*this} + b \cdot \mathtt{*this}\f$.
             - Distributivity of scalar multiplication with respect to vector addition: \f$a \cdot (\mathtt{*this} + x) = a \cdot \mathtt{*this} + a \cdot x\f$.
             - Commutativity of dot (inner) product over the field of reals: \f$(\mathtt{*this}, x)  = (x, \mathtt{*this})\f$.
             - Additivity of dot (inner) product: \f$(\mathtt{*this}, x+y)  = (\mathtt{*this}, x) + (\mathtt{*this}, y)\f$.
             - Consistency of scalar multiplication and norm: \f$\|{\mathtt{*this}} / {\|\mathtt{*this}\|} \| = 1\f$.
             - Reflexivity: \f$\mbox{dual}(\mbox{dual}(\mathtt{*this})) = \mathtt{*this}\f$ .
             - Consistency of apply and dual: \f$\langle \mathtt{*this}, x\rangle = (\mbox{dual}(\mathtt{*this}),x) = (\mathtt{*this},\mbox{dual}(x))\f$.

             The consistency errors are defined as the norms or absolute values of the differences between the left-hand
             side and the right-hand side terms in the above equalities.

             ---
  */
  virtual std::vector<Real> checkVector( const Vector<Real> &x,
                                         const Vector<Real> &y,
                                         const bool printToStream = true,
                                         std::ostream & outStream = std::cout ) const {
    Real zero =  0.0;
    Real one  =  1.0;
    Real a    =  1.234;
    Real b    = -0.4321;
    int width =  94;
    std::vector<Real> vCheck;

    ROL::nullstream bhs; // outputs nothing

    ROL::Ptr<std::ostream> pStream;
    if (printToStream) {
      pStream = ROL::makePtrFromRef(outStream);
    } else {
      pStream = ROL::makePtrFromRef(bhs);
    }

    // Save the format state of the original pStream.
    ROL::nullstream oldFormatState, headerFormatState;
    oldFormatState.copyfmt(*pStream);

    ROL::Ptr<Vector> v    = this->clone();
    ROL::Ptr<Vector> vtmp = this->clone();
    ROL::Ptr<Vector> xtmp = x.clone();
    ROL::Ptr<Vector> ytmp = y.clone();

    *pStream << "\n" << std::setw(width) << std::left << std::setfill('*') << "********** Begin verification of linear algebra. " << "\n\n";
    headerFormatState.copyfmt(*pStream);

    // Commutativity of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->plus(x); xtmp->plus(*this); v->axpy(-one, *xtmp); vCheck.push_back(v->norm());
    *pStream << std::scientific << std::setprecision(12) << std::setfill('>');
    *pStream << std::setw(width) << std::left << "Commutativity of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Associativity of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    ytmp->plus(x); v->plus(*ytmp); xtmp->plus(*this); xtmp->plus(y); v->axpy(-one, *xtmp); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Associativity of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Identity element of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->zero(); v->plus(x); v->axpy(-one, x); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Identity element of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Inverse elements of addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->scale(-one); v->plus(*this); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Inverse elements of addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Identity element of scalar multiplication.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->scale(one); v->axpy(-one, *this); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Identity element of scalar multiplication. Consistency error: " << " " << vCheck.back() << "\n";

    // Consistency of scalar multiplication with field multiplication.
    v->set(*this); vtmp->set(*this);
    v->scale(b); v->scale(a); vtmp->scale(a*b); v->axpy(-one, *vtmp); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Consistency of scalar multiplication with field multiplication. Consistency error: " << " " << vCheck.back() << "\n";

    // Distributivity of scalar multiplication with respect to field addition.
    v->set(*this); vtmp->set(*this);
    v->scale(a+b); vtmp->scale(a); vtmp->axpy(b, *this); v->axpy(-one, *vtmp); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Distributivity of scalar multiplication with respect to field addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Distributivity of scalar multiplication with respect to vector addition.
    v->set(*this); xtmp->set(x); ytmp->set(y);
    v->plus(x); v->scale(a); xtmp->scale(a); xtmp->axpy(a, *this); v->axpy(-one, *xtmp); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Distributivity of scalar multiplication with respect to vector addition. Consistency error: " << " " << vCheck.back() << "\n";

    // Commutativity of dot (inner) product over the field of reals.
    vCheck.push_back(std::abs(this->dot(x) - x.dot(*this)));
    *pStream << std::setw(width) << std::left << "Commutativity of dot (inner) product over the field of reals. Consistency error: " << " " << vCheck.back() << "\n";

    // Additivity of dot (inner) product.
    xtmp->set(x);
    xtmp->plus(y); vCheck.push_back(std::abs(this->dot(*xtmp) - this->dot(x) - this->dot(y))/std::max({static_cast<Real>(std::abs(this->dot(*xtmp))), static_cast<Real>(std::abs(this->dot(x))), static_cast<Real>(std::abs(this->dot(y))), one}));
    *pStream << std::setw(width) << std::left << "Additivity of dot (inner) product. Consistency error: " << " " << vCheck.back() << "\n";

    // Consistency of scalar multiplication and norm.
    v->set(*this);
    Real vnorm = v->norm();
    if (vnorm == zero) {
      v->scale(a);
      vCheck.push_back(std::abs(v->norm() - zero));
    } else {
      v->scale(one/vnorm);
      vCheck.push_back(std::abs(v->norm() - one));
    }
    *pStream << std::setw(width) << std::left << "Consistency of scalar multiplication and norm. Consistency error: " << " " << vCheck.back() << "\n";

    // Reflexivity.
    v->set(*this);
    xtmp = ROL::constPtrCast<Vector>(ROL::makePtrFromRef(this->dual()));
    ytmp = ROL::constPtrCast<Vector>(ROL::makePtrFromRef(xtmp->dual()));
    v->axpy(-one, *ytmp); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Reflexivity. Consistency error: " << " " << vCheck.back() << "\n";

    // Consistency of apply and dual.
    v->set(*this);
    xtmp = x.dual().clone(); xtmp->set(x.dual());
    Real vx  = v->apply(*xtmp);
    Real vxd = v->dot(xtmp->dual());
    Real vdx = xtmp->dot(v->dual());
    if (vx == zero) {
      vCheck.push_back(std::max(std::abs(vx-vxd),std::abs(vx-vdx)));
    }
    else {
      vCheck.push_back(std::max(std::abs(vx-vxd),std::abs(vx-vdx))/std::abs(vx));
    }
    *pStream << std::setw(width) << std::left << "Consistency of apply and dual:" << " " << vCheck.back() << "\n\n";

    //*pStream << "************   End verification of linear algebra.\n\n";

    // Restore format state of pStream used for the header info.
    pStream->copyfmt(headerFormatState);
    *pStream << std::setw(width) << std::left << "********** End verification of linear algebra. " << "\n\n";

    // Restore format state of the original pStream.
    pStream->copyfmt(oldFormatState);

    return vCheck;
  }

}; // class Vector

} // namespace ROL

#endif
