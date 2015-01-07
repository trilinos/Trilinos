
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_VECTOR_H
#define ROL_VECTOR_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_oblackholestream.hpp"

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
class Vector {
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
  virtual Teuchos::RCP<Vector> clone() const = 0;


  /** \brief Compute \f$y \leftarrow \alpha x + y\f$ where \f$y = \mathtt{*this}\f$.

             @param[in]      alpha is the scaling of @b x.
             @param[in]      x     is a vector.

             On return \f$\mathtt{*this} = \mathtt{*this} + \alpha x \f$.
             Uses #clone, #set, #scale and #plus for the computation.
             Please overload if a more efficient implementation is needed.

             ---
  */
  virtual void axpy( const Real alpha, const Vector &x ) {
    Teuchos::RCP<Vector> ax = x.clone();
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
  virtual Teuchos::RCP<Vector> basis( const int i ) const {return Teuchos::null;}


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

             The consistency errors are defined as the norms or absolute values of the differences between the left-hand
             side and the right-hand side terms in the above equalities.

             ---
  */
  virtual std::vector<Real> checkVector( const Vector<Real> &x,
                                         const Vector<Real> &y,
                                         const bool printToStream = true,
                                         std::ostream & outStream = std::cout ) const {
    Real one  =  1.0;
    Real a    =  1.234;
    Real b    = -432.1;
    int width =  94;
    std::vector<Real> vCheck;

    Teuchos::oblackholestream bhs; // outputs nothing

    Teuchos::RCP<std::ostream> pStream;
    if (printToStream) {
      pStream = Teuchos::rcp(&outStream, false);
    } else {
      pStream = Teuchos::rcp(&bhs, false);
    }

    std::ios::fmtflags f( outStream.flags() );

    Teuchos::RCP<Vector> v    = this->clone();
    Teuchos::RCP<Vector> vtmp = this->clone();
    Teuchos::RCP<Vector> xtmp = x.clone();
    Teuchos::RCP<Vector> ytmp = y.clone();

    *pStream << "\n************ Begin verification of linear algebra.\n\n";

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
    xtmp->plus(y); vCheck.push_back(std::abs(this->dot(*xtmp) - x.dot(*this) - y.dot(*this)));
    *pStream << std::setw(width) << std::left << "Additivity of dot (inner) product. Consistency error: " << " " << vCheck.back() << "\n";

    // Consistency of scalar multiplication and norm.
    v->set(*this);
    v->scale(one/v->norm()); vCheck.push_back(std::abs(v->norm() - one));
    *pStream << std::setw(width) << std::left << "Consistency of scalar multiplication and norm. Consistency error: " << " " << vCheck.back() << "\n";

    // Reflexivity.
    v->set(*this);
    xtmp = Teuchos::rcp_const_cast<Vector>(Teuchos::rcpFromRef(this->dual()));
    ytmp = Teuchos::rcp_const_cast<Vector>(Teuchos::rcpFromRef(xtmp->dual()));
    v->axpy(-one, *ytmp); vCheck.push_back(v->norm());
    *pStream << std::setw(width) << std::left << "Reflexivity. Consistency error: " << " " << vCheck.back() << "\n\n";

    *pStream << "************   End verification of linear algebra.\n\n";

    outStream.flags( f );

    return vCheck;
  }

}; // class Vector

} // namespace ROL

#endif
