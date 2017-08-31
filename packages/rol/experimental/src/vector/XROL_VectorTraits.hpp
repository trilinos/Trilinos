
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

#pragma once

#include <tuple>

#include "XROL_Exception.hpp"

namespace XROL {

template<class> struct ElementTraits; 


/** @ingroup la_group
    \class XROL::VectorCreationTraits
    \brief Defines a traits-based approach for creating vectors.
*/

template<class V>
struct VectorCreationTraits {

  using IndexType = typename ElementTraits<V>::IndexType;

    /** \brief Clone to make a new (uninitialized) vector.
             @param[in]      The vector to clone.  
             @return         A smart pointer to the cloned vector.

             Provides the means of allocating temporary memory in ROL.
             ---             
  */
  static auto clone( const V& x ) { return nullptr; }

    /** \brief Return i-th basis vector.

             @param[in] A vector of from the desired vector space
             @param[in] i is the index of the basis function.
             @return A smart pointer to the basis vector with index @b i.

             Overloading the basis is only required if the default gradient 
             implementation is used, which computes a finite-difference 
             approximation.

             These methods must be specialized for concrete vectors
             ---
  */
  static auto basis( const V& x, IndexType i ) { return nullptr; }

}; // VectorCreationTraits



/** @ingroup la_group
    \class XROL::VectorSpaceTraits
    \brief Defines a traits-based approach for accessing vector space properties.
*/


template<class V>
struct VectorSpaceTraits {

  /** \brief Return dimension of the vector space.

             @return The dimension of the vector space, i.e., the total number 
                     of basis vectors.
             ---
  */
  using ElementType = typename ElementTraits<V>::ElementType;
  static auto dimension( const V& v ) { return ElementType(0); }
  
  /** \brief Compute the dual representation of a vector v, for example,
             the result of applying a Riesz map, change of basis, or change
             of memory layout. 

             @param[out] The dual repesentation of the input vector
             @param[in]  The vector for which to compute the dual
             ---
  */
  static void dual( V& vdual, const V& v ) { }

  template<class ...Vs>
  static void CheckDimensions( Vs... vecs ) {
    // TODO: Verify that Vs... consists of type V only
    UniformTypeCheck(values...);
    auto first = std::get<0>(std::make_tuple(vecs...)).dimension();
    auto cond = [first](auto i){ return i.dimension() == first; };
    bool match{true};
    bool _1[] = { ( match &= cond(vecs) )... }; (void)_1;

    if(!match) {
      std::string dims;
      bool _2[] = { ( dims += std::to_string( items.dimension() ) + " ")...}; (void)_2;
      throw IncompatibleDimensions(dims);
    }
  }


};


/** @ingroup la_group
    \class XROL::VectorFunctionTraits
    \brief Defines a traits-based approach for applying elementwise
           functions and reduce operations to vectors 

           These methods must be specialized for concrete vectors

*/


template<class V>
struct VectorFunctionTraits {

    /** Apply a callable function elementwise to a set of input vectors 
        \f$v^1,v^2,\hdots\f$ and write the result to x
  
        /f$ x \leftarrow f(v^1,v^2,\hdots) \f$
        
    */
    template<class Function, class ...Vs>
    static void transform(V& x, const Function& f, const Vs&... vs ) { }

    /** Apply a reduce operation to a vector and return the result 
        (e.g. max, min, sum,...) */
   
    template<class R>
    static ElementType reduce( const Elementwise::Reduce<R> &r, const V& x ) {
      return ElementType(0);
    }    
    

    /** Combines transform and reduce into a single function. 
        Applies a callable function elementwise to a set of input vectors (vs)
        and apply a reduction rule to obtain a scalar quantity

        
        
    */
    template<class Function, class R, class ...Vs>
    static auto transform_and_reduce( const Function &f, const Elementwise::Reduce<R> &r,  const Vs&... vs ) { 
      return ElementType(0);
    }





}; // VectorFunctionTraits


/** @ingroup la_group
    \class XROL::VectorOptionalTraits
    \brief These methods can be either be specialized directly 
           if there is a more efficient implementation for V or they
           have a default implementation based on VectorFunctionTraits
*/
template<class V>
struct VectorOptionalTraits {

  using ElementType   = typename ElementTraits<V>::ElementType;

  using VFT = VectorFunctionTraits<V>;

  /** \brief Compute \f$x \leftarrow x + y\f$, where \f$y = \mathtt{*this}\f$.

             @param[in/out]  x  is the vector to be added to 
             @param[in]      y  is the vector being added

             ---
  */
  static void plus( V& x, const V& y ) { 
    
  }

  /** \brief Compute \f$x \leftarrow \alpha x\f$ 

             @param[in/out]  x is the vector being rescaled
             @param[in]      alpha is the scaling parameter

             ---
  */
  static void scale( V& x, const ElementType alpha ) { }


  /** \brief Compute \f$ \langle x,y \rangle \f$ 

             @param[in]  x is a vector
             @param[in]  y is a vector
             @return     the inner product of x and y

             ---
  */
  static auto dot( const V &x, const V &y ) { }    

  /** \brief Returns \f$ \|x\| \f$
             @param[in]  x is a vector
             @return     The norm of x
  */

  static auto norm( const V &x ) { }

  /** \brief Compute \f$x \leftarrow \alpha x + y\f$ 

             @param[in/out]  The vector to be scaled
             @param[in]      The scaling parameter
             @param[in]      The vector being added

             ---
  */
  static void axpy( V &x, const ElementType alpha, const V &y ) { }


  /** \brief Set vector to zeros 
  */
  static void zero( V &x ) {}

}; // VectorOptionalTraits



template<class V>
struct VectorPrintTraits {
  /** \brief Print vector to stream 
  */
  static void print( const V& x, std::ostream &os ) {}

}; // VectorPrintTraits





} // namespace XROL

