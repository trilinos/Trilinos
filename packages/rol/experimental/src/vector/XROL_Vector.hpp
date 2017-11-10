
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

namespace XROL {

/** \file  XROL_Vector.hpp
    \brief Defines the interfaces for the functions that operate on 
           vectors 
 */

// Forward declare operator
template<class V> struct Operator;

/** \fn         clone
    \brief      Create a new vector from the same vector space as x 
    @param[in]  x  The vector to clone
    \return     a pointer to a new vector
*/
template<class V> 
std::unique_ptr<V> clone( const V& x );// { return nullptr; }

/** \fn         basis 
    \brief      Create the ith canonical vector from the same vector space as x 
     
       \f[ x = \begin{pmatrix} 0 & \hdots & 0 & 1 & 0 & \hdots & 0 \end{pmatrix},
                \quad x_j = \delta_{ij}  \f]

    @param[in,out]  x  Vector from which to create a basis vector
    @param[in]      i  the index of the unit element   
*/ 
template<class V> 
void basis( V& x, index_t<V> i );// { return nullptr; }


/** \fn          dual
    \brief       Given a primal vector \f$x\in X\f$, compute the dual vector \f$x^\ast \in X^\ast\f$
    @param[out]  xdual  The resulting dual vector in \f$X^\ast\f$
    @param[in]   xprim  The given primal vector in \f$X\f$
    
*/  

template<class V> void dual( dual_t<V>& xdual, const V& xprim ); // { } 


/** \fn        dimension 
    \brief     Return the dimensionality of the vector x
    @param[in] x  The given vector
    \return       The dimensionality of the vector x 
*/
template<class V> 
index_t<V> dimension( const V& x ); // { return index_t<V>(0); }


/** \fn             plus 
    \brief          Add vector y to vector x,   \f$ x \leftarrow x + y \f$
    @param[in,out]  x  The vector being added to 
    @param[in]      y  The vector being added
*/
template<class V> void plus( V& x, const V& y ); // {}


/** \fn    scale
    \brief Multiply vector x by alpha,  \f$ x leftarrow \alpha x \f$
*/
template<class V> void scale( V& x, const element_t<V> alpha ); // {}

/** \fn    fill
    \brief Set every element to alpha,  \f$ x leftarrow \alpha e \f$
*/
template<class V> void fill( V& x, const element_t<V> alpha ); // {}

/** \fn    set
    \brief Set vector x equal to y,     \f$ x \leftarrow y
*/
template<class V> void set( V& x, const V& y ); // {}


/** \fn    axpy 
    \brief Scale a vector x and add y to it \f$ x \leftarrow \alpha x + y \f$
*/
template<class V> void axpy( V& y, const element_t<V> alpha, const V& x ); //  {}


/** \fn    zero
    \brief Zero out the elements of x
*/
// template<class V> void zero( V& x ); // {}


/** \fn    dot 
    \brief Compute the inner product of vectors x and y
*/
template<class V> 
element_t<V> dot( const V& x, const V &y ); // { return element_t<V>(0); }


/** \fn     norm 
    \brief Return the norm of the vector x
*/
template<class V> 
magnitude_t<V> norm( const V& x ); // { return magnitude_t<V>(0); }

   
/** \fn     print 
    \brief  Print a vector x to a stream 
*/
template<class V> void print( const V& x, std::ostream& os ); // { }

/** \fn          eval_function
    \brief       Apply a function elementwise to a pack of vectors and write the 
                 result to another vector
    @param[out]  x     The result vector
    @param[in]   f     A function, functor, or lambda taking arbitrarily many scalar
                       arguments of the same type that returns a scalar of the same 
                       type
    @param[in]   vecs  A variadic pack of input vectors of the same type
*/
template<class V, class F, class... Vs> 
void eval_function( V& x, const F& f, const Vs&... vs ); // { }


/** \fn         reduce
    \brief      Apply a reduction operation to the vector
    @param[in]  x  Vector to reduce
    @param[in]  r  Pairwise reduction function r(current,prev). r() should 
                return an initial value, e.g. 0 for sum, 1 for product, 
                numeric_limits::max for min and numeric_limits::lowest for max
    \return     result value
*/
template<class V, class R>
auto reduce( const V& x, const R& r ); // { return r(); }


/** \fn        eval_function_and_reduce
    \brief     Apply a function to a variadic vector pack and reduce the result
    
    @param[in]  r  Pairwise reduction function r(current,prev). r() should 
                return an initial value, e.g. 0 for sum, 1 for product, 
                numeric_limits::max for min and numeric_limits::lowest for max
    @param[in]  f  A function, functor, or lambda taking arbitrarily many scalar
                arguments of the same type that returns a scalar of the same 
                type
    \return     result value
*/
template<class R, class F, class V, class... Vs>
auto eval_function_and_reduce( const R& r, const F& f, const V, const Vs&... vs ); 
// { return r(); }

/** \fn          randomize 
    \brief       Assign pseudorandom values to the elements of a vector given a
                 number generator and probability distribution function 
    @param[out]  x    The vector randomize
    @param[in]   gen  The engine that generates random numbers 
    @param[in]   dist The probability distribution function (functor)
*/
template<class Generator, class Distribution, class V>
void randomize( Generator& g, Distribution& d, V &v ); // {}

/*
template<class V> 
element_t<V> norm( const V& x, const V& y ) {
  auto sum = make_sum(x);
  auto result = eval_function_and_reduce(sum,[](auto a, auto b){ return a*b; },x,y);
  return std::sqrt(result);
} 

template<class V> 
magnitude_t<V> norm( const V& x ) {
  auto sum = make_sum(x);
  auto result = eval_function_and_reduce(sum,[](auto v){ return v*v; },x);
  return std::sqrt(result);
} 
*/



} // namespace XROL
 
