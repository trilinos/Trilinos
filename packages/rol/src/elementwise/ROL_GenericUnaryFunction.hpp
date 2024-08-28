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
#include <utility>

namespace ROL {

/**
 * @brief A template class representing a unary function that can be applied to Vectors.
 * 
 * This class provides a mechanism for applying unary functions to Vectors efficiently,
 * using the Visitor pattern to avoid virtual function calls for each element.
 * 
 * @tparam Real The type used for real numbers in the function and calculations.
 */
template<class Real>
class GenericUnaryFunction {
private:
  template<class> class DerivedVisitor; ///< Forward declaration

protected:
  /**
   * @brief Abstract base class for visitors.
   * 
   * This class is part of the Visitor pattern implementation, allowing for
   * double dispatch without virtual function calls for each vector element.
   */
  struct Visitor { 
    virtual void visit( const GenericUnaryFunction& ) = 0; 
  };

public:
  /** @brief Virtual destructor for proper cleanup of derived classes. */
  virtual ~GenericUnaryFunction() noexcept = default;

  /**
   * @brief Apply the unary function to a single value.
   * 
   * @param x The input value.
   * @return The result of applying the function to x.
   */
  [[nodiscard]] virtual Real operator()( Real x ) const = 0;

  /**
   * @brief Accept a visitor, part of the Visitor pattern implementation.
   * 
   * @param visitor The visitor to accept.
   */
  virtual void accept( Visitor&& visitor ) const = 0;

  /**
   * @brief Apply the unary function to a vector.
   * 
   * This method uses the Visitor pattern to efficiently apply the function
   * to all elements of the vector without virtual function calls per element.
   * 
   * @tparam VecT The type of the vector.
   * @tparam EvalT The type of the evaluator function. 
   * @param vec The vector to apply the function to.
   * @param eval The evaluator function that defines how to apply the unary function to the vector.
   * 
   * Example Usage:
   * @code
   *   template<class Real>
   *   class Vector {   
   *   public:
   *     virtual void applyGenericUnary( const GenericUnaryFunction<Real>& ) {} 
   *   };
   *
   *   template<class Real>
   *   class StdVector : public Vector<Real> {
   *   public:
   *     void applyGenericUnary( const GenericUnaryFunction<Real>& guf ) override {
   *       guf.apply_vectorized(*this,[](StdVector& vec, const auto& f){ vec.applyGenericUnaryImpl(f); });
   *     }
   *
   *     template<class unary_function>
   *     void applyGenericUnaryImpl( const unary_function& f ) { 
   *       for( auto& e : vec_ ) e = f(e);
   *     }
   *   private:
   *     std::vector<Real> vec_;
   *   };
   * @endcode
   */
  template<class VecT, class EvalT> 
  void apply_vectorized( VecT& vec, EvalT&& eval ) const {
    accept(VectorVisitor(vec,std::forward<EvalT>(eval)));
  }

  /**
   * @brief A wrapper class that turns any callable into a GenericUnaryFunction.
   * 
   * This class allows easy creation of GenericUnaryFunction objects from lambdas or other callables.
   * 
   * @tparam Func The type of the callable to wrap.
   * @tparam Base The base class to inherit from, defaults to GenericUnaryFunction.
   */
  template<class Func, class Base=GenericUnaryFunction>
  class Wrapper: public Base {
  public:

    /**
     * @brief Construct a new Wrapper object.
     * 
     * @param f The callable to wrap.
     */
    Wrapper( Func&& f ) : f_{std::forward<Func>(f)} {
      static_assert( std::is_invocable_r_v<Real,std::decay_t<Func>,Real>,
                     "Callable must take and return Real" );
    }

    /**
     * @brief Apply the wrapped function.
     * 
     * @param x The input value.
     * @return The result of applying the wrapped function to x.
     */
    inline Real operator()( Real x ) const override { return f_(x); }

  private:
    /**
     * @brief Accept a visitor, part of the Visitor pattern implementation.
     * 
     * @param visitor The visitor to accept.
     */
    void accept( Visitor&& visitor ) const override { visitor.visit(*this); }

    Func f_; ///< The wrapped callable.
  }; // class Wrapper

  /**
   * @brief Class Template Argument Deduction (CTAD) guide for Wrapper.
   * 
   * This allows the compiler to deduce the template arguments for Wrapper
   * when constructing it from a callable.
   */
  template<class Func> 
  Wrapper( Func&& ) -> Wrapper<std::decay_t<Func>>;

private:
  /**
   * @brief A base class for visitors that implements the visit method.
   * 
   * This class uses the Curiously Recurring Template Pattern (CRTP) to
   * achieve static polymorphism, avoiding virtual function calls.
   * 
   * @tparam Derived The derived visitor class.
   */
  template<class Derived> 
  struct DerivedVisitor : public Visitor {
    /**
     * @brief Visit a GenericUnaryFunction object.
     * 
     * This method casts the visitor to the derived type and calls its visitImpl method.
     * 
     * @param uf The GenericUnaryFunction to visit.
     */
    void visit( const GenericUnaryFunction& uf ) override {
      static_cast<Derived*>(this)->visitImpl(uf);
    }
  }; // struct DerivedVisitor

  /**
   * @brief A visitor that applies a unary function to a vector.
   * 
   * This class implements the actual logic of applying a unary function to a vector.
   * 
   * @tparam VecT The type of the vector.
   * @tparam EvalT The type of the evaluator function.
   */
  template<class VecT, class EvalT>
  class VectorVisitor : public DerivedVisitor<VectorVisitor<VecT, EvalT>> {
  public:
    /**
     * @brief Construct a new VectorVisitor object.
     * 
     * @param vec The vector to apply the function to.
     * @param eval The evaluator function.
     */
    VectorVisitor( VecT& vec, EvalT&& eval ) 
    : vec_{vec}, eval_{std::forward<EvalT>(eval)} {}
  
    /**
     * @brief Apply the unary function to the vector.
     * 
     * This method is called by the visit method of the base DerivedVisitor.
     * 
     * @param uf The GenericUnaryFunction to apply.
     */
    void visitImpl( const GenericUnaryFunction& uf ) {
      eval_(vec_, uf);
    }

  private:
    VecT& vec_;  ///< Reference to the vector.
    EvalT eval_; ///< The evaluator function.
  }; // class VectorVisitor
}; // class GenericUnaryFunction

} // namespace ROL
