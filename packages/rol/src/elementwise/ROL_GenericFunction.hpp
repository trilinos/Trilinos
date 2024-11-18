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

#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace ROL {

constexpr std::size_t MAX_ARITY = 16;

// Function traits to deduce arity
template<typename T>
struct function_traits : public function_traits<decltype(&T::operator())> {}; 

template<typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const> {
  static constexpr std::size_t arity = sizeof...(Args);
  using result_type = ReturnType;
};



/**
 * @brief A template class representing a function that can be applied to Vectors.
 * 
 * This class provides a mechanism for applying functions to Vectors efficiently,
 * using the Visitor pattern to avoid virtual function calls for each element.
 * 
 * @tparam Real The type used for real numbers in the function and calculations.
 */
template<class Real>
class GenericFunction {
protected:
  /**
   * @brief Abstract base class for visitors.
   * 
   * This class is part of the Visitor pattern implementation, allowing for
   * double dispatch without virtual function calls for each vector element.
   */
  struct Visitor { 
    virtual void visit( const GenericFunction& ) = 0; 
  };

  explicit GenericFunction( std::size_t arity ) : arity_{arity} {}

public:

  GenericFunction( const GenericFunction& ) = delete;
  GenericFunction& = operator ( const GenericFunction& ) = delete;
  GenericFunction( GenericFunction&& ) = default;
  GenericFunction& = operator ( GenericFunction&& ) = default;

  /** @brief Virtual destructor for proper cleanup of derived classes. */
  virtual ~GenericFunction() noexcept = default;

  /**
   * @brief Apply the function to a single value.
   * 
   * @param x The input values.
   * @return The result of applying the function to x.
   */
  [[nodiscard]] virtual Real operator()( const std::vector<Real>& x ) const = 0;

  /**
   * @brief Accept a visitor, part of the Visitor pattern implementation.
   * 
   * @param visitor The visitor to accept.
   */
  virtual void accept( Visitor&& visitor ) const = 0;

  /**
   * @brief Apply the function to a vector.
   * 
   * This method uses the Visitor pattern to efficiently apply the function
   * to all elements of the vector without virtual function calls per element.
   * 
   * @tparam VecT The type of the vector.
   * @tparam EvalT The type of the evaluator function. 
   * @param vec The vector to apply the function to.
   * @param eval The evaluator function that defines how to apply the function to the vector.
   * 
   * Example Usage:
   * @code
   *   template<class Real>
   *   class Vector {   
   *   public:
   *     virtual void applyFunction( const GenericFunction<Real>&, const std::vector<const Vector*>& ) {} 
   *     virtual int dimension() const { return 0; }
   *   };
   *
   *   template<class Real>
   *   class StdVector : public Vector<Real> {
   *   public:
   *     void dimension() const override { return vec_.size(); }
   *     void applyGenericUnary( const GenericFunction<Real>& gf, const std::vector<const Vector<Real>*>& vecs ) override {
   *       std::vector<const StdVector*> stdVecs;
   *       stdVecs.reserve(vecs.size());
   *       for(auto vec : vecs) {
   *         stdVecs.push_back(static_cast<const StdVector*>(vec));
   *       }
   *       gf.apply_vectorized(stdVecs,[this](const auto& vecs, const auto& f){ this->applyFunctionImpl(f,vecs); });
   *     }
   *
   *     template<class F>
   *     void applyFunctionImpl( const F& f, const std::vector<const StdVector*>& vecs ) { 
   *       std::vector<Real> args(vecs.size());
   *       for( int i=0; i<vec_.size(); ++i ) {
   *         for( int j=0; j<vecs.size(); ++j ) {
   *           args[j] = (*vecs[j]).vec_[i];
   *         }
   *         vec_[i] = f(args);
   *       }
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
   * @brief A wrapper class that turns any callable into a GenericFunction.
   * 
   * This class allows easy creation of GenericFunction objects from lambdas or other callables.
   * 
   * @tparam Func The type of the callable to wrap.
   * @tparam Base The base class to inherit from, defaults to GenericFunction.
   */
  template<class Func, class Base=GenericFunction>
  class Lambda : public Base {
  public:

    static constexpr std::size_t arity = function_traits<Func>::arity; ///< Number of arguments wrapped callable takes

    /**
     * @brief Construct a new Wrapper object.
     * 
     * @param f The callable to wrap.
     */
    explicit Lambda( Func&& f ) 
    : Base(arity), f_{std::forward<Func>(f)} {
      static_assert( arity > 0, "Callable must take at least one argument" );
      static_assert( arity <= MAX_ARITY, "Callable must not take more than MAX_ARITY arguments" );
    }

    /**
     * @brief Apply the wrapped function.
     * 
     * @param x The input value.
     * @return The result of applying the wrapped function to x.
     */
    Real operator()( const std::vector<Real>& x ) const override {
      if( x.size() != this->arity ) {
        std::stringstream msg;
        msg << "Received vector of " << x.size() << " arguments, but the wrapped callable's arity is " << this->arity << ".";       
        throw std::invalid_argument(msg.str());
      }
      return call_impl(x, std::make_index_sequence<arity>{});
    }

  private:

    /**
     * @brief Expands the elements of a vector into the arguments of a function 
     * @param x The vector of arguments to pass to the callable
     */
    template<std::size_t...I>
    Real call_impl( const std::vector<Real>& x, std::index_sequence<I...> ) const {
      return f_(x[I]...);
    }

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
  Lambda( Func&& ) -> Lambda<std::decay_t<Func>>;

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
     * @brief Visit a GenericFunction object.
     * 
     * This method casts the visitor to the derived type and calls its visitImpl method.
     * 
     * @param uf The GenericFunction to visit.
     */
    void visit( const GenericFunction& gf ) override {
      static_cast<Derived*>(this)->visitImpl(gf);
    }
  }; // struct DerivedVisitor

  /**
   * @brief A visitor that applies a function to a vector.
   * 
   * This class implements the actual logic of applying a function to a vector.
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
    VectorVisitor( const std::vector<const VecT*>& vecs, EvalT&& eval ) 
    : vecs_{vecs}, eval_{std::forward<EvalT>(eval)} {}
  
    /**
     * @brief Apply the function to the vector.
     * 
     * This method is called by the visit method of the base DerivedVisitor.
     * 
     * @param uf The GenericFunction to apply.
     */
    void visitImpl( const GenericFunction& gf ) {
      eval_(vecs_, gf);
    }

  private:
    const std::vector<const VecT*>& vec_;  ///< Reference to the pointers to vectors.
    EvalT eval_; ///< The evaluator function.
  }; // class VectorVisitor

  const std::size_t arity_; ///< Number of Real arguments the function takes

}; // class GenericFunction

} // namespace ROL
