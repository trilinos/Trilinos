#ifndef __TACHO_TASK_VIEW_HPP__
#define __TACHO_TASK_VIEW_HPP__

/// \file Tacho_TaskView.hpp
/// \brief Matrix view with a task handler.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho { 

  template<typename MatViewType>
  class TaskView : public MatViewType {
  public:
    typedef typename MatViewType::value_type   value_type;
    typedef typename MatViewType::ordinal_type ordinal_type;
    typedef typename MatViewType::size_type    size_type;
    typedef typename MatViewType::space_type   space_type;

    typedef Kokkos::Experimental::Future<int,space_type> future_type;

  private:
    future_type _f;

  public:

    /// Interface functions
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)

    KOKKOS_INLINE_FUNCTION
    void setFuture(const future_type &f) { _f = f; }

    KOKKOS_INLINE_FUNCTION
    future_type Future() const { return _f; }

    /// ------------------------------------------------------------------

    /// Constructors
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)

    KOKKOS_INLINE_FUNCTION
    TaskView() 
      : MatViewType(), _f()
    { } 

    /// \brief Copy constructor is deleted to prevent from spreading future too much.
    KOKKOS_INLINE_FUNCTION
    TaskView(const TaskView &b) //= delete ;
      : MatViewType(b), _f(b._f) 
    { }

    /// \brief Wrapping the base object with view
    KOKKOS_INLINE_FUNCTION
    TaskView(typename MatViewType::mat_base_type const & b) 
      : MatViewType(b), _f() 
    { }

    /// \brief Wrapping the base object with view
    KOKKOS_INLINE_FUNCTION
    TaskView(typename MatViewType::mat_base_type const & b,
             const ordinal_type offm, const ordinal_type m,
             const ordinal_type offn, const ordinal_type n) 
      : MatViewType(b, offm, m, offn, n), _f() 
    { }

    /// Destructor
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (o)
    KOKKOS_INLINE_FUNCTION
    ~TaskView() = default;

  };
}

#endif
