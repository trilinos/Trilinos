// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_SCALAR_FLOP_COUNTER_HPP
#define SACADO_SCALAR_FLOP_COUNTER_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_ScalarFlopCounterTraits.hpp"
#include "Sacado_Base.hpp"
#include "Sacado_SFINAE_Macros.hpp"
#include <cmath>
#include <algorithm>    // for std::min and std::max
#include <ostream>      // for std::ostream

namespace Sacado {

  namespace FlopCounterPack {

    /// Class storing flop counts and summary flop counts
    class  FlopCounts {
    public:

      /// Number of total operation supported up till now
      enum { NUM_OPS = 35 };

      /// Enum for operations
      enum EFlopType {
        ASSIGN
        ,PLUS
        ,PLUS_ASSIGN
        ,UNARY_PLUS
        ,MINUS
        ,MINUS_ASSIGN
        ,UNARY_MINUS
        ,MULTIPLY
        ,MULTIPLY_ASSIGN
        ,DIVIDE
        ,DIVIDE_ASSIGN
        ,GREATER_THAN
        ,GREATER_THAN_EQUAL
        ,LESS_THAN
        ,LESS_THAN_EQUAL
        ,EQUAL
        ,EXP
        ,LOG
        ,LOG10
        ,SQRT
        ,CBRT
        ,COS
        ,SIN
        ,TAN
        ,ACOS
        ,ASIN
        ,ATAN
        ,ATAN2
        ,COSH
        ,SINH
        ,TANH
        ,ABS
        ,POW
        ,MAX
        ,MIN
      };

      /// Number of summary operation categories
      enum { NUM_SUMMARY_OPS = 6 };

      /// Enum of summary operation categories
      enum ESummaryFlopType {
        SUMMARY_ASSIGN
        ,SUMMARY_PLUS_MINUS
        ,SUMMARY_MULTIPLY
        ,SUMMARY_DIVIDE
        ,SUMMARY_COMPARISON
        ,SUMMARY_NONLINEAR
      };

      /// Names of individual flops
      static const char* flopCountsNames[NUM_OPS];

      /// Names for summary operation categories
      static const char* summaryFlopCountsNames[NUM_SUMMARY_OPS];

      /*!
       * \brief The number of flops to accumulate as an integer before
       * converting to a double.
       */
      /*!
       * The default value is 100 000 000 and must be less than UINT_MAX-1.
       * Increasing this value may give somewhat better precision for the
       * flop count when counting very large numbers of flops.
       */
      static unsigned int flopGranularity;

      /// Individual flop counts
      double flopCounts[NUM_OPS];

      /// Summary category flop counts
      double summaryFlopCounts[NUM_SUMMARY_OPS];

      /// Total flop count
      double totalFlopCount;

      /// Default constructor
      FlopCounts();

      /// Reset flop counters before starting a block of computations. */
      void reset();

      //// Finalize total flop count after block of computations. */
      void finalize();

      /// Increment an individual flop counter.
      void increment(EFlopType ft);

    private:

      /// Get summary op enum from op enum
      ESummaryFlopType getSummaryType(EFlopType ft);

      /// Partial sum of individual flop counts
      unsigned int partialFlopCounts[NUM_OPS];

      /// Partial sum of summary category flop counts
      unsigned int partialSummaryFlopCounts[NUM_SUMMARY_OPS];
    };

    /** \brief Print a list of flop counts into a single table.
     *
     * \param  n  [in] Number of columns of flop counts
     * \param  names
     *            [in] Array (length <tt>n</tt>) of the names of each
     *            set of flop counts that will be used in the legend.
     * \param  abbr
     *            [in] Array (length <tt>n</tt>) of abbreviated names (less
     *            than 10 chars in length) that will be used as the column
     *            headings.
     * \param  counts
     *            [in] Array (length <tt>n</tt>) of the flop counts themselves.
     * \param  out
     *            [out] Stream that formated table is ouput to.
     */
    std::ostream& printCountersTable(const int n,
                                     const char* names[],
                                     const char* abbr[],
                                     const FlopCounts counts[],
                                     std::ostream &out);

    //
    // Member macros
    //

    ///
#define SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN( OP, OP_NAME )             \
    ScalarFlopCounter<T> operator OP ( const ScalarFlopCounter<T>& s )  \
    {                                                                   \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      val_ OP s.val();                                                  \
      return *this;                                                     \
    }

    /** \brief Templated flop counter class.
     *
     * The purpose of this simple abstract data type is to count flops
     * within a computation.
     */
    template<class T>
    class ScalarFlopCounter : public Base< ScalarFlopCounter<T> > {
    public:

      //! Typename of values
      typedef typename RemoveConst<T>::type value_type;

      //! Typename of scalar's (which may be different from T)
      typedef typename ScalarType<value_type>::type scalar_type;

      //! Turn ScalarFlopCounter into a meta-function class usable with mpl::apply
      template <typename U>
      struct apply {
        typedef ScalarFlopCounter<U> type;
      };

      /** @name Static functions for general clients (apply to all object with type <tt>T</tt> */
      //@{

      /** \brief Reset static flop counters before starting a block of computations. */
      static void resetCounters() { flopCounts_.reset(); }

      /** \brief Finalize total flop count after block of computations. */
      static void finalizeCounters() { flopCounts_.finalize(); }

      /** \brief Get the flop counts after a block of computations. */
      static FlopCounts getCounters() { return flopCounts_; }

      /** \brief Print the current static flop counts to <tt>out</tt>.
       *
       * This function just calls <tt>printCountersTable()</tt>.
       */
      static std::ostream& printCounters( std::ostream &out ) {
        const int n = 1;
        const char* names[n] = { "Current" };
        const char* abbr[n]  = { "count" };
        const FlopCounts counts[n] = { flopCounts_ };
        return printCountersTable( n, &names[0], &abbr[0], &counts[0], out );
      }

      //@{

      /** @name Object functions */
      //@{

      /// Construct to uninitialized
      ScalarFlopCounter() {}

      /// Construct to scalar value
      template <typename S>
      ScalarFlopCounter(const S &v, SACADO_ENABLE_VALUE_CTOR_DECL) : val_(v) {}

      /// Return the current value
      const T& val() const { return val_; }

      /// Set the current value
      void val(const T& a) { val_ = a; }

      SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(=,FlopCounts::ASSIGN);
      SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(+=,FlopCounts::PLUS_ASSIGN);
      SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(-=,FlopCounts::MINUS_ASSIGN);
      SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(*=,FlopCounts::MULTIPLY_ASSIGN);
      SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(/=,FlopCounts::DIVIDE_ASSIGN);

      //@}

    private:

      // Static members
      static FlopCounts flopCounts_;

      // Object members
      T val_;

    public:

      /** Static public functions for non-member functions (not for general clients) */
      //@{

      /** \brief Increment an individual flop counter.
       *
       * Note, this function is ment to be used by nonmember operator
       * functions and not by general clients.
       */
      static void incrCounter( const FlopCounts::EFlopType& ft ) {
        flopCounts_.increment(ft);
      }

      //@}
    };

    // Static data members

    template<class T> FlopCounts ScalarFlopCounter<T>::flopCounts_;

    // ////////////////////////////////////////
    // Nonmember operator function macros

#define SCALAR_FLOP_COUNTER_BINARY_OP( OP, OP_NAME )                    \
    template<class T>                                                   \
    ScalarFlopCounter<T> operator OP (                                  \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>(a.val() OP b.val());                  \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> operator OP (                                  \
      const typename ScalarFlopCounter<T>::value_type& a,               \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>(a OP b.val());                        \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> operator OP (                                  \
      const int& a,                                                     \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>(a OP b.val());                        \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> operator OP (                                  \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const typename ScalarFlopCounter<T>::value_type& b )              \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>(a.val() OP b);                        \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> operator OP (                                  \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const int& b )                                                    \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>(a.val() OP b);                        \
    }

#define SCALAR_FLOP_COUNTER_UNARY_OP( OP, OP_NAME )                     \
    template<class T>                                                   \
    ScalarFlopCounter<T> operator OP (                                  \
      const Base< ScalarFlopCounter<T> >& aa )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>( OP a.val() );                        \
    }

#define SCALAR_FLOP_COUNTER_UNARY_FUNC( OP, OP_NAME )                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> OP(                                            \
      const Base< ScalarFlopCounter<T> >& aa )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>( std::OP( a.val() ) );                \
    }

#define SCALAR_FLOP_COUNTER_BINARY_FUNC( OP, OP_NAME )                  \
    template<class T>                                                   \
    ScalarFlopCounter<T> OP (                                           \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>( std::OP( a.val(), b.val() ) );       \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> OP (                                           \
      const typename ScalarFlopCounter<T>::value_type& a,               \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>( std::OP( a, b.val() ) );             \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> OP (                                           \
      const int& a,                                                     \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>( std::OP( a, b.val() ) );             \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> OP (                                           \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const typename ScalarFlopCounter<T>::value_type& b )              \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>( std::OP( a.val(), b ) );             \
    }                                                                   \
    template<class T>                                                   \
    ScalarFlopCounter<T> OP (                                           \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const int& b )                                                    \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return ScalarFlopCounter<T>( std::OP(a.val(), b ) );              \
    }

#define SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP( OP, OP_NAME )         \
    template<class T>                                                   \
    bool operator OP (                                                  \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return (a.val() OP b.val());                                      \
    }                                                                   \
    template<class T>                                                   \
    bool operator OP (                                                  \
      const typename ScalarFlopCounter<T>::value_type& a,               \
      const Base< ScalarFlopCounter<T> >& bb )                          \
    {                                                                   \
      const ScalarFlopCounter<T>& b = bb.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return (a OP b.val());                                            \
    }                                                                   \
    template<class T>                                                   \
    bool operator OP (                                                  \
      const Base< ScalarFlopCounter<T> >& aa,                           \
      const typename ScalarFlopCounter<T>::value_type& b )              \
    {                                                                   \
      const ScalarFlopCounter<T>& a = aa.derived();                     \
      ScalarFlopCounter<T>::incrCounter(OP_NAME);                       \
      return (a.val() OP b);                                            \
    }

    // ////////////////////////////////////////////
    // Nonmember operator and other functions

    // Binary operations
    SCALAR_FLOP_COUNTER_BINARY_OP(+,FlopCounts::PLUS)
    SCALAR_FLOP_COUNTER_BINARY_OP(-,FlopCounts::MINUS)
    SCALAR_FLOP_COUNTER_BINARY_OP(*,FlopCounts::MULTIPLY)
    SCALAR_FLOP_COUNTER_BINARY_OP(/,FlopCounts::DIVIDE)

    // Unary operations
    SCALAR_FLOP_COUNTER_UNARY_OP(+,FlopCounts::UNARY_PLUS)
    SCALAR_FLOP_COUNTER_UNARY_OP(-,FlopCounts::UNARY_MINUS)

    // Unary functions
    SCALAR_FLOP_COUNTER_UNARY_FUNC(exp,FlopCounts::EXP)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(log,FlopCounts::LOG)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(log10,FlopCounts::LOG10)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(sqrt,FlopCounts::SQRT)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(cbrt,FlopCounts::CBRT)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(cos,FlopCounts::COS)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(sin,FlopCounts::SIN)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(tan,FlopCounts::TAN)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(acos,FlopCounts::ACOS)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(asin,FlopCounts::ASIN)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(atan,FlopCounts::ATAN)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(cosh,FlopCounts::COSH)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(sinh,FlopCounts::SINH)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(tanh,FlopCounts::TANH)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(abs,FlopCounts::ABS)
    SCALAR_FLOP_COUNTER_UNARY_FUNC(fabs,FlopCounts::ABS)

    // Binary functions
    SCALAR_FLOP_COUNTER_BINARY_FUNC(atan2,FlopCounts::ATAN2)
    SCALAR_FLOP_COUNTER_BINARY_FUNC(pow,FlopCounts::POW)
    SCALAR_FLOP_COUNTER_BINARY_FUNC(max,FlopCounts::MAX)
    SCALAR_FLOP_COUNTER_BINARY_FUNC(min,FlopCounts::MIN)

    // Comparison
    SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(>,FlopCounts::GREATER_THAN)
    SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(>=,FlopCounts::GREATER_THAN_EQUAL)
    SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(<,FlopCounts::LESS_THAN)
    SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(<=,FlopCounts::LESS_THAN_EQUAL)
    SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(==,FlopCounts::EQUAL)

  } // namespace FlopCounterPack
} // namespace Sacado

#endif // SACADO_SCALAR_FLOP_COUNTER_HPP
