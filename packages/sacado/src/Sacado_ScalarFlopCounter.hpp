// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_SCALAR_FLOP_COUNTER_HPP
#define SACADO_SCALAR_FLOP_COUNTER_HPP

#include "Sacado_ConfigDefs.h"

namespace Sacado {

namespace FlopCounterPack {

//
// Constants etc.
//

/// Number of total operation supported up till now
enum { NUM_OPS = 18 };
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
  ,SQRT
  ,ABS
};
/// Names for operations
const char* flopCountsNames[NUM_OPS]
= {
  "="
  ,"+"
  ,"+="
  ,"unary +"
  ,"-"
  ,"-="
  ,"unary -"
  ,"*"
  ,"*="
  ,"/"
  ,"/="
  ,">"
  ,">="
  ,"<"
  ,"<="
  ,"=="
  ,"sqrt"
  ,"abs"
};
/// Number of summary operation categories
enum { NUM_SUMMARY_OPS = 7 };
/// Enum of summary operation categories
enum ESummaryFlopType {
  SUMMARY_ASSIGN
  ,SUMMARY_PLUS_MINUS
  ,SUMMARY_MULTIPLY
  ,SUMMARY_DIVIDE
  ,SUMMARY_COMPARISON
  ,SUMMARY_SQRT
  ,SUMMARY_ABS
};
/// Names for summary operation categories
const char* summaryFlopCountsNames[NUM_SUMMARY_OPS]
= {
  "="
  ,"all +-"
  ,"all *"
  ,"all /"
  ,"<,>,=="
  ,"sqrt"
  ,"abs"
};
/// Aggregate struct for flop counts and summary flop counts
struct FlopCounts {
  /// Individual flop counts
  int flopCounts[NUM_OPS];
  /// Summary category flop counts
  int summaryFlopCounts[NUM_SUMMARY_OPS];
};

//
// Member macros
//

///
#define SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN( OP, OP_NAME ) \
ScalarFlopCounter<T> operator OP ( const ScalarFlopCounter<T>& s ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  val_ OP s.val(); \
  return *this; \
}

/** \brief Templated flop counter class.
 *
 * The purpose of this simple abstract data type is to count flops
 * within a computation.
 */
template<class T>
class ScalarFlopCounter {
public:
  
  /** @name Static functions for general clients (apply to all object with type <tt>T</tt> */
  //@{ 

  /** \breif Reset static flop counters before starting a block of computations. */
  static void resetCounters();

  /** \brief Get the flop counts after a block of computations. */
  static FlopCounts getCounters();

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
  static std::ostream& printCountersTable(
    const int n, const char* names[], const char* abbr[]
    ,const FlopCounts counts[], std::ostream &out
    );

  /** \brief Print the current static flop counts to <tt>out</tt>.
   *
   * This function just calls <tt>printCountersTable()</tt>.
   */
  static std::ostream& printCounters( std::ostream &out );
  
  //@{

  /** @name Object functions */
  //@{

  /// Construct to uninitialized
  ScalarFlopCounter();
  /// Construct to scalar value
  ScalarFlopCounter(const T &val);
  /// Return the current value
  const T& val() const;
  /// Set the current value
  void val(const T& a);
  SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(=,ASSIGN);
  SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(+=,PLUS_ASSIGN);
  SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(-=,MINUS_ASSIGN);
  SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(*=,MULTIPLY_ASSIGN);
  SCALAR_FLOP_COUNTER_BINARY_OP_ASSIGN(/=,DIVIDE_ASSIGN);

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
  static void incrCounter( const EFlopType& ft );

  //@}
};

// ////////////////////////////////////////
// Inline members

template<class T>
inline
ScalarFlopCounter<T>::ScalarFlopCounter() {}

template<class T>
inline
ScalarFlopCounter<T>::ScalarFlopCounter(const T &val) : val_(val) {}

template<class T>
inline
const T& ScalarFlopCounter<T>::val() const { return val_; }

template<class T>
inline
void ScalarFlopCounter<T>::val(const T& a) { val_ = a; }

// ///////////////////////////////////////
// Member implementations

// Static functions

template<class T>
void ScalarFlopCounter<T>::resetCounters()
{
  std::fill_n( &flopCounts_.flopCounts[0], int(NUM_OPS), int(0) );
  std::fill_n( &flopCounts_.summaryFlopCounts[0], int(NUM_SUMMARY_OPS), int(0) );
}

template<class T>
FlopCounts ScalarFlopCounter<T>::getCounters()
{
  return flopCounts_;
}

template<class T>
std::ostream& ScalarFlopCounter<T>::printCountersTable(
  const int n, const char* names[], const char* abbr[]
  ,const FlopCounts counts[], std::ostream &out
  )
{
  assert( n >= 1 && names && abbr && counts );
  const int w = 10;
  const char spacer[] = "----------";
  // Print legond
  if(names) {
    out << "\nLegend\n------\n";
    for( int j = 0; j < n; ++j )
      out << "  " << abbr[j] << " = " << names[j] << std::endl;
    out << std::endl;
  }
  // Print table header
  out << std::left << "  " << std::setw(w) << "op\\count";
  for( int j = 0; j < n; ++j ) out << "  " << std::setw(w) << abbr[j];
  out << std::endl;
  out << std::right << "  " << std::setw(w) << spacer;
  for( int j = 0; j < n; ++j ) out << "  " << std::setw(w) << spacer;
  out << std::endl;
  // Print rows of all operation counts
  for( int i = 0; i < NUM_OPS; ++i ) {
    int theseFlops = 0;
    for( int j = 0; j < n; ++j ) theseFlops += counts[j].flopCounts[i];
    if(theseFlops) {
      out << "  " << std::setw(w) << flopCountsNames[i];
      for( int j = 0; j < n; ++j ) out << "  " << std::setw(w) << counts[j].flopCounts[i];
      out << std::endl;
    }
  }
  out << std::right << "  " << std::setw(w) << spacer;
  for( int j = 0; j < n; ++j ) out << "  " << std::setw(w) << spacer;
  out << std::endl;
  // Print summary rows
  std::vector<int> totalFlops(n);
  for( int i = 0; i < NUM_SUMMARY_OPS; ++i ) {
    int theseFlops = 0;
    for( int j = 0; j < n; ++j ) {
      const int flops = counts[j].summaryFlopCounts[i];
      theseFlops += flops;
      totalFlops[j] += flops;
    }
    if(theseFlops) {
      out << "  " << std::setw(w) << summaryFlopCountsNames[i];
      for( int j = 0; j < n; ++j )
        out << "  " << std::setw(w) << counts[j].summaryFlopCounts[i];
      out << std::endl;
    }
  }
  out << std::right << "  " << std::setw(w) << spacer;
  for( int j = 0; j < n; ++j ) out << "  " << std::setw(w) << spacer;
  out << std::endl;
  // Print total flops
  out << "  " << std::setw(w) << "all flops";
  for( int j = 0; j < n; ++j ) out << "  " << std::setw(w) << totalFlops[j];
  out << std::endl;
  //
  return out;
}

template<class T>
std::ostream& ScalarFlopCounter<T>::printCounters( std::ostream &out )
{
  const int n = 1;
  const char* names[n] = { "Current" };
  const char* abbr[n]  = { "count" };
  const FlopCounts counts[n] = { flopCounts_ };
  return printCountersTable( n, &names[0], &abbr[0], &counts[0], out );
} 

// Static data members

template<class T> FlopCounts ScalarFlopCounter<T>::flopCounts_;

// Static public functions for nom-member funtions


template<class T>
void ScalarFlopCounter<T>::incrCounter( const EFlopType& ft )
{
  ++flopCounts_.flopCounts[ft];
  switch(ft) {
    case ASSIGN:
      ++flopCounts_.summaryFlopCounts[SUMMARY_ASSIGN];
      break;
    case PLUS:
    case PLUS_ASSIGN:
    case UNARY_PLUS:
    case MINUS:
    case MINUS_ASSIGN:
    case UNARY_MINUS:
      ++flopCounts_.summaryFlopCounts[SUMMARY_PLUS_MINUS];
      break;
    case MULTIPLY:
    case MULTIPLY_ASSIGN:
      ++flopCounts_.summaryFlopCounts[SUMMARY_MULTIPLY];
      break;
    case DIVIDE:
    case DIVIDE_ASSIGN:
      ++flopCounts_.summaryFlopCounts[SUMMARY_DIVIDE];
      break;
    case SQRT:
      ++flopCounts_.summaryFlopCounts[SUMMARY_SQRT];
      break;
    case ABS:
      ++flopCounts_.summaryFlopCounts[SUMMARY_ABS];
      break;
    case GREATER_THAN:
    case GREATER_THAN_EQUAL:
    case LESS_THAN:
    case LESS_THAN_EQUAL:
    case EQUAL:
      ++flopCounts_.summaryFlopCounts[SUMMARY_COMPARISON];
      break;
    default:
      assert(0);
  }
}

// ////////////////////////////////////////
// Nonmember operator function macros

#define SCALAR_FLOP_COUNTER_BINARY_OP( OP, OP_NAME ) \
template<class T> \
ScalarFlopCounter<T> operator OP ( const ScalarFlopCounter<T>& a, const ScalarFlopCounter<T>& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return ScalarFlopCounter<T>(a.val() OP b.val()); \
} \
template<class T> \
ScalarFlopCounter<T> operator OP ( const T& a, const ScalarFlopCounter<T>& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return ScalarFlopCounter<T>(a OP b.val()); \
} \
template<class T> \
ScalarFlopCounter<T> operator OP ( const int& a, const ScalarFlopCounter<T>& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return ScalarFlopCounter<T>(a OP b.val()); \
} \
template<class T> \
ScalarFlopCounter<T> operator OP ( const ScalarFlopCounter<T>& a, const T& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return ScalarFlopCounter<T>(a.val() OP b); \
} \
template<class T> \
ScalarFlopCounter<T> operator OP ( const ScalarFlopCounter<T>& a, const int& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return ScalarFlopCounter<T>(a.val() OP b); \
}

#define SCALAR_FLOP_COUNTER_UNARY_OP( OP, OP_NAME ) \
template<class T> \
ScalarFlopCounter<T> operator OP ( const ScalarFlopCounter<T>& a ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return ScalarFlopCounter<T>( OP a.val() ); \
}

#define SCALAR_FLOP_COUNTER_UNARY_FUNC( OP, OP_NAME ) \
template<class T> \
ScalarFlopCounter<T> OP( const ScalarFlopCounter<T>& a ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return ScalarFlopCounter<T>( std::OP( a.val() ) ); \
}

#define SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP( OP, OP_NAME ) \
template<class T> \
bool operator OP ( const ScalarFlopCounter<T>& a, const ScalarFlopCounter<T>& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return (a.val() OP b.val()); \
} \
template<class T> \
bool operator OP ( const T& a, const ScalarFlopCounter<T>& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return (a OP b.val()); \
} \
template<class T> \
bool operator OP ( const ScalarFlopCounter<T>& a, const T& b ) \
{ \
  ScalarFlopCounter<T>::incrCounter(OP_NAME); \
  return (a.val() OP b); \
}

// ////////////////////////////////////////////
// Nonmember operator and other functions

// Binary operations
SCALAR_FLOP_COUNTER_BINARY_OP(+,PLUS);
SCALAR_FLOP_COUNTER_BINARY_OP(-,MINUS);
SCALAR_FLOP_COUNTER_BINARY_OP(*,MULTIPLY);
SCALAR_FLOP_COUNTER_BINARY_OP(/,DIVIDE);
// Unary operations
SCALAR_FLOP_COUNTER_UNARY_OP(+,UNARY_PLUS);
SCALAR_FLOP_COUNTER_UNARY_OP(-,UNARY_MINUS);
// Unary functions
SCALAR_FLOP_COUNTER_UNARY_FUNC(sqrt,SQRT);
// Unary functions
SCALAR_FLOP_COUNTER_UNARY_FUNC(abs,ABS);
// Comparison
SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(>,GREATER_THAN);
SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(>=,GREATER_THAN_EQUAL);
SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(<,LESS_THAN);
SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(<=,LESS_THAN_EQUAL);
SCALAR_FLOP_COUNTER_BINARY_COMPARISON_OP(==,EQUAL);

} // namespace FlopCounterPack

} // namespace Sacado

#endif // SACADO_SCALAR_FLOP_COUNTER_HPP
