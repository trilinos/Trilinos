#pragma once
#ifndef __TEST_MACRODEF_HPP__
#define __TEST_MACRODEF_HPP__

/// \file test_macordef.hpp
/// \brief Simple test macros
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#undef __DOT_LINE__
#define __DOT_LINE__ cout << string(60, '=') << endl;

#undef  __TestSuiteDoUnitTests__
#define __TestSuiteDoUnitTests__(VT,OT,ST,SpT,MeT) {                  \
    TestSuite<VT,OT,ST,SpT,MeT>::label =                              \
      "TestSuite<" #VT "," #OT "," #ST "," #SpT "," #MeT ">";         \
    r_val = TestSuite<VT,OT,ST,SpT,MeT>::doUnitTests();               \
  }

#undef  __ASSERT_TRUE__
#define __ASSERT_TRUE__(cond) {                                         \
    if (cond) {                                                         \
    } else {                                                            \
      cout << endl                                                      \
           << string(10, '>') << "Assertion Failed in :" << endl     \
           << string(10, ' ') << __FILE__ << ", " << __LINE__ << endl;  \
      ++r_val;                                                          \
    }                                                                   \
  }

#undef  __EVAL_STRING__
#define __EVAL_STRING__(r_value,eval_string) {          \
    if (r_value)                                        \
      eval_string = "FAILED " + to_string(r_value);     \
    else                                                \
      eval_string = "PASSED";                           \
  }


/*
extern int g_funct_counter;

#define FUNCT_ENTER                                                     \
  cout << setw(10) << right << string(++g_funct_counter, '+')  << ", "  \
  << __FUNCT__ << endl;

#define FUNCT_EXIT                                                      \
  cout << setw(10) << right << string(g_funct_counter--, '-')  << ", "  \
  << __FUNCT__ << endl;
*/


#endif
