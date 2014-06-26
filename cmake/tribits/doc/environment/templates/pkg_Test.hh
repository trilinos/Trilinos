//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<spkg>_test.hh
 * \author <user>
 * \date   <date>
 * \brief  Testing harness for <pkg>.
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: pkg_Test.hh,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#ifndef <spkg>_test_hh
#define <spkg>_test_hh

#include <iostream>
#include <string>

namespace <namespace>_test
{

//===========================================================================//
// PASS/FAILURE LIMIT
//===========================================================================//

// Returns true for pass
// Returns false for fail
// Failure functions also set <namespace>_test::passed to false

// These can be used in any combination in a test to print output messages  
// if no fail functions are called then the test will pass
// (<namespace>_test::passed will have its default value of true)

// Needless to say, these can be used in many different combinations or
// ways.  We do not constrain nemesis tests except that the output must be of
// the form "Test: pass/fail"

bool fail(int line);

bool fail(int line, char *file);

bool pass_msg(const std::string &);

bool fail_msg(const std::string &);

void unit_test(const bool pass, int line, char *file);

//---------------------------------------------------------------------------//
// PASSING CONDITIONALS
//---------------------------------------------------------------------------//

extern bool passed;

} // end namespace <namespace>_test

//===========================================================================//
// TEST MACROS
//
// USAGE:
// if (!condition) ITFAILS;
//
// These are a convenience only
//===========================================================================//

#define ITFAILS      <namespace>_test::fail(__LINE__);
#define FAILURE      <namespace>_test::fail(__LINE__, __FILE__);
#define PASSMSG(a)   <namespace>_test::pass_msg(a);
#define FAILMSG(a)   <namespace>_test::fail_msg(a);
#define UNIT_TEST(x) <namespace>_test::unit_test(x, __LINE__, __FILE__)
    
#endif // <spkg>_test_hh

//---------------------------------------------------------------------------//
//     end of <pkg>/<spkg>_test.hh
//---------------------------------------------------------------------------//
