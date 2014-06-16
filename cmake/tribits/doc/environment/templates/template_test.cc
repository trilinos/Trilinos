//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <tpkg>/<class>.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: template_test.cc,v 1.8 2008/01/04 22:50:12 9te Exp $
//---------------------------------------------------------------------------//

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

#include "harness/DBC.hh"
#include "harness/Scalar_Unit_Test.hh"
#include "harness/Soft_Equivalence.hh"
#include "../Release.hh"

using namespace std;
using namespace nemesis;
using namespace <namespace>;

#define ITFAILS ut.failure(__LINE__);
#define UNIT_TEST(a) if (!(a)) ut.failure(__LINE__);

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//



//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    Scalar_Unit_Test ut(argc, argv, release);
    
    try
    {
        // >>> UNIT TESTS
    }
    catch (std::exception &err)
    {
        std::cout << "ERROR: While testing <class>, " 
                  << err.what()
                  << endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing <class>, " 
                  << "An unknown exception was thrown."
                  << endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
//                        end of <class>.cc
//---------------------------------------------------------------------------//
