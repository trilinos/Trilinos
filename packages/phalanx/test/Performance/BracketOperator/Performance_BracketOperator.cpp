// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

// From test/Utilities directory
#include "Traits.hpp"

/*! \brief Test to check performance of bracket operator

    The Field class uses a bracket operator to access data elements.
    The operation is inlined and compilers will optimize this function
    call away so that it should be as fast as raw access.  This test
    will allow a comparison against raw access to verify that your
    compiler is at an acceptable optimization level.

    This test shows that the fields are a factor of 2 slower if
    running at -O0 instead of -O3 on linux gnu 4.2.4 compilers.  For
    the high optimizaiton level, there is virtually no difference
    in runtimes between the field and the raw data.

*/

struct Point : public shards::ArrayDimTag {
  Point(){};
  const char * name() const ;
  static const Point& tag();
};

const char * Point::name() const 
{ static const char n[] = "Point" ; return n ; }
const Point & Point::tag() 
{ static const Point myself ; return myself ; }

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    RCP<Time> arcp_time = TimeMonitor::getNewTimer("Field<double> Time");
    RCP<Time> mda_time = TimeMonitor::getNewTimer("MDField<double,Point> Time");
    RCP<Time> raw_time = TimeMonitor::getNewTimer("double* Time");
    
    TimeMonitor tm(*total_time);
    
    {
      
      const int num_loops = 1000;
      const int size = 200000;
      double value = 2.0;
      
      double* raw_density = new double[size];

      RCP<DataLayout> dl = rcp(new MDALayout<Point,Point>(1,size));
      Field<double> density("density", dl);
      ArrayRCP<double> a_density = arcp<double>(size);
      density.setFieldData(a_density);
      
      ArrayRCP<double> mda_density = arcp<double>(size);
      RCP<DataLayout> mddl = rcp(new MDALayout<Point,Point>(1,size));
      MDField<double,Point,Point> mddensity("density", mddl);
      mddensity.setFieldData(mda_density);
      
      for (int i=0; i < size; ++i) {
	raw_density[i] = 1.0;
	density[i] = 1.0;
	mddensity[i] = 1.0;
      }

      cout << "Field" << endl;
      {
	TimeMonitor tm(*arcp_time);
	for (int i=0; i < num_loops; ++i)
	  for (int j=0; j < size; ++j)
	    density[j] = value;
      }
      
      cout << "MDField" << endl;
      {
	TimeMonitor tm(*mda_time);
	for (int i=0; i < num_loops; ++i)
	  for (int j=0; j < size; ++j)
	    mddensity[j] = value;
      }
      
      cout << "double*" << endl;
      {
	TimeMonitor tm(*raw_time);
	for (int i=0; i < num_loops; ++i)
	  for (int j=0; j < size; ++j)
	    raw_density[j] = value;
      }
      
      delete [] raw_density;
    }

    // *********************************************************************
    // *********************************************************************
    cout << "\nTest passed!\n" << endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const exception& e) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Exception Caught!" << endl;
    cout << "Error message is below\n " << e.what() << endl;
    cout << "************************************************" << endl;
  }
  catch (...) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Unknown Exception Caught!" << endl;
    cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
