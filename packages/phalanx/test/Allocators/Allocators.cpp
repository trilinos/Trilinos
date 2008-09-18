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


template<typename T>
bool arrayCompare(Teuchos::ArrayRCP<T> a, Teuchos::ArrayRCP<T> b, double tol)
{
  TEST_FOR_EXCEPTION(a.size() != b.size(), std::logic_error, "The arrays in checkTolerance() are not of the same size!");

  double error = 0.0;

  for (typename Teuchos::ArrayRCP<T>::Ordinal i = 0; i < a.size(); ++i)
    error += fabs(a[i] - b[i]);

  if (error > tol)
    return true;
  
  return false;
}

class MyObject {
public:
  static bool destructor_called;
  MyObject() {}
  ~MyObject() 
  {
    destructor_called = true;
    std::cout << "      In Destructor!" << std::endl;
  }
};

bool MyObject::destructor_called = false;

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Allocator Testing
    // *********************************************************************
    {
      cout << "Testing  NewAllocator:" << endl;

      cout << "  Testing ctor...";
      NewAllocator na;
      cout << "passed!" << endl;

      cout << "  Testing allocate()...";
      ArrayRCP<double> vec = na.allocate<double>(4);
      TEST_FOR_EXCEPTION(vec.size() != 4, std::runtime_error, 
			 "Allocator is broken!");
      cout << "passed!" << endl;

    }

    {
      cout << "\nTesting ContiguousAllocator:" << endl;
      cout << "  Testing ctor...";
      ContiguousAllocator<double> ca;
      cout << "passed!" << endl;

      cout << "  Testing addRequiredChunk()...";
      const int size = 10;
      const int num_bytes = size * sizeof(double);
      ca.addRequiredChunk(sizeof(double), size);
      ca.addRequiredChunk(sizeof(double), size);
      ca.addRequiredChunk(sizeof(double), size);
      ca.addRequiredChunk(sizeof(double), 2 * size);
      cout << "passed!" << endl;
      
      cout << "  Testing getTotalBytes()...";
      const int total_bytes = ca.getTotalBytes();
      TEST_FOR_EXCEPTION(total_bytes != 5 * num_bytes, std::logic_error,
			 "addRequiredBytes() is failing!");
      cout << "passed!" << endl;
            
      cout << "  Testing setup()...";
      ca.setup();
      cout << "passed!" << endl;

      cout << "  Testing allocate()...";
      ArrayRCP<double> a = ca.allocate<double>(size);
      ArrayRCP<double> b = ca.allocate<double>(size);
      ArrayRCP<double> c = ca.allocate<double>(size);
      ArrayRCP<double> d = ca.allocate<double>(2 * size);
      TEST_FOR_EXCEPTION(a.size() != size, std::runtime_error, 
			 "Allocator is broken!");
      TEST_FOR_EXCEPTION(d.size() != 2 * size, std::runtime_error, 
			 "Allocator is broken!");
      cout << "passed!" << endl;

      
      cout << "  Testing array integrity...";
      
      for (ArrayRCP<double>::Ordinal i = 0; i < a.size(); ++i)
	a[i] = 1.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < b.size(); ++i)
	b[i] = 2.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < c.size(); ++i)
	c[i] = 3.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < d.size(); ++i)
	d[i] = 4.0;

//       for (ArrayRCP<double>::Ordinal i = 0; i < a.size(); ++i)
// 	cout << "a[" << i << "] = " << a[i] << endl;
//       for (ArrayRCP<double>::Ordinal i = 0; i < b.size(); ++i)
// 	cout << "b[" << i << "] = " << b[i] << endl;
//       for (ArrayRCP<double>::Ordinal i = 0; i < c.size(); ++i)
// 	cout << "c[" << i << "] = " << c[i] << endl;
//       for (ArrayRCP<double>::Ordinal i = 0; i < d.size(); ++i)
// 	cout << "d[" << i << "] = " << d[i] << endl;

      ArrayRCP<double> a_ref = arcp<double>(size);
      ArrayRCP<double> b_ref = arcp<double>(size);
      ArrayRCP<double> c_ref = arcp<double>(size);
      ArrayRCP<double> d_ref = arcp<double>(2 * size);

      for (ArrayRCP<double>::Ordinal i = 0; i < a.size(); ++i)
	a_ref[i] = 1.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < b.size(); ++i)
	b_ref[i] = 2.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < c.size(); ++i)
	c_ref[i] = 3.0;
      for (ArrayRCP<double>::Ordinal i = 0; i < d.size(); ++i)
	d_ref[i] = 4.0;

      TEST_FOR_EXCEPTION(arrayCompare(a, a_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      TEST_FOR_EXCEPTION(arrayCompare(b, b_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      TEST_FOR_EXCEPTION(arrayCompare(c, c_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      TEST_FOR_EXCEPTION(arrayCompare(d, d_ref, 1.0e-8), std::runtime_error,
			 "Array Integrity failed.  Arrays are different!");
      
      cout << "passed!" << endl;


      cout << "  Testing reset()...";
      ca.reset();
      cout << "passed!" << endl;

      
      {
	cout << "  Testing that ArrayRCP embedded object is destroyed properly:" << endl;

	cout << "    building objects...";
	int size = 5;
	ArrayRCP<MyObject> contiguous_block = arcp<MyObject>(size);
	MyObject* ptr = &contiguous_block[3];
	ArrayRCP<MyObject> subview = 
	  arcpWithEmbeddedObjPostDestroy( ptr, 0, 2, contiguous_block, false);
	cout << "finished" << endl;
	
	cout << "    deleting contiguous_block ptr...";
	TEST_FOR_EXCEPTION(MyObject::destructor_called, std::logic_error,
			   "MyObject is not initialized correctly!");
	contiguous_block = null;
	TEST_FOR_EXCEPTION(MyObject::destructor_called, std::logic_error,
			   "MyObject array should not be destroyed yet!");
	cout << "finished" << endl;

	cout << "    deleting subview ptr (should print destructor message "
	     << size << " times):" << endl;
	subview = null;
	TEST_FOR_EXCEPTION(!MyObject::destructor_called, std::logic_error,
			   "Failed to destroy objects correctly!");

	cout << "  Testing that ArrayRCP embedded object is destroyed properly...passed" << endl;
      }
      
    }

    // *********************************************************************
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
