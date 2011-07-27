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

#include <iostream>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Inefficient vector/tensor objects
#include "Special_AlgebraicTypes.hpp"

#include "Shards_Array.hpp"

// TVMET - efficient expression template vector/tensor objects
#ifdef HAVE_PHALANX_TVMET
#include "tvmet/Vector.h"
#include "tvmet/Matrix.h"
#endif

struct Point : public shards::ArrayDimTag {
  Point(){};
  const char * name() const ;
  static const Point& tag();
};

struct Dim : public shards::ArrayDimTag {
  Dim(){};
  const char * name() const ;
  static const Dim& tag();
};

const char * Point::name() const 
{ static const char n[] = "Point" ; return n ; }
const Point & Point::tag() 
{ static const Point myself ; return myself ; }

const char * Dim::name() const 
{ static const char n[] = "Dim" ; return n ; }
const Dim & Dim::tag() 
{ static const Dim myself ; return myself ; }

int main(int argc, char* argv[])
{

  using namespace std;
  using namespace Teuchos;
  using namespace shards;

  GlobalMPISession mpi_session(&argc, &argv);

  const int num_samples = 3;
  const int num_loops = 5000;
  const int size = 1000;
  const int num_vectors = 3;

  TEST_FOR_EXCEPTION(num_loops * size != 5000000, std::logic_error,
		     "Work amount is not constant!");

  // Make all vectors in a contiguous block

  
  // 1. Dumb operator overloaded objects
  MyVector<double>* vector_array = new MyVector<double>[num_vectors * size];
  ArrayRCP< MyVector<double> > a = 
    arcp< MyVector<double> >(vector_array, 0, size, false);
  ArrayRCP< MyVector<double> > b = 
    arcp< MyVector<double> >(&vector_array[size], 0, size, false);
  ArrayRCP< MyVector<double> > c = 
    arcp< MyVector<double> >(&vector_array[2*size], 0, size, false);

  // 2. TVMET
#ifdef HAVE_PHALANX_TVMET
  tvmet::Vector<double, 3>* tvmet_array = 
    new tvmet::Vector<double, 3>[num_vectors * size];
  ArrayRCP< tvmet::Vector<double, 3> > d = 
    arcp< tvmet::Vector<double, 3> >(tvmet_array, 0, size, false);
  ArrayRCP< tvmet::Vector<double, 3> > e = 
    arcp< tvmet::Vector<double, 3> >(&tvmet_array[size], 0, size, false);
  ArrayRCP< tvmet::Vector<double, 3> > f = 
    arcp< tvmet::Vector<double, 3> >(&tvmet_array[2*size], 0, size, false);
#endif

  // 3. MultiDimensional Array Support
  double* mda_array = new double[num_vectors * size * 3];
  shards::Array<double,NaturalOrder,Point,Dim> mda_a(mda_array,size,3);
  shards::Array<double,NaturalOrder,Point,Dim> mda_b(&mda_array[size*3],size,3);
  shards::Array<double,NaturalOrder,Point,Dim> mda_c(&mda_array[2*size*3],size,3);

  // 4. Raw vector support
  double* raw_array = new double[num_vectors * size * 3];

  double* raw_a = raw_array;
  double* raw_b = &raw_array[size*3];
  double* raw_c = &raw_array[2*size*3];

  for (int i=0; i < a.size(); ++i)
    a[i] = 1.0;
  for (int i=0; i < b.size(); ++i)
    b[i] = 2.0;
  for (int i=0; i < c.size(); ++i)
    c[i] = 3.0;
#ifdef HAVE_PHALANX_TVMET
  for (int i=0; i < d.size(); ++i)
    d[i] = 1.0;
  for (int i=0; i < e.size(); ++i)
    e[i] = 2.0;
  for (int i=0; i < f.size(); ++i)
    f[i] = 3.0;
#endif
  for (int i=0; i < size; ++i) {
    int offset = i * 3;
    for (int j=0; j < 3; ++j) {
      raw_a[offset + j] = 1.0;
      raw_b[offset + j] = 2.0;
      raw_c[offset + j] = 3.0;
    }
  }

  RCP<Time> vector_time = TimeMonitor::getNewTimer("Vector Time");
#ifdef HAVE_PHALANX_TVMET
  RCP<Time> tvmet_time = TimeMonitor::getNewTimer("TVMET Time");
#endif
  RCP<Time> mda_time = TimeMonitor::getNewTimer("MultiDimensional Array");
  RCP<Time> raw_time = TimeMonitor::getNewTimer("Raw Time");

  for (int sample = 0; sample < num_samples; ++sample) {
    
    cout << "Vector" << endl;
    {
      TimeMonitor t(*vector_time);
      for (int i=0; i < num_loops; ++i)
	for (int j=0; j < c.size(); ++j)
	  c[j] = -4.0 * a[j] + b[j] * b[j];
    } 
    
#ifdef HAVE_PHALANX_TVMET
    cout << "TVMET" << endl;
    {
      TimeMonitor t(*tvmet_time);
      for (int i=0; i < num_loops; ++i)
	for (int j=0; j < d.size(); ++j)
	  f[j] = -4.0 * d[j] + e[j] *e[j];
    }
#endif
    
    cout << "MultiDimensionalArray" << endl;
    {
      TimeMonitor t(*mda_time);
      for (int i=0; i < num_loops; ++i) {
	for (int j=0; j < size; ++j) {
	  for (int k=0; k < 3; ++k)
	    mda_c(j,k) = -4.0 * mda_a(j,k) + mda_b(j,k) * mda_b(j,k);
	}
      }
    }

    cout << "Raw" << endl;
    {
      TimeMonitor t(*raw_time);
      for (int i=0; i < num_loops; ++i) {
	for (int j=0; j < size; ++j) {
	  int offset = j * 3;
	  for (int k=0; k < 3; ++k)
	    raw_c[offset + k] =
	      -4.0 * raw_a[offset + k] + raw_b[offset + k] * raw_b[offset + k];
	}
      }
    }
    
  } // end loop over samples

  cout << num_samples << " X " << num_loops << " X " << size << endl;
  TimeMonitor::summarize();
  
  double f_vector = vector_time->totalElapsedTime() / raw_time->totalElapsedTime();
#ifdef HAVE_PHALANX_TVMET
  double f_tvmet = tvmet_time->totalElapsedTime() / raw_time->totalElapsedTime();
#endif
  double f_mda = mda_time->totalElapsedTime() / raw_time->totalElapsedTime();
  double f_raw = raw_time->totalElapsedTime() / raw_time->totalElapsedTime();

  std::cout << "vector = " << f_vector << std::endl;
#ifdef HAVE_PHALANX_TVMET
  std::cout << "tvmet  = " << f_tvmet << std::endl;
#endif
  std::cout << "mda    = " << f_mda << std::endl;
  std::cout << "raw    = " << f_raw << std::endl;

  delete [] vector_array;
#ifdef HAVE_PHALANX_TVMET
  delete [] tvmet_array;
#endif
  delete [] mda_array;
  delete [] raw_array;

  std::cout << "\nTest passed!\n" << std::endl; 
    
  return 0;
}

