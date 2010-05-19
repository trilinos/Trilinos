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

#ifndef MOCK_OBJECT_HPP
#define MOCK_OBJECT_HPP

#include <iostream>

template<typename Scalar> class MockObject;

/** \brief Class to mimic a sacado::DFad object

    The DFad object allocates its own internal memory so its dtor need
    to be called to properly delete all data.  This object is used to
    prove that the dtors are called correctly on each chunk of the
    array.

    NOTE: If you really want the benefits of a contiguous block of
    memory in cache by using this allocator, you should never use a
    DFad object since the internal allocation will bring in
    non-contiguous memory.  You should be using something like SFad.
*/
template<typename Scalar>
class MockObject {

private:

  Scalar m_value;
  Scalar* m_derivative_values;
  int m_derivative_size;

public:
  
  static int m_num_times_empty_ctor_called;
  static int m_num_times_sized_ctor_called;
  static int m_num_times_copy_ctor_called;
  static int m_num_times_assignment_op_called;
  static int m_num_times_dtor_called;

public:

  MockObject() : 
    m_derivative_values(NULL), 
    m_derivative_size(0)
  {++m_num_times_empty_ctor_called;}
  
  MockObject(int derivative_size) : 
    m_derivative_values(NULL),
    m_derivative_size(derivative_size) 
  {
    m_derivative_values = new Scalar[m_derivative_size];
    ++m_num_times_sized_ctor_called;
  }

  MockObject(const MockObject& other) 
  {
    m_derivative_values = NULL;
    m_derivative_size = other.m_derivative_size;
    if (m_derivative_size > 0)
      m_derivative_values = new Scalar[m_derivative_size];
    ++m_num_times_copy_ctor_called;
  }
  
  MockObject& operator=(const MockObject& other)
  {
    if (m_derivative_values != NULL) {
      delete[] m_derivative_values;
      m_derivative_values = NULL;
    }
    m_derivative_size = other.m_derivative_size;
    if (m_derivative_size > 0)
      m_derivative_values = new Scalar[m_derivative_size];
    ++m_num_times_assignment_op_called;
    return *this;
  }

  ~MockObject()
  {
    delete[] m_derivative_values;
    ++m_num_times_dtor_called;
  }

private:


};

template<typename Scalar>
void print_num_ctor_dtor_calls(std::ostream& os)
{
  os << "    # times empty ctor called = " 
     << MockObject<Scalar>::m_num_times_empty_ctor_called << std::endl;
  os << "    # times sized ctor called = " 
     << MockObject<Scalar>::m_num_times_sized_ctor_called << std::endl;
  os << "    # times copy ctor called = " 
     << MockObject<Scalar>::m_num_times_copy_ctor_called << std::endl;
  os << "    # times assignment op called = " 
     << MockObject<Scalar>::m_num_times_assignment_op_called << std::endl;
  os << "    # times dtor called = " 
     << MockObject<Scalar>::m_num_times_dtor_called << std::endl;
}

template<>
int MockObject<int>::m_num_times_empty_ctor_called = 0;
template<>
int MockObject<int>::m_num_times_sized_ctor_called = 0;
template<>
int MockObject<int>::m_num_times_copy_ctor_called = 0;
template<>
int MockObject<int>::m_num_times_assignment_op_called = 0;
template<>
int MockObject<int>::m_num_times_dtor_called = 0;



#endif
