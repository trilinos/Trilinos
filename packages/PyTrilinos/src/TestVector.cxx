#include "TestVector.h"
#include <iostream>

#ifndef NumPyArray_h
#  include "NumPyArray.h"
#endif

TestVector::TestVector(int length)
  : m_vector(length)
{
  // Does nothing here yet
}

TestVector::~TestVector()
{
  // Does nothing here yet
}

void TestVector::write() const
{
  int index;
  const int length = m_vector.size();
  for (index = 0; index < length; ++index)
  {
    std::cout << "Vector[" << index << "] = " << m_vector[index] << std::endl;
  }
}

void TestVector::fromNumericArray(PyObject * p_pyObject)
{

  NumPyArrayContiguous numPyArray(p_pyObject);

  const int      numericArrayNumDims      = numPyArray.getNumDims    ();
  const int    * p_numericArrayDimLengths = numPyArray.getDimLengths ();
  const int      numericArrayTotalLength  = numPyArray.getTotalLength();
  double       * p_numericArray           = numPyArray.getDataArray  ();

  int index;

  std::cout << "numericArrayNumDims      = " << numericArrayNumDims     << std::endl;
  std::cout << "p_numericArrayDimLengths = ";
  char * p_comma = "";
  for (int index = 0; index < numericArrayNumDims; ++index) {
    std::cout << p_comma << p_numericArrayDimLengths[index];
    p_comma = ", ";
  }
  std::cout << std::endl;
  std::cout << "numericArrayTotalLength  = " << numericArrayTotalLength << std::endl;

  m_vector.resize(numericArrayTotalLength);


  memcpy (&m_vector[0], p_numericArray, numericArrayTotalLength*sizeof(double));

//   for (int index = 0; index < numericArrayTotalLength; ++index) {
//     m_vector[index] = p_numericArray[index];
//   }
}

void TestVector::testit(PyObject * p_pyObject)
{
#if 0
  // Here we assume p_pyObject is a Python list with a single entry.
  // We unpack the entry to get at the PyArrayObject.  Do this so
  // can make sure object is not mucked with by any SWIG code.
  // For testing, can also send a Python string in the list which will be printed.
  
  std::cout << "In TestVector::testit()" << std::endl;

  import_array();

  PyArrayObject * array     = 0;
  PyObject      * listItem0 = 0;

  if (PyString_Check(p_pyObject)) {
    std::cerr << "p_pyObject is a PyString object." << std::endl;
  } else {
    std::cerr << "p_pyObject is NOT a PyString object." << std::endl;
  }

  if (PyList_Check(p_pyObject)) {
    std::cerr << "p_pyObject is a PyList object." << std::endl;
  } else {
    std::cerr << "p_pyObject is NOT a PyList object." << std::endl;
  }
  
  if (PyList_Check(p_pyObject)) {
    std::cerr << "About to get list item..." << std::endl;
    listItem0 = PyList_GetItem(p_pyObject, 0);
    
    std::cerr << "About to check if list item is a string..." << std::endl;

    if (PyString_Check(listItem0)) {
      char * p_string = PyString_AsString(listItem0);
      std::cerr << "Found the string '" << p_string << "' in the list." << std::endl;
    } else {
      // If here assume it is a PyArrayObject
      std::cerr << "List item is not a string, lets treat it as a PyArrayObject..." << std::endl;
      array = 
        (PyArrayObject *)PyArray_ContiguousFromObject(listItem0,
                                                      PyArray_DOUBLE, 0, 0);
    }
  } else {
    std::cerr << "Error: Object sent in list is NOT a PyList object." << std::endl;
  }

  std::cout << "Leaving TestVector::testit()" << std::endl;
#endif
}
