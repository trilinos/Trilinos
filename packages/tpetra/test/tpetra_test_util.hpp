// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef TPETRA_TEST_UTIL_HPP
#define TPETRA_TEST_UTIL_HPP

#include "Tpetra_ConfigDefs.hpp" // for iostream, etc.
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include "Tpetra_PacketTraits.hpp"
#include "Tpetra_Util.hpp" // for toString
#include "Tpetra_Version.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#endif // TPETRA_MPI

/*! \file tpetra_test_util.hpp
    \brief This file contains utility functions used by the Tpetra Tests.
    
    <ol>
    <li> convenience functions for building up values using only 
    one() and zero() from OrdinalTraits/ScalarTraits.
    <li> output functions that provide style sheet-like functionality
    <li> generator for producing non-contiguous values in an arbitrary type
    </ol>
*/

//! cin used as MPI breakpoint (used for debugging parallel code)
void mpiBreakpoint(int myImageID) {
  if(myImageID == 0) {
    cout << "[TPETRA-DEBUG mpiBreakpoint] ";
    cin >> myImageID;
  }
#ifdef TPETRA_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//======================================================================
// convenience functions for building up values in test suites.
// this way we don't have to use any literals.
//======================================================================

//! safely converts an int value to any defined OrdinalType
template <typename T>
T intToOrdinal(int n) {
  T result = Teuchos::OrdinalTraits<T>::zero();
  while(n > 0) {
    result += Teuchos::OrdinalTraits<T>::one();
    n--;
  }
  while(n < 0) {
    result -= Teuchos::OrdinalTraits<T>::one();
    n++;
  }
  return(result);
};

//! safely converts an int value to any defined ScalarType
template <typename T>
T intToScalar(int n) {
  T result = Teuchos::ScalarTraits<T>::zero();
  while(n > 0) {
    result += Teuchos::ScalarTraits<T>::one();
    n--;
  }
  while(n < 0) {
    result -= Teuchos::ScalarTraits<T>::one();
    n++;
  }
  return(result);
};

//======================================================================
// functions for outputting nicely formatted text
//======================================================================

//! Message outputed at beginning of a test program
void outputStartMessage(std::string const className) {
  cout << "\n************************************************************" << endl;
  cout << "Starting " << className << "Test..." << endl;
  cout << Tpetra::Tpetra_Version() << endl;
  cout << "************************************************************" << endl;
};

//! Message outputed at end of a test program
void outputEndMessage(std::string const className, bool passed) {
  cout << "************************************************************" << endl;
  cout << className << " test ";
  if(passed)
    cout << "passed." << endl;
  else
    cout << "failed." << endl;
  cout << "************************************************************\n" << endl;
};

//! Merssage outputed to seperate sections
void outputHeading(std::string const message) {
  cout << "**************************************************" << endl;
  cout << message << endl;
  cout << "**************************************************" << endl;
};

//! Message outputed to seperate subsections
void outputSubHeading(std::string const message) {
  cout << message << endl;
};

//! Outputs a message from all images in deterministic order
void outputData(int const myImageID, int const numImages, std::string const message) {
#ifdef TPETRA_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  for(int i = 0; i < numImages; i++) {
    if(myImageID == i)
      cout << "[Image " << myImageID << "] " << message << endl;
#ifdef TPETRA_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
};

//======================================================================
// functions for generator
//======================================================================

template <typename T>
T generateValue(T const x, T const y) {
  T const two = intToScalar<T>(2);

  // formula for z(x,y) = 0.5(x^2 + y^2 + 3x + y) + xy
  return(((x*x + y*y + x+x+x + y) / two) + (x*y));
}

// specialization for complex so that both real and imaginary portions get written to
// the real portion gets generateValue(x, 2y), and the imaginary portion gets generateValue(x, 2y+1)
template <typename T>
complex<T> generateValue(complex<T> const x, complex<T> const y) {
  T twoY = y.real() * intToScalar<T>(2);
  T real = generateValue(x.real(), twoY);
  T imag = generateValue(x.real(), (twoY + Teuchos::ScalarTraits<T>::one()));
  return(complex<T>(real, imag));
}

// puts the generated values for (x, 0) ... (x, length-1) into vector
template <typename T>
void generateColumn(std::vector<T>& vector, int const x, int const length) {
  vector.resize(length);
  for(int y = 0; y < length; y++)
    vector[y] = generateValue(intToScalar<T>(x), intToScalar<T>(y));
}

// puts the generated values for (firstx, 0)...(firstx, length-1),(firstx+1, 0),(firstx+1, length-1)...(lastx, length-1) into vector
template <typename T>
void generateMultipleColumns(std::vector<T>& vector, int const firstx, int const lastx, int const length) {
  vector.resize(length * (lastx - firstx + 1));
  typename std::vector<T>::iterator i = vector.begin();
  for(int x = firstx; x <= lastx; x++)
    for(int y = 0; y < length; y++) {
      *i = generateValue(intToScalar<T>(x), intToScalar<T>(y));
      i++;
    }
}

// puts the sum of (firstx, i) to (lastx, i) into vector[i], for i = 0 to length-1
template <typename T>
void generateRowSums(std::vector<T>& vector, int const firstx, int const lastx, int const length) {
  vector.assign(length, Teuchos::ScalarTraits<T>::zero());
  for(int x = firstx; x <= lastx; x++)
    for(int y = 0; y < length; y++)
      vector[y] += generateValue(intToScalar<T>(x), intToScalar<T>(y));
}

template <typename T>
void generateRowMaxs(std::vector<T>& vector, int const firstx, int const lastx, int const length) {
  generateColumn(vector, lastx, length); // for this generator, row max is right-most entry
}

template <typename T>
void generateRowMins(std::vector<T>& vector, int const firstx, int const lastx, int const length) {
  generateColumn(vector, firstx, length); // for this generator, row min is left-most entry
}

template <typename T>
T generateAreaMin(int const firstx, int const firsty, int const lastx, int const lasty) {
  return(generateValue(intToScalar<T>(firstx), intToScalar<T>(firsty))); // for this generator, area min is top left entry
}

template <typename T>
T generateAreaMax(int const firstx, int const firsty, int const lastx, int const lasty) {
  return(generateValue(intToScalar<T>(lastx), intToScalar<T>(lasty))); // for this generator, area max is bottom right entry
}

template <typename T>
void generateXCoords(std::vector<T> const& values, std::vector<T>& xCoords) {
  xCoords.resize(values.size());
  typename std::vector<T>::const_iterator inPtr = values.begin();
  typename std::vector<T>::iterator outPtr = xCoords.begin();
  for(; inPtr != values.end(); inPtr++, outPtr++) {
    T value = Teuchos::ScalarTraits<T>::zero(); // current generated value
    T x = Teuchos::ScalarTraits<T>::zero(); // xcoord for current value
    T max = Teuchos::ScalarTraits<T>::zero(); // right-most xcoord of current diagonal
    while(value < *inPtr) {
      value++;
      x++;
      if(x > max) {
        max++;
        x = Teuchos::ScalarTraits<T>::zero();
      }
    }
    *outPtr = x;
  }
}

template <typename T>
void generateYCoords(std::vector<T> const& values, std::vector<T>& yCoords) {
  yCoords.resize(values.size());
  typename std::vector<T>::const_iterator inPtr = values.begin();
  typename std::vector<T>::iterator outPtr = yCoords.begin();
  for(; inPtr != values.end(); inPtr++, outPtr++) {
    T value = Teuchos::ScalarTraits<T>::zero(); // current generated value
    T y = Teuchos::ScalarTraits<T>::zero(); // ycoord for current value
    T max = Teuchos::ScalarTraits<T>::zero(); // bottom-most ycoord of current diagonal
    while(value < *inPtr) {
      value++;
      if(y == Teuchos::ScalarTraits<T>::zero()) {
        max++;
        y = max;
      }
      else
        y--;
    }
    *outPtr = y;
  }
}

// generator tester function - kept around just in case
// although some of the newer additions aren't tested in it yet
/*
template <typename T>
void testGenerator() {
  int length = 10;
  int length2 = 4;
  
  cout << "\nTesting generateValue..." << endl;
  T tlength = intToScalar<T>(length);
  for(T y = Teuchos::ScalarTraits<T>::zero(); y < tlength; y++) {
    for(T x = Teuchos::ScalarTraits<T>::zero(); x < tlength; x++)
      cout << generateValue(x,y) << "\t";
    cout << endl;
  }
  
  cout << "\nTesting generateColumn..." << endl;
  std::vector< std::vector<T> > columns(length);
  for(int i = 0; i < length; i++) {
    generateColumn(columns[i], i, length);
    cout << "Column " << i << ": " << columns[i] << endl;
  }

  cout << "\nTesting generateMultipleColumns..." << endl;
  std::vector<T> multColumns;
  generateMultipleColumns(multColumns, 0, length2, length2);
  cout << "array bounds = (0,0) to (" << length2 << "," << length2 << ")"  << endl;
  cout << "array contents = " << multColumns << endl;
  
  cout << endl;
}
*/

#endif // TPETRA_TEST_UTIL_HPP
