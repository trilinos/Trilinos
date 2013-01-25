// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_RANDOM_H
#define NOX_RANDOM_H

#include "NOX_Common.H"
#include <string>

namespace NOX {

  /*!
    \brief A class to compute uniformly distributed random numbers in (-1,1).
  */
  /*!
    The %Random class computes pseudo-random (double precision) numbers 
    uniformly distributed between -1.0 and 1.0 using a multiplicative
    congruential generator with modulus 2^31-1 (a Lehmer generator).  
    For a numerical and mathematical treatment of the algorithm, see
    "Random number generators:  good ones are hard to find" by Stephen
    K. Park and Keith W. Miller, Communications of the ACM, Vol. 31 No. 10
    (1988).
  */
  class Random {
  public:

    //! Initialize random number generator with a random seed.
    /*! The random seed is computed using the POSIX rand() function. */
    Random();

    //! Initialize random number generator with the given seed
    /*! The seed should be an integer between 1 and 2147483646 = 2^32-2
        (inclusive).  If the supplied seed is invalid, an error message
	is printed and the seed is replaced by 1.
    */
    Random(int s);

    //! Set the seed to s
    /*! The seed should be an integer between 1 and 2147483646 = 2^32-2
        (inclusive).  If the supplied seed is invalid, an error message
	is printed and the seed is replaced by 1.
    */
    static void setSeed(int s);
      
    //! Compute a random number
    static double number();

  protected:

    //! Checks that the supplied seed value is valid and returns a valid seed.
    /*! The %string argument supplies the calling function name for an error
        message in the case that the supplied seed is invalid.
    */
    static int checkSeed(const std::string& func, int s);

  protected:
      
    //! Seed value
    static double seed;
  };

}
#endif
