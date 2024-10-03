// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
