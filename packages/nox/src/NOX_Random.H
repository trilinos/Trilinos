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
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_RANDOM_H
#define NOX_RANDOM_H

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
    Random(double s);

    //! Set the seed to s
    void setSeed(double s);
      
    //! Compute a random number
    double operator() ();

  protected:
      
    //! Seed value
    double seed;
  };

}
#endif
