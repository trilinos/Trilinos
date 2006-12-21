// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_RANDOM_HPP
#define SACADO_RANDOM_HPP

#include "Sacado_ConfigDefs.h"

namespace Sacado {

  /*! 
   * \brief A random number generator that generates random numbers uniformly
   * distributed in the interval (a,b).
   */
  class Random {
  public:
    
    //! Constructor
    Random(double a_, double b_);

    //! Constructor with seed value \c s
    Random(double a_, double b_, int s);

    //! Destructor
    ~Random();
    
    //! Set seed to \c s
    void setSeed(int s);

    //! Get random number
    double number();

  protected:

    // Check seed is valid
    int checkSeed(const std::string& func, int s);
  
  protected:

    //! Lower bound of interval
    double a;
    
    //! Upper bound of interval
    double b;

    //! %Random number seed
    double seed;

  }; // class Random

} // namespace Sacado

#endif // SACADO_RANDOM_HPP
