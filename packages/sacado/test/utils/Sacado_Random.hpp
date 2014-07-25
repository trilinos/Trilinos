// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_RANDOM_HPP
#define SACADO_RANDOM_HPP

#include "Sacado_ConfigDefs.h"

#include <string>

namespace Sacado {

  /*! 
   * \brief A random number generator that generates random numbers uniformly
   * distributed in the interval (a,b).
   */
  template <typename ScalarT>
  class Random {
  public:

    //! Constructor
    Random();
    
    //! Constructor
    Random(ScalarT a_, ScalarT b_);

    //! Constructor with seed value \c s
    Random(ScalarT a_, ScalarT b_, int s);

    //! Destructor
    ~Random();
    
    //! Set seed to \c s
    void setSeed(int s);

    //! Get random number
    ScalarT number();

  protected:

    // Check seed is valid
    int checkSeed(const std::string& func, int s);
  
  protected:

    //! Lower bound of interval
    ScalarT a;
    
    //! Upper bound of interval
    ScalarT b;

    //! %Random number seed
    ScalarT seed;

  }; // class Random

} // namespace Sacdo

#ifdef HAVE_SACADO_COMPLEX
#include <complex>

namespace Sacado {

  /*! 
   * \brief A partial specialization of Random that generates complex random 
   * numbers uniformly distributed in the box 
   * (a.real(),b.real())x(a.imag(),b.imag()).
   */
  template <typename T>
  class Random< std::complex<T> > {
  public:

    //! Constructor
    Random();
    
    //! Constructor
    Random(const std::complex<T>& a, const std::complex<T>& b);

    //! Constructor with seed value \c s
    Random(const std::complex<T>& a, const std::complex<T>& b, int s);

    //! Destructor
    ~Random();
    
    //! Set seed to \c s
    void setSeed(int s);

    //! Get random number
    std::complex<T> number();
  
  protected:

    //! Random number generator for real component
    Random<T> rand_real;
    
    //! Random number generator for imaginary component
    Random<T> rand_imag;

  }; // class Random

} // namespace Sacado

#endif // HAVE_SACADO_COMPLEX

#include "Sacado_RandomImp.hpp"

#endif // SACADO_RANDOM_HPP
