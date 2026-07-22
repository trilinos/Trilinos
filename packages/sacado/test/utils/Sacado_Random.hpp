// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
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
