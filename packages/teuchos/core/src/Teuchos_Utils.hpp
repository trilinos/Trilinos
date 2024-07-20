// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_UTILS_H
#define TEUCHOS_UTILS_H

/*! \file Teuchos_Utils.hpp
    \brief A utilities class for Teuchos
*/

#include "Teuchos_toString.hpp"

/*! \class Teuchos::Utils
    \brief This class provides some basic std::string and floating-point utilities for Teuchos
*/

namespace Teuchos
{
  using std::string;

  class TEUCHOSCORE_LIB_DLL_EXPORT Utils
    {
    public:

      /** \brief Set a number to zero if it is less than
       * <tt>getChopVal()</tt>. */
      static double chop(const double& x);

      /** \brief Get the chopping value, below which numbers are considered to
       * be zero. */
      static double getChopVal() {return chopVal_;}

      /** \brief Set the chopping value, below which numbers are considered to
       * be zero. */
      static void setChopVal(double chopVal) {chopVal_ = chopVal;}

      /** \brief Determine if a char is whitespace or not. */
      static bool isWhiteSpace( const char c )
        { return ( c==' ' || c =='\t' || c=='\n' ); }

      /** \brief Trim whitespace from beginning and end of std::string. */
      static std::string trimWhiteSpace( const std::string& str );

      /** \brief Write a double as a std::string. */
      static std::string toString(const double& x);

      /** \brief Write an int as a std::string. */
      static std::string toString(const int& x);

      /** \brief Write a long long as a std::string. */
      static std::string toString(const long long& x);

      /** \brief Write an unsigned int as a std::string. */
      static std::string toString(const unsigned int& x);

      /** \brief pi. */
#ifdef M_PI
      static double pi() {return M_PI;}
#else
      static double pi() {return 3.14159265358979323846;}
#endif

      /** \brief Get a parallel file name extention . */
      static std::string getParallelExtension(
        int    procRank = -1
        ,int   numProcs = -1
        );

    private:
      static double chopVal_;
    };


} // end namespace Teuchos

#endif


