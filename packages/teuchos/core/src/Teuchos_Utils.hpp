// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
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

  class TEUCHOS_LIB_DLL_EXPORT Utils
    {
    public:

      /** \brief print a description of the current build. */
      static void aboutBuild();

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


