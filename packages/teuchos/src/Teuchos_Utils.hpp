#ifndef TEUCHOS_UTILS_H
#define TEUCHOS_UTILS_H

#include "Teuchos_ConfigDefs.hpp"

namespace Teuchos
{
  using std::string;

  /**\ingroup Utilities
   *
   * Numerical constants, etc.
   */
  class Utils
    {
    public:
      /**
       * print a description of the current build
       */
      static void aboutBuild();

      /**
       * Set a number to zero if it is less than chopVal.
       */
      static double chop(const double& x);

      /**
       * write a real as a string
       */
      static string toString(const double& x);

      /**
       * write an int as a string
       */
      static string toString(const int& x);

      /**
       * IEEE positive infinity
       */
      static double infinity() {return HUGE_VAL;}

      /**
       * IEEE negative infinity
       */
      static double negativeInfinity() {return -HUGE_VAL;}

      /**
       * pi.
       */
      static double pi() {return M_PI;}

      /**
       * Get the chopping value, below which numbers are considered to be zero
       */
      static double getChopVal() {return chopVal_;}
      /**
       * Set the chopping value, below which numbers are considered to be zero
       */
      static void setChopVal(double chopVal) {chopVal_ = chopVal;}

    private:
      static double chopVal_;
    };

  /** \relates Utils */
  inline string toString(const int& x) {return Utils::toString(x);}

  /** \relates Utils */
  inline string toString(const double& x) {return Utils::toString(x);}

  /** \relates Utils */
  inline string toString(const string& x) {return x;}

}

#endif


