// $Id$ 
// $Source$ 

// NOX: An Object-Oriented Nonlinear Solver Package
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#ifndef _NOX_UTILS_H
#define _NOX_UTILS_H
  
#include "NOX_Parameter_List.H"

namespace NOX {

//! Class to provide parallel communication and printing utilities.
class Utils {

public:

  //! Repeat the given character n times to the given stream.
  static ostream& repeat(ostream& stream, int n, char c = '*');

  //! Set the parameters in the utilities class.
  static void setUtils(Parameter::List& p);

  //! Returns true if the current processor is designated as the processor 
  //! send output to screen and/or file.
  /*  This is used to prevent each processor from printing the same 
   *  information. */
  static bool isPrintProc();

  //! Returns true if printLevel is less than or equal to the
  //! "Output Level" flag AND this is the designated printing process. 
  /*! This is used to determine if an output line should be printed.
   *  The user sets an outputLevel in the ParameterList for the MethodManager 
   *  If the print level is less than the output level, this function 
   *  returns true.  The most information is printed at outputLevel 5.  
   *  The least at level 0. */ 
  static bool doPrint(int printLevel);

  //! Returns true if printLevel is less than or equal to the
  //! "Output Level" flag.  
  /*! This is used to determine if an output line should be printed.
   *  The user sets an outputLevel in the ParameterList for the MethodManager 
   *  If the print level is less than the output level, this function 
   *  returns true.  The most information is printed at outputLevel 5.  
   *  The least at level 0. */ 
  static bool doAllPrint(int printLevel);

  //! Returns the processor ID number for parallel jobs. Returns 0 
  //! for serial jobs.
  /*! Returns the rank of calling process */
  static int getMyPID();

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  class Fill {
  public:
    Fill() {};
    ~Fill() {};
    int n;
    char c;
  };

  class Sci {
  public:
    Sci() {};
    ~Sci() {};
    double d;
  };

  static Fill fillobj;
  static Sci sciobj;

#endif

  /*! \brief We can use it to put in an arbitrary number of stars or
   *  any specified fill character.  cout << Utils::fill(10); cout <<
   *  Utils::fill(10,'-'); */

  static Fill& fill(int filln, char fillc = '*');

  /*! \brief We can use it to put a formated double onto the output
   *  with cout << sci(d); Note: There can only be ONE of these per line. 
   *  \todo Fix so that there can be more than one per line. */

  static Sci& sci(double dval);

public:
  

  //! Output precision for floats and doubles.
  static int precision;
  
  //! Processor ID.  Defaults to 0 for serial.
  static int myPID;
  
  //! Output Level. Defaults to level 2.
  static int outputLevel;

  //! Processor responsible for output. Defaults to processor 0.
  static int printProc;

};

}

ostream& operator<<(ostream&, const NOX::Utils::Fill&);
ostream& operator<<(ostream&, const NOX::Utils::Sci&);


#endif
