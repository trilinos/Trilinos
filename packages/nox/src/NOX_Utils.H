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

#ifndef NOX_UTILS_H
#define NOX_UTILS_H

#include "NOX_Common.H"

// Forward declarations
namespace NOX {
namespace Parameter {
class List;
}
}


namespace NOX {

/*!
  \brief Provides static parallel communication and printing utilities.

  The following parameters are used by this class and should be
  defined in the main solver parameter list.

  - "Output Information" - An integer which is a sum of MsgType's to
    specify how much information to show. Defaults to
    NOX::Utils::Warning + NOX::Utils::OuterIteration +
    NOX::Utils::InnerIteration + NOX::Utils::Parameters = 0xf = 15.

  - "Output Processor" - Specifies the designated print process. Defaults to 0.

  - "MyPID" - Specifies this process's ID. Defaults to 0.

  - "Output Precision" - Specifis the default number of decimal places
    to be used when printing floating point numbers. The default is 4.

  The public variables should never be modified directly.
*/
class Utils {


public:
  /*!
    \brief Message types for printing.

    Note that each message type if a power of two - this is \b very
    important. They can be added together to specify which messages
    you wish to receive.
  */
  enum MsgType { 
    Error = 0,			     // Errors are always printed
    Warning = 0x1,		     // 2^0
    OuterIteration = 0x2,	     // 2^1
    InnerIteration = 0x4,	     // 2^2
    Parameters = 0x8,		     // 2^3
    Details = 0x10,		     // 2^4
    OuterIterationStatusTest = 0x20, // 2^5
    LinearSolverDetails = 0x40        // 2^6
  };

public:

  //! Repeat the given character n times to the given stream.
  static ostream& repeat(ostream& stream, int n, char c = '*');

  /*! 
    \brief Set the parameters in the utilities class.
    
    Valid parameters are:
    <li> "Output Information" - Sum of MsgType flags to set what information
    to output.
    <li> "MyPID" - Processor ID.
    <li> "Output Processor" - Ranges from 0 to MyPID-1.
    <li> "Output Precision" - Precision of values to be printed by output.
  */
  static void setUtils(NOX::Parameter::List& p);

  /*! 
    \brief Returns true if the current processor is designated as the
    processor send output to screen and/or file.

    This is used to prevent each processor from printing the same
    information. 
  */
  static bool isPrintProc();

  /*! 
    \brief This function is deprecated.

    \deprecated Use NOX::Utils::doPrint(NOX::Utils::MsgType type)
    instead. This will print many warnings if you used it and is to be
    phased out.

  */ 
  static bool doPrint(int printLevel);

  /*!
    \brief Indicate whether or not the message should be printed.

    Returns true if the specifed message type should be printed \b and
    this is the print process.

    \note Be careful about putting any computations inside a doPrint
    if statement since \b only the print process will execute that
    code.

  */
  static bool doPrint(MsgType type);

  /*! 
    \brief This function is deprecated.

    \deprecated Use NOX::Utils::doAllPrint(NOX::Utils::MsgType type)
    instead. This will print many warnings if you used it and is to be
    phased out.

  */ 
  static bool doAllPrint(int printLevel);

  /*!
    \brief Indicate whether or not the message should be printed.

    Returns true if the specifed message type should be printed.
  */
  static bool doAllPrint(MsgType type);

  //! Returns the process ID number.
  static int getMyPID();

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  class Fill {
  public:
    Fill(int ntimes, char ch) : n(ntimes), c(ch) {};
    ~Fill() {};
    int n;
    char c;
  };

  class Sci {
  public:
    Sci(double val, int precision = -1) : d(val), p(precision) {};
    ~Sci() {};
    double d;
    int p;
  };

  //static Fill fillobj;
  //  static Sci sciobj;

#endif

  /*!  
    \brief Creates a Fill object which can be used in an output stream
    to insert an arbitrary number of a given character. The default
    character is an asterick.

    For example,
    \code cout << Utils::fill(10) << "or" << Utils::fill(10,'-'); \endcode

    This is modeled after the Form and Bound_form objects in
    Stroustrup, C++ Programming Langauge, 3rd ed., Chapter 21.4.
  */
  static Fill fill(int filln, char fillc = '*');

  /*! 

    \brief Creates a Sci object which can be used in an output stream
    for printing a double precision number in scientific format with
    an arbitrary precision. The default precision is that specificed
    by the Utils static object.

    For example,
    \code cout << Utils::sci(d) << "or" << Utils::sci(d,2); \endcode

    This is modeled after the Form and Bound_form objects in
    Stroustrup, C++ Programming Langauge, 3rd ed., Chapter 21.4.
  */
  static Sci sci(double dval, int prec = -1);

public:

  //! Output precision for floats and doubles.
  static int precision;
  
  //! Processor ID.
  static int myPID;
  
  //! Output test. 
  static int printTest;

  //! Processor responsible for output. 
  static int printProc;

}; // class Utils

} // namespace NOX

ostream& operator<<(ostream&, const NOX::Utils::Fill&);
ostream& operator<<(ostream&, const NOX::Utils::Sci&);


#endif
