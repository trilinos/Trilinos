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
  \brief Provides parallel communication and printing utilities.

  The following parameters are used by this class and should be
  defined in the "Printing" sublist of the solver parameter list.

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

public:				// types

  //! Fill object - used to print the given character the number of times specified.
  class Fill {
  public:
    //! Constructor
    Fill(int ntimes, char ch) : n(ntimes), c(ch) {};
    //! Destructor
    ~Fill() {};
    //! Number of times the character should be printed
    int n;
    //! Character to be printed
    char c;
  };

  //! Sci object - used to print the given value with the specified precision
  class Sci {
  public:
    //! Constructor
    Sci(double val, int precision = -1) : d(val), p(precision) {};
    //! Destructor
    ~Sci() {};
    //! Value to be printed
    double d;
    //! Precision 
    int p;
  };


  /*!
    \brief Message types for printing.

    Note that each message type if a power of two - this is \b very
    important. They can be added together to specify which messages
    you wish to receive.
  */
  enum MsgType { 
    //! Errors are always printed
    Error = 0,			     
    //! 2^0
    Warning = 0x1,		     
    //! 2^1
    OuterIteration = 0x2,	     
    //! 2^2
    InnerIteration = 0x4,	     
    //! 2^3
    Parameters = 0x8,		     
    //! 2^4
    Details = 0x10,		     
    //! 2^5
    OuterIterationStatusTest = 0x20, 
    //! 2^6
    LinearSolverDetails = 0x40        
  };


public:				// static functions

  /*! 
    \brief Deprecated
    \deprecated Use Utils() or reset() on instantiated object instead.
  */
  static void setUtils(NOX::Parameter::List& p);

  /*! 
    \brief Deprecated
    \deprecated Use isPrintProcess() on instantiated object instead.
  */
  static bool isPrintProc();

  /*! 
    \brief Deprecated.
    \deprecated Use isPrintProcess() on instantiated object instead.
  */ 
  static bool doPrint(int printLevel);

  /*!
    \brief Deprecated
    \deprecated Use isPrintProcessAndCode() on instantiated object instead.
  */
  static bool doPrint(MsgType type);

  /*! 
    \brief Deprecated
    \deprecated Use isPrintProcessAndCode() on instantiated object instead.
  */ 
  static bool doAllPrint(int printLevel);

  /*!
    \brief Deprecated
    \deprecated Use isPrintCode() on instantiated object instead.
  */
  static bool doAllPrint(MsgType type);

  /*!
    \brief Deprecated
    \deprecated No replacement.
  */
  static int getMyPID();


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
    \brief Deprecated
    \deprecated Use sciformat() on instantiated object instead.
  */
  static Sci sci(double dval, int prec = -1);

public:				// functions

  //! Empty constructor
  Utils();

  //! Constructor
  explicit Utils(NOX::Parameter::List& p);

  //! Copy constructor
  Utils(const NOX::Utils& u);
  
  //! Destructor
  ~Utils();

  //! Reset the utility parameters
  void reset(NOX::Parameter::List& p);

  //! Returns true if this is the designated print process
  bool isPrintProcess() const;

  //! Returns true if this is a valid print type
  bool isPrintType(NOX::Utils::MsgType type) const;

  //! Returns true if this is the designated print process and this is a valid print type
  bool isPrintProcessAndType(NOX::Utils::MsgType type) const;

  /*! 
    \brief Creates a Sci object which can be used in an output stream
    for printing a double precision number in scientific format with
    an arbitrary precision. The precision is that specificed
    by the Utils object.

    For example,
    \code cout << Utils::sci(d) << "or" << Utils::sci(d,2); \endcode

    This is modeled after the Form and Bound_form objects in
    Stroustrup, C++ Programming Langauge, 3rd ed., Chapter 21.4.
  */
  Sci sciformat(double dval) const;

  /*!
    \brief Static version of sciformat().
  */
  static Sci sciformat(double dval, int prec);

private:

  //! Print out a message regarding being deprecated
  static void deprecated(const string& oldname, const string& newname);

private:

  //! Output precision for floats and doubles.
  int precision;
  
  //! Processor ID.
  int myPID;
  
  //! Output test. 
  int printTest;

  //! Processor responsible for output. 
  int printProc;

}; // class Utils

} // namespace NOX

ostream& operator<<(ostream&, const NOX::Utils::Fill&);
ostream& operator<<(ostream&, const NOX::Utils::Sci&);


#endif
