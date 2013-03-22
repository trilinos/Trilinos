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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef NOX_UTILS_H
#define NOX_UTILS_H

#include "NOX_Common.H"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_RCP.hpp"

// Forward declarations
namespace Teuchos {
  class ParameterList;
}

namespace NOX {

/*!
  \brief Provides printing utilities.

  The following parameters are valid for this class and should be
  defined in the "Printing" sublist of the solver parameter list.

  - "Output Information" - Can be a sublist or an integer.  If it is
  an integer, this value is a sum of MsgType's to specify how much
  information to show. Defaults to NOX::Utils::Warning +
  NOX::Utils::OuterIteration + NOX::Utils::InnerIteration +
  NOX::Utils::Parameters = 0xf = 15.  If this is a sublist, the
  following booleans are valid (set to true to enable print option):
       - "Error"
       - "Warning"
       - "Outer Iteration"
       - "Inner Iteration"
       - "Parameters"
       - "Details"
       - "Outer Iteration StatusTest"
       - "Linear Solver Details"
       - "Test Details"
       - "Stepper Iteration"
       - "Stepper Details"
       - "Stepper Parameters"
       - "Debug"

  - "Output Processor" - Specifies the designated print process. Defaults to 0.

  - "MyPID" - Specifies this process's ID. Defaults to 0.

  - "Output Precision" - Specifis the default number of decimal places
    to be used when printing floating point numbers. The default is 4.

  - "Output Stream" - A Teuchos::RCP<std::ostream> object that 
    will be used for standard output. 

  - "Error Stream" - A Teuchos::RCP<std::ostream> object that 
    will be used for output of error information. 

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

    Note that each message type is a power of two - this is \b very
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
    LinearSolverDetails = 0x40,   
    //! 2^7
    TestDetails = 0x80,
    //! 2^8 -- For LOCA
    StepperIteration = 0x0100,
    //! 2^9 -- For LOCA
    StepperDetails = 0x0200,
    //! 2^10 -- For LOCA
    StepperParameters = 0x0400,
    //! 2^12
    Debug = 0x01000
  };


public:		

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

public:

  //! Constructor.  
  /*!The final two arguments are a reference counted pointers to
    std::ostreams.  This defaults to std::cout and std::cerr if not
    supplied.  Users should only supply this argument if they want
    output directed to a std::ostream other than the defaults.  If so,
    they must wrap the std::ostream in a reference counted pointer for safe
    memory management. */
  Utils(int outputInformation=0xf, int MyPID=0, int outputProcess=0, 
	int outputPrecision=3, 
	const Teuchos::RCP<std::ostream>& outputStream=
	Teuchos::null,
	const Teuchos::RCP<std::ostream>& errorStream=
	Teuchos::null);

  //! Constructor via a parameter list
  Utils(Teuchos::ParameterList& p);

  //! Copy constructor
  Utils(const NOX::Utils& u);
  
  //! Destructor
  ~Utils();

  //! Assignment operator
  NOX::Utils& operator=(const NOX::Utils& source);

  //! Reset the printing parameters
  void reset(Teuchos::ParameterList& p);

  //! Returns true if this is a valid print type
  bool isPrintType(NOX::Utils::MsgType type) const;

  //! Returns the std::ostream for printing if this proces is the print process. Returns a Teuchos::oblackholestream otherwise.
  std::ostream& out() const;

  //! Returns the std::ostream for printing if this process is the print process and the print type is valid. Returns a Teuchos::oblackholestream otherwise.
  std::ostream& out(NOX::Utils::MsgType type) const;

  //! Returns the std::ostream for printing regardless of the print processor.  Only use this call if you want all processes to print to the std::ostream.
  std::ostream& pout() const; 

  //! Returns the std::ostream for printing if the print type matches.  Only use this call if you want all processes to print to the std::ostream for the print type.
  std::ostream& pout(NOX::Utils::MsgType type) const; 

  //! Returns the error stream for printing if this is the print process.
  std::ostream& err() const;

  //! Returns the error stream for printing to all processors.  Only use this call if you want all processes to print to the error stream.
  std::ostream& perr() const;
  
  //! Print this object
  void print(std::ostream& os) const;

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

  //! Output precision for floats and doubles.
  int precision;
  
  //! Processor ID.
  int myPID;
  
  //! Output test. 
  int printTest;

  //! Processor responsible for output. 
  int printProc;

  //! std::ostream derived object that deletes the output.
  /*! This is useful for parallel printing to eliminate errors due to printing to all processors.
   */
  Teuchos::RCP<std::ostream> blackholeStream;

  //! Output stream to be used for printing.
  Teuchos::RCP<std::ostream> printStream;

  //! Output stream used for printing to the print process.
  Teuchos::RCP<std::ostream> myStream;

  //! Output stream used for printing to the print process.
  Teuchos::RCP<std::ostream> errorStream;

}; // class Utils

std::ostream& operator<<(std::ostream&, const NOX::Utils::Fill&);
std::ostream& operator<<(std::ostream&, const NOX::Utils::Sci&);
std::ostream& operator<<(std::ostream&, const NOX::Utils& utils);

} // namespace NOX


#endif
