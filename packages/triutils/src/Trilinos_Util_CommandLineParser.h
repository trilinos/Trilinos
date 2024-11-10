// @HEADER
// ***********************************************************************
//
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2011) Sandia Corporation
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

#ifndef _TRILINOS_UTIL_CLP_
#define _TRILINOS_UTIL_CLP_

#if defined(Triutils_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Triutils package is deprecated"
#endif
#endif

#include <iostream>
#include <string>
#include <map>

class Trilinos_Util_Map {

public:


  Trilinos_Util_Map();

  virtual ~Trilinos_Util_Map()
  {  }

  //@}

  //@{ \name Insertion/Removal methods.

  //! Gets the value of the specified option as an integer. If not found, returns the specified default value.
  virtual int    Get( const std::string input, const int def_value) ;

  //! Gets the value of the specified option as a double. If not found, returns the specified default value.
  virtual double Get( const std::string input, const double def_value) ;

 //! Gets the value of the specified option as a std::string. If not found, returns the specified default value.
  virtual std::string Get( const std::string input, const std::string def_value ) ;

  //! Modify the value of a database entry.
  /*!  This method modifies the value of a database entry. If the entry
    does not exist in the database, return \c false. Otherwise, returns
    \c true.
  */
  virtual bool Set(const std::string input, const char * value);
  virtual bool Set(const std::string input, const std::string value);
  virtual bool Set(const std::string input, const int value);
  virtual bool Set(const std::string input, const double value);

  //! Add an entry to the databse
  /*!  This method add an entry to the databse. First, it checks that
    this entry does not exist. If it exists, the method returns \c
    false. Otherwise, it adds the entry and returns \c true.
  */
  virtual bool   Add( const std::string input, const std::string value );

  inline bool SetLabel(std::string Label) {
    Label_ = Label;
    return true;
  }

  inline std::string GetLabel(std::string /* Label */) {
    return( Label_ );
  }

  //@}

  //@{ \name Query methods.

  /*! \brief Check wheter an option is in the database or not

  This method checks whether option \c input is in the databse or not.
  It returns \c true if it is, \c false otherwise.
  */
  virtual bool Has(const std::string input);

  //@}

  //@{ \name Miscellaneous methods.

  //! Show all the databse entries
  virtual void    ShowAll() const;

  //! Show all the databse entries, including entries beginning with "_"
  virtual void    ShowReallyAll() const;

  virtual void Reset(void);

  //@}

  //@{ \name Friend functions .

  friend std::ostream & operator << (std::ostream & os,
        const Trilinos_Util_Map & S);

  //@}

private:

  std::string Label_;

  // map containing the arguments
  std::map<std::string,std::string> Map_;

};

//! Trilinos_Util::CommandLineParser: A class for managing the input arguments and variables.

/* ======================================================================== */
/*!
 \class Trilinos_Util::CommandLineParser

 Using Trilinos_Util::CommandLineParser, it is easy to handle input line arguments and shell
 variables. For instance, the user can write
 \verbatim
 $ ./a.out -nx 10 -tol 1e-6 -solver=cg
 \endverbatim
 and then easily retrive the value of \c nx, \c tol, and \c solver.

 A simple code using this class is as follows:
 \verbatim
 int main(int argc, char *argv[])
  {

   Trilinos_Util::CommandLineParser CLP(argc,argv);
   int nx = CLP.GetInt("-nx", 123);
   int ny = CLP.GetInt("-ny", 145);
   double tol = CLP.GetDouble("-tol", 1e-12);
   std::string solver = CLP.GetInt("-solver");

   cout << "nx = " << nx << endl;
   cout << "ny = " << ny << " (default value)" << endl;
   cout << "tol = " << tol << endl;
   cout << "solver = " << solver << endl;

   return 0;

 }
 \endverbatim

 Each line option can have a value or not. For options with a value,
 the user can specify this values as follows. Let \c -tolerance be the
 name of the option and \c 1e-12 its value. Both choices are valid:
 - \c -option \c 1e-12 (with one or more spaces)
 - \c -option=1e-12 (an `=' sign and no spaces)

 Options are indentified with one or more dashes (`-'). Each option
 cannot have more than one value.

 Note that the user can specify some values without giving them a name.
 This can be done as follows:
  \verbatim
  $ ./a.out value1 value2 value 3 -nx 10 -tol 1e-6 -solver=cg
  \endverbatim
 Here, \c valueX, (X=1,...,9) is stored in the database entry
 \c ARGV_X.

 To use this class, the user has to build the database using the \c
 argc,argv input arguments. Then, to retrive the option value, the user
 has to use one of the following functions:
 - GetInt
 - GetDouble
 - GetString

 If option name is not found in the database, a value of 0, 0.0 or an
 empty std::string is returned. If needed, the user can also specify a
 default value to return when the option name is not found in the
 database. Method \c HaveOption can be used to query the database for
 an option.

 The user can modify the database as well, using
 - Set
 - Add

 (GetInt, GetDouble, GetString, Set and Add are derived from the base class,
 Trilinos_Util_Map).

 Finally, the user can retrive the integer, double or std::string value of a
 shell environmental variable using:
 - GetIntShellVariable
 - GetDoubleShellVariable
 - GetStringShellVariable

 \date Albuquerque, 19-Jan-04

 \author Marzio Sala, SNL 9214

*/
// ================================================ ====== ==== ==== == =

namespace Trilinos_Util {

class CommandLineParser : public Trilinos_Util_Map
{

 public:

  //@{ \name Constructors/Destructor.
  //! Trilinos_Util_ShellOptions constructor using the options given at the shell line.
  CommandLineParser(int argc, char *argv[] );

  // ============= //
  // query and add //
  // ------------- //

  //@}

  //@{ \name Query methods.

  //!  Returns the name of the program as a C++ std::string.
  virtual std::string GetProgramName(void);

  // =============== //
  // shell variables //
  // --------------- //

  //!Returns the value of the environmenta variable \c str as an integer.
  /*! This methods returns the value of the environmental variable \c
    str. If the variable does not exists, returns \c 0. */
  virtual int    GetIntShellVariable( const char *str );

  //!Returns the value of the environmenta variable \c str as an double.
  /*! This methods returns the value of the environmenta variable \c
    str. If the variable does not exists, returns \c 0.0. */
  virtual double GetDoubleShellVariable( const char *str );

  //!Returns the value of the environmenta variable \c str as a C++ std::string.
  /*! This methods returns the value of the environmenta variable \c
    str. If the variable does not exists, returns \c "". */
  virtual std::string GetStringShellVariable( const char *str );

  //@}

};

// ================================================ ====== ==== ==== == =

class InputFileReader : public Trilinos_Util_Map
{

 public:

  InputFileReader(const char FileName[] );
  ~InputFileReader();

  virtual std::string GetFileName(void) const;

  //  virtual bool SetCommentChar(char c);
  virtual void SetCommentChars(const std::string c);

  //  virtual bool SetSeparationChar(char c);
  virtual void SetSeparationChars(const std::string c);

  virtual int ReadFile();
  virtual int ReadFile(const char FileName[]);

private:
  std::string FileName_;
  std::string CommentChars_;
  std::string SeparationChars_;
  bool FileHasBeenRead_;
};

} // namespace Trilinos_Util
#endif
