//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

/* ======================================================================== */
/*!
 \class ShellOptions

 \brief ShellOptions: a class to manage the input arguments and shell variables.

 With this class, it is easy to handle input line arguments and shell
 varibles. For instance, the user can write
 \verbatim
 $ ./a.out -nx 10 -tol 1e-6 -solver=cg
 \endverbatim
 and then easily retrive the value of \c nx, \c tol, and \c solver.
 
 A simple code using this class is as follows:
 \verbatim
 int main(int argc, char *argv[])
  {

   ShellOptions Args(argc,argv);
   int nx = Args.GetIntOption("-nx", 123);
   int ny = Args.GetIntOption("-ny", 145);
   double tol = Args.GetDoubleOption("-tol", 1e-12);
   string solver = Args.GetIntOption("-solver");

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
 - GetIntOption
 - GetDoubleOption
 - GetStringOption
 
 If option name is not found in the database, a value of 0, 0.0 or an
 empty string is returned. If needed, the user can also specify a
 default value to return when the option name is not found in the
 database. Method \c HaveOption can be used to query the database for
 an option.

 The user can modify the database as well, using
 - SetOption
 - AddOption

 Finally, the user can retrive the integer, double or string value of a
 shell environmental variable using:
 - GetIntShellVariable
 - GetDoubleShellVariable
 - GetStringShellVariable
 
 \date Albuquerque, 01.-Oct-03
 
 \author Marzio Sala, SNL 9214

*/
/* ------------------------------------------------------------------------ */

#include "Trilinos_Util.h"
#include "Trilinos_Util_ShellOptions.h"

/* ======================================================================== */
/*!
 \brief Initialize the database using the options given at the shell line.

*/
/* ------------------------------------------------------------------------ */

Trilinos_Util_ShellOptions::Trilinos_Util_ShellOptions(int argc, char *argv[])
{

#ifndef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  NumOptions = 0;
#endif

  string param;
  string value;
  param = "_PROGRAM_NAME_";  value = argv[0];
  SetOption(param,value);

  param = "_N_ARGS_"; value = argc;
  SetOption(param,value);
  
  // first, manage possible arguments without specifier
  // (e.g., a.out 12 -nx 10 -ny 20)
  // Those arguments are called _ARGV_1_, ... , _ARGV_N_
  // and the value of N is given by _N_UNNAMED_ARGS_

  int N_args = 0;
  char str[80];
  int i=1;
  
  for( i=1 ; i<argc ; ++i ) {
    if( *(argv[i]) == '-' ) break;
    N_args++;
    sprintf( str, "ARGV_%d", N_args );
    param = argv[i];
    SetOption(param,value);
  }

  param = "_N_UNNAMED_ARGS_", value = N_args;
  SetOption(param,value);
  
  // now only arguments with a dash (possibly followed by one
  // other specifier)
  
  for( ; i<argc ; ++i ) {
    // check if the option has a `=' inside.
    // If so, split the string into two substrings
    char * pos = strchr( argv[i], '='); 
    if( pos != NULL ) {
      *pos = '\0';
      param = argv[i], value = pos+1;
      SetOption(param,value);
    } else if( i<argc-1 ) {
      if( *(argv[i+1]) != '-' ) {
	param = argv[i], value = argv[i+1];
	SetOption(param,value);
	++i;
      } else {
	param = argv[i], value = "";
	SetOption(param,value);
      }
    }
    
  }

}

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as an integer

 This method returns the integer value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns 0.

*/
/* ------------------------------------------------------------------------ */

int Trilinos_Util_ShellOptions::GetIntOption( const string input )
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atoi(OptionDatabase[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      return( atoi(OptionValue[i].c_str()) );
    }
  }
#endif
  
  return 0;
  
} /* GetIntOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as an integer

 This method returns the integer value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns the specified default value.

*/
/* ------------------------------------------------------------------------ */

int Trilinos_Util_ShellOptions::GetIntOption( const string input, const int def_value)
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atoi(OptionDatabase[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      return( atoi(OptionValue[i].c_str()) );
    }
  }
#endif
  
  return def_value;
   
} /* GetIntOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a double.

 This method returns the double value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns 0.0.

*/
/* ------------------------------------------------------------------------ */

double Trilinos_Util_ShellOptions::GetDoubleOption( const string input)
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atof(OptionDatabase[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      return( atof(OptionValue[i].c_str()) );
    }
  }
#endif
  
  return 0.0;

} /* GetDoubleOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a double.

 This method returns the double value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns the specified default value.

*/
/* ------------------------------------------------------------------------ */

double Trilinos_Util_ShellOptions::GetDoubleOption( const string input, const double def_value)
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atof(OptionDatabase[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      return( atof(OptionValue[i].c_str()) );
    }
  }
#endif
  
  return def_value;

} /* GetDoubleOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a C++ string.

 This method returns the string value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns an empty string ("").

*/
/* ------------------------------------------------------------------------ */

string Trilinos_Util_ShellOptions::GetStringOption( const string input)
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( OptionDatabase[input] );
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      return( OptionValue[i] );
    }
  }
#endif
    
  return "";
  
} /* GetStringOption */

/* ======================================================================== */
/*!
 \brief Get the value of the specified option as a C++ string.

 This method returns the string value assigned to option \c input.
 If \c input is not in the database, or it cannot be converted to an integer,
 returns the default value.

*/
/* ------------------------------------------------------------------------ */

string Trilinos_Util_ShellOptions::GetStringOption( const string input, const string def_value)
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( OptionDatabase[input] );
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      return( OptionValue[i] );
    }
  }
#endif
  
  return def_value;
  
} /* GetStringOption */

/* ======================================================================== */
/*!
 \brief Check wheter an option is in the database or not

 This method checks whether option \c input is in the databse or not.
 It returns \c true if it is, \c false otherwise.

*/
/* ------------------------------------------------------------------------ */

bool Trilinos_Util_ShellOptions::HaveOption( const string input)
{
  
#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return true;
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      return( true );
    }
  }
#endif
    
  return false;
  
} /* HaveOption */

/* ======================================================================== */
/*!
 \brief Show all the databse entries

*/
/* ------------------------------------------------------------------------ */

void Trilinos_Util_ShellOptions::ShowAll() const 
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  cout << "\nTrilinos_Util_ShellOptions :: \n";
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    if( (*ci).first.at(0) != '_' ) 
      cout << (*ci).first << " = " << (*ci).second << endl;
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i].at(0) != '_' ) 
      cout << OptionName[i] << " = " << OptionValue[i] << endl;
  }
#endif
  
} /* ShowAll */

/* ======================================================================== */
/*!
 \brief Show all the databse entries

*/
/* ------------------------------------------------------------------------ */

void Trilinos_Util_ShellOptions::ShowReallyAll() const 
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  cout << "\nTrilinos_Util_ShellOptions :: \n";
  for( map<string,string>::const_iterator ci = OptionDatabase.begin();
       ci != OptionDatabase.end() ; ++ci ) {
    cout << (*ci).first << " = " << (*ci).second << endl;
  }
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    cout << OptionName[i] << " = " << OptionValue[i] << endl;
  }
#endif
  
} /* ShowReallyAll */

/* ======================================================================== */
/*!
 \brief Add an entry to the databse

 This method add an entry to the databse. First, it checks that this
 entry does not exist. If it exists, the method returns \c
 false. Otherwise, it adds the entry and returns \c true.
 
*/
/* ------------------------------------------------------------------------ */

bool Trilinos_Util_ShellOptions::AddOption( const string input, const string value )
{

  // check that "input" has not been already inserted
  if( this->HaveOption(input) == true )
    return false;

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  OptionDatabase[input] = value;
#else
  OptionName[NumOptions] = input;
  OptionValue[NumOptions] = value;
  NumOptions++;
#endif
  
  return true;

} /* AddOption */

/* ======================================================================== */
/*!
 \brief Modify the value of a database entry.

 This method modifies the value of a database entry. If the entry does
 not exist in the database, return \c false. Otherwise, returns \c true.
 
*/
/* ------------------------------------------------------------------------ */

bool Trilinos_Util_ShellOptions::SetOption( const string input, const string value )
{

#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  OptionDatabase[input] = value;
#else
  bool found = false;
  
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == input ) {
      OptionValue[NumOptions] = value;
      found = true;
      break;
    }
  }
  if( found == false ) {
    OptionName[NumOptions] = input;
    OptionValue[NumOptions] = value;
    NumOptions++;
  }
#endif
  
  return true;

} /* SetOption */

/* ======================================================================== */
/*!
 \brief Returns the name of the program as a C++ string.

*/
/* ------------------------------------------------------------------------ */

string  Trilinos_Util_ShellOptions::GetProgramName( void )
{
#ifdef TRILINOS_UTIL_SHELL_OPTIONS_WITH_STL
  return OptionDatabase["_PROGRAM_NAME_"];
#else
  for( int i=0 ; i<NumOptions ; ++i ) {
    if( OptionName[i] == "_PROGRAM_NAME_" ) 
      return OptionValue[i];
  }
#endif
  
}

/* ======================================================================== */
/*!
  
 \brief Returns the value of the environmenta variable \c str as an integer.

 This methods returns the value of the environmenta variable \c
 str. If the variable does not exists, returns \c 0.

*/
/* ------------------------------------------------------------------------ */

int Trilinos_Util_ShellOptions::GetIntShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0;
  
} /* GetIntShellVariable */

/* ======================================================================== */
/*!
  
 \brief Returns the value of the environmenta variable \c str as a double.
 
 This methods returns the value of the environmenta variable \c
 str. If the variable does not exists, returns \c 0.0.

*/
/* ------------------------------------------------------------------------ */

double Trilinos_Util_ShellOptions::GetDoubleShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0.0;
  
} /* GetDoubleShellVariable */

/* ======================================================================== */
/*!
  
 \brief Returns the value of the environmenta variable \c str as a C++ string.

 This methods returns the value of the environmenta variable \c
 str. If the variable does not exists, returns \c "".

*/
/* ------------------------------------------------------------------------ */

string Trilinos_Util_ShellOptions::GetCharShellVariable( const char *str ) 
{

  char * buffer;

  buffer = getenv( str );
  if( buffer == NULL )
    return( "" );

  return( buffer );
  
} /* GetCharShellVariable */

