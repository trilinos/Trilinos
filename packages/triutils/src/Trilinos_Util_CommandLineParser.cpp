// @HEADER
// ***********************************************************************
// 
//                 TriUtils: Trilinos Utilities Package
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/* ------------------------------------------------------------------------ */

#include "Trilinos_Util.h"
#include "Trilinos_Util_CommandLineParser.h"

Trilinos_Util_Map::Trilinos_Util_Map(void)
{

  SetLabel("Trilinos_Util_Map");
  
  return;
  
}

void Trilinos_Util_Map::Reset(void)
{

  SetLabel("");
  /*FIXME
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    ci->pop_back();
  }
*/
  return;
  
}

// ================================================ ====== ==== ==== == =

Trilinos_Util_CommandLineParser::Trilinos_Util_CommandLineParser(int argc, char *argv[])
{

  SetLabel("Trilinos_Util_CommandLineParser");

  char str[80];
  string value, param;

  Set("PROGRAM_NAME_",argv[0]);

  sprintf(str,"%d",argc);
  Set("_N_ARGS_",str);

  // first, manage possible arguments without specifier
  // (e.g., a.out 12 -nx 10 -ny 20)
  // Those arguments are called _ARGV_1_, ... , _ARGV_N_
  // and the value of N is given by _N_UNNAMED_ARGS_

  int N_args = 0;
  int i=1;
  
  for( i=1 ; i<argc ; ++i ) {
    if( *(argv[i]) == '-' ) break;
    N_args++;
    sprintf( str, "ARGV_%d", N_args );
    string param3;
    param3 = argv[i];
    Set(param3,value);
  }

  sprintf(str,"%d",N_args);
  Set("_N_UNNAMED_ARGS_",str);

  // now only arguments with a dash (possibly followed by one
  // other specifier)
  
  for( ; i<argc ; ++i ) {
    // check if the option has a `=' inside.
    // If so, split the string into two substrings
    char * pos = strchr( argv[i], '='); 
    if( pos != NULL ) {
      *pos = '\0';
      param = argv[i], value = pos+1;
      Set(param,value);
    } else if( i<argc-1 ) {
      if( *(argv[i+1]) != '-' ) {
	param = argv[i], value = argv[i+1];
	Set(param,value);
	++i;
      } else {
	param = argv[i], value = "";
	Set(param,value);
      }
    } else {
      param = argv[i], value = "";
      Set(param,value);
    }
    
  }

}

// ================================================ ====== ==== ==== == =

int Trilinos_Util_Map::Get( const string input, const int def_value)
{
  
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atoi(Map_[input].c_str()) );
  }
  
  return def_value;
   
}

// ================================================ ====== ==== ==== == =

double Trilinos_Util_Map::Get( const string input, const double def_value)
{

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atof(Map_[input].c_str()) );
  }
  
  return def_value;

}

// ================================================ ====== ==== ==== == =

string Trilinos_Util_Map::Get( const string input, const string def_value)
{

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( Map_[input] );
  }

  return def_value;
  
}

// ================================================ ====== ==== ==== == =

bool Trilinos_Util_Map::Has( const string input)
{
  
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return true;
  }
  return false;
  
}

// ================================================ ====== ==== ==== == =

void Trilinos_Util_Map::ShowAll() const 
{

  cout << "\n" << Label_ << " :: \n";
  
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first.at(0) != '_' ) 
      cout << (*ci).first << " = " << (*ci).second << endl;
  }
  
} /* ShowAll */

void Trilinos_Util_Map::ShowReallyAll() const 
{

  cout << "\nTrilinos_Util_CommandLineParser :: \n";

  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    cout << (*ci).first << " = " << (*ci).second << endl;
  }
  
} /* ShowReallyAll */

bool Trilinos_Util_Map::Add( const string input, const string value )
{

  // check that "input" has not been already inserted
  if( this->Has(input) == true )
    return false;

  Map_[input] = value;
  
  return true;

} /* AddOption */

bool Trilinos_Util_Map::Set( const string input, const double value )
{

  char value2[80];
  sprintf( value2, "%lf", value);
  return( Set(input,value2) );
}

bool Trilinos_Util_Map::Set( const string input, const int value )
{

  char value2[80];
  sprintf( value2, "%d", value);
  return( Set(input,value2) );
}

bool Trilinos_Util_Map::Set( const string input, const string value )
{

  Map_[input] = value;

  return true;

} /* Set */

bool Trilinos_Util_Map::Set( const string input, const char * value )
{

  string val(value);
  
  Map_[input] = val;

  return true;

} /* Set */

string Trilinos_Util_CommandLineParser::GetProgramName( void )
{
  return( Get("_PROGRAM_NAME_", "UNDEFINED" ) );
  
}

int Trilinos_Util_CommandLineParser::GetIntShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0;
  
} /* GetIntShellVariable */

double Trilinos_Util_CommandLineParser::GetDoubleShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0.0;
  
} /* GetDoubleShellVariable */

string Trilinos_Util_CommandLineParser::GetStringShellVariable( const char *str ) 
{

  char * buffer;

  buffer = getenv( str );
  if( buffer == NULL )
    return( "" );

  return( buffer );
  
} /* GetCharShellVariable */

// ================================================ ====== ==== ==== == =

ostream & operator << (ostream & os,
		       const Trilinos_Util_Map & S)
{
  S.ShowAll();
  return os;
}

// ================================================ ====== ==== ==== == =

Trilinos_Util_FileOptions::Trilinos_Util_FileOptions(const char FileName[]) :
  FileName_(FileName), CommentChars_("#"), SeparationChars_("="),
  FileHasBeenRead_(false)
{
}

Trilinos_Util_FileOptions::~Trilinos_Util_FileOptions() 
{
  
  FileName_ = "";
  CommentChars_ = "";
  SeparationChars_ = "";
  Reset();
  FileHasBeenRead_ = false;
  
}

string Trilinos_Util_FileOptions::GetFileName() const
{
  return FileName_;
}

void Trilinos_Util_FileOptions::SetCommentChars(const string c)
{
  CommentChars_ = c;
  return;
}

void Trilinos_Util_FileOptions::SetSeparationChars(const string c)
{
  SeparationChars_ = c;
  return;
}

#include <iostream>
#include <fstream>

int Trilinos_Util_FileOptions::ReadFile(const char * FileName) 
{
  FileName_ = FileName;

  return( ReadFile() );
}

int Trilinos_Util_FileOptions::ReadFile()
{
  
  ifstream File(FileName_.c_str());

  if( File.good() == false ) {
    std::cerr << "Error opening file `" << FileName_ << "'\n";
    return -1;
  }

  const int CharMax = 255;
  char line[CharMax];
  string Option, Value;

  while( File.eof() == false ) {
    
    File.getline(line,255);
    string StrLine = line;
    for( int k=0 ; k<(int)CommentChars_.length() ; ++k ) {
      int CommentPos = StrLine.find(CommentChars_.at(k));
      if( CommentPos != -1 ) {
	StrLine = StrLine.substr(0,CommentPos);
      }
    }
    int Length = StrLine.length();
    for( int k=0 ; k< (int) SeparationChars_.length() ; ++k ) {    
      int SepPos = StrLine.find(SeparationChars_.at(k));
      if( SepPos > 0 ) {
	Option = StrLine.substr(0,SepPos);
	Value = StrLine.substr(SepPos+1,Length);
	// ~!@ to erase spaces...
	if( Option.length() > 0 ) Set(Option,Value);
	break;
      }
    }
  }
  
  // close file
  File.close();

  return 0;
  
}


  

  
