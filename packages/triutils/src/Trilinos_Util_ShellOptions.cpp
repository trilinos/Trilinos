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
#include "Trilinos_Util_ShellOptions.h"

Trilinos_Util_Map::Trilinos_Util_Map(void)
{

  SetLabel("Trilinos_Util_Map");
  
#ifndef TRILINOS_UTIL_MAP_WITH_STL
  Allocated_ = 100;
  MapName_ = new string[Allocated_];
  MapValue_ = new string[Allocated_];
  NumEntries_ = 0;
#endif

  return;
  
}

void Trilinos_Util_Map::Reset(void)
{

  SetLabel("");
  
#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    ci->pop_back();
  }
#else
  NumEntries_ = 0;
#endif

  return;
  
}

// ================================================ ====== ==== ==== == =

Trilinos_Util_ShellOptions::Trilinos_Util_ShellOptions(int argc, char *argv[])
{

  SetLabel("Trilinos_Util_ShellOptions");

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

  sprintf(str,"%d%",N_args);
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

int Trilinos_Util_Map::GetInt( const string input )
{

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atoi(Map_[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      return( atoi(MapValue_[i].c_str()) );
    }
  }
#endif
  
  return 0;
  
} /* GetInt */

int Trilinos_Util_Map::GetInt( const string input, const int def_value)
{

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atoi(Map_[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      return( atoi(MapValue_[i].c_str()) );
    }
  }
#endif
  
  return def_value;
   
} /* GetInt */

double Trilinos_Util_Map::GetDouble( const string input)
{

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atof(Map_[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      return( atof(MapValue_[i].c_str()) );
    }
  }
#endif
  
  return 0.0;

} /* GetDoubleOption */

double Trilinos_Util_Map::GetDouble( const string input, const double def_value)
{

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( atof(Map_[input].c_str()) );
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      return( atof(MapValue_[i].c_str()) );
    }
  }
#endif
  
  return def_value;

} /* GetDoubleOption */

string Trilinos_Util_Map::GetString( const string input)
{

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( Map_[input] );
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      return( MapValue_[i] );
    }
  }
#endif
    
  return "";
  
} /* GetString */

string Trilinos_Util_Map::GetString( const string input, const string def_value)
{

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return( Map_[input] );
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      return( MapValue_[i] );
    }
  }
#endif
  
  return def_value;
  
} /* GetString */

bool Trilinos_Util_Map::Have( const string input)
{
  
#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first == input ) 
      return true;
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      return( true );
    }
  }
#endif
    
  return false;
  
} /* HaveOption */

void Trilinos_Util_Map::ShowAll() const 
{

  cout << "\n" << Label_ << " :: \n";
  
#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    if( (*ci).first.at(0) != '_' ) 
      cout << (*ci).first << " = " << (*ci).second << endl;
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i].at(0) != '_' ) 
      cout << MapName_[i] << " = " << MapValue_[i] << endl;
  }
#endif
  
} /* ShowAll */

void Trilinos_Util_Map::ShowReallyAll() const 
{

  cout << "\nTrilinos_Util_ShellOptions :: \n";

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  for( map<string,string>::const_iterator ci = Map_.begin();
       ci != Map_.end() ; ++ci ) {
    cout << (*ci).first << " = " << (*ci).second << endl;
  }
#else
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    cout << MapName_[i] << " = " << MapValue_[i] << endl;
  }
#endif
  
} /* ShowReallyAll */

bool Trilinos_Util_Map::Add( const string input, const string value )
{

  // check that "input" has not been already inserted
  if( this->Have(input) == true )
    return false;

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  Map_[input] = value;
#else
  MapName_[NumEntries_] = input;
  MapValue_[NumEntries_] = value;
  NumEntries_++;
  if( NumEntries_ == Allocated_ ) {
    string * MapName2 = new string[Allocated_*2];
    string * MapValue2 = new string[Allocated_*2];
    for( int i=0 ; i<Allocated_ ; i++ ) {
      MapName2[i] = MapName_[i];
      MapValue2[i] = MapValue_[i];
    }
    Allocated_ *= 2;
    delete [] MapName_;
    delete [] MapValue_;

    MapName_ = MapName2;
    MapValue_ = MapValue2;
  }
#endif
  
  return true;

} /* AddOption */

bool Trilinos_Util_Map::Set( const string input, const int value )
{

  char value2[80];
  sprintf( value2, "%d", value);
  return( Set(input,value2) );
}

bool Trilinos_Util_Map::Set( const string input, const string value )
{

#ifdef TRILINOS_UTIL_MAP_WITH_STL
  Map_[input] = value;
#else
  bool found = false;
  
  for( int i=0 ; i<NumEntries_ ; ++i ) {
    if( MapName_[i] == input ) {
      MapValue_[i] = value;
      found = true;
      break;
    }
  }
  if( found == false ) {
    MapName_[NumEntries_] = input;
    MapValue_[NumEntries_] = value;
    NumEntries_++;
  }
  
  if( NumEntries_ == Allocated_ ) {
    string * MapName2 = new string[Allocated_*2];
    string * MapValue2 = new string[Allocated_*2];
    for( int i=0 ; i<Allocated_ ; i++ ) {
      MapName2[i] = MapName_[i];
      MapValue2[i] = MapValue_[i];
    }
    Allocated_ *= 2;
    delete [] MapName_;
    delete [] MapValue_;

    MapName_ = MapName2;
    MapValue_ = MapValue2;
  }
#endif
  
  return true;

} /* Set */

string Trilinos_Util_ShellOptions::GetProgramName( void )
{
  return( GetString("_PROGRAM_NAME_", "UNDEFINED" ) );
  
}

int Trilinos_Util_ShellOptions::GetIntShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0;
  
} /* GetIntShellVariable */

double Trilinos_Util_ShellOptions::GetDoubleShellVariable( const char *str )
{

  char * buffer;

  buffer = getenv( str );
  if( buffer != NULL )
    return( atoi(buffer) );

  return 0.0;
  
} /* GetDoubleShellVariable */

string Trilinos_Util_ShellOptions::GetStringShellVariable( const char *str ) 
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


namespace Trilinos_Util_ShellOptions_Database
{
  Trilinos_Util_ShellOptions * Data = NULL;
}

void Trilinos_Util_ShellOptions_Set(int argc, char * argv[])
{
  
  if( Trilinos_Util_ShellOptions_Database::Data != NULL ) {
    delete Trilinos_Util_ShellOptions_Database::Data;
  }

  Trilinos_Util_ShellOptions_Database::Data = new Trilinos_Util_ShellOptions(argc,argv);

} /* Trilinos_Util_ShellOptions_Set */

int Trilinos_Util_ShellOptions_GetInt( const string input)
{

  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return 0;
  }

  return( Trilinos_Util_ShellOptions_Database::Data->GetInt(input) );
}

int Trilinos_Util_ShellOptions_GetInt( const string input, const int def_value)
{
  
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return 0;
  }

  return( Trilinos_Util_ShellOptions_Database::Data->GetInt(input,def_value) );
}

double Trilinos_Util_ShellOptions_GetDouble( const string input)
{
  
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return 0.0;
  }

  return( Trilinos_Util_ShellOptions_Database::Data->GetDouble(input) );
}

double Trilinos_Util_ShellOptions_GetDouble( const string input, const double def_value)
{
    
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return 0.0;
  }

  return( Trilinos_Util_ShellOptions_Database::Data->GetDouble(input,def_value) );
}

string Trilinos_Util_ShellOptions_GetString( const string input)
{
    
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return "";
  }

  return( Trilinos_Util_ShellOptions_Database::Data->GetString(input) );
}
  
string Trilinos_Util_ShellOptions_GetString( const string input, const string def_value )
{
    
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return "";
  }

  return( Trilinos_Util_ShellOptions_Database::Data->GetString(input,def_value) );
}
  
bool Trilinos_Util_ShellOptions_Set(const string input, const string value)
{
    
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return false;
  }

  return( Trilinos_Util_ShellOptions_Database::Data->Set(input,value) );
}

bool Trilinos_Util_ShellOptions_Set(const string input, const int value)
{
  
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return false;
  }

  return( Trilinos_Util_ShellOptions_Database::Data->Set(input,value) );
}

bool Trilinos_Util_ShellOptions_Have(const string input)
{
  
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return false;
  }
  
  return( Trilinos_Util_ShellOptions_Database::Data->Have(input) );
}

bool Trilinos_Util_ShellOptions_Add( const string input, const string value )
{
  
  if(Trilinos_Util_ShellOptions_Database:: Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return false;
  }
  
  return( Trilinos_Util_ShellOptions_Database::Data->Add(input,value) );
}

void Trilinos_Util_ShellOptions_ShowAll()
{
  
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return;
  }
  
  Trilinos_Util_ShellOptions_Database::Data->ShowAll();
}

void Trilinos_Util_ShellOptions_ShowReallyAll()
{
  
  if( Trilinos_Util_ShellOptions_Database::Data == NULL ) {
    cerr << "ERROR : Database never set!\n";
    return;
  }
  
  Trilinos_Util_ShellOptions_Database::Data->ShowReallyAll();
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
    for( int k=0 ; k<CommentChars_.length() ; ++k ) {
      int CommentPos = StrLine.find(CommentChars_.at(k));
      if( CommentPos != -1 ) {
	StrLine = StrLine.substr(0,CommentPos);
      }
    }
    int Length = StrLine.length();
    for( int k=0 ; k<SeparationChars_.length() ; ++k ) {    
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


  

  
