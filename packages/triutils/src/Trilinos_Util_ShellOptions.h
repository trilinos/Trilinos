#ifndef __TRILINOS_UTIL_SHELL_OPTIONS
#define __TRILINOS_UTIL_SHELL_OPTIONS

#include <iostream>
#include <string>
#include <map>

using namespace std;

class Trilinos_Util_ShellOptions {

public:

  // ============ //
  // constructors //
  // ------------ //
  
  Trilinos_Util_ShellOptions( void );
  Trilinos_Util_ShellOptions( int argc, char *argv[] );

  // ========== //
  // destructor //
  // ---------- //
  
   ~Trilinos_Util_ShellOptions() {}

  // === //
  // get //
  // --- //
  
  int    GetIntOption( const string input) ;
  int    GetIntOption( const string input, const int def_value) ;

  double GetDoubleOption( const string input) ;
  double GetDoubleOption( const string input, const double def_value) ;

  string GetStringOption( const string input) ;
  string GetStringOption( const string input, const string def_value ) ;

  // =========== //
  // set options //
  // ----------- //

  bool SetOption(const string input, const string value);
  bool SetOption(const string input, const int value);

  // ============= //
  // query and add //
  // ------------- //

  bool   HaveOption(const string input);
  string GetProgramName(void);
  bool   AddOption( const string input, const string value );

  // =============== //
  // shell variables //
  // --------------- //
  
  int    GetIntShellVariable( const char *str );
  double GetDoubleShellVariable( const char *str );
  string GetCharShellVariable( const char *str );

  // ====== //
  // others //
  // ------ //
  
  void    ShowAll() const;
  void    ShowReallyAll() const;

private:
  
  // map containing the arguments
  map<string,string> OptionDatabase;
  
};
#endif
