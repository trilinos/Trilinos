#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Hashtable.hpp"

namespace Teuchos
{

  /**
   * \ingroup System
   * Utilities for extracting information from C++ command lines.
   *
   * The CommandLine class provides methods that can:
   * (*) Check for the presence of a string on the cmdline.
   * (*) Read one or more values (int, double, or string) after a
   * given key string.
   *
   * Ordering of the arguments is irrelevant (except, of course, that a
   * value must
   * always follow its key string).
   *
   * All methods are static, so that they can be accessed from anywhere without
   * the need to pass around a CommandLine object. The CommandLine system
   * must be initialized by calling the method CommandLine::init(argc, argv).
   *
   *
   * For example, suppose you have the command line:
   *
   \begin{verbatim}
   testCmd -flag -int 3 -x 3.14 -names Larry Moe Curly
   \end{verbatim}
   *
   * You want to check if "-flag" is present, and read the values after "-int",
   * "-x", and "-names".
   * Additionally, you want to see if "-flag" is present on the command line.
   * Of course, you will want to assign default values in case the keys
   * aren't set. You would do something like this:
   *
   \begin{verbatim}
   double x;
   int i;
   Array<string> names(3);
   if (CommandLine::find("-flag")) {...}
   if (!CommandLine::findInt("-int", i)) i = defaultInt;
   if (!CommandLine::findDouble("-x", x)) x = defaultX;
   if (!CommandLine::findstring("-names", names, 3)) names = defaultNames;
   \end{verbatim}


   * BUGS/FEATURES:
   *
   * If for some strange reason a key appears twice, e.g.
   *    myCode -pi 3.14 -pi 2.72
   * the second value wins.
   *
   * Everything crashes if you try to call init() twice. I'd rather have
   * a crash than the possibility of inconsistent results.
   *
   *  If you forget to call init() and then call one of the CommandLine
   * functions, everything crashes.
   *
   * HISTORY:
   * In the beginning there was chaos, and its name was getopt(). We saw that
   * it was bad, and created a replacement for it.
   *
   * Global variable based C version 12/90 by Kevin Long, UMass Amherst.
   * Static member based C++ version 9/96
   * Maps for faster lookups added 1/99
   *
   */

  class CommandLine
    {
    public:
      /** initialize with argc and argv. Must be called exactly once. */
      static void init(int argc, void** argv);

      /** check for the presence of a key string. */
      static bool find(const string& str) ;

      static bool find(const string& str, int& position) ;

      /** read the first value after the key string, and interpret as a double */
      static bool findDouble(const string& str, double& val) ;
      /** read the first value after the key string, and interpret as an int */
      static bool findInt(const string& str, int& val) ;
      /** read the first value after the key string, and interpret as a string */
      static bool findString(const string& str, string& val) ;

      /** read a list of values after the key string, and interpret as doubles */
      static bool findDouble(const string& str, Array<double>& val,
                             int count) ;
      /** read a list of values after the key string, and interpret as ints */
      static bool findInt(const string& str, Array<int>& val,
                          int count) ;
      /** read a list of values after the key string, and interpret as strings */
      static bool findString(const string& str, Array<string>& val,
                             int count) ;

      static void print();

      static bool& verbose() {static bool v = false; return v;}
    private:
      static void checkInitialization();

      static Array<string> tokens_;
      static Hashtable<string, int> tokenMap_;

      static bool frozen_; // set to true after setting of (argc, argv) to prevent resets
    };

  // utility: broadcast command-line args to all procs:



  //void bcastArgs(int sender, int& argc, void** argv);


}
#endif







