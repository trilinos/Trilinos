// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef TEUCHOS_COMMAND_LINE_PROCESSOR_HPP
#define TEUCHOS_COMMAND_LINE_PROCESSOR_HPP

/*! \file Teuchos_CommandLineProcessor.hpp
  \brief Basic command line parser for input from <tt>(argc,argv[])</tt>
*/

/** \example CommandLineProcessor/cxx_main.cpp
    This is an example of how to use the Teuchos::CommandLineProcessor class.
*/

#include "Teuchos_map.hpp"
#include "Teuchos_any.hpp"
#include "Teuchos_CompileTimeAssert.hpp"
#include "Teuchos_Ptr.hpp"
#include <vector>

/*! \class Teuchos::CommandLineProcessor
 * \brief Class that helps parse command line input arguments from <tt>(argc,argv[])</tt> and set options.
 *
 * This class will process command-line arguments in the form of <tt>(argc,argv[])</tt>
 * and set user-defined options.  This class can also work in a number of modes.
 * This processor can require that all options be recognized or not.
 *
 * This class object will also setup the behavior of
 * <tt>Teuchos::VerboseObjectBase::getDefaultOStream()</tt> if
 * <tt>this->addOutputSetupOptions()==true</tt>
 *
 * Warning, the option --show-timer-summary is only enabled if the subpackage
 * TeuchosComm is enabled!
 */

namespace Teuchos {

class TEUCHOSCORE_LIB_DLL_EXPORT CommandLineProcessor {
public:

  //! @name Public types
  //@{

  /// Thrown if a parse std::exception occurs and  throwExceptions==true
  class ParseError : public std::logic_error
  {public: ParseError(const std::string& what_arg) : std::logic_error(what_arg) {}};

  /// Thrown if --help was specified and throwExceptions==true
  class HelpPrinted : public ParseError
  {public: HelpPrinted(const std::string& what_arg) : ParseError(what_arg) {}};

  /// Thrown if an unrecognized option was found and throwExceptions==true
  class UnrecognizedOption : public ParseError
  {public: UnrecognizedOption(const std::string& what_arg) : ParseError(what_arg) {}};

  /** \brief Return value for <tt>CommandLineProcessor::parse()</tt>.
      Note: These enums are all given non-negative values since they are designed to
      be returned from main().
   */
  enum EParseCommandLineReturn {
    PARSE_SUCCESSFUL              =  0 /*!< Parsing the command line was successful. */
    ,PARSE_HELP_PRINTED            =  1 /*!< The help statement was printed for the command line parser. */
    ,PARSE_UNRECOGNIZED_OPTION     =  2 /*!< The command line parser encountered an unrecognized option. */
    ,PARSE_ERROR                   =  3 /*!< The command line parser encountered an error. */
  };

  //@}

  //! @name Constructors
  //@{

  /** \brief Default Constructor
   *
   * \param  throwExceptions
   *               [in] If <tt>true</tt> then <tt>this->parse()</tt> will throw
   *               exceptions instead of returning <tt>!=PARSE_SUCCESSFUL</tt>.
   * \param  recogniseAllOptions
   *               [in] If <tt>true</tt> then <tt>this->parse()</tt> will return
   *               the appropriate error for any option it does not recognize.
   *               If <tt>false</tt>, then <tt>this->parse()</tt> will simply
   *               ignore options that it does not recognize.  However, a warning will be
   *               printed for any unrecognized option, but no errors will be returned.
   * \param  addOutputSetupOptions
   *               [in] If <tt>true</tt> then options will be automatically added
   *               to setup <tt>Teuchos::VerboseObjectBase::getDefaultOStream()</tt>.
   */
  CommandLineProcessor(
    bool    throwExceptions       = true
    ,bool   recogniseAllOptions   = true
    ,bool   addOutputSetupOptions = false
    );

  /** \brief Destructor.
   */
  ~CommandLineProcessor();

  //@}

  //! @name Behavior modes
  //@{

  /// Set if an std::exception is thrown, there is a parse error, or help is printed.
  void throwExceptions( const bool & throwExceptions );

  /// Returns true if an std::exception is thrown, there is a parse error, or help is printed.
  bool throwExceptions() const;

  /// Set if all options must be recognized or not.
  void recogniseAllOptions( const bool & recogniseAllOptions );

  /// Returns true if all options must be recognized by the parser.
  bool recogniseAllOptions() const;

  /// Set if options will be automatically added to setup <tt>Teuchos::VerboseObjectBase::getDefaultOStream()</tt>.
  void addOutputSetupOptions( const bool &addOutputSetupOptions );

  /// Returns true options will be automatically added to setup <tt>Teuchos::VerboseObjectBase::getDefaultOStream()</tt>.
  bool addOutputSetupOptions() const;

  //@}

  //! @name Set up options
  //@{

  /** \brief Set a documentation sting for the entire program printed when
   * --help is specified. */
  void setDocString( const char doc_string[] );

  /** \brief Set a boolean option.
   *
   * \param  option_true    [in] (null terminated std::string) If this option is found then
   *                        <tt>*option_val = true</tt> will be set.
   * \param  option_false   [in] (null terminated std::string) If this option is found then
   *                        <tt>*option_val = false</tt> will be set.
   * \param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
   *                        of the option (used for printing in --help).  On output,
   *                        will be set according to <tt>(argc,argv[])</tt>.
   * \param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated std::string
   *                        gives the documentation for the option.
   */
  void setOption(
    const char     option_true[]
    ,const char    option_false[]
    ,bool          *option_val
    ,const char    documentation[] = NULL
    );

  /** \brief Set an integer option.
   *
   * \param  option_name    [in] (null terminated std::string) The name of the option
   *                        (without the leading '--' or trailing '=').
   * \param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
   *                        of the option (used for printing in --help).  On output,
   *                        will be set according to <tt>(argc,argv[])</tt>.
   * \param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated std::string
   *                        gives the documentation for the option.
   */
  void setOption(
    const char     option_name[]
    ,int           *option_val
    ,const char    documentation[] = NULL
    ,const bool    required        = false
    );

  /** \brief Set a long integer option.
   *
   * \param  option_name    [in] (null terminated std::string) The name of the option
   *                        (without the leading '--' or trailing '=').
   * \param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
   *                        of the option (used for printing in --help).  On output,
   *                        will be set according to <tt>(argc,argv[])</tt>.
   * \param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated std::string
   *                        gives the documentation for the option.
   */
  void setOption(
    const char     option_name[]
    ,long int      *option_val
    ,const char    documentation[] = NULL
    ,const bool    required        = false
    );

  /** \brief Set a size_t option.
   *
   * \param  option_name    [in] (null terminated std::string) The name of the option
   *                        (without the leading '--' or trailing '=').
   * \param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
   *                        of the option (used for printing in --help).  On output,
   *                        will be set according to <tt>(argc,argv[])</tt>.
   * \param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated std::string
   *                        gives the documentation for the option.
   */
  void setOption(
    const char     option_name[]
    ,size_t        *option_val
    ,const char    documentation[] = NULL
    ,const bool    required        = false
    );

  /** \brief Set a long long integer option.
   *
   * \param  option_name    [in] (null terminated std::string) The name of the option
   *                        (without the leading '--' or trailing '=').
   * \param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
   *                        of the option (used for printing in --help).  On output,
   *                        will be set according to <tt>(argc,argv[])</tt>.
   * \param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated std::string
   *                        gives the documentation for the option.
   */
  void setOption(
    const char     option_name[]
    ,long long int *option_val
    ,const char    documentation[] = NULL
    ,const bool    required        = false
    );

  /** \brief Set a floating-point option.
   *
   * \param  option_name    [in] (null terminated std::string) The name of the option
   *                        (without the leading '--' or trailing '=').
   * \param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
   *                        of the option (used for printing in --help).  On output,
   *                        will be set according to <tt>(argc,argv[])</tt>.
   * \param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated std::string
   *                        gives the documentation for the option.
   */
  void setOption(
    const char     option_name[]
    ,double        *option_val
    ,const char    documentation[] = NULL
    ,const bool    required        = false
    );

  /** \brief Set a std::string option.
   *
   * \param  option_name    [in] (null terminated std::string) The name of the option
   *                        (without the leading '--' or trailing '=').
   * \param  option_val     [in/out] On input, <tt>*option_val</tt> gives the default value
   *                        of the option (used for printing in --help).  On output,
   *                        will be set according to <tt>(argc,argv[])</tt>.
   * \param  documentation  [in] If <tt>!=NULL</tt>, then this null terminated std::string
   *                        gives the documentation for the option.
   */
  void setOption(
    const char     option_name[]
    ,std::string   *option_val
    ,const char    documentation[] = NULL
    ,const bool    required        = false
    );

  /** \brief Set an enumeration option (templated by enumeration type).
   *
   * \param enum_option_name [in] (null terminated std::string) The name of
   * the option (without the leading '--' or trailing '=').
   *
   * \param enum_option_val [in/out] On input, <tt>*enum_option_val</tt> give
   * the default value of the enumeration (used for printing in --help).
   * After <tt>parse()</tt> finished executing successfully,
   * <tt>*enum_option_val</tt> will contain the user-selected value of the
   * enumeration.
   *
   * \param num_enum_opt_values [in] Gives the number of possible option
   * values to select
   *
   * \param enum_opt_values [in] Array (length <tt>num_enum_opt_values</tt>)
   * that gives the numeric values for each option.  The values in this array
   * are used to set the actual option <tt>*enum_option_val</tt>.
   *
   * \param enum_opt_names [in] Array (length <tt>num_enum_opt_values</tt>)
   * that gives the char string names for each option.  The strings in this
   * function are what is used in the commandline.
   *
   * \param documentation [in] If <tt>!=NULL</tt>, then this array of
   * null-terminated char string gives the documentation for the option.
   *
   * Warning! Only use an <tt>enum</tt> or <tt>int</tt> for <tt>EType</tt>.
   * Using any other type for <tt>EType</tt> could be trouble!
   */
  template <class EType>
  void setOption(
    const char    enum_option_name[]
    ,EType        *enum_option_val
    ,const int    num_enum_opt_values
    ,const EType  enum_opt_values[]
    ,const char*  enum_opt_names[]
    ,const char   documentation[] = NULL
    ,const bool   required        = false
    );

  //@}

  //! @name Parse
  //@{

  /** \brief Parse a command line.
   *
   * \param  argc    [in] number of entries in argv[]
   * \param  argv    [in/out] array (length argc) of command line arguments.
   *                 argv[0] should be the name of the program on the shell as usual.
   * \param  errout  [out] If <tt>!=NULL</tt> then error and help messages are sent here.
   *                 The default is set to <tt>&std::cerr</tt>.
   *
   * Postconditions:
   * <ul>
   * <li>If an unrecognized option is found
   *     <ul>
   *     <li>If <tt>this->recogniseAllOptions()==true</tt>
   *         <ul>
   *         <li>An error message will be printed to <tt>*errout</tt> and parsing will stop as follows:
   *         <li>If <tt>this->throwExceptions()==true</tt>
   *             <ul><li>This method will throw an <tt>UnrecognizedOption</tt> std::exception</ul>
   *         <li>else
   *             <ul><li>This method will return <tt>PARSE_UNRECOGNIZED_OPTION</tt></ul>
   *         <li>endif
   *         </ul>
   *     <li>else
   *         <ul><li>A warning message will be printed to <tt>*errout</tt> and parsing will continue</ul>
   *     <li>endif
   *     </ul>
   * <li>else if the option <tt>--help</tt> is found
   *     <ul>
   *     <li>If <tt>this->throwExceptions()==true</tt>
   *         <ul><li>This method will throw a <tt>HelpPrinted</tt> std::exception</ul>
   *     <li>else
   *         <ul><li>This method will return <tt>PARSE_HELP_PRINTED</tt></ul>
   *     <li>endif
   *     </ul>
   * <li>else
   *     <ul><li>This method will return <tt>PARSE_SUCCESSFUL</tt></ul>
   * <li>endif
   * </ul>
   *
   * Note that if the option <tt>--pause-for-debugging</tt> is
   * present, then std::string <tt>Type 0 and press enter to continue
   * :</tt> will be printed to standard error (<tt>std::cerr</tt>) and
   * execution will be suspended until the user enters any non-null
   * character(s).  This option is designed to make it easier to
   * attach a debugger, especially in a parallel MPI program.  If
   * HAVE_MPI is defined, then output/input is only performed with the
   * process with rank 0 and then MPI calls insure that all processors
   * wait (using <tt>MPI_Barrier(MPI_COMM_WORLD)</tt>) until the user
   * has entered something.  This allows the user to attach a debugger
   * to one or more parallel MPI processes and set breakpoints before
   * execution resumes.  Note that the stream <tt>*errout</tt> is not
   * used for this output/input but instead <tt>std::cerr</tt> is
   * directly used.
   *
   * If <tt>Teuchos::VerboseObjectBase::getDefaultOStream().get()!=NULL</tt>
   * and <tt>this->addOutputSetupOptions()</tt>, then any of the default setup
   * options for <tt>Teuchos::VerboseObjectBase::getDefaultOStream()</tt> that
   * are set on the commandline will be set on
   * <tt>Teuchos::VerboseObjectBase::getDefaultOStream()</tt>.
   */
  EParseCommandLineReturn  parse(
    int             argc
    ,char*          argv[]
    ,std::ostream   *errout = &std::cerr
    ) const;

  //@}

  //! @name Miscellaneous
  //@{

  /** \brief Print the help message.
   *
   * \param  out  [in/out] The stream the documentation will be printed to.
   *
   * This will print a formatted set of documentation that shows what
   * options are set, what their default values are and any
   * user-supplied documentation about each option.
   */
   void printHelpMessage( const char program_name[], std::ostream &out ) const;

  /** \brief Call to print timers so that they don't get printed in the
   * destructor.
   *
   * Calling this function after the first call has no effect.
   */
  void printFinalTimerSummary(const Ptr<std::ostream> &out = null);

  //@}

public:
  //
  enum EOptType { OPT_NONE, OPT_BOOL_TRUE, OPT_BOOL_FALSE, OPT_INT, OPT_LONG_INT, OPT_SIZE_T,
  OPT_LONG_LONG_INT,
  OPT_DOUBLE, OPT_STRING, OPT_ENUM_INT };

  // RAB: 2003/10/10: Note: I had to move this out of the private section since
  // the sun compiler (version 7) complained (rightly it now appears after looking
  // up what the ISO/ANSI C++ standard says) about the declaration for opt_val_val_t
  // not being able to access a private member of CommandLineProcessor.

private:

  // /////////////////////////////////
  // Private types

  // ToDo: RAB: 2004/05/25: Clean up these data structures and add
  // support for a templated enum type.  This will clean up usage
  // quite a bit.

  //
  struct opt_val_val_t {
    opt_val_val_t():
      opt_type(OPT_NONE),
      required(false),
      was_read(false)
      {}
    opt_val_val_t( EOptType opt_type_in, const any& opt_val_in, bool required_in )
      :opt_type(opt_type_in),opt_val(opt_val_in),required(required_in),was_read(false)
      {}
    EOptType     opt_type;
    any          opt_val; // Will be bool*, int*, double*, std::string* or a small int (for OPT_ENUM_INT)
    bool         required;
    bool         was_read;
  };

  //
  typedef Teuchos::map<std::string,opt_val_val_t>   options_list_t;

  //
  struct opt_doc_t {
    opt_doc_t()
      :opt_type(OPT_NONE)
      {}
    opt_doc_t(EOptType opt_type_in, const std::string& opt_name_in, const std::string& opt_name_false_in
            ,const std::string &documentation_in, const any &default_val_in )
      :opt_type(opt_type_in),opt_name(opt_name_in),opt_name_false(opt_name_false_in)
      ,documentation(documentation_in),default_val(default_val_in)
      {}
    EOptType     opt_type;
    std::string  opt_name;
    std::string  opt_name_false; // only for bool
    std::string  documentation;
    any          default_val;
  };

  //
  typedef std::vector<opt_doc_t>   options_documentation_list_t;

  //
  struct enum_opt_data_t {
    enum_opt_data_t()
      :enum_option_val(NULL), num_enum_opt_values(0)
      {}
    enum_opt_data_t(
      int          *_enum_option_val
      ,const int    _num_enum_opt_values
      ,const int    _enum_opt_values[]
      ,const char*  _enum_opt_names[]
      )
      :enum_option_val(_enum_option_val)
      ,num_enum_opt_values(_num_enum_opt_values)
      ,enum_opt_values(_enum_opt_values,_enum_opt_values+_num_enum_opt_values)
      {
        for( int k = 0; k < num_enum_opt_values; ++k )
          enum_opt_names.push_back(std::string(_enum_opt_names[k]));
      }
    int                  *enum_option_val;
    int                  num_enum_opt_values;
    std::vector<int>     enum_opt_values;
    std::vector<std::string>  enum_opt_names;
  };

  //
  typedef std::vector<enum_opt_data_t> enum_opt_data_list_t;

  // /////////////////////////////////
  // Private data members

  bool                             throwExceptions_;
  bool                             recogniseAllOptions_;
  bool                             addOutputSetupOptions_;
  std::string                      doc_string_;

  //use pragmas to disable some false positive warnings in windows sharedlib exports
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4251)
#endif
  mutable options_list_t           options_list_;
  options_documentation_list_t     options_documentation_list_;
  enum_opt_data_list_t             enum_opt_data_list_;
#ifdef _MSC_VER
#pragma warning(pop)
#endif

  bool  output_all_front_matter_;
  bool  output_show_line_prefix_;
  bool  output_show_tab_count_;
  bool  output_show_proc_rank_;
  int   output_to_root_rank_only_;
  bool  print_rcpnode_statistics_on_exit_;
  bool  show_timer_summary_on_exit_;

  bool printed_timer_summary_;

  bool  added_extra_output_setup_options_;
  bool  in_add_extra_output_setup_options_;

  static const bool  output_all_front_matter_default_;
  static const bool  output_show_line_prefix_default_;
  static const bool  output_show_tab_count_default_;
  static const bool  output_show_proc_rank_default_;
  static const int   output_to_root_rank_only_default_;
  static const bool  print_rcpnode_statistics_on_exit_default_;
  static const bool  show_timer_summary_on_exit_default_;

  // /////////////////////////////////
  // Private member functions

  // Set the extra output setup options
  void add_extra_output_setup_options() const;

  // Set an integer enumeration option
  void setEnumOption(
    const char    enum_option_name[]
    ,int          *enum_option_val
    ,const int    num_enum_opt_values
    ,const int    enum_opt_values[]
    ,const char*  enum_opt_names[]
    ,const char   documentation[]
    ,const bool   required
    );

  // Set an enum int option
  bool set_enum_value(
    int                  argv_i
    ,char*               argv[]
    ,const std::string   &enum_opt_name
    ,const int           enum_id
    ,const std::string   &enum_str_val
    ,std::ostream        *errout
    ) const;

  // Print the valid enum values
  void print_enum_opt_names(
    const int            enum_id
    ,std::ostream        &out
    ) const;

  // Return the name of the default value for an enum
  std::string enum_opt_default_val_name(
    const std::string    &enum_name
    ,const int           enum_id
    ,std::ostream        *errout
    ) const;

  // Return the index given and option value
  int find_enum_opt_index(
    const std::string           &enum_opt_name
    ,const int                  opt_value
    ,const enum_opt_data_t      &enum_data
    ,std::ostream               *errout
    ) const;

  // Get the option and the value from an entry in argv[].
  // Will return false if entry is not formated properly.
  bool get_opt_val(
    const char     str[]
    ,std::string   *opt_name
    ,std::string   *opt_val_str // May be empty on return
    ) const;

  // String for option type
  std::string opt_type_str( EOptType ) const;

  // Print bad option
  void print_bad_opt(
    int             argv_i
    ,char*          argv[]
    ,std::ostream   *errout
    ) const;

public: // Hidden implementation stuff that clients should never see

  /// \class TimeMonitorSurrogate
  /// \brief Interface by which CommandLineProcessor may use TimeMonitor.
  /// \warning Users should not use this class or rely on it in any way.
  ///   It is an implementation detail.
  ///
  /// \section Teuchos_TimeMonitorSurrogate_Summary Summary
  ///
  /// This class provides an interface by which CommandLineProcessor
  /// may optionally call TimeMonitor::summarize(), without needing to
  /// know that the TimeMonitor class exists.  This allows Teuchos to
  /// put CommandLineProcessor in a separate package from TimeMonitor.
  /// We want to do this because TimeMonitor depends on Comm, and is
  /// therefore in the TeuchosComm subpackage, but
  /// CommandLineProcessor does not depend on Comm, and is therefore
  /// in a different subpackage.  This design lets
  /// CommandLineProcessor automatically support showing summary
  /// timings just by having the TeuchosComm subpackage enabled and
  /// having its libaries linked in.
  ///
  /// \section Teuchos_TimeMonitorSurrogate_Note Note to Teuchos developers
  ///
  /// The TimeMonitorSurrogateImplInserter class in the TeuchosComm
  /// subpackage will ensure that CommandLineProcessor gets informed
  /// about TimeMonitor even before the program starts executing
  /// main().  This happens automatically; you don't need to change
  /// the main() function.
  ///
  /// This is an instance of the <a
  /// href="http://en.wikipedia.org/wiki/Dependency_injection">Dependency
  /// injection</a> design pattern.  CommandLineProcessor is not
  /// supposed to know about TimeMonitor, because
  /// CommandLineProcessor's subpackage does not depend on
  /// TimeMonitor's subpackage.  Thus, CommandLineProcessor interacts
  /// with TimeMonitor through the TimeMonitorSurrogate interface.
  class TimeMonitorSurrogate {
  public:
    ///! brief.
    virtual ~TimeMonitorSurrogate() {}
    //! Summarize timings over all process(es) to the given output stream.
    virtual void summarize(std::ostream &out=std::cout) = 0;
  };

  static void setTimeMonitorSurrogate(const RCP<TimeMonitorSurrogate> &timeMonitorSurrogate);

  static RCP<TimeMonitorSurrogate> getTimeMonitorSurrogate();

private:

  static RCP<TimeMonitorSurrogate>& getRawTimeMonitorSurrogate();

}; // end class CommandLineProcessor


// /////////////////////////
// Inline members


// Behavior modes


inline
void CommandLineProcessor::throwExceptions( const bool & throwExceptions_in )
{ throwExceptions_ = throwExceptions_in; }


inline
bool CommandLineProcessor::throwExceptions() const
{ return throwExceptions_; }


inline
void CommandLineProcessor::recogniseAllOptions( const bool & recogniseAllOptions_in )
{ recogniseAllOptions_ = recogniseAllOptions_in; }


inline
bool CommandLineProcessor::recogniseAllOptions() const
{ return recogniseAllOptions_; }


inline
void CommandLineProcessor::addOutputSetupOptions( const bool &addOutputSetupOptions_in )
{ addOutputSetupOptions_ = addOutputSetupOptions_in; }


inline
bool CommandLineProcessor::addOutputSetupOptions() const
{ return addOutputSetupOptions_; }


template <class EType>
inline
void CommandLineProcessor::setOption(
  const char    enum_option_name[]
  ,EType       *enum_option_val
  ,const int    num_enum_opt_values
  ,const EType  enum_opt_values[]
  ,const char*  enum_opt_names[]
  ,const char   documentation[]
  ,const bool   required
  )
{
  // RAB: 2004/05/25: Every C++ implementation that I know of just
  // represents enumerations as int's and therefore this will compile
  // just fine.  However, the ISO/ANSI C++ standard says that
  // compilers are allowed to use a smaller storage type for an enum
  // but must not require storage any larger than an 'int'.  If the
  // below compile-time assertion does not compile then we need to do
  // something different but it will be a lot of work!
  CompileTimeAssert<sizeof(int)-sizeof(EType)>();
  //CompileTimeAssert<sizeof(int)-sizeof(EType)-1>(); // Uncomment to see compilation error
  setEnumOption(
    enum_option_name
    ,reinterpret_cast<int*>(enum_option_val)
    ,num_enum_opt_values
    ,reinterpret_cast<const int*>(enum_opt_values)
    ,enum_opt_names
    ,documentation
    ,required
    );
}


inline
std::string CommandLineProcessor::opt_type_str( EOptType opt_type ) const
{
  std::string str;
  switch( opt_type ) {
    case OPT_BOOL_TRUE:
      str = "bool";
      break;
    case OPT_INT:
      str = "int";
      break;
    case OPT_LONG_INT:
      str = "long int";
      break;
    case OPT_SIZE_T:
      str = "size_t";
      break;
    case OPT_LONG_LONG_INT:
      str = "long long int";
      break;
    case OPT_DOUBLE:
      str = "double";
      break;
    case OPT_STRING:
      str = "string";
      break;
    case OPT_ENUM_INT:
      str = "enum";
      break;
    default:
      assert(0); // Local programming error only
  }
  return str;
}


} // end namespace Teuchos


#endif // TEUCHOS_COMMAND_LINE_PROCESSOR_HPP
