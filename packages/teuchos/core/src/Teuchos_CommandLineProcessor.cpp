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

// //////////////////////////////////////////////////
// Teuchos_CommandLineProcessor.cpp


#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
//#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_Assert.hpp"

#ifdef HAVE_MPI
#  include "mpi.h"
#endif


namespace {


inline int my_max( int a, int b ) { return a > b ? a : b; }


std::string remove_quotes( const std::string& str )
{
  if(str[0] != '\"')
    return str;
  return str.substr(1,str.size()-2);
}


std::string add_quotes( const std::string& str )
{
  if(str[0] == '\"')
    return str;
  return "\"" + str + "\"";
}


} // end namespace


namespace Teuchos {


const bool  CommandLineProcessor::output_all_front_matter_default_(false);
const bool  CommandLineProcessor::output_show_line_prefix_default_(false);
const bool  CommandLineProcessor::output_show_tab_count_default_(false);
const bool  CommandLineProcessor::output_show_proc_rank_default_(false);
const int   CommandLineProcessor::output_to_root_rank_only_default_(0);
const bool  CommandLineProcessor::print_rcpnode_statistics_on_exit_default_(false);
const bool  CommandLineProcessor::show_timer_summary_on_exit_default_(false);


CommandLineProcessor::CommandLineProcessor(
  bool   throwExceptions_in
  ,bool  recogniseAllOptions_in
  ,bool  addOutputSetupOptions_in
  )
  :throwExceptions_(throwExceptions_in)
  ,recogniseAllOptions_(recogniseAllOptions_in)
  ,addOutputSetupOptions_(addOutputSetupOptions_in)
  ,output_all_front_matter_(output_all_front_matter_default_)
  ,output_show_line_prefix_(output_show_line_prefix_default_)
  ,output_show_tab_count_(output_show_tab_count_default_)
  ,output_show_proc_rank_(output_show_proc_rank_default_)
  ,output_to_root_rank_only_(output_to_root_rank_only_default_)
  ,print_rcpnode_statistics_on_exit_(print_rcpnode_statistics_on_exit_default_)
  ,show_timer_summary_on_exit_(show_timer_summary_on_exit_default_)
  ,printed_timer_summary_(false)
  ,added_extra_output_setup_options_(false)
  ,in_add_extra_output_setup_options_(false)
{}


CommandLineProcessor::~CommandLineProcessor()
{
  printFinalTimerSummary();
}


// Set up options


void CommandLineProcessor::setDocString( const char doc_string[] )
{
  doc_string_ = doc_string;
}


void CommandLineProcessor::setOption(
  const char     option_true[]
  ,const char    option_false[]
  ,bool          *option_val
  ,const char    documentation[]
  )
{
  add_extra_output_setup_options();
  TEUCHOS_TEST_FOR_EXCEPT(!(option_val!=NULL));
  options_list_[std::string(option_true)]
    = opt_val_val_t(OPT_BOOL_TRUE,any(option_val),false);
  options_list_[std::string(option_false)]
    = opt_val_val_t(OPT_BOOL_FALSE,any(option_val),false);
  options_documentation_list_.push_back(
    opt_doc_t(OPT_BOOL_TRUE, option_true, option_false,
      std::string(documentation?documentation:""), any(option_val)) 
    );
}


void CommandLineProcessor::setOption(
  const char     option_name[]
  ,int           *option_val
  ,const char    documentation[]
  ,const bool    required
  )
{
  add_extra_output_setup_options();
  TEUCHOS_TEST_FOR_EXCEPT(!(option_val!=NULL));
  options_list_[std::string(option_name)]
    = opt_val_val_t(OPT_INT,any(option_val),required);
  options_documentation_list_.push_back(
    opt_doc_t(OPT_INT, option_name, "", std::string(documentation?documentation:""),
      any(option_val))
    );
}


void CommandLineProcessor::setOption(
  const char     option_name[]
  ,double        *option_val
  ,const char    documentation[]
  ,const bool    required
  )
{
  add_extra_output_setup_options();
  TEUCHOS_TEST_FOR_EXCEPT(!(option_val!=NULL));
  options_list_[std::string(option_name)]
    = opt_val_val_t(OPT_DOUBLE,any(option_val),required);
  options_documentation_list_.push_back(
    opt_doc_t(OPT_DOUBLE, option_name, "", std::string(documentation?documentation:""),
      any(option_val))
    );
}


void CommandLineProcessor::setOption(
  const char     option_name[]
  ,std::string   *option_val
  ,const char    documentation[]
  ,const bool    required
  )
{
  add_extra_output_setup_options();
  TEUCHOS_TEST_FOR_EXCEPT(!(option_val!=NULL));
  options_list_[std::string(option_name)]
    = opt_val_val_t(OPT_STRING,any(option_val),required);
  options_documentation_list_.push_back(
    opt_doc_t(OPT_STRING, option_name, "", std::string(documentation?documentation:""),
      any(option_val))
    );
}


// Parse command line


CommandLineProcessor::EParseCommandLineReturn
CommandLineProcessor::parse(
  int             argc
  ,char*          argv[]
  ,std::ostream   *errout
  ) const
{
  add_extra_output_setup_options();
  std::string        opt_name;
  std::string        opt_val_str;
  const std::string  echo_cl_opt = "echo-command-line";
  const std::string  help_opt = "help";
  const std::string  pause_opt = "pause-for-debugging";
  int procRank = GlobalMPISession::getRank();
  for( int i = 1; i < argc; ++i ) {
    bool gov_return = get_opt_val( argv[i], &opt_name, &opt_val_str );
    if( !gov_return ) {
      if(procRank == 0) 
        print_bad_opt(i,argv,errout);
      if( recogniseAllOptions() ) 
        return PARSE_UNRECOGNIZED_OPTION;
      else {
        continue;
      }
    }
    if( opt_name == echo_cl_opt ) {
      if(errout && procRank == 0) {
        *errout << "\nEchoing the command-line:\n\n";
        for( int j = 0; j < argc; ++j )
          *errout << argv[j] << " ";
        *errout << "\n\n";
      }
      continue;
    }
    if( opt_name == help_opt ) {
      if(errout) printHelpMessage( argv[0], *errout );
      return PARSE_HELP_PRINTED;
    }
    if( opt_name == pause_opt ) {
      if(procRank == 0) {
        std::cerr << "\nType 0 and press enter to continue : ";
        int dummy_int = 0;
        std::cin >> dummy_int;
      }
#ifdef HAVE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      continue;
    }
    // Lookup the option (we had better find it!)
    options_list_t::iterator  itr = options_list_.find(opt_name);
    if( itr == options_list_.end() ) {
      if(procRank == 0)
        print_bad_opt(i,argv,errout);
      if( recogniseAllOptions() ) 
        return PARSE_UNRECOGNIZED_OPTION;
      else
        continue;
    }
    // Changed access to second value of std::map to not use overloaded arrow operator, 
    // otherwise this code will not compile on Janus (HKT, 12/01/2003) 
    opt_val_val_t &opt_val_val = (*itr).second;
    opt_val_val.was_read = true;
    switch( opt_val_val.opt_type ) {
      case OPT_BOOL_TRUE:
        *(any_cast<bool*>(opt_val_val.opt_val)) = true;
        break;
      case OPT_BOOL_FALSE:
        *(any_cast<bool*>(opt_val_val.opt_val)) = false;
        break;
      case OPT_INT:
        *(any_cast<int*>(opt_val_val.opt_val)) = std::atoi(opt_val_str.c_str());
        break;
      case OPT_DOUBLE:
        *(any_cast<double*>(opt_val_val.opt_val)) = std::atof(opt_val_str.c_str());
        break;
      case OPT_STRING:
        *(any_cast<std::string*>(opt_val_val.opt_val)) = remove_quotes(opt_val_str);
        break;
      case OPT_ENUM_INT:
        if( !set_enum_value( i, argv, opt_name, any_cast<int>(opt_val_val.opt_val),
            remove_quotes(opt_val_str), errout ) )
        {
          return PARSE_UNRECOGNIZED_OPTION;
        }
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT(true); // Local programming error only
    }
  }
  // Look for options that were required but were not set
  for(
    options_list_t::const_iterator itr = options_list_.begin();
    itr != options_list_.end();
    ++itr
    )
  {
    const opt_val_val_t   &opt_val_val  = (*itr).second;
    if( opt_val_val.required && !opt_val_val.was_read ) {
      const std::string     &opt_val_name = (*itr).first;
#define CLP_ERR_MSG \
      "Error, the option --"<<opt_val_name<<" was required but was not set!"
      if(errout)
        *errout << std::endl << argv[0] << " : " << CLP_ERR_MSG << std::endl;
      if( throwExceptions() ) {
        TEUCHOS_TEST_FOR_EXCEPTION( true, ParseError, CLP_ERR_MSG );
      }
      return PARSE_ERROR;
#undef CLP_ERR_MSG
    }
  }
  // Set the options of a default stream exists and if we are asked to
  RCP<FancyOStream>
    defaultOut = VerboseObjectBase::getDefaultOStream();
  if (defaultOut.get() && addOutputSetupOptions_) {
    if (output_all_front_matter_ != output_all_front_matter_default_)
      defaultOut->setShowAllFrontMatter(output_all_front_matter_);
    if (output_show_line_prefix_ != output_show_line_prefix_default_)
      defaultOut->setShowLinePrefix(output_show_line_prefix_);
    if (output_show_tab_count_ != output_show_tab_count_default_)
      defaultOut->setShowTabCount(output_show_tab_count_);
    if (output_show_proc_rank_ != output_show_proc_rank_default_)
      defaultOut->setShowProcRank(output_show_proc_rank_);
    if (output_to_root_rank_only_ != output_to_root_rank_only_default_)
      defaultOut->setOutputToRootOnly(output_to_root_rank_only_);
    RCPNodeTracer::setPrintRCPNodeStatisticsOnExit(print_rcpnode_statistics_on_exit_);
  }
  return PARSE_SUCCESSFUL;
}


void CommandLineProcessor::printHelpMessage( const char program_name[],
  std::ostream &out ) const
{
  add_extra_output_setup_options();
  int procRank = GlobalMPISession::getRank();
  if (procRank == 0) {
    using std::setw;
    using std::endl;
    
    const int opt_type_w = 8;
    const char spc_chars[] = "  ";
    
    // Get the maximum length of an option name
    int opt_name_w = 19; // For the 'pause-for-debugging' option
    options_documentation_list_t::const_iterator itr;
    for (
      itr = options_documentation_list_.begin();
      itr != options_documentation_list_.end();
      ++itr
      )
    {
      opt_name_w = my_max(opt_name_w,itr->opt_name.length());
      if( itr->opt_type )
        opt_name_w = my_max(opt_name_w,itr->opt_name_false.length());
    }
    opt_name_w += 2;
    
    // Some built-in options
    out
      << "Usage: " << program_name << " [options]\n"
      << spc_chars << "options:\n"
      << spc_chars
      << "--"
#ifdef HAVE_STD_IOS_BASE_FMTFLAGS
      << std::left << setw(opt_name_w) << "help"
      << std::left << setw(opt_type_w) << " "
#else
      << std::setiosflags(std::ios::left) << setw(opt_name_w) << "help"
      << std::setiosflags(std::ios::left) << setw(opt_type_w) << " "
#endif
      << "Prints this help message"
      << std::endl
      << spc_chars
      << "--"
#ifdef HAVE_STD_IOS_BASE_FMTFLAGS
      << std::left << setw(opt_name_w) << "pause-for-debugging"
      << std::left << setw(opt_type_w) << " "
#else
      << std::setiosflags(std::ios::left) << setw(opt_name_w) << "pause-for-debugging"
      << std::setiosflags(std::ios::left) << setw(opt_type_w) << " "
#endif
      << "Pauses for user input to allow attaching a debugger"
      << std::endl
      << spc_chars
      << "--"
#ifdef HAVE_STD_IOS_BASE_FMTFLAGS
      << std::left << setw(opt_name_w) << "echo-command-line"
      << std::left << setw(opt_type_w) << " "
#else
      << std::setiosflags(std::ios::left) << setw(opt_name_w) << "echo-command-line"
      << std::setiosflags(std::ios::left) << setw(opt_type_w) << " "
#endif
      << "Echo the command-line but continue as normal"
      << std::endl;
    for(
      itr = options_documentation_list_.begin();
      itr != options_documentation_list_.end();
      ++itr )
    {
      // print top line with option name, type and short documentation string
      out
        << spc_chars
        << "--"
#ifdef HAVE_STD_IOS_BASE_FMTFLAGS
        << std::left << setw(opt_name_w) << itr->opt_name
        << std::left << setw(opt_type_w) << opt_type_str(itr->opt_type)
#else
        << std::setiosflags(std::ios::left) << setw(opt_name_w) << itr->opt_name
        << std::setiosflags(std::ios::left) << setw(opt_type_w) << opt_type_str(itr->opt_type)
#endif
        << ( itr->documentation.length() ? itr->documentation.c_str() : "No documentation" )
        << std::endl;
      // If an enumeration option then the next line is the value options
      if( itr->opt_type == OPT_ENUM_INT ) {
        out
          << spc_chars
          << "  "
          << setw(opt_name_w) << ""
          << setw(opt_type_w) << "";
        print_enum_opt_names( any_cast<int>(itr->default_val), out );
        out
          << std::endl;
      }
      // Now print the line that contains the default values
      if( itr->opt_type == OPT_BOOL_TRUE ) {
        out
          << spc_chars
          << "--"
          << setw(opt_name_w) << itr->opt_name_false;
      }
      else {
        out
          << spc_chars
          << "  "
          << setw(opt_name_w) << " ";
      }
      out
        << setw(opt_type_w) << " "
        << "(default: ";
      switch( itr->opt_type ) {
        case OPT_BOOL_TRUE:
          out << "--" << ( (*(any_cast<bool*>(itr->default_val))) ?
            itr->opt_name : itr->opt_name_false );
          break;
        case OPT_INT:
        case OPT_DOUBLE:
        case OPT_STRING:
        case OPT_ENUM_INT:
          out << "--" << itr->opt_name;
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true); // Local programming error only
      }
      switch( itr->opt_type ) {
        case OPT_BOOL_TRUE:
          break;
        case OPT_INT:
          out << "=" << (*(any_cast<int*>(itr->default_val)));
          break;
        case OPT_DOUBLE:
          out <<  "=" << (*(any_cast<double*>(itr->default_val)));
          break;
        case OPT_STRING:
          out <<  "=" << add_quotes(*(any_cast<std::string*>(itr->default_val)));
          break;
        case OPT_ENUM_INT:
          out <<  "=" << add_quotes(
            enum_opt_default_val_name(itr->opt_name,any_cast<int>(itr->default_val),&out));
          break;
        default:
          TEUCHOS_TEST_FOR_EXCEPT(true); // Local programming error only
      }
      out << ")\n";
    }
    if(doc_string_.length()) {
      out << "\nDETAILED DOCUMENTATION:\n\n" << doc_string_ << std::endl << std::endl;
    }
    if(throwExceptions_)
      TEUCHOS_TEST_FOR_EXCEPTION( true, HelpPrinted, "Help message was printed" );
  }
}


void CommandLineProcessor::printFinalTimerSummary(
  const Ptr<std::ostream> &out_inout
  )
{
  if (!printed_timer_summary_ && show_timer_summary_on_exit_) {
    RCP<std::ostream> out;
    if (nonnull(out_inout)) {
      out = rcpFromPtr(out_inout);
    }
    else {
      out = VerboseObjectBase::getDefaultOStream();
    }
    getTimeMonitorSurrogate()->summarize(*out << "\n");
    printed_timer_summary_ = true;
  }
}


// private


void CommandLineProcessor::add_extra_output_setup_options() const
{
  if(
    // Are we in this function already and calling it recursively?
    in_add_extra_output_setup_options_
    ||
    // Have we already setup these options?
    added_extra_output_setup_options_
    ||
    // Are we not supposed to setup these options?
    !addOutputSetupOptions_
    )
  {
    return; // If any of the above is true, we need to return right away!
  }
  // Set the commandline options for this ...
  CommandLineProcessor
    *clp = const_cast<CommandLineProcessor*>(this);
  clp->in_add_extra_output_setup_options_ = true;
  clp->setOption(
    "output-all-front-matter","output-no-front-matter",&clp->output_all_front_matter_
    ,"Set if all front matter is printed to the default FancyOStream or not"
    );
  clp->setOption(
    "output-show-line-prefix","output-no-show-line-prefix",&clp->output_show_line_prefix_
    ,"Set if the line prefix matter is printed to the default FancyOStream or not"
    );
  clp->setOption(
    "output-show-tab-count","output-no-show-tab-count",&clp->output_show_tab_count_
    ,"Set if the tab count is printed to the default FancyOStream or not"
    );
  clp->setOption(
    "output-show-proc-rank","output-no-show-proc-rank",&clp->output_show_proc_rank_
    ,"Set if the processor rank is printed to the default FancyOStream or not"
    );
  clp->setOption(
    "output-to-root-rank-only",&clp->output_to_root_rank_only_
    ,"Set which processor (the root) gets the output.  If < 0, then all processors get output."
    );
  clp->setOption(
    "print-rcpnode-statistics-on-exit", "no-print-rcpnode-statistics-on-exit",
    &clp->print_rcpnode_statistics_on_exit_,
    "Set if the RCPNode usage statistics will be printed on exit or not.  Warning,"
    " this prints to std::cerr or every process so do not turn this on for very large"
    " parallel runs."
    );
  if (nonnull(getTimeMonitorSurrogate())) {
    clp->setOption(
      "show-timer-summary", "no-show-timer-sumary", &clp->show_timer_summary_on_exit_,
      "If true, then Teuchos::TimeMonitor::summarize() is called in"
      " CommandLineProcessor's destructor (usually at the end of main)."
      );
  }

  clp->added_extra_output_setup_options_ = true;
  clp->in_add_extra_output_setup_options_ = false;
}


void CommandLineProcessor::setEnumOption(
  const char    enum_option_name[]
  ,int          *enum_option_val
  ,const int    num_enum_opt_values
  ,const int    enum_opt_values[]
  ,const char*  enum_opt_names[]
  ,const char   documentation[]
  ,const bool   required
  )
{
  add_extra_output_setup_options();

  TEUCHOS_TEST_FOR_EXCEPT(enum_option_val==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(num_enum_opt_values<=0);
  TEUCHOS_TEST_FOR_EXCEPT(enum_opt_values==NULL);
  TEUCHOS_TEST_FOR_EXCEPT(enum_opt_names==NULL);

  enum_opt_data_list_.push_back(
    enum_opt_data_t(enum_option_val,num_enum_opt_values,enum_opt_values,enum_opt_names)
    );
  const int opt_id = enum_opt_data_list_.size()-1;
  options_list_[std::string(enum_option_name)]
    = opt_val_val_t(OPT_ENUM_INT,any(opt_id),required);
  options_documentation_list_.push_back(
    opt_doc_t(OPT_ENUM_INT,enum_option_name, "",
      std::string(documentation?documentation:""), any(opt_id))
    );
}


bool CommandLineProcessor::set_enum_value(
  int                  argv_i
  ,char*               argv[]
  ,const std::string   &enum_opt_name
  ,const int           enum_id
  ,const std::string   &enum_str_val
  ,std::ostream        *errout
  ) const
{
  const enum_opt_data_t
    &enum_opt_data = enum_opt_data_list_.at(enum_id);
  std::vector<std::string>::const_iterator
    itr_begin = enum_opt_data.enum_opt_names.begin(),
    itr_end   = enum_opt_data.enum_opt_names.end(),
    itr       =  std::find( itr_begin, itr_end, enum_str_val );
  if( itr == itr_end ) {
    const int j = argv_i;
#define CLP_ERR_MSG \
      "Error, the value \"" << enum_str_val << "\" for the " \
      << j<<(j==1?"st":(j==2?"nd":(j==3?"rd":"th"))) << " option --" \
      << enum_opt_name << " was not recognized (use --help)!"
    if(errout)
      *errout << std::endl << argv[0] << " : " << CLP_ERR_MSG << std::endl;
    if( throwExceptions() ) {
      TEUCHOS_TEST_FOR_EXCEPTION( true, UnrecognizedOption, CLP_ERR_MSG );
    }
    else {
      return false;
    }
#undef CLP_ERR_MSG
  }
  const int enum_opt_val_index = itr - itr_begin;
  *enum_opt_data.enum_option_val = enum_opt_data.enum_opt_values.at(enum_opt_val_index);
  return true;
}


void CommandLineProcessor::print_enum_opt_names(
  const int            enum_id
  ,std::ostream        &out
  ) const
{
  const enum_opt_data_t
    &enum_opt_data = enum_opt_data_list_.at(enum_id);
  typedef std::vector<std::string>::const_iterator itr_t;
  out << "Valid options:";
  for(
    itr_t itr = enum_opt_data.enum_opt_names.begin();
    itr != enum_opt_data.enum_opt_names.end();
    ++itr
    )
  {
    if( itr != enum_opt_data.enum_opt_names.begin() ) out << ",";
    out << " " << add_quotes(*itr);
  }
}


std::string
CommandLineProcessor::enum_opt_default_val_name(
  const std::string    &enum_name
  ,const int           enum_id
  ,std::ostream        *errout
  ) const
{
  const enum_opt_data_t
    &enum_opt_data = enum_opt_data_list_.at(enum_id);
  return enum_opt_data.enum_opt_names.at(
    find_enum_opt_index(
      enum_name,*enum_opt_data.enum_option_val,enum_opt_data,errout
      )
    );
}


int CommandLineProcessor::find_enum_opt_index(
  const std::string           &enum_opt_name
  ,const int                  opt_value
  ,const enum_opt_data_t      &enum_data
  ,std::ostream               *errout
  ) const
{
  std::vector<int>::const_iterator
    itr_begin = enum_data.enum_opt_values.begin(),
    itr_end   = enum_data.enum_opt_values.end(),
    itr       =  std::find( itr_begin, itr_end, opt_value );
  if( itr == itr_end ) {
#define CLP_ERR_MSG \
      ( recogniseAllOptions() ? "Error" : "Warning" ) \
      << ", option --" << enum_opt_name << " was given an invalid " \
      "initial option value of " << opt_value << "!"
    if(errout)
      *errout << CLP_ERR_MSG << std::endl;
    if( throwExceptions() )
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument, CLP_ERR_MSG );
#undef CLP_ERR_MSG
  }
  return itr - itr_begin;
}


bool CommandLineProcessor::get_opt_val(
  const char     str[]
  ,std::string   *opt_name
  ,std::string   *opt_val_str
  ) const
{
  const int len = std::strlen(str);
  if( len < 3 )
    return false; // Can't be an option with '--' followed by at least one char
  if( str[0] != '-' || str[1] != '-' )
    return false; // Not a recognised option
  // Find the '='
  int equ_i;
  for( equ_i = 2; equ_i < len && str[equ_i] != '='; ++equ_i );
  // Set opt_name
  opt_name->assign( str + 2, equ_i-2 );
  // Set opt_val_str
  if( equ_i == len ) {
    *opt_val_str = "";
  }
  else {
    opt_val_str->assign( str + equ_i + 1, len - equ_i - 1 );
  }
  return true;
}

void CommandLineProcessor::print_bad_opt(
  int             argv_i
  ,char*          argv[]
  ,std::ostream   *errout
  ) const
{
  const int j = argv_i;
#define CLP_ERR_MSG \
    ( recogniseAllOptions() ? "Error" : "Warning" ) \
    << ", the " << j<<(j==1?"st":(j==2?"nd":(j==3?"rd":"th"))) \
    << " option \'" << argv[argv_i] << "\' was not recognized (use --help)!"
  if(errout)
    *errout << std::endl << argv[0] << " : " << CLP_ERR_MSG << std::endl;
  if( recogniseAllOptions() && throwExceptions() )
    TEUCHOS_TEST_FOR_EXCEPTION( true, UnrecognizedOption, CLP_ERR_MSG );
#undef CLP_ERR_MSG
}


// Hidden stuff


void CommandLineProcessor::setTimeMonitorSurrogate(
  const RCP<CommandLineProcessor::TimeMonitorSurrogate> &timeMonitorSurrogate)
{
  getRawTimeMonitorSurrogate() = timeMonitorSurrogate;
}


RCP<CommandLineProcessor::TimeMonitorSurrogate>
CommandLineProcessor::getTimeMonitorSurrogate()
{
  return getRawTimeMonitorSurrogate();
}


RCP<CommandLineProcessor::TimeMonitorSurrogate>&
CommandLineProcessor::getRawTimeMonitorSurrogate()
{
  static RCP<TimeMonitorSurrogate> timeMonitorSurrogate;
  return timeMonitorSurrogate;
} 


} // end namespace Teuchos
