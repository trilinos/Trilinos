// Copyright (c) 2014, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <fstream>
#include <cstring>
#include <vector>
#include <sstream>
#include <climits>
#include <iostream>
#include <iomanip>
#include <ctime>

#include "aprepro.h"
#include "aprepro_parser.h"
#include "apr_scanner.h"
#include "apr_stats.h"
#include "apr_util.h"


namespace {
  const unsigned int HASHSIZE = 5939;
  const char* version_string = "4.20 (2014/09/24)";
  
  unsigned hash_symbol (const char *symbol)
  {
    unsigned hashval;
    for (hashval = 0; *symbol != '\0'; symbol++)
      hashval = *symbol + 65599 * hashval;
    return (hashval % HASHSIZE);
  }
}

namespace SEAMS {
  Aprepro *aprepro;  // A global for use in the library.  Clean this up...
  int   echo = true;
  
  Aprepro::Aprepro()
    : lexer(NULL), sym_table(HASHSIZE),
      stringInteractive(false), stringScanner(NULL),
      errorStream(&std::cerr), warningStream(&std::cerr), infoStream(&std::cout),
      stateImmutable(false), doLoopSubstitution(true), doIncludeSubstitution(true),
      isCollectingLoop(false)
  {
    ap_file_list.push(file_rec());
    init_table("#");
    aprepro = this;

    // See the random number generator...
    time_t time_val = std::time ((time_t*)NULL);
    srand((unsigned)time_val);
  }

  Aprepro::~Aprepro()
  {
    if(stringScanner && stringScanner != lexer)
      delete stringScanner;

    delete lexer;
    
    for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
      for (symrec *ptr = sym_table[hashval]; ptr != NULL; ) {
	symrec *save = ptr;
	ptr = ptr->next;
	delete save;
      }
    }
    aprepro = NULL;
    cleanup_memory();
  }

  std::string Aprepro::version() const {return version_string;}

  void Aprepro::clear_results()
  {
    parsingResults.str("");
    parsingResults.clear();
  }
  
  bool Aprepro::parse_stream(std::istream& in, const std::string& in_name)
  {
    ap_file_list.top().name = in_name;

    if (!ap_options.include_file.empty()) {
      file_rec include_file(ap_options.include_file.c_str(), 0, false, 0);
      ap_file_list.push(include_file);
      // File included on command line will be processed as immutable and no-echo
      // Will revert to global settings at end of file.
      stateImmutable = true;
      echo = false;
    }

    Scanner *scanner = new Scanner(*this, &in, &parsingResults);
    this->lexer = scanner;

    Parser parser(*this);
    parser.set_debug_level(ap_options.trace_parsing);
    return (parser.parse() == 0);
  }

  bool Aprepro::parse_file(const std::string &filename)
  {
    std::ifstream in(filename.c_str());
    if (!in.good()) return false;
    return parse_stream(in, filename);
  }

  bool Aprepro::parse_string(const std::string &input, const std::string& sname)
  {
    std::istringstream iss(input);
    return parse_stream(iss, sname);
  }

  bool Aprepro::parse_strings(const std::vector<std::string> &input, const std::string& sname)
  {
    std::stringstream iss;
    for (size_t i=0; i < input.size(); i++) {
      iss << input[i] << '\n';
    }
    return parse_stream(iss, sname);
  }

  bool Aprepro::parse_string_interactive(const std::string &input)
  {
    stringInteractive = true;

    if(ap_file_list.size() == 1)
      ap_file_list.top().name = "interactive_input";

    if (!ap_options.include_file.empty()) {
      file_rec include_file(ap_options.include_file.c_str(), 0, false, 0);
      ap_file_list.push(include_file);
      // File included on command line will be processed as immutable and no-echo
      // Will revert to global settings at end of file.
      stateImmutable = true;
      echo = false;
    }

    if(!stringScanner)
      stringScanner = new Scanner(*this, &stringInput, &parsingResults);

    this->lexer = stringScanner;

    stringInput.str(input);
    stringInput.clear();

    Parser parser(*this);
    parser.set_debug_level(ap_options.trace_parsing);
    bool result = parser.parse() == 0;

    stringInteractive = false;
    return result;
  }

  void Aprepro::error(const std::string& msg, bool line_info, bool prefix) const
  {
    std::stringstream ss;
    if (prefix) {
      (*errorStream) << "Aprepro: ERROR: ";
    }

    ss << msg;

    if(line_info) {
      ss << " (" << ap_file_list.top().name <<
	", line " << ap_file_list.top().lineno + 1 << ")";
    }
    ss << "\n";

    // Send it to the user defined stream
    (*errorStream) << ss.str();
  }

  void Aprepro::warning(const std::string &msg, bool line_info, bool prefix) const
  {
    if(!ap_options.warning_msg)
      return;

    std::stringstream ss;
    if (prefix) {
      (*warningStream) << "Aprepro: WARNING: ";
    }

    ss << msg;

    if(line_info) {
      ss << " (" << ap_file_list.top().name <<
            ", line " << ap_file_list.top().lineno + 1 << ")";
    }
    ss << "\n";

    // Send it to the user defined stream
    (*warningStream) << ss.str();
  }

  void Aprepro::info(const std::string &msg, bool line_info, bool prefix) const
  {
    if(!ap_options.info_msg)
      return;

    std::stringstream ss;
    if(prefix) {
      (*infoStream) << "Aprepro: INFO: ";
    }
    ss << msg;

    if(line_info) {
      ss << " (" << ap_file_list.top().name <<
            ", line " << ap_file_list.top().lineno + 1 << ")";
    }
    ss << "\n";

    // Send it to the user defined stream
    (*infoStream) << ss.str();
  }

  void Aprepro::set_error_streams(std::ostream *error,
                                  std::ostream *warning, std::ostream *info)
  {
    errorStream = error;
    warningStream = warning;
    infoStream = info;
  }

  /* Two methods for opening files. In OPEN_FILE, the file must exist
     or else the code will exit. If CHECK_OPEN_FILE, it is OK if the
     file does not exist. A Null 'pointer' will be returned.
  */
  std::fstream* Aprepro::open_file(const std::string &file, const char *mode)
  {
    std::ios::openmode smode = std::ios::in;
    if (mode[0] == 'w')
      smode = std::ios::out;
  
    /* See if file exists in current directory (or as specified) */
    std::fstream *pointer = new std::fstream(file.c_str(), smode);
    if ((pointer == NULL || pointer->bad() || !pointer->good()) && !ap_options.include_path.empty()) {
      /* If there is an include path specified, try opening file there */
      std::string file_path(ap_options.include_path);
      file_path += "/";
      file_path += file;
      pointer = new std::fstream(file_path.c_str(), smode);
    }

    /* If pointer still null, print error message */
    if (pointer == NULL || pointer->bad() || !pointer->good()) {
      error("Can't open " + file, false);
      exit(EXIT_FAILURE);
    }

    return pointer;
  }

  std::fstream *Aprepro::check_open_file(const std::string &file, const char *mode)
  {
    std::ios::openmode smode = std::ios::in;
    if (mode[0] == 'w')
      smode = std::ios::out;

    std::fstream *pointer = new std::fstream(file.c_str(), smode);

    if ((pointer == NULL || pointer->bad() || !pointer->good()) && !ap_options.include_path.empty()) {
      /* If there is an include path specified, try opening file there */
      std::string file_path(ap_options.include_path);
      file_path += "/";
      file_path += file;
      pointer = new std::fstream(file_path.c_str(), smode);
    }
    return pointer;
  }

  symrec *Aprepro::putsym (const std::string &sym_name, SYMBOL_TYPE sym_type, bool is_internal)
  {
    int parser_type = 0;
    bool is_function = false;
    switch (sym_type)
      {
      case VARIABLE:
	parser_type = Parser::token::VAR;
	break;
      case STRING_VARIABLE:
	parser_type = Parser::token::SVAR;
	break;
      case ARRAY_VARIABLE:
	parser_type = Parser::token::AVAR;
	break;
      case IMMUTABLE_VARIABLE:
	parser_type = Parser::token::IMMVAR;
	break;
      case IMMUTABLE_STRING_VARIABLE:
	parser_type = Parser::token::IMMSVAR;
	break;
      case UNDEFINED_VARIABLE:
	parser_type = Parser::token::UNDVAR;
	break;
      case FUNCTION:
	parser_type = Parser::token::FNCT;
	is_function = true;
	break;
      case STRING_FUNCTION:
	parser_type = Parser::token::SFNCT;
	is_function = true;
	break;
      case ARRAY_FUNCTION:
	parser_type = Parser::token::AFNCT;
	is_function = true;
	break;
      }

    // If the type is a function type, it can be overloaded as long as
    // it returns the same type which means that the "parser_type" is
    // the same.  If we have a function, see if it has already been
    // defined and if so, check that the parser_type matches and then
    // retrn that pointer instead of creating a new symrec.

    if (is_function) {
      symrec *ptr = getsym(sym_name.c_str());
      if (ptr != NULL) {
	if (ptr->type != parser_type) {
    error("Overloaded function " + sym_name + "does not return same type", false);
	  exit(EXIT_FAILURE);
	}
	// Function with this name already exists; return that
	// pointer.
	// Note that the info and syntax fields will contain the
	// latest values, not the firstt...
	return ptr;
      }
    }
    
    symrec *ptr = new symrec(sym_name, parser_type, is_internal);
    if (ptr == NULL)
      return NULL;
  
    unsigned hashval = hash_symbol (ptr->name.c_str());
    ptr->next = sym_table[hashval];
    sym_table[hashval] = ptr;
    return ptr;
  }

  bool Aprepro::set_option(const std::string &option)
  {
    // Option should be of the form "--option" or "-O"
    // I.e., double dash for long option, single dash for short.
    // Option can be followed by "=value" if the option requires
    // a value.

    // Some options (--include) 
    
    if (option == "--debug" || option == "-d") {
      ap_options.debugging = true;
    }
    else if (option == "--version" || option == "-v") {
      std::cerr << "Aprepro version " << version() << "\n";
    }
    else if (option == "--nowarning" || option == "-W") {
      ap_options.warning_msg = false;
    }
    else if (option == "--copyright" || option == "-C") {
      copyright();
    }
    else if (option == "--message" || option == "-M") {
      ap_options.info_msg = true;
    }
    else if (option == "--immutable" || option == "-X") {
      ap_options.immutable = true;
      stateImmutable = true;
    }
    else if (option == "--trace" || option == "-t") {
      ap_options.trace_parsing = true;
    }
    else if (option == "--interactive" || option == "-i") {
      ap_options.interactive = true;
    }
    else if (option == "--one_based_index" || option == "-1") {
      ap_options.one_based_index = true;
    }
    else if (option == "--exit_on" || option == "-e") {
      ap_options.end_on_exit = true;
    }
    else if (option.find("--include=") != std::string::npos || (option[1] == 'I' && option[2] == '=')) {
      size_t index = option.find_first_of('=');
      std::string file_or_path = option.substr(index+1);
      if (is_directory(file_or_path)) {
	ap_options.include_path = file_or_path;
      } else {
	ap_options.include_file = file_or_path;
      }
    }
    else if (option == "--keep_history" || option == "-k") {
      ap_options.keep_history = true;
    }
    else if (option == "--help" || option == "-h") {
      std::cerr << "\nAPREPRO PREPROCESSOR OPTIONS:\n"
		<< "          --debug or -d: Dump all variables, debug loops/if/endif\n"
		<< "        --version or -v: Print version number to stderr          \n"
		<< "      --immutable or -X: All variables are immutable--cannot be modified\n"
	        << "--one_based_index or -1: Array indexing is one-based (default = zero-based)\n"
		<< "    --interactive or -i: Interactive use, no buffering           \n"
		<< "    --include=P or -I=P: Include file or include path            \n"
		<< "                       : If P is path, then optionally prepended to all include filenames\n"
		<< "                       : If P is file, then processed before processing input file\n"
		<< "        --exit_on or -e: End when 'Exit|EXIT|exit' entered       \n"
		<< "           --help or -h: Print this list                         \n"
		<< "        --message or -M: Print INFO messages                     \n"
		<< "      --nowarning or -W: Do not print WARN messages              \n"
    << "      --copyright or -C: Print copyright message                 \n"
    << "   --keep_history or -k: Keep a history of aprepro substitutions.\n\n"
	        << "\tUnits Systems: si, cgs, cgs-ev, shock, swap, ft-lbf-s, ft-lbm-s, in-lbf-s\n"
		<< "\tEnter {DUMP_FUNC()} for list of functions recognized by aprepro\n"
		<< "\tEnter {DUMP_PREVAR()} for list of predefined variables in aprepro\n\n"
		<< "\t->->-> Send email to gdsjaar@sandia.gov for aprepro support.\n\n";
    }
    return true;
  }
  
  void Aprepro::add_variable(const std::string &sym_name, const std::string &sym_value, bool immutable)
  {
    if (check_valid_var(sym_name.c_str())) {
    SYMBOL_TYPE type = immutable ? IMMUTABLE_STRING_VARIABLE : STRING_VARIABLE;
    symrec *var = putsym(sym_name, type, false);
    char *tmp = NULL;
    new_string(sym_value.c_str(), &tmp);
    var->value.svar = tmp;
    }
    else {
      warning("Invalid variable name syntax '" + sym_name + "'. Variable not defined.\n", false);
    }
  }

  void Aprepro::add_variable(const std::string &sym_name, double sym_value, bool immutable)
  {
    if (check_valid_var(sym_name.c_str())) {
      SYMBOL_TYPE type = immutable ? IMMUTABLE_VARIABLE : VARIABLE;
      symrec *var = putsym(sym_name, type, false);
      var->value.var = sym_value;
    }
    else {
      warning("Invalid variable name syntax '" + sym_name + "'. Variable not defined.\n", false);
    }
  }

  std::vector<std::string> Aprepro::get_variable_names(bool doInternal)
  {
    std::vector<std::string> names;

    for(unsigned int hashval = 0; hashval < HASHSIZE; hashval++)
    {
      for(symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next)
      {
        if(ptr->isInternal != doInternal)
          continue;

        switch(ptr->type)
        {
          case Parser::token::VAR:
          case Parser::token::IMMVAR:
          case Parser::token::SVAR:
          case Parser::token::IMMSVAR:
          case Parser::token::AVAR:
            // Add to our vector
            names.push_back(ptr->name);
            break;

          default:
            // Do nothing
            break;
        }
      }
    }

    return names;
  }

  void Aprepro::remove_variable(const std::string &sym_name)
  {
    symrec *ptr = getsym(sym_name.c_str());
    bool is_valid_variable =
        (ptr != NULL) && (!ptr->isInternal) &&
        ((ptr->type == Parser::token::VAR) ||
         (ptr->type == Parser::token::SVAR) ||
         (ptr->type == Parser::token::AVAR) ||
         (ptr->type == Parser::token::IMMVAR) ||
         (ptr->type == Parser::token::IMMSVAR) ||
         (ptr->type == Parser::token::UNDVAR));

    if(is_valid_variable)
    {
      int hashval = hash_symbol(sym_name.c_str());
      for (symrec *ptr = sym_table[hashval]; ptr != NULL; ) {
        symrec *save = ptr;
        ptr = ptr->next;
        delete save;
      }

      sym_table[hashval] = NULL;
    }
    else
    {
      warning("Variable '" + sym_name + "' not defined.\n", false);
    }
  }

  symrec *Aprepro::getsym (const char *sym_name) const
  {
    symrec *ptr = NULL;
    for (ptr = sym_table[hash_symbol (sym_name)]; ptr != NULL; ptr = ptr->next)
      if (strcmp (ptr->name.c_str(), sym_name) == 0)
	return ptr;
    return NULL;
  }

  void Aprepro::dumpsym (int type, bool doInternal) const
  {
    const char *comment = getsym("_C_")->value.svar;
    int width = 10; // controls spacing/padding for the variable names
  
    if (type == Parser::token::VAR || type == Parser::token::SVAR || type == Parser::token::AVAR) {
      (*infoStream) << "\n" << comment << "   Variable    = Value" << std::endl;

      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
	for (symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next) {
	  if ((doInternal && ptr->isInternal) || (!doInternal && !ptr->isInternal)) {
	    if (ptr->type == Parser::token::VAR)
        (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name <<
                         "\t= " << std::setprecision(10) <<  ptr->value.var << "}" << std::endl;
	    else if (ptr->type == Parser::token::IMMVAR)
        (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name <<
                         "\t= " << std::setprecision(10) << ptr->value.var << "}\t(immutable)" << std::endl;
      else if (ptr->type == Parser::token::SVAR)
        (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name <<
                         "\t= \"" << ptr->value.svar << "\"}" << std::endl;
	    else if (ptr->type == Parser::token::IMMSVAR)
        (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name <<
                         "\t= \"" << ptr->value.svar << "\"}\t(immutable)" << std::endl;
	    else if (ptr->type == Parser::token::AVAR) {
	      array *arr = ptr->value.avar;
        (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name <<
                         "\t (array) rows = " << arr->rows << ", cols = " << arr->cols <<
                         "} " << std::endl;
	    }
	  }
	}
      }
    }
    else if (type == Parser::token::FNCT || type == Parser::token::SFNCT || type == Parser::token::AFNCT) {
      (*infoStream) << "\nFunctions returning double:" << std::endl;
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
	for (symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next) {
	  if (ptr->type == Parser::token::FNCT) {
      (*infoStream) << std::left << std::setw(2*width) << ptr->syntax <<
                       ":  " << ptr->info << std::endl;
	  }
	}
      }

      (*infoStream) << "\nFunctions returning string:" << std::endl;
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
	for (symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next) {
	  if (ptr->type == Parser::token::SFNCT) {
      (*infoStream) << std::left << std::setw(2*width) << ptr->syntax <<
                       ":  " << ptr->info << std::endl;
	  }
	}
      }
      
      (*infoStream) << "\nFunctions returning array:" << std::endl;
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
	for (symrec *ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next) {
	  if (ptr->type == Parser::token::AFNCT) {
      (*infoStream) << std::left << std::setw(2*width) << ptr->syntax <<
                       ":  " << ptr->info << std::endl;
	  }
	}
      }
    }
  }

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

#define MAXLEN 16
  void Aprepro::statistics(std::ostream *out) const
  {
    std::ostream *output = out;
    if (output == NULL)
      output = &std::cout;
    
    symrec *ptr;
    unsigned entries = 0;
    int maxlen = 0;
    int minlen = INT_MAX;
    int lengths[MAXLEN];
    int longer = 0;
    double hash_ratio = 0.0;

    Stats stats;
    
    memset ((void *) lengths, 0, sizeof (lengths));

    for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
      int chain_len = 0;
      for (ptr = sym_table[hashval]; ptr != NULL; ptr = ptr->next)
	chain_len++;

      hash_ratio += chain_len * (chain_len + 1.0);
      entries += chain_len;
      if (chain_len >= MAXLEN)
	++longer;
      else
	++lengths[chain_len];

      minlen = min (minlen, chain_len);
      maxlen = max (maxlen, chain_len);

      if (chain_len > 0)
	stats.newsample (chain_len);
    }

    hash_ratio = hash_ratio / ((float) entries / HASHSIZE *
			       (float) (entries + 2.0 * HASHSIZE - 1.0));
    (*output) << entries << " entries in " << HASHSIZE << " element hash table, "
	      << lengths[0] << " (" << ((double) lengths[0] / HASHSIZE) * 100.0 << "%) empty.\n"
	      << "Hash ratio = " << hash_ratio << "\n"
	      << "Mean (nonempty) chain length = " << stats.mean()
	      << ", max = " << maxlen
	      << ", min = " << minlen
	      << ", deviation = " << stats.deviation() << "\n";

    for (int i = 0; i < MAXLEN; i++)
      if (lengths[i])
	(*output) << lengths[i] << " chain(s) of length " << i << "\n";

    if (longer)
      (*output) << longer << " chain(s) of length " << MAXLEN << " or longer\n";
  }

  void Aprepro::copyright(std::ostream *out) const
  {
    std::ostream *output = out;
    if (output == NULL)
      output = &std::cout;
    
    (*output) << "COPYRIGHT NOTICE\n";
  }

  void Aprepro::add_history(const std::string& original, const std::string& substitution)
  {
    if(!ap_options.keep_history)
      return;

    if(!original.empty())
    {
      history_data hist;
      hist.original = original;
      hist.substitution = substitution;
      hist.index = outputStream.top()->tellp();

      history.push_back(hist);
    }
  }

  const std::vector<history_data>& Aprepro::get_history()
  {
    return history;
  }

  void Aprepro::clear_history()
  {
    if(ap_options.keep_history)
      history.clear();
  }
}
