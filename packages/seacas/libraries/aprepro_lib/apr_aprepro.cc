// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#include "apr_scanner.h"    // for Scanner
#include "apr_stats.h"      // for Stats
#include "aprepro.h"        // for Aprepro, symrec, file_rec, etc
#include "aprepro_parser.h" // for Parser, Parser::token, etc
#include "terminal_color.h"
#include <climits>  // for INT_MAX
#include <cstddef>  // for size_t
#include <cstdlib>  // for exit, EXIT_SUCCESS, etc
#include <cstring>  // for memset, strcmp
#include <fstream>  // for operator<<, basic_ostream, etc
#include <iomanip>  // for operator<<, setw, etc
#include <iostream> // for left, cerr, cout, streampos
#include <stack>    // for stack
#include <stdexcept>
#include <string> // for string, operator==, etc
#include <unistd.h>
#include <vector> // for allocator, vector

namespace {
  const unsigned int HASHSIZE       = 5939;
  const char *       version_string = "5.17 (2020/07/27)";

  void output_copyright();

  unsigned hash_symbol(const char *symbol)
  {
    unsigned hashval;
    for (hashval = 0; *symbol != '\0'; symbol++) {
      hashval = *symbol + 65599 * hashval;
    }
    return (hashval % HASHSIZE);
  }
} // namespace

namespace SEAMS {
  Aprepro *aprepro = nullptr; // A global for use in the library.  Clean this up...
  bool     echo    = true;

  Aprepro::Aprepro()
  {
    sym_table.resize(HASHSIZE);
    ap_file_list.push(file_rec());
    init_table("$");
    aprepro = this;
  }

  Aprepro::~Aprepro()
  {
    if (!outputStream.empty()) {
      outputStream.top()->flush();
    }

    if ((stringScanner != nullptr) && stringScanner != lexer) {
      delete stringScanner;
    }

    delete lexer;

    for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
      for (symrec *ptr = sym_table[hashval]; ptr != nullptr;) {
        symrec *save = ptr;
        ptr          = ptr->next;
        if (save->type == Parser::token::AVAR) {
          delete save->value.avar;
        }
        delete save;
      }
    }
    aprepro = nullptr;
    cleanup_memory();
  }

  std::string Aprepro::version() const { return version_string; }

  void Aprepro::clear_results()
  {
    parsingResults.str("");
    parsingResults.clear();
  }

  bool Aprepro::parse_stream(std::istream &in, const std::string &in_name)
  {
    std::istream *in_cpy    = &in;
    ap_file_list.top().name = in_name;

    auto scanner = new Scanner(*this, in_cpy, &parsingResults);
    this->lexer  = scanner;

    if (!ap_options.include_file.empty()) {
      stateImmutable = true;
      echo           = 0;
      scanner->add_include_file(ap_options.include_file, true);
    }

    Parser parser(*this);
    parser.set_debug_level(ap_options.trace_parsing);
    return (parser.parse() == 0);
  }

  bool Aprepro::parse_file(const std::string &filename)
  {
    std::ifstream in(filename);
    if (!in.good()) {
      return false;
    }
    return parse_stream(in, filename);
  }

  bool Aprepro::parse_string(const std::string &input, const std::string &sname)
  {
    std::istringstream iss(input);
    return parse_stream(iss, sname);
  }

  bool Aprepro::parse_strings(const std::vector<std::string> &input, const std::string &sname)
  {
    std::stringstream iss;
    for (auto &elem : input) {
      iss << elem << '\n';
    }
    return parse_stream(iss, sname);
  }

  bool Aprepro::parse_string_interactive(const std::string &input)
  {
    stringInteractive = true;

    if (ap_file_list.size() == 1) {
      ap_file_list.top().name = "interactive_input";
    }
    if (stringScanner == nullptr) {
      stringScanner = new Scanner(*this, &stringInput, &parsingResults);
    }

    if (!ap_options.include_file.empty()) {
      stateImmutable = true;
      echo           = 0;
      stringScanner->add_include_file(ap_options.include_file, true);
    }

    this->lexer = stringScanner;

    stringInput.str(input);
    stringInput.clear();

    Parser parser(*this);
    parser.set_debug_level(ap_options.trace_parsing);
    bool result = parser.parse() == 0;

    stringInteractive = false;
    return result;
  }

  void Aprepro::error(const std::string &msg, bool line_info, bool prefix) const
  {
    bool              colorize = (errorStream == &std::cerr) && isatty(fileno(stderr));
    std::stringstream ss;
    if (prefix) {
      if (colorize) {
        (*errorStream) << trmclr::red;
      }
      (*errorStream) << "Aprepro: ERROR: ";
    }

    ss << msg;
    if (prefix && colorize) {
      ss << trmclr::normal;
    }

    if (line_info) {
      ss << " (" << ap_file_list.top().name << ", line " << ap_file_list.top().lineno + 1 << ")";
    }
    ss << "\n";

    // Send it to the user defined stream
    (*errorStream) << ss.str();
    parseErrorCount++;
  }

  void Aprepro::warning(const std::string &msg, bool line_info, bool prefix) const
  {
    if (!ap_options.warning_msg) {
      return;
    }

    bool              colorize = (warningStream == &std::cerr) && isatty(fileno(stderr));
    std::stringstream ss;
    if (prefix) {
      if (colorize) {
        (*warningStream) << trmclr::yellow;
      }
      (*warningStream) << "Aprepro: WARNING: ";
    }

    ss << msg;

    if (prefix && colorize) {
      ss << trmclr::normal;
    }

    if (line_info) {
      ss << " (" << ap_file_list.top().name << ", line " << ap_file_list.top().lineno + 1 << ")";
    }
    ss << "\n";

    // Send it to the user defined stream
    (*warningStream) << ss.str();
    parseWarningCount++;
  }

  void Aprepro::info(const std::string &msg, bool line_info, bool prefix) const
  {
    if (!ap_options.info_msg) {
      return;
    }

    bool              colorize = (infoStream == &std::cerr) && isatty(fileno(stderr));
    std::stringstream ss;
    if (prefix) {
      if (colorize) {
        (*infoStream) << trmclr::blue;
      }
      (*infoStream) << "Aprepro: INFO: ";
    }
    ss << msg;

    if (prefix && colorize) {
      ss << trmclr::normal;
    }

    if (line_info) {
      ss << " (" << ap_file_list.top().name << ", line " << ap_file_list.top().lineno + 1 << ")";
    }
    ss << "\n";

    // Send it to the user defined stream
    (*infoStream) << ss.str();
  }

  void Aprepro::set_error_streams(std::ostream *c_error, std::ostream *c_warning,
                                  std::ostream *c_info)
  {
    if (c_error != nullptr) {
      errorStream = c_error;
    }
    if (c_warning != nullptr) {
      warningStream = c_warning;
    }
    if (c_info != nullptr) {
      infoStream = c_info;
    }
  }

  /* Two methods for opening files:

     In OPEN_FILE, the file must exist or else the code will exit in
     batch mode or return null pointer if interactive.

     In CHECK_OPEN_FILE, it is OK if the file does not exist. A Null
     'pointer' will be returned.
  */
  std::fstream *Aprepro::open_file(const std::string &file, const char *mode)
  {
    std::ios::openmode smode = std::ios::in;
    if (mode[0] == 'w') {
      smode = std::ios::out;
    }

    /* See if file exists in current directory (or as specified) */
    auto pointer = new std::fstream(file, smode);
    if ((pointer == nullptr || pointer->bad() || !pointer->good()) &&
        !ap_options.include_path.empty()) {
      /* If there is an include path specified, try opening file there */
      std::string file_path(ap_options.include_path);
      file_path += "/";
      file_path += file;
      delete pointer;
      pointer = new std::fstream(file_path, smode);
    }

    /* If pointer still null, print error message */
    if (pointer == nullptr || pointer->bad() || !pointer->good()) {
      std::string err = "Can't open " + file;
      error(err, false);
      delete pointer;
      pointer = nullptr;
      if (!stringInteractive) {
        throw std::runtime_error(err);
      }
    }

    return pointer;
  }

  std::fstream *Aprepro::check_open_file(const std::string &file, const char *mode)
  {
    std::ios::openmode smode = std::ios::in;
    if (mode[0] == 'w') {
      smode = std::ios::out;
    }

    auto pointer = new std::fstream(file, smode);

    if ((pointer == nullptr || pointer->bad() || !pointer->good()) &&
        !ap_options.include_path.empty()) {
      /* If there is an include path specified, try opening file there */
      std::string file_path(ap_options.include_path);
      file_path += "/";
      file_path += file;
      delete pointer;
      pointer = new std::fstream(file_path, smode);
    }
    return pointer;
  }

  symrec *Aprepro::putsym(const std::string &sym_name, SYMBOL_TYPE sym_type, bool is_internal)
  {
    int  parser_type = 0;
    bool is_function = false;
    switch (sym_type) {
    case SYMBOL_TYPE::VARIABLE: parser_type = Parser::token::VAR; break;
    case SYMBOL_TYPE::STRING_VARIABLE: parser_type = Parser::token::SVAR; break;
    case SYMBOL_TYPE::ARRAY_VARIABLE: parser_type = Parser::token::AVAR; break;
    case SYMBOL_TYPE::IMMUTABLE_VARIABLE: parser_type = Parser::token::IMMVAR; break;
    case SYMBOL_TYPE::IMMUTABLE_STRING_VARIABLE: parser_type = Parser::token::IMMSVAR; break;
    case SYMBOL_TYPE::UNDEFINED_VARIABLE: parser_type = Parser::token::UNDVAR; break;
    case SYMBOL_TYPE::FUNCTION:
      parser_type = Parser::token::FNCT;
      is_function = true;
      break;
    case SYMBOL_TYPE::STRING_FUNCTION:
      parser_type = Parser::token::SFNCT;
      is_function = true;
      break;
    case SYMBOL_TYPE::ARRAY_FUNCTION:
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
      symrec *ptr = getsym(sym_name);
      if (ptr != nullptr) {
        if (ptr->type != parser_type) {
          std::string err = "Overloaded function '" + sym_name + "' does not return same type";
          error(err, false);
          throw std::runtime_error(err);
        }
        // Function with this name already exists; return that
        // pointer.
        // Note that the info and syntax fields will contain the
        // latest values, not the firstt...
        return ptr;
      }
    }
    else {
      symrec *ptr = getsym(sym_name);
      if (ptr != nullptr) {
        std::string err = "Internal Error: Variable '" + sym_name + "' is already defined.";
        error(err, false);
        throw std::runtime_error(err);
      }
    }

    auto ptr = new symrec(sym_name, parser_type, is_internal);

    unsigned hashval   = hash_symbol(ptr->name.c_str());
    ptr->next          = sym_table[hashval];
    sym_table[hashval] = ptr;
    return ptr;
  }

  int Aprepro::set_option(const std::string &option, const std::string &optional_value)
  {
    // Option should be of the form "--option" or "-O"
    // I.e., double dash for long option, single dash for short.
    // Option can be followed by "=value" if the option requires
    // a value.  The value will either be part of 'option' if it contains
    // an '=', or it will be in 'optional_value'.  if 'optional_value' used,
    // return 1; else return 0;

    // Some options (--include)
    int ret_value = 0;

    if (option == "--debug" || option == "-d") {
      ap_options.debugging = true;
    }
    if (option == "--dumpvars" || option == "-D") {
      ap_options.dumpvars = true;
    }
    else if (option == "--version" || option == "-v") {
      std::cerr << "Algebraic Preprocessor (Aprepro) version " << version() << "\n";
      exit(EXIT_SUCCESS);
    }
    else if (option == "--nowarning" || option == "-W") {
      ap_options.warning_msg = false;
    }
    else if (option == "--copyright" || option == "-C") {
      output_copyright();
      exit(EXIT_SUCCESS);
    }
    else if (option == "--message" || option == "-M") {
      ap_options.info_msg = true;
    }
    else if (option == "--immutable" || option == "-X") {
      ap_options.immutable = true;
      stateImmutable       = true;
    }
    else if (option == "--errors_fatal" || option == "-f") {
      ap_options.errors_fatal = true;
    }
    else if (option == "--errors_and_warnings_fatal" || option == "-F") {
      ap_options.errors_and_warnings_fatal = true;
      ap_options.errors_fatal              = true;
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
    else if (option.find("--include") != std::string::npos || (option[1] == 'I')) {
      std::string value;

      size_t index = option.find_first_of('=');
      if (index != std::string::npos) {
        value = option.substr(index + 1);
      }
      else {
        value     = optional_value;
        ret_value = 1;
      }

      if (is_directory(value)) {
        ap_options.include_path = value;
      }
      else {
        ap_options.include_file = value;
      }
    }
    else if (option == "--keep_history" || option == "-k") {
      ap_options.keep_history = true;
    }

    else if (option.find("--comment") != std::string::npos || (option[1] == 'c')) {
      std::string comment;
      // In short version, do not require equal sign (-c# -c// )
      if (option[1] == 'c' && option.length() > 2 && option[2] != '=') {
        comment = option.substr(2);
      }
      else {
        size_t index = option.find_first_of('=');
        if (index != std::string::npos) {
          comment = option.substr(index + 1);
        }
        else {
          comment   = optional_value;
          ret_value = 1;
        }
      }
      symrec *ptr = getsym("_C_");
      if (ptr != nullptr) {
        ptr->value.svar = comment;
      }
    }

    else if (option == "--help" || option == "-h") {
      std::cerr
          << "\nAprepro version " << version() << "\n"
          << "\nUsage: aprepro [options] [-I path] [-c char] [var=val] [filein] [fileout]\n"
          << "          --debug or -d: Dump all variables, debug loops/if/endif\n"
          << "       --dumpvars or -D: Dump all variables at end of run        \n"
          << "        --version or -v: Print version number to stderr          \n"
          << "      --immutable or -X: All variables are immutable--cannot be modified\n"
          << "   --errors_fatal or -f: Exit program with nonzero status if errors are "
             "encountered\n"
          << " --errors_and_warnings_fatal or -F: Exit program with nonzero status if "
             "warnings are encountered\n"
          << "--one_based_index or -1: Array indexing is one-based (default = zero-based)\n"
          << "    --interactive or -i: Interactive use, no buffering           \n"
          << "    --include=P or -I=P: Include file or include path            \n"
          << "                       : If P is path, then optionally prepended to all include "
             "filenames\n"
          << "                       : If P is file, then processed before processing input file\n"
          << "                       : variables defined in P will be immutable.\n"
          << "        --exit_on or -e: End when 'Exit|EXIT|exit' entered       \n"
          << "           --help or -h: Print this list                         \n"
          << "        --message or -M: Print INFO messages                     \n"
          << "      --nowarning or -W: Do not print WARN messages              \n"
          << "  --comment=char or -c=char: Change comment character to 'char'      \n"
          << "      --copyright or -C: Print copyright message                 \n"
          << "   --keep_history or -k: Keep a history of aprepro substitutions.\n"
          << "                         (not for general interactive use)       \n"
          << "          --quiet or -q: (deprecated, option is ignored)         \n"
          << "                var=val: Assign value 'val' to variable 'var'    \n"
          << "                         Use var=\\\"sval\\\" for a string variable\n\n"
          << "\tUnits Systems: si, cgs, cgs-ev, shock, swap, ft-lbf-s, ft-lbm-s, in-lbf-s\n"
          << "\tEnter {DUMP()} for list of user-defined variables\n"
          << "\tEnter {DUMP_FUNC()} for list of functions recognized by aprepro\n"
          << "\tEnter {DUMP_PREVAR()} for list of predefined variables in aprepro\n\n"
          << "\t->->-> Send email to gdsjaar@sandia.gov for aprepro support.\n\n";
      exit(EXIT_SUCCESS);
    }
    return ret_value;
  }

  void Aprepro::add_variable(const std::string &sym_name, const std::string &sym_value,
                             bool immutable, bool internal)
  {
    if (check_valid_var(sym_name.c_str())) {
      SYMBOL_TYPE type =
          immutable ? SYMBOL_TYPE::IMMUTABLE_STRING_VARIABLE : SYMBOL_TYPE::STRING_VARIABLE;
      symrec *var = getsym(sym_name);
      if (var == nullptr) {
        var = putsym(sym_name, type, internal);
      }
      else {
        var->type = immutable ? Parser::token::IMMSVAR : Parser::token::SVAR;
      }
      var->value.svar = sym_value;
    }
    else {
      warning("Invalid variable name syntax '" + sym_name + "'. Variable not defined.\n", false);
    }
  }

  void Aprepro::add_variable(const std::string &sym_name, double sym_value, bool immutable,
                             bool internal)
  {
    if (check_valid_var(sym_name.c_str())) {
      SYMBOL_TYPE type = immutable ? SYMBOL_TYPE::IMMUTABLE_VARIABLE : SYMBOL_TYPE::VARIABLE;
      symrec *    var  = getsym(sym_name);
      if (var == nullptr) {
        var = putsym(sym_name, type, internal);
      }
      else {
        var->type = immutable ? Parser::token::IMMVAR : Parser::token::VAR;
      }
      var->value.var = sym_value;
    }
    else {
      warning("Invalid variable name syntax '" + sym_name + "'. Variable not defined.\n", false);
    }
  }

  void Aprepro::add_variable(const std::string &sym_name, array *value)
  {
    if (check_valid_var(sym_name.c_str())) {
      SYMBOL_TYPE type = SYMBOL_TYPE::ARRAY_VARIABLE;
      symrec *    var  = getsym(sym_name);
      if (var == nullptr) {
        var = putsym(sym_name, type, false);
      }
      else {
        var->type = Parser::token::AVAR;
      }
      var->value.avar = value;
    }
    else {
      warning("Invalid variable name syntax '" + sym_name + "'. Variable not defined.\n", false);
    }
  }

  std::vector<std::string> Aprepro::get_variable_names(bool doInternal)
  {
    std::vector<std::string> names;

    for (unsigned int hashval = 0; hashval < HASHSIZE; hashval++) {
      for (symrec *ptr = sym_table[hashval]; ptr != nullptr; ptr = ptr->next) {
        if (ptr->isInternal != doInternal) {
          continue;
        }

        switch (ptr->type) {
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
    symrec *ptr = getsym(sym_name);
    bool    is_valid_variable =
        (ptr != nullptr) && (!ptr->isInternal) &&
        ((ptr->type == Parser::token::VAR) || (ptr->type == Parser::token::SVAR) ||
         (ptr->type == Parser::token::AVAR) || (ptr->type == Parser::token::IMMVAR) ||
         (ptr->type == Parser::token::IMMSVAR) || (ptr->type == Parser::token::UNDVAR));

    if (is_valid_variable) {
      int     hashval  = hash_symbol(sym_name.c_str());
      symrec *hash_ptr = sym_table[hashval];

      // Handle the case if the variable we want to delete is first in the
      // linked list.
      if (ptr == hash_ptr) {
        // NOTE: If ptr is the only thing in the linked list, ptr->next will be
        // nullptr, which is what we want in sym_table when we delete ptr.
        sym_table[hashval] = ptr->next;
        delete ptr;
      }

      // Handle the case where the variable we want to delete is somewhere
      // in the middle or at the end of the linked list.
      else {
        // Find the preceding ptr (singly linked list).
        // NOTE: We don't have a check for nullptr here because the fact that
        // ptr != hash_ptr tells us that we must have more than one item in our
        // linked list, in which case hash_ptr->next will not be nullptr until we
        // reach the end of the list. hash_ptr->next should be equal to ptr
        // before that happens.
        while (hash_ptr->next != ptr) {
          hash_ptr = hash_ptr->next;
        }

        // NOTE: If ptr is at the end of the list ptr->next will be nullptr, in
        // which case this will change hash_ptr to be the end of the list.
        hash_ptr->next = ptr->next;
        delete ptr;
      }
    }
    else {
      warning("Variable '" + sym_name + "' not defined.\n", false);
    }
  }

  symrec *Aprepro::getsym(const char *sym_name) const
  {
    symrec *ptr = nullptr;
    for (ptr = sym_table[hash_symbol(sym_name)]; ptr != nullptr; ptr = ptr->next) {
      if (strcmp(ptr->name.c_str(), sym_name) == 0) {
        return ptr;
      }
    }
    return nullptr;
  }

  symrec *Aprepro::getsym(const std::string &sym_name) const { return getsym(sym_name.c_str()); }

  void Aprepro::dumpsym(const char *type, bool doInternal) const
  {
    if (type[0] == 'v') {
      dumpsym(Parser::token::VAR, doInternal);
    }
    else if (type[0] == 'f') {
      dumpsym(Parser::token::FNCT, doInternal);
    }
  }

  void Aprepro::dumpsym(int type, bool doInternal) const { dumpsym(type, nullptr, doInternal); }

  void Aprepro::dumpsym(int type, const char *pre, bool doInternal) const
  {
    std::string comment = getsym("_C_")->value.svar;
    int         width   = 10; // controls spacing/padding for the variable names
    int         fwidth  = 20; // controls spacing/padding for the function names
    std::string spre;

    if (pre) {
      spre = pre;
    }

    if (type == Parser::token::VAR || type == Parser::token::SVAR || type == Parser::token::AVAR) {
      (*infoStream) << "\n" << comment << "   Variable    = Value" << '\n';

      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
        for (symrec *ptr = sym_table[hashval]; ptr != nullptr; ptr = ptr->next) {
          if (pre == nullptr || ptr->name.find(spre) != std::string::npos) {
            if (doInternal == ptr->isInternal) {
              if (ptr->type == Parser::token::VAR) {
                (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name
                              << "\t= " << std::setprecision(10) << ptr->value.var << "}" << '\n';
              }
              else if (ptr->type == Parser::token::IMMVAR) {
                (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name
                              << "\t= " << std::setprecision(10) << ptr->value.var
                              << "}\t(immutable)" << '\n';
              }
              else if (ptr->type == Parser::token::SVAR) {
                if (strchr(ptr->value.svar.c_str(), '\n') != nullptr ||
                    strchr(ptr->value.svar.c_str(), '"') != nullptr) {
                  (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name
                                << "\t= '" << ptr->value.svar << "'}" << '\n';
                }
                else {
                  (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name
                                << "\t= \"" << ptr->value.svar << "\"}" << '\n';
                }
              }
              else if (ptr->type == Parser::token::IMMSVAR) {
                if (strchr(ptr->value.svar.c_str(), '\n') != nullptr ||
                    strchr(ptr->value.svar.c_str(), '"') != nullptr) {
                  (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name
                                << "\t= '" << ptr->value.svar << "'}\t(immutable)" << '\n';
                }
                else {
                  (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name
                                << "\t= \"" << ptr->value.svar << "\"}\t(immutable)" << '\n';
                }
              }
              else if (ptr->type == Parser::token::AVAR) {
                array *arr = ptr->value.avar;
                (*infoStream) << comment << "  {" << std::left << std::setw(width) << ptr->name
                              << "\t (array) rows = " << arr->rows << ", cols = " << arr->cols
                              << "} " << '\n';
              }
            }
          }
        }
      }
    }
    else if (type == Parser::token::FNCT || type == Parser::token::SFNCT ||
             type == Parser::token::AFNCT) {
      (*infoStream) << trmclr::blue << "\nFunctions returning double:" << trmclr::normal << '\n';
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
        for (symrec *ptr = sym_table[hashval]; ptr != nullptr; ptr = ptr->next) {
          if (pre == nullptr || ptr->name.find(spre) != std::string::npos) {
            if (ptr->type == Parser::token::FNCT) {
              (*infoStream) << std::left << trmclr::green << std::setw(fwidth) << ptr->syntax
                            << trmclr::normal << ":  " << ptr->info << '\n';
            }
          }
        }
      }

      (*infoStream) << trmclr::blue << trmclr::blue
                    << "\nFunctions returning string:" << trmclr::normal << '\n';
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
        for (symrec *ptr = sym_table[hashval]; ptr != nullptr; ptr = ptr->next) {
          if (pre == nullptr || ptr->name.find(spre) != std::string::npos) {
            if (ptr->type == Parser::token::SFNCT) {
              (*infoStream) << std::left << trmclr::green << std::setw(fwidth) << ptr->syntax
                            << trmclr::normal << ":  " << ptr->info << '\n';
            }
          }
        }
      }

      (*infoStream) << trmclr::blue << "\nFunctions returning array:" << trmclr::normal << '\n';
      for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
        for (symrec *ptr = sym_table[hashval]; ptr != nullptr; ptr = ptr->next) {
          if (pre == nullptr || ptr->name.find(spre) != std::string::npos) {
            if (ptr->type == Parser::token::AFNCT) {
              (*infoStream) << std::left << trmclr::green << std::setw(fwidth) << ptr->syntax
                            << trmclr::normal << ":  " << ptr->info << '\n';
            }
          }
        }
      }
    }
  }

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))

#define MAXLEN 16
  void Aprepro::statistics(std::ostream *out) const
  {
    std::ostream *output = out;
    if (output == nullptr) {
      output = &std::cout;
    }

    symrec * ptr;
    unsigned entries = 0;
    int      maxlen  = 0;
    int      minlen  = INT_MAX;
    int      lengths[MAXLEN];
    int      longer = 0;

    Stats stats;

    memset((void *)lengths, 0, sizeof(lengths));

    for (unsigned hashval = 0; hashval < HASHSIZE; hashval++) {
      int chain_len = 0;
      for (ptr = sym_table[hashval]; ptr != nullptr; ptr = ptr->next) {
        chain_len++;
      }

      entries += chain_len;
      if (chain_len >= MAXLEN) {
        ++longer;
      }
      else {
        ++lengths[chain_len];
      }

      minlen = min(minlen, chain_len);
      maxlen = max(maxlen, chain_len);

      if (chain_len > 0) {
        stats.newsample(chain_len);
      }
    }

    (*output) << entries << " entries in " << HASHSIZE << " element hash table, " << lengths[0]
              << " (" << (static_cast<double>(lengths[0]) / HASHSIZE) * 100.0 << "%) empty.\n"
              << "Mean (nonempty) chain length = " << stats.mean() << ", max = " << maxlen
              << ", min = " << minlen << ", deviation = " << stats.deviation() << "\n";

    for (int i = 0; i < MAXLEN; i++) {
      if (lengths[i] != 0) {
        (*output) << lengths[i] << " chain(s) of length " << i << "\n";
      }
      if (longer != 0) {
        (*output) << longer << " chain(s) of length " << MAXLEN << " or longer\n";
      }
    }
  }

  void Aprepro::add_history(const std::string &original, const std::string &substitution)
  {
    if (!ap_options.keep_history) {
      return;
    }

    if (!original.empty()) {
      history_data hist{};
      hist.original     = original;
      hist.substitution = substitution;
      hist.index        = outputStream.empty() ? std::streampos(0) : outputStream.top()->tellp();

      history.push_back(hist);
    }
  }

  const std::vector<history_data> &Aprepro::get_history() { return history; }

  void Aprepro::clear_history()
  {
    if (ap_options.keep_history) {
      history.clear();
    }
  }
} // namespace SEAMS

namespace {
  void output_copyright()
  {
    std::cerr << "\n\tCopyright (c) 2014-2017 National Technology & Engineering Solutions\n"
              << "of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with\n"
              << "NTESS, the U.S. Government retains certain rights in this software.\n"
              << "\n"
              << "\tRedistribution and use in source and binary forms, with or without\n"
              << "\tmodification, are permitted provided that the following conditions\n"
              << "\tare met:\n"
              << "\n"
              << "\t   * Redistributions of source code must retain the above copyright\n"
              << "\t     notice, this list of conditions and the following disclaimer.\n"
              << "\t   * Redistributions in binary form must reproduce the above\n"
              << "\t     copyright notice, this list of conditions and the following\n"
              << "\t     disclaimer in the documentation and/or other materials provided\n"
              << "\t     with the distribution.\n"
              << "\t   * Neither the name of NTESS nor the names of its\n"
              << "\t     contributors may be used to endorse or promote products derived\n"
              << "\t     from this software without specific prior written permission.\n"
              << "\n"
              << "\tTHIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n"
              << "\t'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n"
              << "\tLIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n"
              << "\tA PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT\n"
              << "\tOWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,\n"
              << "\tSPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT\n"
              << "\tLIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n"
              << "\tDATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY\n"
              << "\tTHEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
              << "\t(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE\n"
              << "\tOF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n";
  }
} // namespace
