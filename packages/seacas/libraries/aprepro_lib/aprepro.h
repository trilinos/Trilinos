// Copyright(C) 1999-2023,  National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

// Might be good to add a callback function which would be called
// when there was output -- In LexerOutput for example.  Default
// could be to just write to std::cout or to resultsOutput stringstream...
#pragma once

#include <cstdlib>
#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include "apr_symrec.h"

#if defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) ||                \
    defined(__MINGW32__) || defined(_WIN64) || defined(__MINGW64__)
#include <io.h>
#define isatty _isatty
#endif

/** The SEAMS namespace is used to encapsulate the three parser classes
 * SEAMS::Parser, SEAMS::Scanner and SEAMS::Aprepro */
namespace SEAMS {

  struct Symtable;

  /* Global options */
  struct aprepro_options
  {
    std::string include_path{};
    std::string include_file{};
    bool        end_on_exit{false};
    bool        errors_fatal{false};
    bool        errors_and_warnings_fatal{false};
    bool        require_defined{false}; // flag to treat undefined vars as errors
    bool        warning_msg{true};
    bool        info_msg{false};
    bool        debugging{false};
    bool        dumpvars{false};
    bool        dumpvars_json{false};
    bool        interactive{false};
    bool        immutable{false};
    bool        trace_parsing{false}; // enable debug output in the bison parser
    bool        one_based_index{false};
    bool        keep_history{false}; // Flag to keep a history of Aprepro substitutions
    aprepro_options() = default;
  };

  /* Structure for holding file names and line counters */
  struct file_rec
  {
    std::string name{"STDIN"};
    symrec     *loop_index{nullptr};
    double      loop_increment{1};
    int         lineno{1};
    int         loop_count{0};
    int         loop_level{0};
    bool        tmp_file{false};

    file_rec(const char *my_name, int line_num, bool is_temp, int loop_cnt)
        : name(my_name), lineno(line_num), loop_count(loop_cnt), tmp_file(is_temp)
    {
    }
    file_rec() = default;
  };

  /* Structure for holding aprepro substitution info */
  struct history_data
  {
    std::string    original;
    std::string    substitution;
    std::streampos index; // Character index in the output where the substitution begins.
  };

  /** The Aprepro class brings together all components. It creates an instance of
   * the Parser and Scanner classes and connects them. Then the input stream is
   * fed into the scanner object and the parser gets it's token
   * sequence. Furthermore the aprepro object is available in the grammar rules as
   * a parameter. Therefore the aprepro class contains a reference to the
   * structure into which the parsed data is saved. */
  class Aprepro
  {
  public:
    /// construct a new parser aprepro context
    Aprepro();
    ~Aprepro();
    Aprepro(const Aprepro &)            = delete;
    Aprepro &operator=(const Aprepro &) = delete;

    enum class SYMBOL_TYPE {
      INTERNAL                  = 0,
      VARIABLE                  = 1,
      STRING_VARIABLE           = 2,
      UNDEFINED_VARIABLE        = 5,
      FUNCTION                  = 3,
      STRING_FUNCTION           = 4,
      ARRAY_FUNCTION            = 6,
      ARRAY_VARIABLE            = 7,
      IMMUTABLE_VARIABLE        = 11,
      IMMUTABLE_STRING_VARIABLE = 12
    };

    bool state_is_immutable() const { return stateImmutable; }

    /** Return an  std::ostringstream reference to get the results of
        the parse_* call (* = stream, file, or string).
    */
    const std::ostringstream &parsing_results() const { return parsingResults; }
    void                      clear_results();

    /** Return string representation of current version of aprepro + commit date.  */
    static const std::string &version();

    /** Return string representation of current version of aprepro.  */
    static const std::string &short_version();

    /** Return long version: `# Algebraic Preprocessor (Aprepro) version X.X (date)` */
    std::string long_version() const;

    /** Invoke the scanner and parser for a stream.
     * @param in        input stream
     * @param in_name   stream name for error messages
     * @return          true if successfully parsed
     */
    bool parse_stream(std::istream &in, const std::string &in_name = "stream input");

    /** Invoke the scanner and parser on an input string.
     * @param input     input string
     * @param sname     stream name for error messages
     * @return          true if successfully parsed
     */
    bool parse_string(const std::string &input, const std::string &sname = "string stream");

    /** Invoke the scanner and parser on a vector of strings.
     * @param input     vector of input strings
     * @param sname     stream name for error messages
     * @return          true if successfully parsed
     */
    bool parse_strings(const std::vector<std::string> &input, const std::string &sname);

    /** Invoke the scanner and parser on an input string in an
     * interactive manner.
     * @param input input stringInput
     * @return true if successfully parsed
     */
    bool parse_string_interactive(const std::string &input);

    /** Get the string interactive flag, which indicates if we are in
     * the middle of parsing a string in an interactive manner.
     */
    bool string_interactive() const { return stringInteractive; }

    /** Invoke the scanner and parser on a file. Use parse_stream with a
     * std::ifstream if detection of file reading errors is required.
     * @param filename  input file name
     * @return          true if successfully parsed
     */
    bool parse_file(const std::string &filename);

    void statistics(std::ostream *out = nullptr) const;

    aprepro_options      ap_options;
    std::stack<file_rec> ap_file_list{};

    std::stack<std::ostream *> outputStream{};

    SYMBOL_TYPE    get_symbol_type(const SEAMS::symrec *symbol) const;
    SEAMS::symrec *getsym(const char *sym_name) const;
    SEAMS::symrec *getsym(const std::string &sym_name) const;
    SEAMS::symrec *putsym(const std::string &sym_name, SYMBOL_TYPE sym_type, bool is_internal);

    void add_variable(const std::string &sym_name, const std::string &sym_value,
                      bool immutable = false, bool internal = false);
    void add_variable(const std::string &sym_name, double sym_value, bool immutable = false,
                      bool internal = false);
    void add_variable(const std::string &sym_name, array *value);

    std::vector<std::string> get_variable_names(bool doInternal = false);
    void                     remove_variable(const std::string &sym_name);

    int set_option(const std::string &option, const std::string &optional_value = std::string(""));

    std::fstream *open_file(const std::string &file, const char *mode);
    std::fstream *check_open_file(const std::string &file, const char *mode);

    /** Pointer to the current lexer instance, this is used to connect the
     * parser to the scanner. It is used in the yylex macro. */
    class Scanner *lexer{nullptr};

    /** Error handling. */
    int get_error_count() const
    {
      return parseErrorCount;
    } /** Return number of errors reported during parse */
    int get_warning_count() const
    {
      return parseWarningCount;
    } /** Return number of warnings reported during parse */
    void error(const std::string &msg, bool line_info = true, bool prefix = true) const;
    void warning(const std::string &msg, bool line_info = true, bool prefix = true) const;
    void info(const std::string &msg, bool line_info = false, bool prefix = true) const;

    // The info stream. To only print out info messages if the -M option was
    // specified, use info(...) instead.
    bool          closeInfo{false};
    std::ostream *infoStream{&std::cout};

    void set_error_streams(std::ostream *c_error, std::ostream *c_warning, std::ostream *c_info);
    void set_error_streams(std::ostream *c_error, std::ostream *c_warning, std::ostream *c_info,
                           bool close_error, bool close_warning, bool close_info);

    void dumpsym(const char *type, bool doInternal) const;
    void dumpsym(int type, bool doInternal) const;
    void dumpsym_json() const;
    void dumpsym(int type, const char *pre, bool doInternal) const;

    array *make_array(int r, int c);
    array *make_array(const array &from);
    void   redefine_array(array *data);

  private:
    std::unique_ptr<Symtable> sym_table;

    void                 init_table(const char *comment);
    std::vector<array *> array_allocations{};
    std::ostringstream   parsingResults{};

    // Input stream used with parse_string_interactive
    std::istringstream stringInput{};

    bool           stringInteractive{false};
    class Scanner *stringScanner{nullptr};

    // For error handling
    std::ostream *errorStream{&std::cerr};
    std::ostream *warningStream{&std::cerr};
    bool          closeError{false};
    bool          closeWarning{false};

    // For substitution history.
    std::vector<history_data> history{};

    // For repeatble and user-friendly help/dump output.
    std::vector<SEAMS::symrec *> get_sorted_sym_table() const;

    mutable int parseErrorCount{0};
    mutable int parseWarningCount{0};

  public:
    bool stateImmutable{false};

    // Flag to do Aprepro substitutions within loops. Default value is true. If set to
    // false, content within the loop will be treated as verbatim text.
    bool doLoopSubstitution{true};

    // Flag to do Aprepro substitutions when including a file. Default value is true.
    // If set to false, content within the file will be treated as verbatim text that
    // needs to be sent through Aprepro again later.
    bool doIncludeSubstitution{true};

    // Flag to indicate whether Aprepro is in the middle of collecting lines for a
    // loop.
    bool isCollectingLoop{false};

    // Record the substitution of the current Aprepro statement. This function will also
    // reset the historyString and add an entry to the substitution map.
    void add_history(const std::string &original, const std::string &substitution);

    // Used to avoid undefined variable warnings in old ifdef/ifndef construct
    mutable bool inIfdefGetvar{false};

    const std::vector<history_data> &get_history();
    void                             clear_history();
  };

} // namespace SEAMS
