// Copyright(C) 1999-2020 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

// Might be good to add a callback function which would be called
// when there was output -- In LexerOutput for example.  Default
// could be to just write to std::cout or to resultsOutput stringstream...

#ifndef SEAMS_DRIVER_H
#define SEAMS_DRIVER_H

#include <cstdlib>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include <cmath>
#ifndef math_errhandling
#define math_errhandling MATH_ERRNO
#endif

#if defined(_MSC_VER)
#include <io.h>
#define isatty _isatty
#endif

/** The SEAMS namespace is used to encapsulate the three parser classes
 * SEAMS::Parser, SEAMS::Scanner and SEAMS::Aprepro */
namespace SEAMS {

  struct array
  {
    std::vector<double> data{};
    int                 rows{0};
    int                 cols{0};

    array(int r, int c) : rows(r), cols(c) { data.resize(r * c); }
    array()  = default;
    ~array() = default;
  };

  struct symrec
  {
    std::string name{};
    std::string syntax{};
    std::string info{};
    int         type;
    bool        isInternal;
    struct value
    {
      double var{0};
      double (*fnctptr)(){nullptr};
      double (*fnctptr_d)(double){nullptr};
      double (*fnctptr_c)(char *){nullptr};
      double (*fnctptr_dc)(double, char *){nullptr};
      double (*fnctptr_cd)(char *, double){nullptr};
      double (*fnctptr_cc)(char *, char *){nullptr};
      double (*fnctptr_dd)(double, double){nullptr};
      double (*fnctptr_ddd)(double, double, double){nullptr};
      double (*fnctptr_ccc)(char *, char *, char *){nullptr};
      double (*fnctptr_ccd)(char *, char *, double){nullptr};
      double (*fnctptr_dddd)(double, double, double, double){nullptr};
      double (*fnctptr_ddddc)(double, double, double, double, char *){nullptr};
      double (*fnctptr_dddddd)(double, double, double, double, double, double){nullptr};
      double (*fnctptr_a)(const array *){nullptr};
      std::string svar{};
      const char *(*strfnct)(){nullptr};
      const char *(*strfnct_c)(char *){nullptr};
      const char *(*strfnct_d)(double){nullptr};
      const char *(*strfnct_a)(const array *){nullptr};
      const char *(*strfnct_dd)(double, double){nullptr};
      const char *(*strfnct_cc)(char *, char *){nullptr};
      const char *(*strfnct_ccc)(char *, char *, char *){nullptr};
      const char *(*strfnct_dcc)(double, char *, char *){nullptr};
      const char *(*strfnct_dcccc)(double, char *, char *, char *, char *){nullptr};
      array *avar{nullptr}; /* Array Variable */
      array *(*arrfnct_c)(const char *){nullptr};
      array *(*arrfnct_cc)(const char *, const char *){nullptr};
      array *(*arrfnct_cd)(const char *, double){nullptr};
      array *(*arrfnct_ddd)(double, double, double){nullptr};
      array *(*arrfnct_dd)(double, double){nullptr};
      array *(*arrfnct_d)(double){nullptr};
      array *(*arrfnct_a)(const array *){nullptr};

      value() = default;
    } value;
    symrec *next;

    symrec(const char *my_name, int my_type, bool is_internal = false)
        : name(my_name), syntax(my_name), info("UNDEFINED"), type(my_type), isInternal(is_internal),
          next(nullptr)
    {
      value.var = 0;
    }

    symrec(const std::string &my_name, int my_type, bool is_internal = false)
        : name(my_name), syntax(my_name), info("UNDEFINED"), type(my_type), isInternal(is_internal),
          next(nullptr)
    {
      value.var = 0;
    }
  };

  /* Global options */
  struct aprepro_options
  {
    std::string include_path{};
    std::string include_file{};
    bool        end_on_exit{false};
    bool        errors_fatal{false};
    bool        errors_and_warnings_fatal{false};
    bool        warning_msg{true};
    bool        info_msg{false};
    bool        debugging{false};
    bool        dumpvars{false};
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
    int         lineno{0};
    int         loop_count{0};
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

    enum class SYMBOL_TYPE {
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

    /** Return string representation of current version of aprepro.  */
    std::string version() const;

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
    bool string_interactive() { return stringInteractive; }

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
    std::ostream *infoStream{&std::cout};

    void set_error_streams(std::ostream *c_error, std::ostream *c_warning, std::ostream *c_info);

    void dumpsym(const char *type, bool doInternal) const;
    void dumpsym(int type, bool doInternal) const;
    void dumpsym(int type, const char *pre, bool doInternal) const;

  private:
    void                  init_table(const char *comment);
    std::vector<symrec *> sym_table{};
    std::ostringstream    parsingResults{};

    // Input stream used with parse_string_interactive
    std::istringstream stringInput{};

    bool           stringInteractive{false};
    class Scanner *stringScanner{nullptr};

    // For error handling
    std::ostream *errorStream{&std::cerr};
    std::ostream *warningStream{&std::cerr};

    // For substitution history.
    std::vector<history_data> history{};

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

    // Flag to inidicate whether Aprepro is in the middle of collecting lines for a
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

#endif // SEAMS_DRIVER_H
