/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Governement retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
/** \file aprepro.h Declaration of the SEAMS::Aprepro class. */

/* Might be good to add a callback function which would be called
   when there was output -- In LexerOutput for example.  Default
   could be to just write to std::cout or to resultsOutput stringstream...
*/

#ifndef SEAMS_DRIVER_H
#define SEAMS_DRIVER_H

#include <string>
#include <cstdlib>
#include <vector>
#include <stack>
#include <ostream>
#include <sstream>

/** The SEAMS namespace is used to encapsulate the three parser classes
 * SEAMS::Parser, SEAMS::Scanner and SEAMS::Aprepro */
namespace SEAMS {

  /* Array structure */
  struct array {
    std::vector<double> data;
    int rows;
    int cols;
    
    array(int r, int c) : rows(r), cols(c) {data.resize(r*c);}
    array() : rows(0), cols(0) {}
    ~array() {}
  };

  struct symrec
  {
    std::string name;
    std::string syntax;
    std::string info;
    int   type;
    bool  isInternal;  
    struct value {
      double var;
      double (*fnctptr)();
      double (*fnctptr_d)(double);
      double (*fnctptr_c)(char*);
      double (*fnctptr_dc)(double,char*);
      double (*fnctptr_cd)(char*,double);
      double (*fnctptr_cc)(char*,char*);
      double (*fnctptr_dd)(double,double);
      double (*fnctptr_ddd)(double,double,double);
      double (*fnctptr_ccd)(char*, char*,double);
      double (*fnctptr_dddd)(double,double,double,double);
      double (*fnctptr_ddddc)(double,double,double,double,char*);
      double (*fnctptr_dddddd)(double,double,double,double,double,double);
      double (*fnctptr_a)(const array*);
      const char *svar;
      const char *(*strfnct)();
      const char *(*strfnct_c)(char*);
      const char *(*strfnct_d)(double);
      const char *(*strfnct_a)(const array*);
      const char *(*strfnct_dd)(double,double);
      const char *(*strfnct_ccc)(char*,char*,char*);
      const char *(*strfnct_dcc)(double,char*,char*);
      const char *(*strfnct_dcccc)(double, char*, char*, char*, char*);
      array *avar; /* Array Variable */
      array *(*arrfnct_c)(const char*);
      array *(*arrfnct_cc)(const char*,const char*);
      array *(*arrfnct_cd)(const char*,double);
      array *(*arrfnct_dd)(double,double);
      array *(*arrfnct_d)(double);
      array *(*arrfnct_a)(const array*);

      value() :
	var(0),
	fnctptr(NULL),
	fnctptr_d(NULL),
	fnctptr_c(NULL),
	fnctptr_dc(NULL),
	fnctptr_cd(NULL),
	fnctptr_cc(NULL),
	fnctptr_dd(NULL),
	fnctptr_ddd(NULL),
	fnctptr_ccd(NULL),
	fnctptr_dddd(NULL),
	fnctptr_ddddc(NULL),
	fnctptr_dddddd(NULL),
	fnctptr_a(NULL),
	svar(NULL),
	strfnct(NULL),
	strfnct_c(NULL),
	strfnct_d(NULL),
	strfnct_a(NULL),
	strfnct_dd(NULL),
	strfnct_ccc(NULL),
	strfnct_dcc(NULL),
	strfnct_dcccc(NULL),
	avar(NULL),
	arrfnct_c(NULL),
	arrfnct_cc(NULL),
	arrfnct_cd(NULL),
	arrfnct_dd(NULL),
	arrfnct_d(NULL),
	arrfnct_a(NULL) {}
    } value;
    symrec *next;

    symrec(const char *my_name, int my_type, bool is_internal = false)
      : name(my_name), syntax(my_name), info("UNDEFINED"), type(my_type),
	isInternal(is_internal), next(NULL)
    {
      value.var = 0;
    }

    symrec(const std::string &my_name, int my_type, bool is_internal = false)
      : name(my_name), syntax(my_name), info("UNDEFINED"), type(my_type),
	isInternal(is_internal), next(NULL)
    {
      value.var = 0;
    }
  };

  /* Global options */
  struct aprepro_options
  {
    std::string include_path;
    std::string include_file;
    bool end_on_exit;
    bool warning_msg;
    bool info_msg;
    bool debugging;
    bool interactive;
    bool immutable;
    bool trace_parsing;    // enable debug output in the bison parser
    bool one_based_index;
    bool keep_history;  // Flag to keep a history of Aprepro substitutions
    aprepro_options() :
      end_on_exit(false),
      warning_msg(true),
      info_msg(false),
      debugging(false),
      interactive(false),
      immutable(false),
      trace_parsing(false),
      one_based_index(false),
      keep_history(false)
    {}
  };

  /* Structure for holding file names and line counters */
  struct file_rec
  {
    std::string name;
    int	  lineno;
    int	  loop_count;
    bool  tmp_file;

    file_rec(const char *my_name, int line_num, bool is_temp, int loop_cnt)
      : name(my_name), lineno(line_num), loop_count(loop_cnt), tmp_file(is_temp) {}
    file_rec()
      : name("STDIN"), lineno(0), loop_count(0), tmp_file(false) {}
  };

  /* Structure for holding aprepro substitution info */
  struct history_data
  {
    std::string original;
    std::string substitution;
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

    enum SYMBOL_TYPE {
      VARIABLE = 1,
      STRING_VARIABLE = 2,
      UNDEFINED_VARIABLE = 5,
      FUNCTION = 3,
      STRING_FUNCTION = 4,
      ARRAY_FUNCTION = 6,
      ARRAY_VARIABLE = 7,
      IMMUTABLE_VARIABLE = 11,
      IMMUTABLE_STRING_VARIABLE = 12
    };
    
    bool state_is_immutable() const {return stateImmutable;}
    
    /** Return an  std::ostringstream reference to get the results of
	the parse_* call (* = stream, file, or string).
    */
    const std::ostringstream &parsing_results() const {return parsingResults;}
    void clear_results();
    
    /** Return string representation of current version of aprepro.  */
    std::string version() const;
    
    /** Invoke the scanner and parser for a stream.
     * @param in	input stream
     * @param sname	stream name for error messages
     * @return		true if successfully parsed
     */
    bool parse_stream(std::istream& in,
		      const std::string& sname = "stream input");

    /** Invoke the scanner and parser on an input string.
     * @param input	input string
     * @param sname	stream name for error messages
     * @return		true if successfully parsed
     */
    bool parse_string(const std::string& input,
		      const std::string& sname = "string stream");

    /** Invoke the scanner and parser on a vector of strings.
     * @param input	vector of input strings
     * @param sname	stream name for error messages
     * @return		true if successfully parsed
     */
    bool parse_strings(const std::vector<std::string> &input, const std::string& sname);

    /** Invoke the scanner and parser on an input string in an
     * interactive manner.
     * @param input input stringInput
     * @return true if successfully parsed
     */
    bool parse_string_interactive(const std::string &input);

    /** Get the string interactive flag, which indicates if we are in
     * the middle of parsing a string in an interactive manner.
     */
    bool string_interactive() {return stringInteractive;}

    /** Invoke the scanner and parser on a file. Use parse_stream with a
     * std::ifstream if detection of file reading errors is required.
     * @param filename	input file name
     * @return		true if successfully parsed
     */
    bool parse_file(const std::string& filename);

    void statistics(std::ostream *out = NULL) const;
    void copyright(std::ostream *out = NULL) const;
    
    aprepro_options ap_options;
    std::stack<file_rec> ap_file_list;

    std::stack<std::ostream*> outputStream;
    
    SEAMS::symrec *getsym(const char *) const;
    SEAMS::symrec *putsym(const std::string &sym_name, SYMBOL_TYPE sym_type, bool is_internal);

    void add_variable(const std::string &sym_name, const std::string &sym_value, bool is_immutable=false);
    void add_variable(const std::string &sym_name, double sym_value, bool is_immutable=false);
    std::vector<std::string> get_variable_names(bool doInternal = false);
    void remove_variable(const std::string &sym_name);

    bool set_option(const std::string &option);
    
    std::fstream *open_file(const std::string &file, const char *mode);
    std::fstream *check_open_file(const std::string &file, const char *mode);
    
    /** Pointer to the current lexer instance, this is used to connect the
     * parser to the scanner. It is used in the yylex macro. */
    class Scanner* lexer;

    /** Error handling. */
    void error(const std::string& msg,
               bool line_info=true, bool prefix=true) const;
    void warning(const std::string& msg,
                 bool line_info=true, bool prefix=true) const;
    void info(const std::string& msg,
              bool line_info=false, bool prefix=true) const;

    // The info stream. To only print out info messages if the -M option was
    // specified, use info(...) instead.
    std::ostream *infoStream;

    void set_error_streams(std::ostream* error, std::ostream* warning,
                           std::ostream* info);

    void dumpsym (int type, bool doInternal) const;
  private:
    
    void init_table(const char *comment);
    std::vector<symrec*> sym_table;
    std::ostringstream parsingResults;

    // Input stream used with parse_string_interactive
    std::istringstream stringInput;

    bool stringInteractive;
    class Scanner* stringScanner;

    // For error handling
    std::ostream *errorStream;
    std::ostream *warningStream;

    // For substitution history.
    std::vector<history_data> history;

  public:
    bool stateImmutable;

    // Flag to do Aprepro substitutions within loops. Default value is true. If set to
    // false, content within the loop will be treated as verbatim text.
    bool doLoopSubstitution;

    // Flag to do Aprepro substitutions when including a file. Default value is true.
    // If set to false, content within the file will be treated as verbatim text that
    // needs to be sent through Aprepro again later.
    bool doIncludeSubstitution;

    // Flag to inidicate whether Aprepro is in the middle of collecting lines for a
    // loop.
    bool isCollectingLoop;

    // Record the substitution of the current Aprepro statement. This function will also
    // reset the historyString and add an entry to the substitution map.
    void add_history(const std::string& original, const std::string& substitution);

    const std::vector<history_data> &get_history();
    void clear_history();
  };

} // namespace SEAMS

#endif // SEAMS_DRIVER_H
