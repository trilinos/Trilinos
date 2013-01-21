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

  struct symrec
  {
    std::string name;
    std::string syntax;
    std::string info;
    int   type;
    bool  isInternal;  
    union {
      double var;
      double (*fnctptr)();
      double (*fnctptr_d)(double);
      double (*fnctptr_c)(char*);
      double (*fnctptr_dc)(double,char*);
      double (*fnctptr_cd)(char*,double);
      double (*fnctptr_cc)(char*,char*);
      double (*fnctptr_dd)(double,double);
      double (*fnctptr_ddd)(double,double,double);
      double (*fnctptr_dddd)(double,double,double,double);
      double (*fnctptr_dddddd)(double,double,double,double,double,double);
      const char *svar;
      const char *(*strfnct)();
      const char *(*strfnct_c)(char*);
      const char *(*strfnct_d)(double);
      const char *(*strfnct_ccc)(char*,char*,char*);
      const char *(*strfnct_dcc)(double,char*,char*);
      const char *(*strfnct_dcccc)(double, char*, char*, char*, char*);
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
    char comment;
    bool end_on_exit;
    bool warning_msg;
    bool info_msg;
    bool debugging;
    bool interactive;
    bool immutable;
    bool trace_parsing;    // enable debug output in the bison parser

    aprepro_options() :
      end_on_exit(false),
      warning_msg(true),
      info_msg(false),
      debugging(false),
      interactive(false),
      immutable(false),
      trace_parsing(false)
    {}
  };

  /* Structure for holding file names and line counters */
  struct file_rec
  {
    std::string name;
    int	  lineno;
    bool  tmp_file;
    int	  loop_count;
    file_rec(const char *my_name, int line_num, bool is_temp, int loop_cnt)
      : name(my_name), lineno(line_num), tmp_file(is_temp), loop_count(loop_cnt) {}
    file_rec()
      : name("UNKNOWN_FILE_NAME"), lineno(0), tmp_file(false), loop_count(0) {}
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
      IMMUTABLE_VARIABLE = 11,
      IMMUTABLE_STRING_VARIABLE = 12
    };
    
    bool state_is_immutable() const {return stateImmutable;}
    
    /** Return an  std::ostringstream reference to get the results of
	the parse_* call (* = stream, file, or string).
    */
    const std::ostringstream &parsing_results() const {return parsingResults;}
    void clear_results() {parsingResults.str("");}
    
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

    void add_variable(const std::string &sym_name, const std::string &sym_value);
    void add_variable(const std::string &sym_name, double sym_value);
    
    std::fstream *open_file(const std::string &file, const char *mode);
    std::fstream *check_open_file(const std::string &file, const char *mode);
    
    /** Pointer to the current lexer instance, this is used to connect the
     * parser to the scanner. It is used in the yylex macro. */
    class Scanner* lexer;

    /** Error handling. */
    void error(const std::string& m) const;

    void dumpsym (int type, bool doInternal) const;
  private:
    
    void init_table(char comment);
    std::vector<symrec*> sym_table;
    std::ostringstream parsingResults;
  public:
    bool stateImmutable;
  };

} // namespace SEAMS

#endif // SEAMS_DRIVER_H
