#ifndef STK_STK_BALANCE_STK_BALANCE_INTERNAL_BALANCECOMMANDLINE_HPP_
#define STK_STK_BALANCE_STK_BALANCE_INTERNAL_BALANCECOMMANDLINE_HPP_
#include <mpi.h>
#include <stk_util/command_line/CommandLineParser.hpp>
#include <string>
#include <stk_balance/internal/balanceDefaults.hpp>
#include <stk_balance/balanceUtils.hpp>

namespace stk { class CommandLineParserParallel; }

namespace stk {
namespace balance {

struct CommandLineOptions
{
    stk::CommandLineOption infile{"infile", "i",
                                  "undecomposed serial input mesh file"};
    stk::CommandLineOption outputDirectory{"outputDirectory", "o",
                                           "output directory for decomposition"};

    stk::CommandLineOption smDefaults{"sm", "",
                                      "Use settings suitable for solving Solid Mechanics problems. "
                                      "This flag implies:\n"
                                      "    --face-search-rel-tol=0.15\n"
                                      "    Face search graph vertex weight multiplier = 10\n"
                                      "    Face search graph edge weight = 3"};
    stk::CommandLineOption sdDefaults{"sd", "",
                                      "Use settings suitable for solving Structural Dynamics problems. "
                                      "This flag implies:\n"
                                      "    --face-search-abs-tol=0.0001\n"
                                      "    Face search graph vertex weight multiplier = 5\n"
                                      "    Face search graph edge weight = 15\n"
                                      "    Handle spider elements"};

    stk::CommandLineOption contactSearch{"contact-search", "",
                                         "Use proximity search for contact [on|off]"};
    stk::CommandLineOption faceSearchAbsTol{"face-search-abs-tol", "",
                                            "Use an absolute tolerance for face contact search. "
                                            "Optionally provide a numeric tolerance value."};
    stk::CommandLineOption faceSearchRelTol{"face-search-rel-tol", "",
                                            "Use a tolerance relative to the face size for face contact search. "
                                            "Optionally provide a numeric tolerance value."};
    stk::CommandLineOption decompMethod{"decomp-method", "",
                                        "Use this geometric decomposition method [rcb|rib|multijagged] "
                                        "or graph-based decomposition method [parmetis].\n"
                                        "Note that geometric methods do not use contact search and "
                                        "ignore all search-related options, as well as ignoring spider elements."};
};

struct ParsedOptions
{
    std::string m_inFile;
    std::string outputDirectory;
    enum AppTypeDefaults appTypeDefaults;
    double faceSearchAbsTol;
    double faceSearchRelTol;
    bool useContactSearch;
    std::string decompMethod;

    enum ParserFlags {
      APP_TYPE              = 1,
      FACE_SEARCH_ABS_TOL   = 2,
      FACE_SEARCH_REL_TOL   = 4,
      CONTACT_SEARCH        = 8,
      DECOMP_METHOD         = 16
    };

    ParsedOptions()
      : ParsedOptions("", "")
    {
    }

    ParsedOptions(const std::string& inFile_, const std::string& outputDirectory_)
      : m_inFile(inFile_),
        outputDirectory(outputDirectory_),
        appTypeDefaults(NO_DEFAULTS),
        faceSearchAbsTol(defaultFaceSearchTolerance),
        faceSearchRelTol(0.15),
        useContactSearch(true),
        decompMethod("parmetis"),
        m_parsedState(0)
    {
    }

    void set_contact_search(bool useSearch)
    {
      useContactSearch = useSearch;
      m_parsedState |= CONTACT_SEARCH;
    }

    void set_face_search_abs_tol(double tolerance)
    {
      faceSearchAbsTol = tolerance;
      m_parsedState |= FACE_SEARCH_ABS_TOL;
    }

    void set_face_search_rel_tol(double tolerance)
    {
      faceSearchRelTol = tolerance;
      m_parsedState |= FACE_SEARCH_REL_TOL;
    }

    void set_app_type_default(AppTypeDefaults appTypeDefault_)
    {
      appTypeDefaults = appTypeDefault_;
      m_parsedState |= APP_TYPE;
    }

    void set_decomp_method(std::string method)
    {
      decompMethod = method;
      m_parsedState |= DECOMP_METHOD;
    }

    bool is_option_provided(ParserFlags parserFlags) const
    {
      return (m_parsedState & parserFlags);
    }

private:
    unsigned m_parsedState;
};

std::string get_quick_example(const std::string &execName,
                              const std::string &infileName,
                              const std::string &outputDirectoryName,
                              MPI_Comm comm);

ParsedOptions parse_balance_command_line(int argc,
                                 const char**argv,
                                 const std::string &execName,
                                 MPI_Comm comm);

void print_running_msg(const std::string &execName, const ParsedOptions &balanceOptions, MPI_Comm comm);

std::string construct_output_file_name(const std::string& outputDirectory, const std::string& inputFile);

enum SearchToleranceType {
  ABSOLUTE,
  RELATIVE,
  UNDEFINED
};

bool get_use_contact_search(std::string useContactSearch);
std::string validate_decomp_method(std::string decompMethod);

}}

#endif /* STK_STK_BALANCE_STK_BALANCE_INTERNAL_BALANCECOMMANDLINE_HPP_ */

