#ifndef _SystemInterface_h
#define _SystemInterface_h

#include <CodeTypes.h>
#include <GetLongOpt.h>

#include <iosfwd>
class SystemInterface
{
 public:
  SystemInterface();
  ~SystemInterface();

  bool parse_options(int argc, char **argv);
  
  char field_suffix() const {return fieldSuffix_;}

  StringIdVector global_var_names() const {return globalVarNames_;}
  StringIdVector node_var_names() const {return nodeVarNames_;}
  StringIdVector elem_var_names() const {return elemVarNames_;}
  StringIdVector nset_var_names() const {return nsetVarNames_;}
  StringIdVector sset_var_names() const {return ssetVarNames_;}
  StringIdVector vars_to_list()   const {return varsToList_;}
  bool list_vars() const {return listVars_;}
  
  std::string input_file() const  {return inputFile_;}
  std::string output_file() const {return outputFile_;}
      
  double minimum_time() const {return minimumTime_;}
  double maximum_time() const {return maximumTime_;}
  
  //! Dumps representation of data in this class to cerr
  void dump(std::ostream &str) const;
  
  static void show_version();
  
 private:
  void enroll_options();
  void parse_exclude(const char *list);

  double minimumTime_;
  double maximumTime_;
  
  GetLongOpt options_; //!< Options parsing
  
  std::string inputFile_;
  std::string outputFile_;

  StringIdVector globalVarNames_;
  StringIdVector nodeVarNames_;
  StringIdVector elemVarNames_;
  StringIdVector nsetVarNames_;
  StringIdVector ssetVarNames_;
  StringIdVector varsToList_;

  bool listVars_;
  char fieldSuffix_;
};
#endif
