// Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
// * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//           
// * Redistributions in binary form must reproduce the above
//   copyright notice, this list of conditions and the following
//   disclaimer in the documentation and/or other materials provided
//   with the distribution.
//                         
// * Neither the name of Sandia Corporation nor the names of its
//   contributors may be used to endorse or promote products derived
//   from this software without specific prior written permission.
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

#ifndef _SystemInterface_h
#define _SystemInterface_h

#include <EML_CodeTypes.h>              // for StringIdVector
#include <GetLongOpt.h>                 // for GetLongOption
#include <iosfwd>                       // for ostream
#include <string>                       // for string

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
  
  GetLongOption options_; //!< Options parsing
  
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
