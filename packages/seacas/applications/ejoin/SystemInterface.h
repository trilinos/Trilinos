// Copyright(C) 2010 Sandia Corporation.
// 
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
#ifndef Sierra_SystemInterface_h
#define Sierra_SystemInterface_h

#include <CodeTypes.h>
#include <GetLongOpt.h>
#include "Vector3.h"

#include <iosfwd>
class SystemInterface
{
 public:
  SystemInterface();
  ~SystemInterface();

  bool parse_options(int argc, char **argv);
  
  int debug() const {return debugLevel_;}
  int screen_width() const {return screenWidth_;}
      
  double tolerance() const {return tolerance_;}
  bool match_node_ids() const {return matchNodeIds_;}
  bool match_node_xyz() const {return matchNodeXYZ_;}
  bool match_elem_ids() const {return matchElemIds_;}
  bool match_elem_xyz() const {return matchElemXYZ_;}
  bool omit_nodesets() const {return omitNodesets_;}
  bool omit_sidesets() const {return omitSidesets_;}
  bool convert_nodes_to_nodesets(int part_number) const;
  bool disable_field_recognition() const {return disableFieldRecognition_;}
  
  int step_min() const {return stepMin_;}
  int step_max() const {return stepMax_;}
  int step_interval() const {return stepInterval_;}

  Vector3 offset() const {return offset_;}
  std::vector<int> node_nodeset_conversions() const {return nodesetConvertParts_;}
  StringIdVector global_var_names() const {return globalVarNames_;}
  StringIdVector node_var_names() const {return nodeVarNames_;}
  StringIdVector elem_var_names() const {return elemVarNames_;}
  StringIdVector nset_var_names() const {return nsetVarNames_;}
  StringIdVector sset_var_names() const {return ssetVarNames_;}

  const Omissions &block_omissions() const {return blockOmissions_;}
  const Omissions &nset_omissions() const {return nsetOmissions_;}
  const Omissions &sset_omissions() const {return ssetOmissions_;}
  
  const std::string &element_status_variable() const {return elementStatusVariable_;}
  const std::string &nodal_status_variable()   const {return nodalStatusVariable_;}

  //! Dumps representation of data in this class to cerr
  void dump(std::ostream &str) const;
  
  static void show_version();
  
  // Make this private eventually...
  StringVector   inputFiles_;
  std::string outputName_;

 private:
  void enroll_options();
  bool decompose_filename(const std::string &cs);
  void parse_exclude(const char *list);

  /*! The defined formats for the count attribute are:<br>
    - <missing> -- default -- 1 <= count <= oo  (all steps)
    - "X"                  -- X <= count <= X  (just step X) (if X == -1, do last step only)
    - "X:Y"                -- X to Y by 1
    - "X:"                 -- X to oo by 1
    - ":Y"                 -- 1 to Y by 1
    - "::Z"                -- 1 to oo by Z

    The count and step must always be >= 0
  */
  void parse_step_option(const char *token);

  GetLongOpt options_; //!< Options parsing
  

  int debugLevel_;
  int screenWidth_;
  int stepMin_;
  int stepMax_;
  int stepInterval_;
  bool omitNodesets_;
  bool omitSidesets_;
  bool matchNodeIds_;
  bool matchNodeXYZ_;
  bool matchElemIds_;
  bool matchElemXYZ_;
  bool disableFieldRecognition_;
  
  std::string elementStatusVariable_;
  std::string nodalStatusVariable_;

  Vector3 offset_;
  double tolerance_;
      
  Omissions blockOmissions_;
  Omissions nsetOmissions_;
  Omissions ssetOmissions_;

  std::vector<int> nodesetConvertParts_;

  StringIdVector globalVarNames_;
  StringIdVector nodeVarNames_;
  StringIdVector elemVarNames_;
  StringIdVector nsetVarNames_;
  StringIdVector ssetVarNames_;
};
#endif
