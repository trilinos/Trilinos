/*
 * Copyright(C) 2012 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
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
#ifndef info_SystemInterface_h
#define info_SystemInterface_h

#include "Ioss_GetLongOpt.h"

#include <string>
#include <iosfwd>
#include <vector>

namespace Info {
  class Interface
    {
    public:
      Interface();
      ~Interface();

      bool parse_options(int argc, char **argv);
  
      int summary() const {return summary_;}
      bool check_node_status() const {return checkNodeStatus_;}
      bool compute_volume()  const {return computeVolume_;}
      bool compute_bbox()  const {return computeBBox_;}
      bool adjacencies() const {return adjacencies_;}
      bool ints_64_bit() const {return ints64Bit_;}
      
      double maximum_time() const {return maximumTime_;}
      double minimum_time() const {return minimumTime_;}
      int surface_split_scheme() const {return surfaceSplitScheme_;}
      char field_suffix_separator() const {return fieldSuffixSeparator_;}

      std::string cwd() const {return cwd_;}
      std::string filename() const {return filename_;}
      std::string type() const {return filetype_;}
      
  
      //! Dumps representation of data in this class to cerr
      void dump(std::ostream &str) const;
  
      static void show_version();
  
    private:
      void enroll_options();

      Ioss::GetLongOption options_;
      
      bool checkNodeStatus_;
      bool computeVolume_;
      bool adjacencies_;
      bool ints64Bit_;
      bool computeBBox_;
 
      char fieldSuffixSeparator_;
      
      int summary_;
      int surfaceSplitScheme_;
      
      double minimumTime_;
      double maximumTime_;
      std::string cwd_;
      std::string filetype_;
      std::string filename_;
    };
}
#endif
