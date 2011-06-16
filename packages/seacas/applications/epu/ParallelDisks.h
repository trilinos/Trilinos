/*
 * Copyright(C) 2010 Sandia Corporation.  Under the terms of Contract
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

#ifndef ParallelDisksH
#define ParallelDisksH

#include <CodeTypes.h>
#include <string>

/*****************************************************************************/
namespace Excn {
  class ParallelDisks 
    {
      //: This class declares the unique, global object that represents the
      //: disk farm on a parallel machine.  Its function is to map
      //: processors onto disk files.

    public:

      ParallelDisks();
      ~ParallelDisks();

      void Number_of_Raids(int);
      void Raid_Offset(int);

      int Number_of_Raids() const;
      int Raid_Offset() const;

      static void Create_IO_Filename(std::string &, int num, int total_num);
      int  IO_Node_Number(int) const;
      int  Logical_Node_Number(int) const;

      void rename_single_file_for_mp(const std::string& dir, std::string& name) const;
      void rename_file_for_mp(const std::string& rootdir, const std::string& subdir, 
			      std::string& name, int num, int total_num) const;
      void rename_file_for_no_overwrite(std::string& name);

    private:

      std::string* disk_names;

      int    number_of_raids;
      int    raid_offset;


      void create_disk_names();

      // not defined (the parallel disks object is unique!)
      ParallelDisks(const ParallelDisks&);
      ParallelDisks& operator=(const ParallelDisks&);

    };
}
#endif
