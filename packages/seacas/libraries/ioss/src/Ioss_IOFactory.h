// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
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

#ifndef IOSS_Ioss_IOFactory_h
#define IOSS_Ioss_IOFactory_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <Ioss_DBUsage.h>

#include <map>
#include <vector>

namespace Ioss {

class IOFactory;

typedef std::vector<std::string> NameList;
typedef std::map<std::string, IOFactory*, std::less<std::string> > IOFactoryMap;
typedef IOFactoryMap::value_type IOFactoryValuePair;

class DatabaseIO;
class IOFactory {
public:
  virtual ~IOFactory() {};
  static DatabaseIO* create(const std::string& type,
                            const std::string& filename,
                            DatabaseUsage db_usage,
                            MPI_Comm communicator = MPI_COMM_WORLD);

  static int describe(NameList *names);
  static void clean();
    
protected:
  explicit IOFactory(const std::string& type);

  virtual DatabaseIO* make_IO(const std::string& filename,
                              DatabaseUsage db_usage,
                              MPI_Comm communicator) const = 0;

  static void alias(const std::string& base, const std::string& syn);

private:
  //    static IOFactoryMap* registry_;
  static IOFactoryMap* registry();
};

}
#endif
