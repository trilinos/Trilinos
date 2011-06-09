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

#ifndef IOSS_Ioss_Transform_h
#define IOSS_Ioss_Transform_h

#include <string>
#include <vector>
#include <map>

namespace Ioss {
  class Field;
  class VariableType;
}
namespace Ioss {
  class Transform
    {
    public:
      virtual ~Transform();
      virtual const
	Ioss::VariableType *output_storage(const Ioss::VariableType *in) const = 0;
      virtual int output_count(int in) const = 0;

      bool execute(const Ioss::Field &field, void *data);

      virtual void set_property(const std::string &name, int value);
      virtual void set_property(const std::string &name, double value);
      virtual void set_properties(const std::string &name,
				  const std::vector<int> &values);
      virtual void set_properties(const std::string &name,
				  const std::vector<double> &values);
    protected:
      Transform();

      virtual bool internal_execute(const Ioss::Field &field, void *data) = 0;
    };
}

namespace Iotr {

  class Factory;

  typedef std::vector<std::string> NameList;
  typedef std::map<std::string, Factory*, std::less<std::string> > FactoryMap;
  typedef FactoryMap::value_type FactoryValuePair;

  class Factory {
  public:
    virtual ~Factory() {};
    static Ioss::Transform* create(const std::string& type);

    static int describe(NameList *names);

  protected:
    explicit Factory(const std::string& type);
    virtual Ioss::Transform* make(const std::string&) const = 0;
    static void alias(const std::string& base, const std::string& syn);

  private:
    static FactoryMap* registry();
  };
}

#endif // IOSS_Ioss_Transform_h
