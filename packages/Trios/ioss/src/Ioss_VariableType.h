/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_VariableType_h
#define SIERRA_Ioss_VariableType_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <cstring>
#include <vector>
#include <map>

namespace Ioss {
  typedef std::vector<std::string> NameList;

  class VariableType;
  
  typedef std::map<std::string, VariableType*, std::less<std::string> > VariableTypeMap;
  typedef VariableTypeMap::value_type VTM_ValuePair;

  class Registry {
  public:
    void insert(const Ioss::VTM_ValuePair &value, bool delete_me);
    VariableTypeMap::iterator begin() {return m_registry.begin();}
    VariableTypeMap::iterator end()   {return m_registry.end();}
    VariableTypeMap::iterator find(const std::string &type) {return m_registry.find(type);}

    ~Registry();
  private:
    Ioss::VariableTypeMap m_registry;
    std::vector<Ioss::VariableType*> m_deleteThese;
  };

#define MAX_SUFFIX 8
  struct Suffix {
    Suffix(const char new_data[MAX_SUFFIX]) {std::strncpy(data, new_data,         MAX_SUFFIX); data[MAX_SUFFIX]='\0';}
    Suffix(const std::string &new_data)     {std::strncpy(data, new_data.c_str(), MAX_SUFFIX); data[MAX_SUFFIX]='\0';}
    bool operator==(const std::string &str) const
    {return std::strncmp(data, str.c_str(), MAX_SUFFIX) == 0;}
    bool operator!=(const std::string &str) const
    {return std::strncmp(data, str.c_str(), MAX_SUFFIX) != 0;}
    char data[MAX_SUFFIX+1];
  };

  class VariableType {
  public:

    static void alias(const std::string& base, const std::string& syn);
    static int describe(NameList *names);

    virtual ~VariableType();
    int component_count() const;

    // Override this function if the derived class has no suffices
    // For example, a 'vector_2d' has suffices "x" and "y"
    // A 'quad4' has no suffices...
    virtual int suffix_count() const;
    std::string name() const;

    static std::string numeric_label(int which, int ncomp,
					      const std::string &name);
    virtual std::string label(int which, const char suffix_sep='_') const = 0;
    virtual std::string label_name(const std::string& base, int which,
				      const char suffix_sep='_') const;
    virtual bool match(const std::vector<Suffix> &suffices) const;

    static const VariableType* factory(const std::string& name, int copies = 1);
    static const VariableType* factory(const std::vector<Suffix> &suffices);

  protected:
    VariableType(const std::string& type, int comp_count, bool delete_me = false);
    static Registry& registry();

  private:

    const std::string name_;
    int componentCount;

    VariableType(const VariableType&); // Do not implement...
    VariableType& operator=(const VariableType&); // Do not implement...

    static bool build_variable_type(const std::string& type);
  };
}
inline std::string Ioss::VariableType::name() const
{return name_;}

inline int Ioss::VariableType::component_count() const
{return componentCount;}

inline int Ioss::VariableType::suffix_count() const
{return componentCount;}
#endif
