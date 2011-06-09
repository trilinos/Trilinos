#ifndef IOSS_Ioss_NamedSuffixVariableType_h
#define IOSS_Ioss_NamedSuffixVariableType_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <Ioss_VariableType.h>

namespace Ioss {
  class NamedSuffixVariableType : public VariableType {
  public:
    
    //  'which' is 1-based
    std::string label(int which, const char suffix_sep='_') const
      {
	return suffixList[which-1];
      }
    
    NamedSuffixVariableType(const std::string& my_name, int number_components, bool delete_me)
      : Ioss::VariableType(my_name, number_components, delete_me)
      {
	suffixList.resize(number_components);
	suffixList.assign(number_components, "UNSET");
      }

      //! Define the suffix list for this field.
      //  'which' is 1-based to conform to the 'label' function usage.
      // If user doesn't add suffices, then 'label' will return "UNSET"
      void add_suffix(size_t which, const std::string &suffix)
      {
	suffixList[which-1] = suffix;
      }

  private:
      NamedSuffixVariableType(const NamedSuffixVariableType&);
      std::vector<std::string> suffixList;
  };
}

#endif
