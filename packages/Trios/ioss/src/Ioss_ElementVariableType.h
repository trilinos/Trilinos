/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_ElementVariableType_h
#define SIERRA_Ioss_ElementVariableType_h

#include <Ioss_CodeTypes.h>
#include <Ioss_VariableType.h>
#include <string>

namespace Ioss {
  class ElementVariableType : public Ioss::VariableType
    {
    public:
      std::string label(int, const char) const {return "";}
      std::string label_name(const std::string& base, int,const char) const {return base;}
      int suffix_count() const {return 0;}
    protected:
      ElementVariableType(const std::string& type, int comp_count);
    };

  inline ElementVariableType::ElementVariableType(const std::string& type,
						  int comp_count)
    : Ioss::VariableType(type, comp_count, false) {}
}
#endif
