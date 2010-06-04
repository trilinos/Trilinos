/*--------------------------------------------------------------------*/
/*    Copyright 2006 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_CompositeVariableType_h
#define SIERRA_Ioss_CompositeVariableType_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <Ioss_VariableType.h>

namespace Ioss {
  class CompositeVariableType : public VariableType {
  public:
    static std::string composite_name(const std::string &base, int copies);
    static VariableType* composite_variable_type(const VariableType *inst, int copies);

    std::string label(int which, const char suffix_sep='_') const;
    CompositeVariableType(const std::string& name, int number_components, bool delete_me);
    CompositeVariableType(const VariableType *base_type, int copies, bool delete_me);

  private:
    CompositeVariableType(const CompositeVariableType&);
    const VariableType *baseType;
    int copies_;
  };
}
#endif
