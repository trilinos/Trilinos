/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_ConstructedVariableType_h
#define IOSS_Ioss_ConstructedVariableType_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <Ioss_VariableType.h>

namespace Ioss {
  class ConstructedVariableType : public VariableType {
  public:
    std::string label(int which, const char suffix_sep='_') const;
    ConstructedVariableType(const std::string& name, int number_components, bool delete_me);
    explicit ConstructedVariableType(int number_components, bool delete_me);

  private:
    ConstructedVariableType(const ConstructedVariableType&);
  };
}
#endif
