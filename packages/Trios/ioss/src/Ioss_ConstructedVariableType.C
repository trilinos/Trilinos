/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_ConstructedVariableType.h>
#include <assert.h>
#include <cstdio>
#include <Ioss_VariableType.h>
#include <Ioss_Utils.h>
#include <string>

Ioss::ConstructedVariableType::ConstructedVariableType(const std::string& my_name,
						       int number_components, bool delete_me)
  : Ioss::VariableType(my_name, number_components, delete_me)
{}

Ioss::ConstructedVariableType::ConstructedVariableType(int number_components, bool delete_me)
  : Ioss::VariableType(std::string("Real[")+Utils::to_string(number_components)+std::string("]"),
		       number_components, delete_me)
{}

std::string Ioss::ConstructedVariableType::label(int which, const char) const
{
  assert(which > 0 && which <= component_count());
  if (component_count() == 1)
    return "";
  return VariableType::numeric_label(which, component_count(), name());
}
