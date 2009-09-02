/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_common.h"
#if defined(HAVE_ML_MLAPI)

#include <iostream>
#include "MLAPI_BaseObject.h"

namespace MLAPI {

int BaseObject::count_ = 0;

std::ostream& operator<< (std::ostream& os, const BaseObject& obj)
{
  return(obj.Print(os));
}

} // namespace MLAPI

#endif
