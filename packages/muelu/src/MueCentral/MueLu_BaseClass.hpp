#ifndef MUELU_BASECLASS_HPP
#define MUELU_BASECLASS_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_VerboseObject.hpp"
#include "MueLu_Describable.hpp"

namespace MueLu {

  //! Base class for MueLu classes
  class BaseClass
    : public VerboseObject, public Describable
  {
    
  public:
    
    //! Destructor.
    virtual ~BaseClass() {}

    //@}

  }; // class BaseClass

} // namespace MueLu

//! Helper macro for implementing Describable::describe() for BaseClass objects.
//  This macro defines ostream out0 that print only on root node. It print description() and indent the ostream.
//  Note: Runtime1 displays basic parameter information when Parameters0 is not enabled.
#define MUELU_DESCRIBE                                                  \
  using std::endl;                                                      \
  Teuchos::FancyOStream& out0 = (VerboseObject::GetProcRankVerbose() == 0) ? out : VerboseObject::GetBlackHole(); \
                                                                        \
  if ((verbLevel & Runtime1) && (!verbLevel & Parameters0))             \
    out << description() << endl;                                       \
  else if (verbLevel & Runtime0)                                        \
    out << BaseClass::description() << endl;                            \
                                                                        \
  Teuchos::OSTab tab1(out);                                             \
  //

#define MUELU_BASECLASS_SHORT
#endif // ifndef MUELU_BASECLASS_HPP
