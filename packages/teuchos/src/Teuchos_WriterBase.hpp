#ifndef TEUCHOS_WRITERBASE_H
#define TEUCHOS_WRITERBASE_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RefCountPtrDecl.hpp"
#include <string>
#include <stdexcept>

namespace Teuchos
{
  using std::string;

  /** 
   * Teuchos::WriterBase is the base class for user-defined writers
   */
  class WriterBase
    {
    public:
      /** Empty ctor */
      WriterBase(){;}

      /** virtual dtor */
      virtual ~WriterBase(){;}

      /** print a string */
      virtual void print(const string& msg) = 0 ;

      /** print a string followed by a newline */
      virtual void println(const string& msg) = 0 ;

      /** Return a RefCountPtr containing self. 
          This is easy to implement: {return rcp(this);} */
      virtual RefCountPtr<WriterBase> getRcp() = 0 ;

    private:
    };
}
#endif
