#ifndef TEUCHOS_DEFAULTWRITER_H
#define TEUCHOS_DEFAULTWRITER_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_WriterBase.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include <iostream>

namespace Teuchos
{
  using std::string;
  using std::ostream;

  /** \ingroup IO
   *
   */
  class DefaultWriter : public WriterBase
    {
    public:

      /** */
      DefaultWriter();

      /** */
      DefaultWriter(ostream& os);

      /** */
      DefaultWriter(const RefCountPtr<ostream>& os_ptr);

      /** */
      virtual ~DefaultWriter(){;}

      /** */
      virtual void print(const string& msg);

      /** */
      virtual void println(const string& msg);

      /** Return a RefCountPtr containing self. */
      virtual RefCountPtr<WriterBase> getRcp() {return rcp(this, true);}
    private:
      RefCountPtr<ostream>   os_ptr_;
      std::ostream           &os_;
      static string& header();
    };
}
#endif
