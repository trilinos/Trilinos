#ifndef Teuchos_STRINGINPUTSOURCE_H
#define Teuchos_STRINGINPUTSOURCE_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputSource.hpp"


namespace Teuchos
{

  using std::string;

  /** \ingroup XML
   * StringInputSource reads XML from a String
   */

  class StringInputSource : public XMLInputSource
    {
    public:
      /** ctor */
      StringInputSource(const string& text);
      /** virtual dtor */
      virtual ~StringInputSource(){;}

      /** create a StringInputStream */
      virtual RefCountPtr<XMLInputStream> stream() const;

    private:
      string text_;
    };


}
#endif

