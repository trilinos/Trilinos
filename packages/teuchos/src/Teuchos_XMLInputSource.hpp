#ifndef TEUCHOS_XMLINPUTSOURCE_H
#define TEUCHOS_XMLINPUTSOURCE_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLInputStream.hpp"

namespace Teuchos
{
  /** 
   * XMLInputSource represents a source of XML input that can be parsed
   * to produce an XMLObject. The source might be a file, a socket, a
   * string. The XMLObject is created with a call to the getObject() method.
   *
   * The source gets its data from a XMLInputStream object that is
   * created (internally) to work with this source.
   *
   * getObject() is implemented with EXPAT if HAVE_EXPAT is 1, or with
   * Xerces if HAVE_XERCES is 1.
   */
  class XMLInputSource
    {
    public:
      XMLInputSource(){;}

      /** virtual dtor */
      virtual ~XMLInputSource(){;}

      /**  Virtual input source interface */
      virtual RefCountPtr<XMLInputStream> stream() const = 0 ;

      /** get an object by invoking the TreeBuildingXMLHandler on the
       * input data */
      XMLObject getObject() const ;

    };

}
#endif

