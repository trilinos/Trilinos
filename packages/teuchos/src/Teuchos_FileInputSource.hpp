#ifndef Teuchos_FILEINPUTSOURCE_H
#define Teuchos_FILEINPUTSOURCE_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputSource.hpp"


namespace Teuchos
{

  /** \ingroup XML
   * Input source that reads XML from a file.
   */

  class FileInputSource : public XMLInputSource
    {
    public:
      /** ctor */
      FileInputSource(const string& filename);

      /** virtual dtor */
      virtual ~FileInputSource(){;}

      /** create a FileInputStream */
      virtual RefCountPtr<XMLInputStream> stream() const;

    private:
      string filename_;
    };

}
#endif

