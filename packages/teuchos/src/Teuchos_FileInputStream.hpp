#ifndef TEUCHOS_FILEINPUTSTREAM_H
#define TEUCHOS_FILEINPUTSTREAM_H

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLInputStream.hpp"

namespace Teuchos
{
  using std::string;

  /** \ingroup XML
   * Input stream for reading an entire document from a file. This
   * is a low-level object and should not be needed at the user level.
   * FileInputSource is the user-level object.
   */

  class FileInputStream : public XMLInputStream
    {
    public:
      /** construct with a filename */
      FileInputStream(const string& filename);
      /** TUVD */
      virtual ~FileInputStream() {;}

      /** read up to maxToRead bytes */
      virtual unsigned int readBytes(unsigned char* const toFill,
                                     const unsigned int maxToRead);

    private:
      FILE* file_;
    };
}
#endif

