#ifndef TEUCHOS_OUT_H
#define TEUCHOS_OUT_H

#include "Teuchos_ConfigDefs.hpp"
#include <string>
#include <stdexcept>
#include "Teuchos_WriterBase.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos
{

  /**
   * Class Out allows customizable control of output in a parallel environment.
   * Out writes text using a user-selected writer. The default writer is
   * class DefaultWriter which prepends a processor identifier to every line; this 
   * is helpful in debugging.
   * 
   * All methods of Out are static. Methods such as println() and printf() work
   * on every processor. The "root" methods such as rootPrinln() and rootPrintf()
   * write if called on the root processor, but are ignored on all non-root processors.
   */
  class Out
    {
    public:
      /** Print a string using the selected writer */
      static void print(const std::string& msg);

      /** Print a string plus a newline using the selected writer */
      static void println(const std::string& msg);

      /** Print a string on the root processor, do a no-op on all other processors. */
      static void rootPrintln(const std::string& msg);

      /** Print using standard C formatting */
      static void printf(const char* format, ...);

      /** Print on the root processor using standard C formatting */
      static void rootPrintf(const char* format, ...);

      /** Set the writer to used by all methods of Out */
      static void setWriter(const RefCountPtr<WriterBase>& writer);


    private:
      static RefCountPtr<WriterBase> writer_;

      /** */
      static void vprintf(const char* format, va_list args);

			/** */
			static int hack_vsnprintf(char* str, size_t size, const char* format, va_list args);
    };
}

#endif
