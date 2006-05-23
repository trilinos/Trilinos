#ifndef PHX_CORE_UTILS_H
#define PHX_CORE_UTILS_H

#include "Epetra_Map.h"

#include "Teuchos_Hashtable.hpp"

namespace phx {
namespace core {

class Utils 
{
  public:
    Utils() {}
    ~Utils() {}

    static int getNumDimensions()
    {
      return(numDimensions_);
    }

    static void setNumDimensions(const int numDimensions)
    {
      numDimensions_ = numDimensions;
    }

    static int numDimensions_;

  private:

}; // class Utils

} // namespace core
} // namespace phx

#endif
