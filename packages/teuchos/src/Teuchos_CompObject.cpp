// Kris
// 07.08.03 -- Move into Teuchos package/namespace

// Constructor

#include "Teuchos_CompObject.hpp"

namespace Teuchos
{

CompObject::CompObject() : flopCounter_(0)
{
}

// Copy Constructor

CompObject::CompObject(const CompObject& source) : flopCounter_(source.flopCounter_)
{
}

// Destructor

CompObject::~CompObject()
{
  flopCounter_ = 0;
}

} // namespace Teuchos
