#ifndef __Krylov_ProductIsZero_hpp
#define __Krylov_ProductIsZero_hpp

#include <exception>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Krylov {

  class ProductIsZero : public std::exception 
  {
  public:
    ProductIsZero (const int index, const bool modified);
    const char* what() const;

    int where_;
    bool modified_;

  private:
    std::string errMsg_;
  };

  ProductIsZero::ProductIsZero (const int index, const bool modified)
    : where_ (index), modified_ (modified)
  {}

  const char* ProductIsZero::what() const 
  { 
    if (errMsg_.size() == 0)
      {
	// This method really shouldn't throw, but ostringstream
	// operations might throw due to out-of-memory.  However, if we
	// run out of memory with this little work, the program is
	// completely trashed, so I'm not going to worry about it.
	std::ostringstream os;
	const std::string modifiedString = modified_ ? "(modified) " : "";
	os << "Product to maximize in " << modifiedString << "Leja ordering is ze"
	  "ro at index j = " << where_ << ".  This may be due either to underfl"
	  "ow, or repeated shifts not counted in the multiplicies vector.";
	errMsg_ = os.str();
      }
    return errMsg_.c_str();
  }

} // namespace Krylov
#endif // __Krylov_ProductIsZero_hpp
