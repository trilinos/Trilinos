// Kris
// 07.08.03 -- Move into Teuchos package/namespace

// This is the Arbitrary class straight out of NOX with very little modification. Only change is the removal of inclusion of NOX_Common.H; instead, <string> is included directly.

// The Arbitrary class is included in ParameterList for backward compatibility with NOX and any other packages that use the NOX-style parameter list.

// #include "NOX_Common.H"		// class data element (string)
#include <string>
#include <iostream>

namespace Teuchos {

class Arbitrary {
  
public:
  
  Arbitrary() {};
  virtual ~Arbitrary() {};
  virtual Arbitrary* clone() const = 0;
  virtual const std::string& getType() const = 0;
  virtual std::ostream& print(std::ostream& stream, int indent = 0) const {return stream;} ;

};

}
