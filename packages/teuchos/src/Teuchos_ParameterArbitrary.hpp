// Kris
// 07.08.03 -- Move into Teuchos package/namespace

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
