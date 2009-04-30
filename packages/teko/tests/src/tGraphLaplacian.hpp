#ifndef __tGraphLaplacian_hpp__
#define __tGraphLaplacian_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"

#include <string>

#include "Test_Utils.hpp"

namespace PB {
namespace Test {

class tGraphLaplacian : public UnitTest {
public:
   virtual ~tGraphLaplacian() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_single_array(int verbosity,std::ostream & os);
   bool test_multi_array(int verbosity,std::ostream & os);

protected:
   bool compareMatrix(const Epetra_CrsMatrix & gl,const std::string & name,int verbosity,std::ostream & os) const;
   double tolerance_;
};

} // end namespace Tests
} // end namespace PB

#endif
