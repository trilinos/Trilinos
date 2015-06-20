
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "ROL_Distribution.hpp"

int main(int argc, char* argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  Teuchos::RCP<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing
  if (iprint > 0)
    outStream = Teuchos::rcp(&std::cout, false);
  else
    outStream = Teuchos::rcp(&bhs, false);

  int errorFlag  = 0;

  try {
    Teuchos::RCP<ROL::Distribution<double> > dist;
    ROL::EDistribution ed;
    std::vector<double> data;
   
    // Dirac 
    std::cout << "\nDIRAC DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_DIRAC;
    data.resize(1,0.0);
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data));
    dist->test(*outStream);

    // Gaussian 
    std::cout << "\nGAUSSIAN DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_GAUSSIAN;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Gaussian 
    std::cout << "\nTRUNCATED GAUSSIAN DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_TRUNCATEDGAUSSIAN;
    data.resize(4,0.0);
    data[0] = -1.0; data[1] = 1.0; data[2] = 0.0; data[3] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Uniform
    std::cout << "\nUNIFORM DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_UNIFORM;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Logistic
    std::cout << "\nLOGISTIC DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_LOGISTIC;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Triangle
    std::cout << "\nTRIANGLE DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_TRIANGLE;
    data.resize(3,0.0);
    data[0] = 0.0; data[1] = 0.5; data[2] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Parabolic
    std::cout << "\nPARABOLIC DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_PARABOLIC;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Raised Cosine
    std::cout << "\nRAISED COSINE DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_RAISEDCOSINE;
    data.resize(2,0.0);
    data[0] = 1.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Laplace
    std::cout << "\nLAPLACE DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_LAPLACE;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Cauchy
    std::cout << "\nCAUCHY DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_CAUCHY;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Smale
    std::cout << "\nSMALE DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_SMALE;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Arcsine
    std::cout << "\nARCSINE DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_ARCSINE;
    data.resize(2,0.0);
    data[0] = 0.0; data[1] = 1.0;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);

    // Kumaraswamy
    std::cout << "\nKUMARASWAMY DISTRIBUTION\n\n";
    ed = ROL::DISTRIBUTION_KUMARASWAMY;
    data.resize(4,0.0);
    data[0] = 0.0; data[1] = 1.0; data[2] = 0.5; data[3] = 0.5;
    dist = Teuchos::rcp(new ROL::Distribution<double>(ed,data) );
    dist->test(*outStream);
  } 
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try
    
  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
