// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   
//
// ***********************************************************************
// @HEADER
//
// This is another test where you need to look at the output to
// know if it's right.  This should be fixed.


#include <Zoltan2_StridedInput.hpp>   
#include <ErrorHandlingForTests.hpp>   

#include <Teuchos_RCP.hpp>   
#include <Teuchos_ArrayRCP.hpp>   
#include <Teuchos_Comm.hpp>   
#include <Teuchos_DefaultComm.hpp>   

using Zoltan2::StridedInput;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using namespace std;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = 
    Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();

  if (rank > 0)
    return 0;

  Teuchos::RCP<const Zoltan2::Environment> envPtr = 
    Teuchos::rcp(new Zoltan2::Environment);

  // Strided input with stride 1

  Array<float> input1(12);
  for (int i=0; i < 12; i++)
    input1[i] = (i+1) * 5; 

  RCP<StridedInput<int, float> > s1;

  try{
    s1 = rcp(new StridedInput<int, float>(envPtr, input1.view(0, 12), 1));
  }
  catch{
    TEST_FAIL_AND_EXIT(comm, 0, "Error in constructor 1", 1);
  }

  std::cout << std::endl;
  std::cout << "Test 1, input: " << input1 << std::endl;
  std::cout << "[] test: ";
  for (int i=0; i < 12; i++)
    std::cout << (*s1)[i] << " ";
  std::cout << std::endl;

  ArrayRCP<const float> fromS1;
  s1->getInputArray(fromS1);
  std::cout << "getInputArray test: ";
  for (int i=0; i < 12; i++)
    std::cout << fromS1[i] << " ";
  std::cout << std::endl;

  StridedInput s1Copy<int, float>;
  s1Copy = *s1;

  std::cout << "assignment operator test: ";
  for (int i=0; i < 12; i++)
    std::cout << s1Copy[i] << " ";
  std::cout << std::endl;

  // Strided input with stride 3

  Array<long> input2(12, -1);
  for (int i=0; i < 12; i+=3)
    input2[i] = (i+1) * 5; 

  RCP<StridedInput<int, long> > s2;

  try{
    s2 = rcp(new StridedInput<int, long>(envPtr, input2.view(0, 12), 3));
  }
  catch{
    TEST_FAIL_AND_EXIT(comm, 0, "Error in constructor 2", 2);
  }

  std::cout << std::endl;
  std::cout << "Test 2, input: " << input2 << std::endl;
  std::cout << "[] test: ";
  for (int i=0; i < 12; i++)
    std::cout << (*s2)[i] << " ";
  std::cout << std::endl;

  ArrayRCP<const long> fromS2;
  s2->getInputArray(fromS2);
  std::cout << "getInputArray test: ";
  for (int i=0; i < 4; i++)
    std::cout << fromS2[i] << " ";
  std::cout << std::endl;

  StridedInput s2Copy<int, long>;
  s2Copy = *s2;

  std::cout << "assignment operator test: ";
  for (int i=0; i < 4; i++)
    std::cout << s2Copy[i] << " ";
  std::cout << std::endl;

  return 0;
}
