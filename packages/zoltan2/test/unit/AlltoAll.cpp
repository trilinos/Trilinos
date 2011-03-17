// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER

// TODO: doxygen comments

#include <iostream>
#include <Zoltan2_AlltoAll.hpp>
#include <Teuchos_GlobalMPISession.hpp>   
#include <Teuchos_RCP.hpp>   
#include <Teuchos_ArrayRCP.hpp>   
#include <Teuchos_Comm.hpp>   
#include <Teuchos_DefaultComm.hpp>   

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  int rank = session.getRank();
  int nprocs = session.getNProc();

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // test of Z2::AlltoAll using a ordinal type

  int *outBuf = new int [2*nprocs];
  Teuchos::ArrayRCP<int> sendBuf = Teuchos::arcp(outBuf, 0, 2*nprocs);

  for (int i=0, j=1; i < 2*nprocs ; i+=2,j++){
    outBuf[i] = j*10;
    outBuf[i+1] = j*10 + 1;
  } 

  Teuchos::ArrayRCP<int> recvBuf;
  
  Z2::AlltoAll<int, long, int>(*comm, 
                sendBuf,    // ints to send from this process to all the others
                2,          // two ints per process
                recvBuf);   // will be allocated and filled in AlltoAll

  int *inBuf = recvBuf.get();

  bool pass = true;

  int myvals[2] = {(rank+1) * 10, (rank+1) * 10 + 1};

  for (int i=0; i < 2*nprocs; i+=2){
    if (inBuf[i] != myvals[0] && inBuf[i+1] != myvals[1]){
      pass = false;
      break;
    }
  }

  if (pass)
    std::cout << "PASS" << std::endl;
  else
    std::cout << "FAIL" << std::endl;
}
