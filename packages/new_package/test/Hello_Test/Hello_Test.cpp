#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Newp_Hello.h"
#include <iostream>
#include <sstream>
#include <istream>
#include <ostream>
#include <string>

int main(){
  using namespace std;
  Epetra_SerialComm esc;
  Epetra_Comm * ec;
  ec = dynamic_cast<Epetra_Comm*>(&esc);
  stringbuf stringb;
  streambuf * sb;
  sb = dynamic_cast<streambuf*>(&stringb);
  iostream io(sb);
  ostream * os;
  os = dynamic_cast<ostream*>(&io);
  Newp_Hello nph(*ec);
  nph.Print(*os);
  char results[100];
  io.readsome(results, 99);
  char * expected = "This will print out one line for each of the 1 processes \n\nHello.  I am process 0\n";
  if(strcmp(results, expected) != 0){
    cout << "Test failed!" << endl;
    cout << "Expected:" << endl << expected << endl;
    cout << "Got:" << endl << results << endl;
    return 1;
  }
  cout << "Test passed!" << endl;
  return 0;
}
