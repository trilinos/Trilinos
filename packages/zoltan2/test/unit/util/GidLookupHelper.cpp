// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
// @HEADER
// ***********************************************************************
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
// ***********************************************************************
// @HEADER
//
// Test the GidLookupHelper class.

#include <Zoltan2_GidLookupHelper.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using Teuchos::Comm;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();
  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);
  ArrayRCP<int> idx;

  if (rank == 0){
    int *buf = new int [10];
    for (int i=0; i < 10; i++)
      buf[i] = (i+1)*10;

    // ten distinct gids
    cout << endl << "10 distinct int gids" << endl;
    ArrayRCP<const int> intGids = Teuchos::arcp(buf, 0, 10, true);
    Zoltan2::GidLookupHelper<int, int> glhInt(env, intGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 10; i++){
      cout << "Gid: " << buf[i] << " idx: " << glhInt.lookup(buf[i]) << endl;
    }

    // ten gids with duplicates
    cout << endl << "10 int gids with duplicates" << endl;
    for (int i=0; i < 5; i++) 
      buf[i] = buf[i+5] = i;
    Zoltan2::GidLookupHelper<int, int> glhIntDup(env, intGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 10; i++){
      cout << "Gid: " << buf[i] << " idx: " << glhIntDup.lookup(buf[i]) << endl;
    }

    // string gids
    cout << endl << "4 string gids" << endl;
    string *sbuf = new string [4];
    sbuf[0] = string("gid number 0");
    sbuf[1] = string("gid number 1");
    sbuf[2] = string("gid number 2");
    sbuf[3] = string("gid number 3");
    ArrayRCP<const string> sGids = Teuchos::arcp(sbuf, 0, 4, true);
    Zoltan2::GidLookupHelper<string, int> glhString(env, sGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 4; i++){
      cout << "Gid: " << sbuf[i] << " idx: " << glhString.lookup(sbuf[i]) << endl;
    }
    
    // string gids with duplicates
    cout << endl << "4 string gids with a duplicate" << endl;
    sbuf[3] = string("gid number 0");
    Zoltan2::GidLookupHelper<string, int> glhStringDup(env, sGids);
    cout << "Lookup:" << endl;
    for (int i=0; i < 4; i++){
      cout << "Gid: " << sbuf[i] << " idx: " << glhStringDup.lookup(sbuf[i]) << endl;
    }

    cout << "PASS" << endl;
  }
}
