// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_Hashtable.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_ErrorPolling.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_Version.hpp"

using namespace Teuchos;
using std::string;

/* Test of Teuchos container classes */

int main( int argc, char* argv[] )
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try
    {

      if (MPIComm::world().getRank() == 0)
        {
          /*------- do several tests of the Array class --------------- */

          /* You can create an empty array and append to it. */
          Array<double> x;
          x.append(1.0);
          x.append(4.0);
          x.append(9.0);

          /* there is a toString() method of Array that writes a list bounded by
           * curly braces */
          std::cerr << "x = " << x.toString() << std::endl;


          /* You can create an array of a specified size */
          Array<double> y(3);

          /* Array elements can be set using the [] indexing operator */
          y[0] = 3.14;
          y[1] = 2.72;
          y[2] = 1.42;

          /* Array elements can be read using the const [] indexing operator */
          for (int i=0; i<y.length(); i++)
            {
              std::fprintf(stderr, "%d %g\n", i, y[i]);
            }

          /* If compiled with boundschecking, Array will catch bounds violations. */
          if (Array<double>::hasBoundsChecking())
            {
              bool caughtBoundsError = false;
              try
                {
                  double z = y[10];
                  (void)z;
                }
              catch(std::exception& eb)
                {
                  caughtBoundsError = true;
                  std::cerr << "caught [expected] bounds error: \n" <<  eb.what() << std::endl;
                }
              if (!caughtBoundsError)
                {
                  std::cerr << "FAILED TO CATCH BOUNDS ERROR" << std::endl;
                }
            }
          else
            {
              std::cerr << "Teuchos compiled w/o array boundschecking" << std::endl;
            }


          /*------------ test packing of Arrays ------------------- */

          std::cerr << "testing packing of arrays..." << std::endl;
          Array<std::string> schools;
          schools.append("Cornell");
          schools.append("Princeton");
          schools.append("Carnegie Mellon");
          schools.append("North Carolina");
          schools.append("Maryland");
          schools.append("Caltech");
          schools.append("Chicago");

          Array<char> packed;
          MPIContainerComm<std::string>::pack(schools, packed);

          Array<std::string> unpacked;
          MPIContainerComm<std::string>::unpack(packed, unpacked);

          bool ok = true;
          if (unpacked.size() != schools.size())
            {
              ok = false;
            }
          else
            {
              for (unsigned int i=0; i<schools.size(); i++)
                {
                  if (unpacked[i] != schools[i])
                    {
                      ok = false;
                    }
                }
            }
      
          if (!ok)
            {
              std::cerr << "pack/unpack FAILED!" << std::endl;
              std::cerr << "original: " << schools << std::endl;
              std::cerr << "unpacked: " << unpacked << std::endl;
            }
          else
            {
              std::cerr << "pack/unpack OK" << std::endl;
            }
        }

      /*------------ test gathering of Arrays ------------------- */
      
      Array<std::string> teams;
      int rank = MPIComm::world().getRank();
      if (rank==0)
        {
          teams.append("Orioles");
          teams.append("Yankees");
          teams.append("Pirates");
        }
      else if (rank==1)
        {
          teams.append("Ravens");
          teams.append("Jets");
          teams.append("Steelers");
          teams.append("Colts");
          teams.append("Raiders");
        }
      else 
        {
          teams.append("Bruins");
          teams.append("Rangers");
        }

      Array<Array<std::string> > allTeams;
      MPIContainerComm<std::string>::gatherv(teams, allTeams, 0, MPIComm::world());

      if (rank==0)
        {
          std::cout << "all teams = " << allTeams << std::endl;
        }
      

      /*------------ test polling of exceptions across procs --------- */
      
      try
        {
          try
            {
              TEST_FOR_EXCEPTION(MPIComm::world().getRank()==1,
                              std::runtime_error, 
                              "std::exception [expected]");
            }
          catch(std::exception& ex1)
            {
              std::cerr << "successful detection of std::exception on proc="
                   << MPIComm::world().getRank() << std::endl;
              ErrorPolling::reportFailure(MPIComm::world());
              TEUCHOS_TRACE(ex1);
            }
          TEUCHOS_POLL_FOR_FAILURES(MPIComm::world());
        }
      catch(std::exception& ex)
        {
          std::cerr << ex.what() << std::endl;
        }
      std::cerr << "p=" << MPIComm::world().getRank() 
           << ": std::exception polling successful" << std::endl;

      if (MPIComm::world().getRank()==0)
        {
      

          /*-------- do several tests of the HashSet class ------------- */
      
          /* Add entries to a set using the put method */
          HashSet<std::string> trilinosPackages;
          trilinosPackages.put("epetra");
          trilinosPackages.put("ml");
          trilinosPackages.put("TSF");
          trilinosPackages.put("nox");
          trilinosPackages.put("meros");
      
          /* count entries using the size() method */
          std::fprintf(stderr, "trilinos has %d packages\n", trilinosPackages.size());

          /* write to a std::string using the toString() method */
          std::cerr << "trilinos packages are: " << trilinosPackages.toString() << std::endl;

          /* test for the presence of a member using the containsKey() method */
      
          if (trilinosPackages.containsKey("epetra"))
            {
              std::cerr << "epetra is in the list of trilinos packages" << std::endl;
            }
          else
            {
              std::cerr << "epetra is not in the list of trilinos packages" << std::endl;
            }

          if (trilinosPackages.containsKey("Space Invaders"))
            {
              std::cerr << "Space Invaders is in the list of trilinos packages" << std::endl;
            }
          else
            {
              std::cerr << "Space Invaders is not in the list of trilinos packages" << std::endl;
            }


          /*-------------- do several tests of the Hashtable class -------- */

          /* add entries using the put() method */
          Hashtable<std::string, int> battles;
      
          battles.put("hastings",    1066);
          battles.put("waterloo",    1815);
          battles.put("gettysburg",  1863);
          battles.put("verdun",      1916);
          battles.put("midway",      1942);
          battles.put("normandy",    1944);

          /* write to a std::string using the toString() method */
          std::cerr << "hashtable is: " << battles.toString() << std::endl;
      
          /* test for the presence of a key using the containsKey() method */
          if (battles.containsKey("cannae"))
            {
              std::fprintf(stderr, "the battle of cannae occured in %d\n", battles.get("cannae"));
            }
          else
            {
              std::cerr << "cannae is not in our hashtable" << std::endl;
            }

          /* test for the presence of a key using the containsKey() method */
          if (battles.containsKey("verdun"))
            {
              std::fprintf(stderr, "the battle of verdun occured in %d\n", battles.get("verdun"));
            }
          else
            {
              std::cerr << "verdun is not in our hashtable" << std::endl;
            }

          /* remove a member of the hashtable (bug# 2983)*/
          battles.remove( "waterloo" );

          /* write to a std::string using the toString() method */
          std::cerr << "hashtable is (after removal of waterloo): " << battles.toString() << std::endl;
   

          /*-------------- do several tests of the StrUtils class --------- */

          /* stringTokenizer() splits a std::string into whitespace-delimited tokens */
          std::string test = "Sandia National Laboratories";

          Array<std::string> tokens = StrUtils::stringTokenizer(test);

          std::cerr << "tokens = " << tokens.toString() << std::endl;

          /* atof() converts a std::string to its double value */
          double pi = StrUtils::atof("3.14159265358");
          std::fprintf(stderr, "pi = %g, tan(pi/4)=%g\n", pi, std::tan(pi/4.0));

          /* atoi() converts a std::string to its integer value */
          int a = StrUtils::atoi("-101");
          std::fprintf(stderr, "a = %d\n", a);

          /* allCaps() converts to upper case */
          std::cerr << "all caps: " << StrUtils::allCaps(test) << std::endl;
        }      
    }
  catch(std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }

}
