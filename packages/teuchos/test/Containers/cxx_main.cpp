// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Hashtable.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_ErrorPolling.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_Version.hpp"

using namespace Teuchos;
using std::string;

/* Test of Teuchos container classes */

int main(int argc, char** argv)
{
  cout << Teuchos::Teuchos_Version() << endl << endl;

  try
    {
      MPISession::init(&argc, (void***) &argv);


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
          cerr << "x = " << x.toString() << endl;


          /* You can create an array of a specified size */
          Array<double> y(3);

          /* Array elements can be set using the [] indexing operator */
          y[0] = 3.14;
          y[1] = 2.72;
          y[2] = 1.42;

          /* Array elements can be read using the const [] indexing operator */
          for (int i=0; i<y.length(); i++)
            {
              fprintf(stderr, "%d %g\n", i, y[i]);
            }

          /* If compiled with boundschecking, Array will catch bounds violations. */
          if (Array<double>::hasBoundsChecking())
            {
              bool caughtBoundsError = false;
              try
                {
                  double z = y[10];
                }
              catch(std::exception& eb)
                {
                  caughtBoundsError = true;
                  cerr << "caught [expected] bounds error: \n" <<  eb.what() << endl;
                }
              if (!caughtBoundsError)
                {
                  cerr << "FAILED TO CATCH BOUNDS ERROR" << endl;
                }
            }
          else
            {
              cerr << "Teuchos compiled w/o array boundschecking" << endl;
            }


          /*------------ test packing of Arrays ------------------- */

          cerr << "testing packing of arrays..." << endl;
          Array<string> schools;
          schools.append("Cornell");
          schools.append("Princeton");
          schools.append("Carnegie Mellon");
          schools.append("North Carolina");
          schools.append("Maryland");
          schools.append("Caltech");
          schools.append("Chicago");

          Array<char> packed;
          MPIContainerComm<string>::pack(schools, packed);

          Array<string> unpacked;
          MPIContainerComm<string>::unpack(packed, unpacked);

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
              cerr << "pack/unpack FAILED!" << endl;
              cerr << "original: " << schools << endl;
              cerr << "unpacked: " << unpacked << endl;
            }
          else
            {
              cerr << "pack/unpack OK" << endl;
            }
        }

      /*------------ test gathering of Arrays ------------------- */
      
      Array<string> teams;
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

      Array<Array<string> > allTeams;
      MPIContainerComm<string>::gatherv(teams, allTeams, 0, MPIComm::world());

      if (rank==0)
        {
          cout << "all teams = " << allTeams << endl;
        }
      

      /*------------ test polling of exceptions across procs --------- */
      
      try
        {
          try
            {
              TEST_FOR_EXCEPTION(MPIComm::world().getRank()==1,
                              runtime_error, 
                              "exception [expected]");
            }
          catch(std::exception& ex1)
            {
              cerr << "successful detection of exception on proc="
                   << MPIComm::world().getRank() << endl;
              ErrorPolling::reportFailure(MPIComm::world());
              TEUCHOS_TRACE(ex1);
            }
          TEUCHOS_POLL_FOR_FAILURES(MPIComm::world());
        }
      catch(std::exception& ex)
        {
          cerr << ex.what() << endl;
        }
      cerr << "p=" << MPIComm::world().getRank() 
           << ": exception polling successful" << endl;

      if (MPIComm::world().getRank()==0)
        {
      

          /*-------- do several tests of the HashSet class ------------- */
      
          /* Add entries to a set using the put method */
          HashSet<string> trilinosPackages;
          trilinosPackages.put("epetra");
          trilinosPackages.put("ml");
          trilinosPackages.put("TSF");
          trilinosPackages.put("nox");
          trilinosPackages.put("meros");
      
          /* count entries using the size() method */
          fprintf(stderr, "trilinos has %d packages\n", trilinosPackages.size());

          /* write to a string using the toString() method */
          cerr << "trilinos packages are: " << trilinosPackages.toString() << endl;

          /* test for the presence of a member using the containsKey() method */
      
          if (trilinosPackages.containsKey("epetra"))
            {
              cerr << "epetra is in the list of trilinos packages" << endl;
            }
          else
            {
              cerr << "epetra is not in the list of trilinos packages" << endl;
            }

          if (trilinosPackages.containsKey("Space Invaders"))
            {
              cerr << "Space Invaders is in the list of trilinos packages" << endl;
            }
          else
            {
              cerr << "Space Invaders is not in the list of trilinos packages" << endl;
            }


          /*-------------- do several tests of the Hashtable class -------- */

          /* add entries using the put() method */
          Hashtable<string, int> battles;
      
          battles.put("hastings",    1066);
          battles.put("waterloo",    1815);
          battles.put("gettysburg",  1863);
          battles.put("verdun",      1916);
          battles.put("midway",      1942);
          battles.put("normandy",    1944);

          /* write to a string using the toString() method */
          cerr << "hashtable is: " << battles.toString() << endl;
      
          /* test for the presence of a key using the containsKey() method */
          if (battles.containsKey("cannae"))
            {
              fprintf(stderr, "the battle of cannae occured in %d\n", battles.get("cannae"));
            }
          else
            {
              cerr << "cannae is not in our hashtable" << endl;
            }

          /* test for the presence of a key using the containsKey() method */
          if (battles.containsKey("verdun"))
            {
              fprintf(stderr, "the battle of verdun occured in %d\n", battles.get("verdun"));
            }
          else
            {
              cerr << "verdun is not in our hashtable" << endl;
            }

          /* remove a member of the hashtable (bug# 2983)*/
          battles.remove( "waterloo" );

          /* write to a string using the toString() method */
          cerr << "hashtable is (after removal of waterloo): " << battles.toString() << endl;
   

          /*-------------- do several tests of the StrUtils class --------- */

          /* stringTokenizer() splits a string into whitespace-delimited tokens */
          string test = "Sandia National Laboratories";

          Array<string> tokens = StrUtils::stringTokenizer(test);

          cerr << "tokens = " << tokens.toString() << endl;

          /* atof() converts a string to its double value */
          double pi = StrUtils::atof("3.14159265358");
          fprintf(stderr, "pi = %g, tan(pi/4)=%g\n", pi, tan(pi/4.0));

          /* atoi() converts a string to its integer value */
          int a = StrUtils::atoi("-101");
          fprintf(stderr, "a = %d\n", a);

          /* allCaps() converts to upper case */
          cerr << "all caps: " << StrUtils::allCaps(test) << endl;
        }      
    }
  catch(std::exception& e)
    {
      cerr << e.what() << endl;
    }
  MPISession::finalize();
}
