#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Hashtable.hpp"
#include "Teuchos_HashSet.hpp"
#include "Teuchos_Out.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;
using std::string;

/* Test of Teuchos container classes */

int main(int argc, char** argv)
{
  try
    {

      /*-------------- do several tests of the Array class ------------------- */

      /* You can create an empty array and append to it. */
      Array<double> x;
      x.append(1.0);
      x.append(4.0);
      x.append(9.0);

      /* there is a toString() method of Array that writes a list bounded by
       * curly braces */
      Out::println("x = " + x.toString());


      /* You can create an array of a specified size */
      Array<double> y(3);

      /* Array elements can be set using the [] indexing operator */
      y[0] = 3.14;
      y[1] = 2.72;
      y[2] = 1.42;

      /* Array elements can be read using the const [] indexing operator */
      for (int i=0; i<y.length(); i++)
        {
          Out::printf("%d %g\n", i, y[i]);
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
              Out::println(string("caught bounds error: \n") + eb.what());
            }
          if (!caughtBoundsError)
            {
              Out::println("FAILED TO CATCH BOUNDS ERROR");
            }
        }
      else
        {
          Out::println("Teuchos compiled w/o array boundschecking");
        }
      

      /*-------------- do several tests of the HashSet class ------------------- */
      
      /* Add entries to a set using the put method */
      HashSet<string> trilinosPackages;
      trilinosPackages.put("epetra");
      trilinosPackages.put("ml");
      trilinosPackages.put("TSF");
      trilinosPackages.put("nox");
      trilinosPackages.put("meros");
      
      /* count entries using the size() method */
      Out::printf("trilinos has %d packages\n", trilinosPackages.size());

      /* write to a string using the toString() method */
      Out::println("trilinos packages are: " + trilinosPackages.toString());

      /* test for the presence of a member using the containsKey() method */
      
      if (trilinosPackages.containsKey("epetra"))
        {
          Out::println("epetra is in the list of trilinos packages");
        }
      else
        {
          Out::println("epetra is not in the list of trilinos packages");
        }

      if (trilinosPackages.containsKey("Space Invaders"))
        {
          Out::println("Space Invaders is in the list of trilinos packages");
        }
      else
        {
          Out::println("Space Invaders is not in the list of trilinos packages");
        }

      /*-------------- do several tests of the Hashtable class ------------------- */

      /* add entries using the put() method */
      Hashtable<string, int> battles;
      
      battles.put("hastings",    1066);
      battles.put("waterloo",    1815);
      battles.put("gettysburg",  1863);
      battles.put("verdun",      1916);
      battles.put("midway",      1942);
      battles.put("normandy",    1944);

      /* write to a string using the toString() method */
      Out::println("hashtable is: " + battles.toString());
      
      /* test for the presence of a key using the containsKey() method */
      if (battles.containsKey("cannae"))
        {
          Out::printf("the battle of cannae occured in %d\n", battles.get("cannae"));
        }
      else
        {
          Out::println("cannae is not in our hashtable");
        }

      /* test for the presence of a key using the containsKey() method */
      if (battles.containsKey("verdun"))
        {
          Out::printf("the battle of verdun occured in %d\n", battles.get("verdun"));
        }
      else
        {
          Out::println("verdun is not in our hashtable");
        }

      

      /*-------------- do several tests of the StrUtils class ------------------- */

      /* stringTokenizer() splits a string into whitespace-delimited tokens */
      string test = "Sandia National Laboratories";

      Array<string> tokens = StrUtils::stringTokenizer(test);

      Out::println("tokens = " + tokens.toString());

      /* atof() converts a string to its double value */
      double pi = StrUtils::atof("3.14159265358");
      Out::printf("pi = %g, tan(pi/4)=%g\n", pi, tan(pi/4.0));

      /* atoi() converts a string to its integer value */
      int a = StrUtils::atoi("-101");
      Out::printf("a = %d\n", a);

      /* allCaps() converts to upper case */
      Out::println("all caps: " + StrUtils::allCaps(test));
      
      return 0;
    }
  catch(std::exception& e)
    {
      cerr << e.what() << endl;
    }
}
