#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Hashtable.hpp"
#include "Teuchos_HashSet.hpp"
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
              cerr << "caught bounds error: \n" <<  eb.what() << endl;
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
      

      /*-------------- do several tests of the HashSet class ------------------- */
      
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

      

      /*-------------- do several tests of the StrUtils class ------------------- */

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
      
      return 0;
    }
  catch(std::exception& e)
    {
      cerr << e.what() << endl;
    }
}
