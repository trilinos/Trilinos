#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_FileInputSource.hpp"



using namespace Teuchos;
using std::string;

/* Test of Teuchos XML handling classes */

int main(int argc, void** argv)
{
  try
    {
      MPISession::init(&argc, &argv);

      /* create an XML object */
      XMLObject problem("Problem");
      XMLObject solver("Solver");
      XMLObject prec("Preconditioner");

      solver.addAttribute("type", "gmres");
      solver.addInt("maxiters", 1000);
      solver.addInt("restarts", 100);
      solver.addDouble("tol", 1.0e-10);

      solver.addChild(prec);

      prec.addAttribute("type", "ILUk");
      prec.addInt("k", 2);

      problem.addChild(solver);

      string str = problem.toString();
      cerr << str << endl;

#ifdef HAVE_EXPAT

      /* parse XML in a string */
      StringInputSource src(str);
      XMLObject reread = src.getObject();
      
      cerr << reread << endl;

      /* write to a file, and then read and parse the file */
      ofstream of("tmp.xml");
      of << reread << endl;
      
      FileInputSource fileSrc("tmp.xml");
      XMLObject fileXML = fileSrc.getObject();
      
      cerr << fileXML << endl;
      
#endif

      return 0;
    }
  catch(std::exception& e)
    {
      cerr << e.what() << endl;
    }
  MPISession::finalize();
}
