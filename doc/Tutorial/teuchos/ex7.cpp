#include "Teuchos_CommandLineProcessor.hpp"

int main(int argc, char* argv[])
{
  // Creating an empty command line processor looks like:
  Teuchos::CommandLineProcessor My_CLP;

  /* To set and option, it must be given a name and default value.  Additionally,
     each option can be given a help string.  Although it is not necessary, a help
     string aids a users comprehension of the acceptable command line arguments.
     Some examples of setting command line options are:
  */
  // Set an integer command line option.
  int NumIters = 1550;
  My_CLP.setOption("iterations", &NumIters, "Number of iterations");
  // Set a double-precision command line option.
  double Tolerance = 1e-10;    
  My_CLP.setOption("tolerance", &Tolerance, "Tolerance");
  // Set a string command line option.
  string Solver = "GMRES";
  My_CLP.setOption("solver", &Solver, "Linear solver");
  // Set a boolean command line option.    
  bool Precondition;
  My_CLP.setOption("precondition","no-precondition",
		   &Precondition,"Preconditioning flag");

  /* There are also two methods that control the strictness of the command line processor.
     For a command line processor to be sensitive to any bad command line option that it 
     does not recognize use:
  */
  My_CLP.recogniseAllOptions(false);
  
  /* Then, if the parser finds a command line option it doesn't recognize, it will
     throw an exception.  To prevent a command line processor from throwing an exception 
     when it encounters a unrecognized option or help is printed, use:
  */
  My_CLP.throwExceptions(false);
  
  //Finally, to parse the command line, argc and argv are passed to the parse method:
  My_CLP.parse( argc, argv );

  return 0;
}
