#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Version.hpp"

int main(int argc, char* argv[])
{
  cout << Teuchos::Teuchos_Version() << endl << endl;

  // Creating an empty parameter list looks like:
  Teuchos::ParameterList My_List;

  // Setting parameters in this list can be easily done:
  My_List.set("Max Iters", 1550);
  My_List.set("Tolerance", 1e-10);
  My_List.set("Solver", "GMRES");

  /* The templated ``set'' method should cast the input {\it value} to the
     correct data type.  However, in the case where the compiler is not casting the input
     value to the expected data type, an explicit cast can be used with the ``set'' method:
  */
  My_List.set("Tolerance", (float)(1e-10));

  /* A hierarchy of parameter lists can be constructed using {\tt Teuchos::ParameterList}.  This 
     means another parameter list is a valid {\it value} in any parameter list.  To create a sublist
     in a parameter list and obtain a reference to it:
  */
  Teuchos::ParameterList& Prec_List = My_List.sublist("Preconditioner");

  // Now this parameter list can be filled with values:
  Prec_List.set("Type", "ILU");
  Prec_List.set("Drop Tolerance", 1e-3);

  // The parameter list can be queried about the existance of a parameter, sublist, or type:
  // Has a solver been chosen?
  bool solver_defined = My_List.isParameter("Solver");
  // Has a preconditioner been chosen?
  bool prec_defined = My_List.isSublist("Preconditioner"); 
#ifdef HAVE_TEMPLATE_QUALIFIER
  // Has a tolerance been chosen and is it a double-precision number?
  bool tol_double = My_List.template isType<double>("Tolerance");
#endif
  // Has a drop tolerance been chosen and is it a double-precision number?
  bool dtol_double = Teuchos::isParameterType<double>(Prec_List, "Drop Tolerance"); 

  /* The last two methods for checking the parameter type are equivalent.
     There is some question as to whether the syntax of the first type-checking
     method is acceptable to older compilers.  Thus, the second type-checking method
     is offered as a portable alternative.
  */
  // Parameters can be retrieved from the parameter list in quite a few ways:
  // Get method that creates and sets the parameter if it doesn't exist.
  int its = My_List.get("Max Iters", 1200);
  float tol;
#ifdef HAVE_TEMPLATE_QUALIFIER
  // Get method that retrieves a parameter of a particular type.
  tol = My_List.template get<float>("Tolerance");

  /* In the above example, the first ``get'' method is a safe way of
     obtaining a parameter when its existence is indefinite but required.
     The second ``get'' method should be used when the existense of the parameter
     is definite.  This method will throw an exception if the parameter doesn't exist. 
     The safest way to use the second ``get'' method
     is in a try/catch block:
  */
  try {
    tol = My_List.template get<float>("Tolerance");
  }
  catch ( std::exception& e) {
    tol = 1e-6;
  }
#endif

  /* The second ``get'' method uses a syntax that may not be
     acceptable to older compilers.  Optionally, there is another portable templated 
     ``get'' function that can be used in the place of the second ``get'' method:
  */
  try {
    tol = Teuchos::getParameter<float>(My_List, "Tolerance");
  }
  catch ( std::exception& e) {
    tol = 1e-6;
  }

  // A parameter list can be sent to the output stream:
  cout<< My_List << endl;

  /* It is important to note that mispelled parameters 
     (with additional space characters, capitalizations, etc.) may be ignored.  
     Therefore, it is important to be aware that a given parameter has not been used. 
     Unused parameters can be printed with method:
  */ 
  My_List.unused( cout );

  return 0;
}
