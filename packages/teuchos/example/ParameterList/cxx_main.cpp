/*
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
*/

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_Array.hpp"	
#include "Teuchos_Version.hpp"

int main(int argc, char* argv[])
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  // Creating an empty parameter list looks like:
  Teuchos::ParameterList My_List;

  // Setting parameters in this list can be easily done:
  My_List.set("Max Iters", 1550, "Determines the maximum number of iterations in the solver");
  My_List.set("Tolerance", 1e-10, "The tolerance used for the convergence check");
  
  // For the "Solver" option, create a validator that will automatically
  // create documentation for this parameter but will also help in validation.
  Teuchos::RCP<Teuchos::StringToIntegralParameterEntryValidator<int> >
    solverValidator = Teuchos::rcp(
      new Teuchos::StringToIntegralParameterEntryValidator<int>(
        Teuchos::tuple<std::string>( "GMRES", "CG", "TFQMR" )
        ,"Solver"
        )
      );
  My_List.set(
    "Solver"
    ,"GMRES" // This will be validated by solverValidator right here!
    ,"The type of solver to use."
    ,solverValidator
    );

  /* The templated ``set'' method should cast the input {\it value} to the
     correct data type.  However, in the case where the compiler is not casting the input
     value to the expected data type, an explicit cast can be used with the ``set'' method:
  */
  My_List.set("Tolerance", (float)(1e-10), "The tolerance used for the convergence check");

  /* Reference-counted pointers can also be passed through a Teuchos::ParameterList.
     To illustrate this we will use the Teuchos::Array class to create an array of 10 doubles
     representing an initial guess for a linear solver, whose memory is being managed by a 
     Teuchos::RCP.
   */

  Teuchos::RCP<Teuchos::Array<double> > rcp_Array = 
    Teuchos::rcp( new Teuchos::Array<double>( 10, 0.0 ) );
  
  My_List.set("Initial Guess", rcp_Array, "The initial guess as a RCP to an array object.");

  /* A hierarchy of parameter lists can be constructed using {\tt Teuchos::ParameterList}.  This 
     means another parameter list is a valid {\it value} in any parameter list.  To create a sublist
     in a parameter list and obtain a reference to it:
  */
  Teuchos::ParameterList&
    Prec_List = My_List.sublist("Preconditioner",false,"Sublist that defines the preconditioner.");

  // Now this parameter list can be filled with values:
  Prec_List.set("Type", "ILU", "The tpye of preconditioner to use");
  Prec_List.set("Drop Tolerance", 1e-3
                ,"The tolerance below which entries from the\n""factorization are left out of the factors.");

  // The parameter list can be queried about the existance of a parameter, sublist, or type:
  // Has a solver been chosen?
  bool solver_defined = false, prec_defined = false, dtol_double = false;
  solver_defined = My_List.isParameter("Solver");
  // Has a preconditioner been chosen?
  prec_defined = My_List.isSublist("Preconditioner"); 
  // Has a tolerance been chosen and is it a double-precision number?
  bool tol_double = false;
  tol_double = My_List.INVALID_TEMPLATE_QUALIFIER isType<double>("Tolerance");
  // Has a drop tolerance been chosen and is it a double-precision number?
  dtol_double = Teuchos::isParameterType<double>(Prec_List, "Drop Tolerance"); 

  /* The last two methods for checking the parameter type are equivalent.
     There is some question as to whether the syntax of the first type-checking
     method is acceptable to older compilers.  Thus, the second type-checking method
     is offered as a portable alternative.
  */
  // Parameters can be retrieved from the parameter list in quite a few ways:
  // Get method that creates and sets the parameter if it doesn't exist.
  int its = 0;
  its = My_List.get("Max Iters", 1200);
  float tol;
  // Get method that retrieves a parameter of a particular type.
  tol = My_List.INVALID_TEMPLATE_QUALIFIER get<float>("Tolerance");
  // Get the "Solver" value and validate!
  std::string
    solver = solverValidator->validateString(
      Teuchos::getParameter<std::string>(My_List,"Solver")
      );

  /* In the above example, the first ``get'' method is a safe way of
     obtaining a parameter when its existence is indefinite but required.
     The second ``get'' method should be used when the existense of the parameter
     is definite.  This method will throw an std::exception if the parameter doesn't exist. 
     The safest way to use the second ``get'' method
     is in a try/catch block:
  */
  try {
    tol = My_List.INVALID_TEMPLATE_QUALIFIER get<float>("Tolerance");
  }
  catch ( const std::exception& e) {
    tol = 1e-6;
  }

  /* The second ``get'' method uses a syntax that may not be
     acceptable to older compilers.  Optionally, there is another portable templated 
     ``get'' function that can be used in the place of the second ``get'' method:
  */
  try {
    tol = Teuchos::getParameter<float>(My_List, "Tolerance");
  }
  catch ( const std::exception& e) {
    tol = 1e-6;
  }

  /* We can use this same syntax to get reference-counted objects out of the parameter
     list, like the initial guess.
  */
  try {
    Teuchos::RCP<Teuchos::Array<double> > init_guess = 
      Teuchos::getParameter<Teuchos::RCP<Teuchos::Array<double> > >(My_List, "Initial Guess");
  }
  catch ( std::exception& e) {
    std::cout << e.what() << std::endl;
  }

  std::cout << "\n# Printing this parameter list using opeator<<(...) ...\n\n";
  std::cout << My_List << std::endl;

  std::cout << "\n# Printing the parameter list only showing documentation fields ...\n\n";
  My_List.print(std::cout,Teuchos::ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true)); 

  /* It is important to note that mispelled parameters 
     (with additional space characters, capitalizations, etc.) may be ignored.  
     Therefore, it is important to be aware that a given parameter has not been used. 
     Unused parameters can be printed with method:
  */ 
  std::cout << "\n# Showing unused parameters ...\n\n";
  My_List.unused( std::cout );

  return 0;
}
