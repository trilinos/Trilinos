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

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_Version.hpp"


using std::string;
using Teuchos::XMLObject;
using Teuchos::StringInputSource;
using Teuchos::FileInputSource;

/* Test of Teuchos XML handling classes */

int main(int argc, char** argv)
{
  std::cout << Teuchos::Teuchos_Version() << std::endl << std::endl;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  try
   {
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

      int foundIndex  = problem.findFirstChild("Solver");
      if(foundIndex == -1)
      {
        std::cerr << "Find child didn't find the child!"
          <<std::endl << std::endl;
        return -1;
      }

      const XMLObject foundChild = problem.getChild(foundIndex);
      if(foundChild.getTag() != solver.getTag())
      {
        std::cerr << "Find child found the wrong tag!" << std::endl <<
          "Found index was: " << foundIndex << std::endl <<
          std::endl << std::endl;
        return -1;
      }

      if(problem.findFirstChild("NON EXSISTENT CHILD") != -1){
        std::cerr << "First first child didn't return -1 when it was "
          "suppose to!" <<std::endl << std::endl;
        return -1;
      }

      std::string str = problem.toString();
      std::cerr << str << std::endl;

      /* parse XML in a std::string */
      StringInputSource src(str);
      XMLObject reread = src.getObject();
      
      std::cerr << reread << std::endl;

      /* write to a file, and then read and parse the file */
      std::ofstream of("tmp.xml");
      of << reread << std::endl;
      
      FileInputSource fileSrc("tmp.xml");
      XMLObject fileXML = fileSrc.getObject();
      
      std::cerr << fileXML << std::endl;

      return 0;
    }
  catch(std::exception& e)
    {
      std::cerr << e.what() << std::endl;
    }
}
