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
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_Version.hpp"


using namespace Teuchos;
using std::string;

/* Test of Teuchos XML handling classes */

int main(int argc, void** argv)
{
  cout << Teuchos::Teuchos_Version() << endl << endl;

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
