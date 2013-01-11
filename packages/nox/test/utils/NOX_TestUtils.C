//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
#include "NOX_TestUtils.H"
#include <sstream>
#include <fstream>

bool NOX::getNextQuotedString(const std::string& line, std::string::size_type& pos, std::string& value)
{
  // Initialize value
  value = "";

  // Compute the length of the line
  std::string::size_type linelength = line.length();

  // Find the location of the first quote
  std::string::size_type pos1 = line.find('"', pos); 
  
  // Check that the operation succeeded and that we're not at the end of the line
  if ((pos1 == std::string::npos) || (pos1 == linelength - 1))
  {
    pos = std::string::npos;
    return false;
  }

  // Advance to first character after the first quote
  pos1 = pos1 + 1;

  // Find the location of the second quote
  std::string::size_type pos2 = line.find('"', pos1); 

  // Check that the operation was successful
  if (pos2 == std::string::npos)
  {
    pos = std::string::npos;
    return false;
  }

  // Compute the length of the expression
  std::string::size_type length = pos2 - pos1;

  // Compute the final position
  pos = (pos2 == (linelength - 1)) ? std::string::npos : pos2 + 1;

  // Extract the substring
  value = line.substr(pos1, length);

  // Return true
  return true;
}

bool NOX::getNextString(const std::string& line, std::string::size_type& pos, std::string& value)
{
  // Initialize value
  value = "";

  // Compute the length of the line
  std::string::size_type linelength = line.length();

  // Find the location of the first non-space character
  std::string::size_type pos1 = line.find_first_not_of(' ', pos); 
  
  // Check that the operation succeeded 
  if (pos1 == std::string::npos) 
  {
    pos = std::string::npos;
    return false;
  }

  // Find the location of the next space, if any
  std::string::size_type pos2 = line.find(' ', pos1); 

  // Compute the length of the expression
  std::string::size_type length = (pos2 == std::string::npos) ? linelength - pos1 : pos2 - pos1;

  // Compute the final position
  pos = (pos2 == (linelength - 1)) ? std::string::npos : pos2 + 1;

  // Extract the substring
  value = line.substr(pos1, length);

  // Return true
  return true;
}

bool NOX::getNextInt(const std::string& line, std::string::size_type& pos, int& value)
{
  std::string field;
  if ((!NOX::getNextString(line,pos,field)) || (field.size() == 0))
    return false;

  return ((sscanf(field.c_str(), "%d", &value)) == 1);
}

bool NOX::getNextDouble(const std::string& line, std::string::size_type& pos, double& value)
{
  std::string field;
  if ((!NOX::getNextString(line,pos,field)) || (field.size() == 0))
    return false;

  return ((sscanf(field.c_str(), "%le", &value)) == 1);
}

bool NOX::processTextInputFileLine(const std::string& line, Teuchos::ParameterList& params, 
					Teuchos::ParameterList*& subPtr)
{
  
  std::string::size_type pos;	// current position inside line

  std::string field;			// name of parameter

  //bool tmpbool;			// for reading in a boolean
  int tmpint;			// for reading in an int
  double tmpdouble;		// for reading in a double
  std::string tmpstring;		// for reading in a std::string
  //Value tmpvalue;		// for reading in an Value
  //NOX::Entry tmpvalue;          // darin: use NOX::entry instead of Value
  //Vector tmpvector;		// for reading in an Vector

  std::string type;			// parameter type (used in new style)

  static int n = -1;			// size of vector (used in old style)
  char c;			// used to read in '=' character in old style

  if (line.size() == 0) 	// empty line
  {
    return false;
  }
  else if (line[0] == '#')	// comment
  {
    return false;
  }
  else if (line[0] == '@')	// sublist command
  {
    subPtr = &params;
    
    pos = 0;
    
    while (pos != std::string::npos)
    {
      if ((NOX::getNextQuotedString(line, pos, field)) && (field.size() > 0))
	subPtr = &(subPtr->sublist(field));
    }

    return true;
  }
  else if (line[0] == '"')	// new style
  {

    // Get the name
    pos = 0;
    if ((!NOX::getNextQuotedString(line, pos, field)) || (field.empty()))
      return false;
    
    // Read in the type
    if (!NOX::getNextString(line, pos, type))
      return false;

    if (type == "int")
    {
      if (!NOX::getNextInt(line, pos, tmpint))
      {
	return false;
      }
      else
      {
	subPtr->set(field, tmpint);
	return true;
      }
    }
    else if (type == "bool")
    {
      if (!NOX::getNextString(line, pos, tmpstring))
      {
	return false;
      }
      else
      {
	subPtr->set(field, (tmpstring == "true") );
	return true;
      }
    }
    else if (type == "double")
    {
      if (!NOX::getNextDouble(line, pos, tmpdouble))
      {
	return false;
      }
      else
      {
	subPtr->set(field, tmpdouble);
	return true;
      }
    }
    else if (type == "string")
    {
      if ((!NOX::getNextQuotedString(line, pos, tmpstring)) || (tmpstring.empty()))
      {
	return false;
      }
      else
      {
	subPtr->set(field, tmpstring);
	return true;
      }
    }
    /*else if (type == "vector")
    {
      // get the size
      if (!getNextInt(line, pos, tmpint))
	return false;

      if (tmpint < 0)
	return false;

      tmpvector.resize(tmpint);
      for (int i = 0; i < tmpint; i ++)
      {
	if (!getNextDouble(line, pos, tmpvector[i]))
	  return false;
      }
      
      subPtr->set(field, tmpvector);	
      return true;
      }*/
    else
      return false;

  } // end new style

  else				// old style
  {
    std::istringstream linein(line);
    linein >> field;
    
    if ((field == "n_parameters") || (field == "continuous_design"))
    {
      linein >> c;
      linein >> n;
    }    
    
    /*else if ((field == "initial_point") || (field == "cdv_initial_point"))
    {
      if (n == -1)
      {
	cout << "Error: initial_point specified before n_parameters " << std::endl;
	return false;
      }
      
      linein >> c;
      tmpvector.resize(n);
      for (int i = 0; i < n; i ++)
	linein >> tmpvector[i];
      params.sublist("Solver").set("Initial X", tmpvector);
    }
    */

    /*else if (field == "initial_f") 
    {
      linein >> c;
      linein >> tmpdouble;
      tmpvalue.setValueTo(true, tmpdouble);
      params.sublist("Solver").set("Initial F", tmpvalue);
    }
    */

    /*else if ((field == "lower_bounds") || (field == "cdv_lower_bounds"))
    {
      if (n == -1)
      {
	cout << "Error:lower_bounds specified before n_parameters " << std::endl;
	return false;
      }
      
      linein >> c;
      tmpvector.resize(n);
      for (int i = 0; i < n; i ++)
	linein >> tmpvector[i];
      params.sublist("Bounds").set("Lower", tmpvector);
      
    }
    */

    /*else if ((field == "upper_bounds") || (field == "cdv_upper_bounds"))
    {
      if (n == -1)
      {
	cout << "Error:upper_bounds specified before n_parameters " << std::endl;
	return false;
      }
      
      linein >> c;
      tmpvector.resize(n);
      for (int i = 0; i < n; i ++)
	linein >> tmpvector[i];
      params.sublist("Bounds").set("Upper", tmpvector);
      
    }
    */

    else if (field == "executable") 
    {
      linein >> c;
      linein >> tmpstring;
      params.sublist("Evaluator").set("Executable Name", tmpstring);
    }
    
    else if (field == "params_prefix") 
    {
      linein >> c;
      linein >> tmpstring;
      params.sublist("Evaluator").set("Input Prefix", tmpstring);
    }
    
    else if (field == "result_prefix") 
    {
      linein >> c;
      linein >> tmpstring;
      params.sublist("Evaluator").set("Output Prefix", tmpstring);
    }
    
    else 
    {
      std::cout << "Ignoring unrecognized field." << field << "\n";
    }
    
  } /* old style */
  return true;
}

bool NOX::parseTextInputFile(const std::string filename, Teuchos::ParameterList& params)
{
  // Open the input file
  std::ifstream fin;
  fin.open(filename.c_str());

  if (!fin) 
  {
    std::cout << "Error: Cannot find input file " << filename << std::endl;
    return false;
  }

  std::string line;			// one line from input file
  Teuchos::ParameterList* subPtr;	// sublist pointer

  subPtr = &params;

  while (!fin.eof())
  {
    getline(fin, line);
      
    if (!NOX::processTextInputFileLine(line, params, subPtr))
    {
      std::cout << "Ignoring input file line: " << line << std::endl;
    }
  }

  fin.close();

  
  return true;
}
		    
