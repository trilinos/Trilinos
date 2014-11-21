// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER
#include "TrilinosCouplings_Pamgen_Utils.hpp"

#include <stdexcept>

// Pamgen includes
#include "create_inline_mesh.h"
#include "pamgen_extras.h"

void TrilinosCouplings::pamgen_error_check(std::ostream & os, long long cr_result){
  if (cr_result == ERROR_PARSING_DEFINITION){
    int essz = getPamgenEchoStreamSize();
    char * echo_char_array = new char[essz+1];
    os<<"PARSE ERROR\n";
    echo_char_array[essz] = '\0';
    echo_char_array = getPamgenEchoStream(echo_char_array);
    if(echo_char_array) os<<echo_char_array;
    if(cr_result == ERROR_CREATING_IMD) os<<"ERROR Failure to create Inline_Mesh_Desc creation\n";
    delete [] echo_char_array;
  }
    
  if(cr_result == ERROR_CREATING_MS){
    int essz = getPamgenErrorStreamSize();
    char * error_char_array = new char[essz+1];
    error_char_array[essz] = '\0';
    error_char_array = getPamgenErrorStream(error_char_array);
    if(error_char_array) os<<error_char_array;
    os<<"\nERROR Failure to create Mesh_Specification\n";
    delete [] error_char_array;  
  }
        
  int wssz = getPamgenWarningStreamSize();
  if(wssz){
    char * warning_char_array = new char[wssz+1];
    warning_char_array[wssz] = '\0';
    warning_char_array = getPamgenWarningStream(warning_char_array);
    os<<"WARNING Records\n";
    os<<warning_char_array;
    delete [] warning_char_array;  
  }
   
  if(cr_result!=ERROR_FREE_CREATION) throw std::runtime_error("Error in pamgen mesh creation.");
}
