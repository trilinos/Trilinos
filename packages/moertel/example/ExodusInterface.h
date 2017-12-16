/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/

#ifndef MOERTEL_EXODUS_H
#define MOERTEL_EXODUS_H

#include <iostream>
#include <iomanip>
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Galeri_AbstractGrid.h"

// Exodus stuff

#include "exodusII.h"

#define log_ex_err(o, s) o << "Line No: " << __LINE__ << \
                   " Function: " << __FUNCTION__ << \
                   " Error: " << s // note I leave ; out

class ExodusInterface
{

public:

  ExodusInterface(Epetra_Comm& Comm) : Comm_(Comm) {}

  ~ExodusInterface() {}

  const Epetra_Comm& Comm() const
  {
    return(Comm_);
  }

  void Write(const Galeri::FiniteElements::AbstractGrid& data, const string& BaseName,
             const Epetra_MultiVector& Field){

	int comp_ws = sizeof(double); // = 8
	int io_ws = sizeof(double); // = 8

	string FileName = BaseName + ".exo";
	int ex_id = ex_create(FileName.c_str(), EX_CLOBBER, &comp_ws, &io_ws);

	int num_dim = data.NumDimensions();
	int num_nodes = data.NumGlobalVertices(); 
	int num_elem = data.NumGlobalElements(); 
	int num_elem_blk = 1;
	int num_node_sets = 0;
	int num_side_sets = 0;

	int ex_err = ex_put_init(ex_id, FileName.c_str(), num_dim, 
	   num_nodes, num_elem, num_elem_blk, 
	   num_node_sets, num_side_sets);

    if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");

    std::vector<double> coord(3);
    std::vector<int>    vertices(data.NumVerticesPerElement());

    const Epetra_Map& RowMap = data.RowMap();
//    const Epetra_Map& VertexMap = data.VertexMap();

    std::vector<double> x(data.NumMyVertices());
    std::vector<double> y(data.NumMyVertices());
    std::vector<double> z(data.NumMyVertices());

    for (int i = 0 ; i < data.NumMyVertices() ; ++i)
    {
      data.VertexCoord(i, &coord[0]);
      x[i] = coord[0];
      y[i] = coord[1];
      z[i] = coord[2];
    }

    int n = 0;
    if (Field.Comm().MyPID() == 0)
      n = RowMap.NumGlobalElements();

    Epetra_Map SingleProcMap(-1, n, 0, Field.Comm());
    Epetra_MultiVector SingleProcField(SingleProcMap, 1);

    Epetra_Import FieldImporter(SingleProcMap, RowMap);
    SingleProcField.Import(Field, FieldImporter, Insert);

    if (Comm().MyPID() == 0)
    {
      switch (data.NumDimensions()) {
      case 2:{
	       const char* coord_names[] = {"x", "y"};
	       ex_err = ex_put_coord_names(ex_id, (char**)coord_names);
           if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	       ex_err = ex_put_coord(ex_id, &x[0], &y[0], NULL);
           if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	     }
        break;
      case 3: {
	       const char* coord_names[] = {"x", "y", "z"};
	       ex_err = ex_put_coord_names(ex_id, (char**)coord_names);
           if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	       ex_err = ex_put_coord(ex_id, &x[0], &y[0], &z[0]);
           if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	      }
        break;
      default:
        throw(-1);
      }

    }
    Comm().Barrier();

    for (int ProcID = 0 ; ProcID < Comm().NumProc() ; ++ProcID) {

      if (Comm().MyPID() == ProcID) {

        if (ProcID == 0) {
          string type = data.ElementType();

          if (type == "GALERI_TRIANGLE"){
		       const char * elem_type = "TRIANGLE";
		       ex_err = ex_put_block(ex_id, EX_ELEM_BLOCK, 1, (char *)elem_type, data.NumGlobalElements(), 3, 0, 0, 0);
               if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	  }
          else if (type == "GALERI_QUAD"){
		       const char * elem_type = "QUAD4";
		       ex_err = ex_put_block(ex_id, EX_ELEM_BLOCK, 1, (char *)elem_type, data.NumGlobalElements(), 4, 0, 0, 0);
               if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	  }
          else if (type == "GALERI_TET"){
		       const char * elem_type = "TETRA";
		       ex_err = ex_put_block(ex_id, EX_ELEM_BLOCK, 1, (char *)elem_type, data.NumGlobalElements(), 4, 0, 0, 0);
               if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	  }
          else if (type == "GALERI_HEX"){
		       const char * elem_type = "HEX";
		       ex_err = ex_put_block(ex_id, EX_ELEM_BLOCK, 1, (char *)elem_type, data.NumGlobalElements(), 8, 0, 0, 0);
               if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	  }
          else
          {
            cerr << "Incorrect element type (" << type << ")" << endl;
            throw(-1);
          }
        }

	    std::vector<int> connect_tmp(data.NumMyElements() * data.NumVerticesPerElement());
	    int cnt = 0;

        for (int i = 0 ; i < data.NumMyElements() ; ++i) {
          data.ElementVertices(i, &vertices[0]);
          for (int j = 0 ; j < data.NumVerticesPerElement() ; ++j)
		  connect_tmp[cnt++] = data.VertexMap().GID(vertices[j]) + 1;
        }

	    ex_err = ex_put_conn(ex_id, EX_ELEM_BLOCK, 1, &connect_tmp[0], 0, 0);
        if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");
	

      }

	  int num_nodal_fields = 1;
      std::vector<double> field(data.NumMyVertices());

	  ex_err = ex_put_variable_param (ex_id, EX_NODAL, num_nodal_fields);
      if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");

	  const char* var_names[] = {"u"};
      ex_err = ex_put_variable_names (ex_id, EX_NODAL, num_nodal_fields, (char**)var_names);
      if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");

	  for(int i = 0; i < data.NumMyVertices(); i++)
			  field[i] = Field[0][i];

	  ex_err = ex_put_var (ex_id, 1, EX_NODAL, 1, 1, data.NumMyVertices(), &field[0]);
      if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");

	  ex_err = ex_close(ex_id);
      if(ex_err) log_ex_err(std::cout, "ex_error is non-zero");

      Comm().Barrier();

    } // for Procs, write elements

  }

private:
  const Epetra_Comm& Comm_;

};

#endif
