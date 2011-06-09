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
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/

#ifndef MOERTEL_EXODUS
#define MOERTEL_EXODUS

#include <iostream>
#include <iomanip>
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Galeri_AbstractGrid.h"

// Exodus stuff

#include "exodusII.h"

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

    vector<double> coord(3);
    vector<int>    vertices(data.NumVerticesPerElement());

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
	       ex_err = ex_put_coord(ex_id, &x[0], &y[0], NULL);
	     }
        break;
      case 3: {
	       const char* coord_names[] = {"x", "y", "z"};
	       ex_err = ex_put_coord_names(ex_id, (char**)coord_names);
	       ex_err = ex_put_coord(ex_id, &x[0], &y[0], &z[0]);
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
		       ex_err = ex_put_elem_block(ex_id, 1, (char *)elem_type, data.NumGlobalElements(), 3, 0);
	  }
          else if (type == "GALERI_QUAD"){
		       const char * elem_type = "QUAD4";
		       ex_err = ex_put_elem_block(ex_id, 1, (char *)elem_type, data.NumGlobalElements(), 4, 0);
	  }
          else if (type == "GALERI_TET"){
		       const char * elem_type = "TETRA";
		       ex_err = ex_put_elem_block(ex_id, 1, (char *)elem_type, data.NumGlobalElements(), 4, 0);
	  }
          else if (type == "GALERI_HEX"){
		       const char * elem_type = "HEX";
		       ex_err = ex_put_elem_block(ex_id, 1, (char *)elem_type, data.NumGlobalElements(), 8, 0);
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

	    ex_err = ex_put_elem_conn(ex_id, 1, &connect_tmp[0]);
	

      }

 	/* Write the field data out
	
		ex_err = ex_put_nodal_var(int exodus_file_id, int time_step, 
		 int nodal_var_index, // which field set is being written
		 int num_nodes, // number of nodes worth of data
		 void *nodal_var_vals  // the data
		);
	*/


	  int num_nodal_fields = 1;
      std::vector<double> field(data.NumMyVertices());

	  ex_err = ex_put_var_param (ex_id, "N", num_nodal_fields);

	  const char* var_names[] = {"u"};
      ex_err = ex_put_var_names (ex_id, "N", num_nodal_fields, (char**)var_names);

//	  for(int i = 0; i < data.NumMyVertices(); i++)
//			  field[i] = SingleProcField[0][i];
	  for(int i = 0; i < data.NumMyVertices(); i++)
			  field[i] = Field[0][i];

	  ex_err = ex_put_nodal_var (ex_id, 1, 1, data.NumMyVertices(), &field[0]);

	  ex_err = ex_close(ex_id);

      Comm().Barrier();

    } // for Procs, write elements

  }

private:
  const Epetra_Comm& Comm_;

};

#endif
