// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_EigenVerify_hpp
#define percept_EigenVerify_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <iostream>

#include <stk_mesh/base/Field.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

namespace stk
{
  namespace io
  {
    class StkMeshIoBroker;
  }
}

  namespace percept
  {

    class EigenVerify
    {
    public:
      EigenVerify(stk::ParallelMachine comm_in) : 
    comm(comm_in),
	num_meshes(2),
	meshes_in(num_meshes,""),
	mesh_out(),
	file_out(),
	field_name(),
	mesh_data(),
	clp(),
	num_time_steps_all(num_meshes),
	time_steps_all(num_meshes),
	inputField(num_meshes),
	errorField(NULL),
	fieldAll(num_meshes),
	xferFieldAll(num_meshes)
	{}

      void run(int argc,  char** argv);

    private:
      void process_options();

      void create_mesh_data(stk::io::StkMeshIoBroker * mesh_data, 
			    const std::string &filename);

      void load_time_data(const int m);

      void create_fields(const int num_time_steps);

      void load_field_data(const int m);

      int get_num_time_steps();

      stk::ParallelMachine comm;
      const int num_meshes;

      std::vector<std::string> meshes_in;
      std::string mesh_out;
      std::string file_out;

      std::string field_name;

      std::vector<stk::io::StkMeshIoBroker *> mesh_data;

      Teuchos::CommandLineProcessor clp;

      std::vector<int> num_time_steps_all;
      std::vector<std::vector<double> > time_steps_all;

      std::vector<stk::mesh::Field<double> *> inputField;
      stk::mesh::Field<double> * errorField;
      std::vector<stk::mesh::Field<double> *> fieldAll;
      std::vector<stk::mesh::Field<double> *> xferFieldAll;
    };

      
  }//namespace percept

#endif
