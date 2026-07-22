// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MiniEM_AddFieldsToMesh.hpp"
#include "Panzer_String_Utilities.hpp"
#include "Panzer_STK_Interface.hpp"
#include <string>
#include <vector>

void mini_em::addFieldsToMesh(panzer_stk::STK_Interface & mesh,
                               const Teuchos::ParameterList & output_list)
{
  // register cell averaged scalar fields
  const Teuchos::ParameterList & cellAvgQuants = output_list.sublist("Cell Average Quantities");
  for(Teuchos::ParameterList::ConstIterator itr=cellAvgQuants.begin();
      itr!=cellAvgQuants.end();++itr) {
    const std::string & blockId = itr->first;
    const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
    std::vector<std::string> tokens;

    // break up comma seperated fields
    panzer::StringTokenizer(tokens,fields,",",true);

    for(std::size_t i=0;i<tokens.size();i++)
      mesh.addCellField(tokens[i],blockId);
  }

  // register cell averaged components of vector fields
  // just allocate space for the fields here. The actual calculation and writing
  // are done by panzer_stk::ScatterCellAvgVector.
  const Teuchos::ParameterList & cellAvgVectors = output_list.sublist("Cell Average Vectors");
  for(Teuchos::ParameterList::ConstIterator itr = cellAvgVectors.begin();
      itr != cellAvgVectors.end(); ++itr) {
    const std::string & blockId = itr->first;
    const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
    std::vector<std::string> tokens;

    // break up comma seperated fields
    panzer::StringTokenizer(tokens,fields,",",true);

    for(std::size_t i = 0; i < tokens.size(); i++) {
      std::string d_mod[3] = {"X","Y","Z"};
      for(std::size_t d = 0; d < mesh.getDimension(); d++)
        mesh.addCellField(tokens[i]+d_mod[d],blockId);
    }
  }

  // register cell quantities
  const Teuchos::ParameterList & cellQuants = output_list.sublist("Cell Quantities");
  for(Teuchos::ParameterList::ConstIterator itr=cellQuants.begin();
      itr!=cellQuants.end();++itr) {
    const std::string & blockId = itr->first;
    const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
    std::vector<std::string> tokens;

    // break up comma seperated fields
    panzer::StringTokenizer(tokens,fields,",",true);

    for(std::size_t i=0;i<tokens.size();i++)
      mesh.addCellField(tokens[i],blockId);
  }

  // register ndoal quantities
  const Teuchos::ParameterList & nodalQuants = output_list.sublist("Nodal Quantities");
  for(Teuchos::ParameterList::ConstIterator itr=nodalQuants.begin();
      itr!=nodalQuants.end();++itr) {
    const std::string & blockId = itr->first;
    const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
    std::vector<std::string> tokens;

    // break up comma seperated fields
    panzer::StringTokenizer(tokens,fields,",",true);

    for(std::size_t i=0;i<tokens.size();i++)
      mesh.addSolutionField(tokens[i],blockId);
  }

  const Teuchos::ParameterList & allocNodalQuants = output_list.sublist("Allocate Nodal Quantities");
  for(Teuchos::ParameterList::ConstIterator itr=allocNodalQuants.begin();
      itr!=allocNodalQuants.end();++itr) {
    const std::string & blockId = itr->first;
    const std::string & fields = Teuchos::any_cast<std::string>(itr->second.getAny());
    std::vector<std::string> tokens;

    // break up comma seperated fields
    panzer::StringTokenizer(tokens,fields,",",true);

    for(std::size_t i=0;i<tokens.size();i++)
      mesh.addSolutionField(tokens[i],blockId);
  }

}

