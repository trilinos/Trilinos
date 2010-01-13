//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA

************************************************************************
*/
//@HEADER

#ifndef _Isorropia_EpetraZoltanLib_hpp_
#define _Isorropia_EpetraZoltanLib_hpp_

#include <Isorropia_ConfigDefs.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Isorropia_EpetraCostDescriber.hpp>
#include <Isorropia_EpetraLibrary.hpp>

#include <QueryObject.hpp>
#include <zoltan_cpp.h>

#ifdef HAVE_EPETRA
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_CrsGraph;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

namespace Isorropia {

namespace Epetra {
  class CostDescriber;


class ZoltanLibClass : public Library {
public:

  ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph, int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph, 
		 Teuchos::RCP<const Epetra_MultiVector> input_coords, int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
	         Teuchos::RCP<CostDescriber> costs, int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_CrsGraph> input_graph, Teuchos::RCP<CostDescriber> costs, 
                 Teuchos::RCP<const Epetra_MultiVector> input_coords, Teuchos::RCP<const Epetra_MultiVector> weights, 
                 int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, 
                 Teuchos::RCP<const Epetra_MultiVector> input_coords, int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix,
	         Teuchos::RCP<CostDescriber> costs, int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_RowMatrix> input_matrix, Teuchos::RCP<CostDescriber> costs, 
		 Teuchos::RCP<const Epetra_MultiVector> input_coords, Teuchos::RCP<const Epetra_MultiVector> weights,
                 int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_MultiVector> input_coords, int inputType=unspecified_input_);
  ZoltanLibClass(Teuchos::RCP<const Epetra_MultiVector> input_coords,
            Teuchos::RCP<const Epetra_MultiVector> weights, int inputType=unspecified_input_);




  /** Method to partition the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern partitioning. 

      \param[out]  newPartitions The new partition for each of my objects, in
                                   local ID order.  The objects may be rows or
                               non-zeroes (for
                               CrsGraph and RowMatrix input) or coordinates (for
                               MultiVector input).  Partition numbers can range from
                               zero to numProcs-1.
      \param[out]  exportsSize  The number of my objects that will be exported to
                              another process under the new partitioning.  This is
                             also the number of elements in newPartitions that are
                             not equal to my process rank.
      \param[out]  imports   A list of the global IDs of the objects that will be
                            imported to my process under the new partitioning
   */

  virtual int
  repartition(Teuchos::ParameterList& paramlist,
	      std::vector<int>& newPartitions,
	      int& exportsSize,
	      std::vector<int>& imports);

  /** Method to color the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern coloring. 

      \param[out]  colorAssignment A list of integers indicating the coloring of
                              the object, in local ID order.
  */
  virtual int
  color(Teuchos::ParameterList& paramlist,
	std::vector<int>& colorAssignment);

  /** Method to order the object that the ZoltanLibClass was contructed with.

      \param[in] paramlist  Parameters to govern ordering . 

      \param[out]  orderAssignment A list of integers indicating the ordering of
                              the object, in local ID order.
  */
  virtual int
  order(Teuchos::ParameterList& paramlist,
	std::vector<int>& orderAssignment);

protected:
  virtual int precompute();
  virtual int postcompute();
  void computeCost();
  void preCheckPartition();

  void setParameterList(Teuchos::ParameterList& zoltanParamList);

private:
  Teuchos::ParameterList zoltanParamList_;
  std::string partMethod_; // stores partitioning method used, perhaps should be in EpetraLibrary?
  Zoltan *zz_;
  Teuchos::RCP<ZoltanLib::QueryObject> queryObject_;
  int num_obj_;

};//class ZoltanLibClass

}//namespace Epetra
}//namespace Isorropia

#endif //HAVE_EPETRA

#endif

