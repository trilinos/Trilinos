//@HEADER
//************************************************************************
//
//              Isorropia: Partitioning and Load Balancing Package
//                Copyright (2006) Sandia Corporation
//
//Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
//license for use of this work by or on behalf of the U.S. Government.
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
//************************************************************************
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
  ZoltanLibClass(Teuchos::RCP<const Epetra_BlockMap> input_map, int inputType=unspecified_input_);




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


#if defined(Isorropia_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Isorropia package is deprecated"
#endif
#endif

