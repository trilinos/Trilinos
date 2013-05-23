/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK2_TOMMHDPARTITIONER_DEF_HPP
#define IFPACK2_TOMMHDPARTITIONER_DEF_HPP
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_TomMHDPartitioner_decl.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include <vector>
#include <string>

namespace Ifpack2 {
//==============================================================================
template<class GraphType>
TomMHDPartitioner<GraphType>::TomMHDPartitioner(const Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> >& Graph) :
  NumLocalParts_(1),
  Graph_(Graph),
  OverlappingLevel_(0),
  IsComputed_(false),
  verbose_(false),
  UserPart_(Graph)
{
}

//==============================================================================
template<class GraphType>
TomMHDPartitioner<GraphType>::~TomMHDPartitioner()
{
}

//==============================================================================
template<class GraphType>
typename GraphType::local_ordinal_type TomMHDPartitioner<GraphType>::numLocalParts() const 
{
  return( NumLocalParts_);
}

//==============================================================================
template<class GraphType>
size_t TomMHDPartitioner<GraphType>::overlappingLevel() const
{
  return(OverlappingLevel_);
}

//==============================================================================
template<class GraphType>
typename GraphType::local_ordinal_type TomMHDPartitioner<GraphType>::operator() (LocalOrdinal MyRow) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(((MyRow < 0) || ((size_t)MyRow > Graph_->getNodeNumRows())), std::runtime_error, "Ifpack2::TomMHDPartitioner::operator() invalid row.");
   
  return(Partition_[MyRow]);
}


//==============================================================================
template<class GraphType>
typename GraphType::local_ordinal_type TomMHDPartitioner<GraphType>::operator() (LocalOrdinal i, LocalOrdinal j) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( ((i < 0) || ((i > NumLocalParts_) || (j < 0) || (j > Parts_[i].size()))),
			      std::runtime_error, "Ifpack2::TomMHDPartitioner::operator() invalid row or node.");
  return(Parts_[i][j]);
}

//==============================================================================
template<class GraphType>
size_t TomMHDPartitioner<GraphType>::numRowsInPart(const LocalOrdinal Part) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( ((Part < 0) || (Part > NumLocalParts_)),
			      std::runtime_error, "Ifpack2::TomMHDPartitioner::numRowsInPart() invalid partition.");
  return(Parts_[Part].size());
}

//==============================================================================
template<class GraphType>
void TomMHDPartitioner<GraphType>::rowsInPart(const LocalOrdinal Part,  Teuchos::ArrayRCP<LocalOrdinal> & List) const
{
  // Let numRowsInPart do the sanity checking...
  size_t numrows= numRowsInPart(Part); 
  for (size_t i = 0 ; i < numrows; i++)
    List[i] = Parts_[Part][i];
}

//==============================================================================
template<class GraphType>
Teuchos::ArrayView<const typename GraphType::local_ordinal_type>  TomMHDPartitioner<GraphType>::nonOverlappingPartition() const
{
  return(Partition_.view(0, Graph_->getNodeNumRows()));
}

//==============================================================================
template<class GraphType>
void TomMHDPartitioner<GraphType>::setParameters(Teuchos::ParameterList& List)
{
  NumLocalParts_    = List.get("partitioner: local parts", NumLocalParts_);
  OverlappingLevel_ = List.get("partitioner: overlap", OverlappingLevel_);
  verbose_          = List.get("partitioner: print level", verbose_);

  if (NumLocalParts_ < 0)
    NumLocalParts_ = Graph_->getNodeNumRows() / (-NumLocalParts_);
  if (NumLocalParts_ == 0)
    NumLocalParts_ = 1;
  
  // Sanity checking
  TEUCHOS_TEST_FOR_EXCEPTION( ((NumLocalParts_ < 0) || (size_t)NumLocalParts_ > Graph_->getNodeNumRows()),
			      std::runtime_error, "Ifpack2::TomMHDPartitioner::setParameters() invalid NumLocalParts_");
  TEUCHOS_TEST_FOR_EXCEPTION( (OverlappingLevel_ < 0),
			      std::runtime_error, "Ifpack2::TomMHDPartitioner::setParameters() invalid OverlappingLevel_");

  setPartitionParameters(List);


  // Pass the parameters to the overlapping partition:
  UserPart_.setParameters(List);
}


template<class GraphType>
void TomMHDPartitioner<GraphType>::setPartitionParameters(Teuchos::ParameterList& List)
{
  Map_ = List.get("partitioner: map",Map_);
  if (Map_ == Teuchos::null)
    {
      throw std::runtime_error("AAAAARRRGGGGGHHHh!");
    }
  NV_ = List.get("partitioner: nv",NV_);
  NP_ = List.get("partitioner: np",NP_);
  
}


//==============================================================================
template<class GraphType>
void TomMHDPartitioner<GraphType>::compute()
{
  TEUCHOS_TEST_FOR_EXCEPTION( ((NumLocalParts_ < 1) ||  (OverlappingLevel_ < 0)),
			      std::runtime_error, "Ifpack2::TomMHDPartitioner::compute() invalid NumLocalParts_ or OverlappingLevel_");

  std::string PrintMsg_("TomMHDPartitioner: ");

  // some output
  if (verbose_ && (Graph_->getComm()->getRank() == 0)) {
    std::cout << PrintMsg_ << "Number of local parts          = " << NumLocalParts_ << std::endl;
    std::cout << PrintMsg_ << "Approx. Number of global parts = " 
	      << NumLocalParts_ * Graph_->getComm()->getSize() << std::endl;
    std::cout << PrintMsg_ << "Amount of overlap              = " << OverlappingLevel_ << std::endl;
  }

  // 1.- Compute the user partition
  UserPart_.compute();

  // 2.- allocate memory
  // This is ok because the nonoverlapping partition and the number of partitions don't change
  Partition_.resize(Graph_->getNodeNumRows());
  Parts_.resize(NumLocalParts_);

  // 3.- sanity checks on input graph
  // Sure! Why not?
  TEUCHOS_TEST_FOR_EXCEPTION( (!Graph_->isFillComplete() || Graph_->getGlobalNumRows() != Graph_->getGlobalNumCols()),
			      std::runtime_error, "Ifpack2::TomMHDPartitioner::compute() input graph error");
 
  // 4.- perform non-overlapping partition 
  computePartitions();

  // 4.- compute the partitions with overlapping
  computeOverlappingPartitions();

  // 5.- mark as computed
  IsComputed_ = true;
}

//==============================================================================

template<class GraphType>
void TomMHDPartitioner<GraphType>::TomMHDPartitioner::computePartitions()
{

  if (!UserPart_.isComputed())
    {
      // Throw some exception
      throw std::runtime_error("AAAAARRRGGGGGHHHh!");
    }

  //Copy the nonoverlapping partition from the UserPartition
  for (size_t ii=0 ; ii < Graph_->getNodeNumRows() ; ++ii)
    {
      Partition_[ii] = UserPart_(ii);
    }


}

template<class GraphType>
void TomMHDPartitioner<GraphType>::computeOverlappingPartitions()
{

  // Most of this work will be done by UserPart_.compute(). All we
  // need to do is copy the UserPart_ partition and expand it to
  // include our magnetics



  //THIS WHOLE FUNCTION ONLY WORKS IN SERIAL!!
  using std::vector;

  //create storage for our temporary partition
  vector<vector<size_t> > tmp;
  tmp.resize(NumLocalParts_);


  // loop over all current partitions.
  for (int part=0 ; part < NumLocalParts_ ; ++part)
    {
      
      // The first entry in each partition is PRESSURE!
      //tmp[part].push_back(UserPart_(part,0));
      
      for(size_t col=1; col < (size_t)UserPart_.numRowsInPart(part) ; ++col) 
        {
          LocalOrdinal LRID = UserPart_(part,col);
          
          if ( (col==1) && (LRID % 2) ) 
            { //First entry (no previous), odd)
              tmp[part].push_back(LRID-1);
              tmp[part].push_back(LRID);
              continue;
            }
          else if ( (LRID % 2) && (LRID != UserPart_(part,col-1) + 1) )
            {
              // ODD and no matching previous
              tmp[part].push_back(LRID-1);
              tmp[part].push_back(LRID);
              continue;
            }
          else if ( !(LRID % 2) && (LRID+1 != UserPart_(part,col+1)) )
            {
              // EVEN and no matching follower
              tmp[part].push_back(LRID);
              tmp[part].push_back(LRID+1);
              continue;
            }
          else 
            {
              tmp[part].push_back(LRID);
            }
          
        }//for(size_t col=0...
      
      tmp[part].push_back(UserPart_(part,0));

      //ADD MAGNETICS COUPLING -- THIS ONLY WORKS IN SERIAL!!
      size_t velocity_size = tmp[part].size();
      for(size_t kk=1; kk < velocity_size ; kk=kk+2)
        {
          tmp[part].push_back(NV_+NP_+tmp[part][kk]/2);
        }
      
      
      //Now convert STL vectors to ArrayRCPs
      Parts_[part].resize(tmp[part].size());
      for (size_t ll=0; ll < tmp[part].size() ; ++ll)
        {
          Parts_[part][ll] = tmp[part][ll];
        }
      
      
    }// for (int part=0...


}

//==============================================================================
template<class GraphType>
bool TomMHDPartitioner<GraphType>::isComputed() const
{
  return(IsComputed_);
}

//==============================================================================
/*
template<class GraphType>
std::ostream& TomMHDPartitioner<GraphType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}
*/

  //! Prints the partitions to an output stream
template<class GraphType>
std::ostream& TomMHDPartitioner<GraphType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  
  for(int ii=0 ; ii < this->NumLocalParts_ ; ++ii){
    fos << std::setw(3) << ii << ": ";
    for (int jj=0 ; jj < this->Parts_[ii].size() ; ++jj)
        {
          fos <<  this->Parts_[ii][jj] << " ";
        }
    fos << std::endl;
    
  }

  return(os);
}


//==============================================================================
template<class GraphType>
std::string TomMHDPartitioner<GraphType>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isComputed()) {
    oss << "{status = computed";
  }
  else {
    oss << "{status = is not computed";
  }
  oss <<"}";
  return oss.str();
}

//==============================================================================
template<class GraphType>
void  TomMHDPartitioner<GraphType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  if(verbLevel==Teuchos::VERB_NONE) return;

  os << "================================================================================" << endl;
  os << "Ifpack2::TomMHDPartitioner" << endl;
  os << "Number of local rows  = " << Graph_->getNodeNumRows() << endl;
  os << "Number of global rows = " << Graph_->getGlobalNumRows() << endl;
  os << "Number of local parts = " << NumLocalParts_ << endl;
  os << "Overlapping level     = " << OverlappingLevel_ << endl;
  os << "Is computed           = " << IsComputed_ << endl;
  os << "================================================================================" << endl;
}


}// namespace Ifpack2

#endif // IFPACK2_TOMMHDPARTITIONER_DEF_HPP
