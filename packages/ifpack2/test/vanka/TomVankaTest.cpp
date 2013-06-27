/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

// This will be a driver to test the development and implementation of
// the Vanka Smoother in Ifpack2. The tasks to be completed are as
// follows:
//
// 1. Translate Ifpack_UserPartition interface into Ifpack2 -- DONE
// 2. Feed it pressure LIDs -- SORTA DONE (ACTUALLY GIDs IN MY EXAMPLE)
// 3. Specify overlap = 1 -- DONE
// 4. Verify that this spits out the correct partitioning for Stokes Vanka -- DONE
// 5. Write a unit test -- TODO
// 6. Figure out how to efficiently add magnetics coupling -- EH KINDA SORTA DONE
// 7. Write a unit test -- TODO
//

#include <Teuchos_ConfigDefs.hpp>
#include <Ifpack2_ConfigDefs.hpp>

#include <Teuchos_DefaultComm.hpp>
#ifdef HAVE_MPI
#include <Teuchos_DefaultMpiComm.hpp>
#else
#include <Teuchos_DefaultSerialComm.hpp>
#endif
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <Teuchos_FancyOStream.hpp>

#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Partitioner.hpp>
#include <Ifpack2_OverlappingPartitioner.hpp>
#include <Ifpack2_UserPartitioner_def.hpp>
#include <Ifpack2_LinearPartitioner_def.hpp>
#include <Ifpack2_TomBlockRelaxation_def.hpp>
#include <Ifpack2_ILUT.hpp>
#include <Ifpack2_SparseContainer.hpp>
#include <Ifpack2_AdditiveSchwarz_def.hpp>

#include "Ifpack2_Parameters.hpp"

//#include <Kokkos_ConfigDefs.hpp>
//#include <Kokkos_SerialNode.hpp>

#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>

#include <Ifpack2_TomMHDPartitioner_def.hpp>

int main(int argc, char *argv[]){
  
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef double Scalar;
  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;

  
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> sparse_matrix_type;
  typedef Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> crsgraph;
  typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
  //typedef Kokkos::SerialNode node_type;
  using Teuchos::RCP;
  

  LocalOrdinal invalid = Teuchos::OrdinalTraits<LocalOrdinal>::invalid();


  // Initialize a "FancyOStream" to output to standard out (cout)
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);

  //Define communicator:
  Tpetra::DefaultPlatform::DefaultPlatformType& platform = Tpetra::DefaultPlatform::getDefaultPlatform();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = platform.getComm();



  // **************************************** //
  // MHD TEST                                 //
  // **************************************** //


  std::string mhdFileName = "MHDQ2Q1Q2Matrix_16x16.mm";

  RCP<sparse_matrix_type> myMHDA = reader_type::readSparseFile(mhdFileName, comm, platform.getNode(), true);

  std::cout << "The size of A is: " << myMHDA->getGlobalNumRows() << " x "
            << myMHDA->getGlobalNumCols() << std::endl;
  

  // START THE PARTITIONER TEST!
  Teuchos::ArrayRCP<LocalOrdinal> myMHDMap((size_t) 3556,invalid);
  for (LocalOrdinal jj=2178; jj<2178+289; jj++)
    {
      myMHDMap[jj] = jj - 2178;
    }

  Teuchos::ParameterList MHDList;
  MHDList.set("partitioner: type"       ,"TomMHD");
  MHDList.set("partitioner: local parts",(LocalOrdinal) 289);
  MHDList.set("partitioner: map"        ,myMHDMap);
  MHDList.set("partitioner: overlap"    ,(size_t) 1);
  MHDList.set("partitioner: nv"         ,(GlobalOrdinal) 2178);
  MHDList.set("partitioner: np"         ,(GlobalOrdinal) 289);

  Ifpack2::TomMHDPartitioner<crsgraph> MyMHDPart(myMHDA->getGraph());
  MyMHDPart.setParameters(MHDList);
  MyMHDPart.compute();


  MHDList.set("relaxation: type","Symmetric Gauss-Seidel");
  MHDList.set("relaxation: sweeps", (size_t) 1);
  MHDList.set("relaxation: damping factor", 1.0);
  //  MHDList.set("relaxation: 

  std::cout << std::endl << std::endl << std::endl;

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> CRS;
  typedef Ifpack2::ILUT< Tpetra::CrsMatrix<Scalar,LocalOrdinal,LocalOrdinal,Node>    > ILUTlo;
  typedef Ifpack2::ILUT< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>   > ILUTgo;

  /*
  Ifpack2::TomBlockRelaxation<CRS, Ifpack2::SparseContainer<CRS,Ifpack2::ILUT<CRS> > > prec(myMHDA);
  
  prec.setParameters(MHDList);
  prec.initialize();
  prec.compute();
  */

  ///////////////////////////////////////

  std::cout << "Beginning construction of AdditiveSchwarz stuff...\n\n";
  
  int overlapLevel=0;
  Ifpack2::AdditiveSchwarz<CRS,Ifpack2::TomBlockRelaxation<CRS, Ifpack2::SparseContainer<CRS,Ifpack2::ILUT<CRS> > > > myprec(myMHDA,overlapLevel);
  myprec.setParameters(MHDList);

  std::cout << "Beginning AdditiveSchwarz::initialize()...\n\n";
  myprec.initialize();

  std::cout << "Beginning AdditiveSchwarz::compute()...\n\n";
  myprec.compute();  

  //  Ifpack2::AdditiveSchwarz< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>, Ifpack2::TomBlockRelaxation< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>, Ifpack2::ILUT< Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > > > addschw(myMHDA,0);



  //std::ofstream MyOutStream("FullMHDIfpackVerify.txt");
  //MyMHDPart.print(MyOutStream);
  //MyMHDPart.print(std::cout);
  
  //MyOutStream.close();

  //prec.describe(*out,Teuchos::VERB_HIGH);

  //typedef Ifpack2::TomBlockRelaxation<CRS, Ifpack2::SparseContainer<CRS,ILUTlo> > TOM_Block_Relax;

  //Let's try AdditiveSchwarz
  //Ifpack2::AdditiveSchwarz<CRS, TOM_Block_Relax> addSchwarz(myMHDA,0); 

  //  addSchwarz.setParameters(MHDList);
  //  addSchwarz.initialize();
  //  addSchwarz.compute();

  /*  for (int kk=0; kk<289; kk++)
    {
      Teuchos::ArrayRCP<LocalOrdinal> myMHDList(76,0);
      MyMHDPart.rowsInPart(kk,myMHDList);
      MHDFile << std::setw(3) << std::right << kk << ": ";
      for (int ll=0; ll < myMHDList.size(); ll++)
        {
          MHDFile << myMHDList[ll] << " ";
        }
      MHDFile << std::endl;

    }
  */


}//int main
