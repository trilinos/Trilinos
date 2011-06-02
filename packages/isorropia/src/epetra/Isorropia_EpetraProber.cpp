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
#include <Isorropia_EpetraProber.hpp>
#ifdef HAVE_EPETRA
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_Map.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Comm.h>

namespace Isorropia{
namespace Epetra{

  
Prober::Prober():input_graph_(0),colorer_(0),has_colored(false){}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Prober::Prober(Teuchos::RCP<const Epetra_CrsGraph> input_graph,
               const Teuchos::ParameterList& paramlist,
               bool compute_now)
:input_graph_(input_graph),colorer_(0),has_colored(false)
{
  setList(paramlist);
  if(compute_now) color();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Prober::Prober(const Epetra_CrsGraph * input_graph,
               const Teuchos::ParameterList& paramlist, bool compute_now)
	       :input_graph_(Teuchos::rcp<const Epetra_CrsGraph>(input_graph,false)),
                colorer_(0),has_colored(false)
{
  setList(paramlist);
  if(compute_now) color();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Prober::Prober(Teuchos::RCP<const Epetra_CrsMatrix> input_matrix,
               const Teuchos::ParameterList & paramlist, bool compute_now)
:input_graph_(Teuchos::rcp<const Epetra_CrsGraph>(&input_matrix->Graph(),false)),
 colorer_(0),List_(paramlist),has_colored(false)
{
  setList(paramlist);
  if(compute_now) color();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Prober::Prober(const Epetra_CrsMatrix * input_matrix,
               const Teuchos::ParameterList & paramlist, bool compute_now)
:input_graph_(Teuchos::rcp<const Epetra_CrsGraph>(&input_matrix->Graph(),false)),
 colorer_(0),List_(paramlist),has_colored(false)
{
  setList(paramlist);
  if(compute_now) color();
}
////////////////////////////////////////////////////////////////////////////////
  


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Prober::setList(const Teuchos::ParameterList& paramlist)
{
  List_=paramlist;

  /* Ensure Distance 2 */
  Teuchos::ParameterList dummy;
  string two,s_dummy;
  Teuchos::ParameterList zoltan=List_.get("ZOLTAN",dummy);
//   two=zoltan.get("DISTANCE",s_dummy);
//   if(strcmp(two.c_str(),"2")){
//     if(!input_graph_->Comm().MyPID())
//       cerr<<"WARNING: Prober requires distance 2 coloring.  Resetting coloring type for you."<<endl;
//     zoltan.set("DISTANCE","2");
//     List_.set("ZOLTAN",zoltan);
//   }

}
////////////////////////////////////////////////////////////////////////////////
  
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void Prober::color()
{
  if(!has_colored) 
  {
    delete colorer_;
    colorer_=new Isorropia::Epetra::Colorer(input_graph_,List_,false);
  }

  colorer_->color(true);
  has_colored=true;
}
////////////////////////////////////////////////////////////////////////////////

int Prober::getNumOrthogonalVectors()
{
    if (has_colored)
        return colorer_->numColors();
    else
        return -1;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int Prober::probe(const Epetra_Operator & op, Epetra_CrsMatrix & out_matrix)
{
  /* Sanity Checks */
  if(input_graph_.is_null()) return -1;
  if(input_graph_->DataPtr() != out_matrix.Graph().DataPtr()) return -1;

  if(!has_colored) color();
  int Ncolors=colorer_->numColors();
  int N=out_matrix.NumMyRows();

  if(Ncolors==0) return -1;

  /* Allocs */
  Epetra_MultiVector temp1(out_matrix.DomainMap(),Ncolors,true);
  Epetra_MultiVector temp2(out_matrix.RangeMap(),Ncolors,false);

  /* Grab the color data */
  Epetra_MapColoring * col;
  Teuchos::RCP<Epetra_MapColoring> col_=colorer_->generateRowMapColoring();
  if(out_matrix.Comm().NumProc() == 1 ||
    input_graph_->Importer() == 0 )
    col=&(*col_);
  else {
    col= new Epetra_MapColoring(input_graph_->ColMap());
    col->Import(*col_,*input_graph_->Importer(),Insert);
  }
   
  /* Probing */
  for(int i=0;i<N;i++)
    temp1[(*col)[i]-1][i]=1.0;
  op.Apply(temp1,temp2);
  
  /* Matrix Fill*/
  out_matrix.FillComplete();
  int entries, *indices;
  double *values;
  for(int i=0;i<N;i++){
    out_matrix.ExtractMyRowView(i,entries,values,indices);
    for(int j=0;j<entries;j++)
      values[j]=temp2[(*col)[indices[j]]-1][i];
  }  

  /* Cleanup */
  if(out_matrix.Comm().NumProc() != 1 && input_graph_->Importer() != 0)
    delete col;
      
  return 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
Teuchos::RCP<Epetra_CrsMatrix> Prober::probe(const Epetra_Operator & op)
{
  Teuchos::RCP<Epetra_CrsMatrix> out_matrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*input_graph_));
  Teuchos::RCP<Epetra_CrsMatrix> null;
  int rv=probe(op,*out_matrix); 
  if(rv==0) return out_matrix;
  else return null;
}
////////////////////////////////////////////////////////////////////////////////

  
}//epetra
}//isorropia

  
#endif
