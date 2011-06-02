//============================================================================
// Name        : parallel_HK.cpp
// Author      : Arif Khan
// Email       : khan58@purdue.edu
// Copyright   : Purdue University
// Description : Parallel Hopkroft_Karp in C++/openMP, Ansi-style
//============================================================================

#include"Isorropia_EpetraMatcher.hpp"
#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#ifdef HAVE_EPETRA
//#ifdef HAVE_MPI
//#include <Epetra_MpiComm.h>
//#else
#include <Epetra_SerialComm.h>
//#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#ifdef HAVE_EPETRAEXT
#include <EpetraExt_CrsMatrixIn.h>
#endif

#endif

#include <Teuchos_CommandLineProcessor.hpp>

#include <unistd.h>

using namespace std;
//__sync_fetch_and_add(&a,1);

matcher::matcher(char * infile)
{
	int i;
	process_mtx_compressed(infile);
	finish=false;

	Node t;
	t.layer_num=-1;
	t.scanned=1;
	t.edgelist.clear();
	
	for(i=0;i<U;i++)
	{	
		LU.push_back(t);
		mateU.push_back(-1);
	}
	for(i=0;i<V;i++)
	{	
		LV.push_back(t);
		mateV.push_back(-1);
		vlayered.push_back(0);
	}
}

matcher::~matcher() {
	mateU.clear();
	mateV.clear();
	del_m.clear();
	LV.clear();
	LU.clear();
	vlayered.clear();
	graph.clear();
}

void matcher::process_mtx_compressed(char *fname)
{
	/*****************************************/
	int rc=0;
	#ifdef HAVE_EPETRAEXT
	 // bool verbose = false;
	  int localProc = 0;
	//   std::string *fstr;

	//#ifdef HAVE_MPI
	  /*int numProcs;
	  MPI_Init(&argc, &argv);
	  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
	  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	  const Epetra_MpiComm Comm(MPI_COMM_WORLD);
	  const Epetra_MpiComm Comm;*/
	  
	//#else
	  const Epetra_SerialComm Comm;
	//#endif

	  //string *inputFile = new string("a.mtx");
	  //const char *fname1 = inputFile->c_str();

	  // Read in the matrix market file and distribute its rows across the
	  // processes.
	  //
	  // This reader uses the default Epetra_Map for number of rows for the
	  // RowMap() and for the RangeMap().  For non-square matrices it uses
	  // the default Epetra_Map for the number of columns for the DomainMap(),
	  // otherwise it uses the RowMap().
	  //
	  // The maps can be specified with other versions of MMFtoCrsMatrix().

	  Epetra_CrsMatrix *matrixPtr;
	  rc = EpetraExt::MatrixMarketFileToCrsMatrix(fname, Comm, matrixPtr);
	  if (rc < 0){
		 if (localProc==0){
		   cout << "error reading input file" << std::endl << "FAIL" << std::endl;
		 }
		 exit(1);
	  }
	  else
	  		cout<<"Crs Matrix Created!!!"<<endl;
	#else
	  fail = 0;
	  if (localProc == 0){
		 std::cout << "Test not run because it requires EPETRA_EXT" << std::endl;
	  }
	#endif

	//#ifdef HAVE_MPI
	  //MPI_Finalize();
	//#endif

	rc=matrixPtr->ExtractCrsDataPointers(CRS_Pointers,CRS_Indices,CRS_Vals);
	if(rc==0)
		cout<<"Input Processing Done"<<endl;
	else
		cout<<"Input Processing Failed"<<endl;
		
	/*for(int i=0;i<matrixPtr->NumGlobalRows();i++)
	{
		for(int j=*(CRS_Pointers+i);j<*(CRS_Pointers+i+1);j++)
			cout<<*(CRS_Indices+j)<<" ";
		cout<<endl;
	}*/
	
	U=matrixPtr->NumGlobalRows();
	V=matrixPtr->NumGlobalCols();
	E=matrixPtr->NumGlobalNonzeros();
	cout<<"(U,V,E):"<<U<<","<<V<<","<<E<<endl;
}


int matcher::vlayer_clear()
{
	int i;
	
	#pragma omp parallel for
	for(i=0;i<V;i++)
		vlayered[i]=0;
		
	return 1;
}

int matcher::is_intersect(int k)
{
	unsigned int i;
	int flag=0;

	#pragma omp parallel for
	for(i=0;i<vlist.size();i++)
		if(mateV[vlist[i]]==-1)
			flag=1;
	
	if(flag==1)
		return 1;
	return 0;
}

void matcher::delete_matched_v()
{
	int i;
	
	#pragma omp parallel for
	for(i=0;i<V;i++)
		if(LV[i].layer_num==k_star)
			if(mateV[i]!=-1) // !=-1 means already has a mate,i.e., this node is matched...so remove it
				LV[i].scanned=1; // removing it and make it scanned=1 actually equivalent
}

void matcher::update_vlayered(int k)
{
	int i;
	if(vlist.size()>0)
		vlist.clear();
	#pragma omp parallel for schedule(dynamic, 126)
	for(i=0;i<V;i++)
	{		
		if(LV[i].layer_num==k)
		{	
			vlayered[i]=1;     // updating_valayered...
			#pragma omp critical
			vlist.push_back(i); //critical section
		}
	}
}

int matcher::construct_layered_graph()
{
	int k,t,flag,i,j;
	k=k_star=0;
	
	vlayer_clear();
	
	Node tmp;
	tmp.layer_num=-1;
	tmp.scanned=1;
	
	#pragma omp parallel for
	for(i=0;i<V;i++)
	{	
		LV[i].layer_num=-1;
		LV[i].scanned=1;
	}
	
	// Creating L0
	if(icm>1)
	{
		#pragma omp parallel for
		for(i=0;i<U;i++)     // if mateU[i]==-1 it means that it is not matched
		{
			if(mateU[i]==-1)
			{
				LU[i].layer_num=0;
				LU[i].scanned=0;
			}
			else
			{
				LU[i].layer_num=-1;
				LU[i].scanned=1;
			}
		}
	}
	
	//BFS layer construction
	k=k_star=0;
	int counting;
	while(true)  //change condition
	{
		nvlist.clear();
		counting=0;
		flag=0;
		//omp_set_num_threads(8);			
		if(k==0)                
		{
				if(icm==1)
				{
					#pragma omp parallel for private(j,t)
					for(i=0;i<U;i++)
					{
						LU[i].layer_num=0;
						LU[i].scanned=0;
						LU[i].edgelist.clear();
						
						for(t=*(CRS_Pointers+i);t<*(CRS_Pointers+i+1);t++)
						//for(t=graph[1][i];t<graph[1][i+1];t++)
						{
							//j=graph[0][t];
							j=*(CRS_Indices+t);						
							if(mateU[i]!=j && vlayered[j]==0) //there is a edge & not in Matching & not in Layers
							{	
								LV[j].layer_num=k+1;
								LV[j].scanned=0;
								LU[i].edgelist.push_back(j);
								flag=1;
							}	
						}
					}
				}
				else
				{
					#pragma omp parallel for private(j,t)
					for(i=0;i<U;i++)
					{
						if(LU[i].layer_num==0)
						{
							LU[i].edgelist.clear();
							for(t=*(CRS_Pointers+i);t<*(CRS_Pointers+i+1);t++)
							//for(t=graph[1][i];t<graph[1][i+1];t++)
							{
								//j=graph[0][t];
								j=*(CRS_Indices+t);		
								
								if(mateU[i]!=j && vlayered[j]==0) //there is a edge & not in Matching & not in Layers
								{	
									LV[j].layer_num=k+1;
									LV[j].scanned=0;
									LU[i].edgelist.push_back(j);
									flag=1;
								}	
							}
						}
					}
				}
		}
		else
		{
			flag=0;
			#pragma omp parallel for private(j,t)
			for(i=0;i<(signed)vlist.size();i++)
			{
				int id=mateV[vlist[i]];
				LU[id].edgelist.clear();
				for(t=*(CRS_Pointers+id);t<*(CRS_Pointers+id+1);t++)
				//for(t=graph[1][id];t<graph[1][id+1];t++)
				{
					//j=graph[0][t];
					j=*(CRS_Indices+t);
						
					if(mateU[id]!=j && vlayered[j]==0) //there is a edge & not in Matching & not in Layers
					{	
						LV[j].layer_num=k+1;
						LV[j].scanned=0;
						LU[id].edgelist.push_back(j);
						flag=1;
					}	
				}
			}
		}
		
		vlist.clear();
				
		//omp_set_num_threads(24);			
		update_vlayered(k+1);  // updating new inserted V nodes.*/
						
		if(flag==0)
		{
			finish=true;
			return k;
		}

		if(is_intersect(k+1)==1) //Lk+1 intersect V0 not empty
		{	
			k_star=k+1;
			break;
		}
		else  //build k+2.....
		{
			
			//omp_set_num_threads(8);				
			#pragma omp parallel for private(j)
			for(j=0;j<(signed)vlist.size();j++)
			{
				int id=mateV[vlist[j]];
				//cout<<"Construction..."<<endl;
				LU[id].layer_num=k+2;
				LU[id].scanned=0;
				LV[vlist[j]].edgelist.clear();
				LV[vlist[j]].edgelist.push_back(id);
			}
		}
		k=k+2;
	}
	vlayer_clear();
	return k_star;
}

int matcher::recursive_path_finder(int k, int p, vector<int>* path)
{
	int i,ind,res=0;
	path->push_back(p);
	
	if(k==k_star)
		return 1;

	if(k%2==0) // Layers where vertices are from set U
	{
		for(i=0;(unsigned)i<LU[p].edgelist.size();i++)
		{
			ind=LU[p].edgelist[i];
			
			res=0;
			#pragma omp critical
			if(LV[ind].scanned==0)              /// I am locking the whole L 
			{	
				LV[ind].scanned=1;
				res=1;
			}
				
			if(res==1 && recursive_path_finder(k+1,ind,path)==1)
					return 1;	
			
		}
	}
	else
		if(recursive_path_finder(k+1,LV[p].edgelist[0],path)==1)
			return 1;
				
	path->pop_back();
	return 0;
}

int matcher::iterative_path_finder(int k, int p, vector<int>* path)
{
	int i,ind,res=0,cur_p,cur_k,last_k;
	stack <vector<int> > st;
	vector<int> temp;
		
	temp.clear();
	temp.push_back(k);
	temp.push_back(p);
	st.push(temp);
	last_k=-1;        
	
	while(!st.empty())
	{
		temp=st.top();
		st.pop();
		cur_k=temp[0];
		cur_p=temp[1];
		
		if(cur_k<last_k)
		{	
			path->pop_back(); 
			path->pop_back();
		}
		last_k=cur_k;
		
		if(cur_k==k_star)
		{	
			res=0;
			#pragma omp critical
			if(LV[cur_p].scanned==0)
			{
				LV[cur_p].scanned=1;
				res=1;
			}
			if(res==1)
			{
				path->push_back(cur_p);
				return 1;
			}
			
		}
		
		if((cur_k%2)==0)
		{	
			path->push_back(cur_p);

			for(i=0;(unsigned)i<LU[cur_p].edgelist.size();i++)
			{
				ind=LU[cur_p].edgelist[i];
				vector<int> t;
				t.push_back(cur_k+1);
				t.push_back(ind);
				st.push(t);
				
			}
		}
		else
		{
			res=0;
			#pragma omp critical
			if(LV[cur_p].scanned==0)
			{
				LV[cur_p].scanned=1;
				res=1;
			}
			if(res==1)
			{
				path->push_back(cur_p);
				ind=LV[cur_p].edgelist[0];
				vector<int> t;
				t.push_back(cur_k+1);
				t.push_back(ind);
				st.push(t);
			}
		}	
	}
	
	return 0;
}

void matcher::find_set_del_M()
{
	int i;
	delete_matched_v();
	del_m.clear();
	
	//omp_set_num_threads(2);
	#pragma omp parallel for
	for(i=0;i<U;i++)
	{
		if(LU[i].layer_num==0)
		{	
			LU[i].scanned=1;
			vector<int>* path=new vector<int>;
			if(k_star>2000)
			{	
				if(iterative_path_finder(0,i,path)==1)
					#pragma omp critical
					del_m.push_back(*path);
			}
			else
			{	
				if(recursive_path_finder(0,i,path)==1)
					#pragma omp critical
					del_m.push_back(*path);
			}
			delete path;
		}
	}
}

bool matcher::remove_edge(int s, int t)
{
	mateV[s]=-1;
	mateU[t]=-1;
	return true;
}

void matcher::add_edge(int s, int t)
{
	mateU[s]=t;
	mateV[t]=s;
}

int matcher::augment_matching()
{
	int i,j,count;

	count=del_m.size();
	#pragma omp parallel for private(j)
	for(i=0;i<count;i++)
	{
		for(j=1;(unsigned)j<del_m[i].size()-2;j=j+2)
		{
			if(!remove_edge(del_m[i][j],del_m[i][j+1]))
				cout<<"Error: Remove Edge...not found"<<endl;
		}
			
		for(j=0;(unsigned)j<del_m[i].size()-1;j=j+2)
			add_edge(del_m[i][j],del_m[i][j+1]);
	}
	del_m.clear();
	return count;
}

bool matcher::U0_empty()
{
	int i,flag=0;
	
	#pragma omp parallel for
	for(i=0;i<U;i++)
		if(mateU[i]==-1)
			flag=1;
	if(flag==1)
		return false;
	else return true;
}

vector<int> matcher::get_matching()
{
	int totc=0,count;
	time_t t1,t2,t3,t4,t_st,t_end;
	time(&t_st);
	icm=0;
	while(true)
	{
		icm++;
		time(&t1);
		construct_layered_graph();
		//cout<<"bfs"<<endl;
		time(&t2);
		if(finish || U0_empty())
			break;
		find_set_del_M();
		//cout<<"dfs"<<endl;
		time(&t3);
		count=augment_matching();
		time(&t4);
		totc+=count;
		cout<<"["<<icm<<"] Layers="<<k_star+1<<" BFS="<<difftime(t2,t1)<<" DFS="<<difftime(t3,t2)<<" Time="<<difftime(t4,t1)<<" matched="<<count<<" size="<<totc<<endl;
	}
	time(&t_end);
	
	cout<<"Total time is less than "<<(t_end-t_st)+1<<" seconds"<<endl;
	/// Returning the matched edges in int*
	return mateU;
}

