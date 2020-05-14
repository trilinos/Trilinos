#ifndef __Zoltan2_ice_util_h__
#define __Zoltan2_ice_util_h__

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

//function to read the boundary information from testing files
//File format is:
//two-word description
//#edges
//v1 v2
//v3 v4
//...
//
//where v1, v2, v3 and v4 are vertex identifiers.
template<typename gno_t>
void read_boundary_file(const char *filename, size_t& num_edges, 
                        Teuchos::Array<gno_t> & boundary_edges){
  std::ifstream fin(filename);
  if(!fin){
    std::cout<<"Unable to open file "<<filename<<"\n";
    exit(0);
  }
  //the first two words are a description of the test,
  //so we need to skip them
  std::string throwaway;
  fin>>throwaway>>throwaway;
  //this line gives the number of boundary edges
  int arrlength;
  fin>>arrlength;
  
  num_edges = 2*arrlength;
  gno_t a, b;
  //get the list of boundary edges instead of flags representing articulation points,
  //the integration code takes a list of boundary edges.
  //int count = 0;
  while(fin>>a>>b){
    boundary_edges.push_back(a-1);
    boundary_edges.push_back(b-1);
  }
}

//function to read the grounding information from the testing files
//File format is:
//#vertices
//<grounding for vtx 0>
//<grounding for vtx 1>
//...
//<grounding for vtx n>
//
//grounding can be a float or an integer, a zero value means floating,
//nonzero means touching the ground.
void read_grounded_file(const char* filename, size_t& nVtx, 
                        Teuchos::Array<int>& grounded_flags){
  std::ifstream fin(filename);
  if(!fin){
    std::cout<<"Unable to open "<<filename<<"\n";
    exit(0);
  }
  //the first number is the number of vertices
  fin>>nVtx;
  //the rest of the numbers are basal friction data
  for(size_t i = 0; i < nVtx; i++){
    float gnd;
    fin>>gnd;
    grounded_flags.push_back((gnd > 0.0));
  }
}

#endif
