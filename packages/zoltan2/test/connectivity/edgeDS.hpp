#ifndef __edgeDS_hpp__
#define __edgeDS_hpp__
#include<unordered_map>
#include<assert.h>
#include<unordered_set>
#include<iostream>
#include<stdlib.h>

class edge {
public:
  //ensure the smaller node goes into u,
  //so all edges are represented the same way
  edge(int a, int b){
    if(a < b){
      u = a;
      v = b;
    } else {
      v = a;
      u = b;
    }
  }

  int getU() const {return u;}
  int getV() const {return v;}
  void validate(void){
    if(u > v){
      std::cout<<u<<" is not less than or equal to "<<v<<"\n";
    }
    assert(u<=v);
  }
  
  bool operator==(const edge& rhs) const{
    return (u == rhs.u) && (v == rhs.v);
  }
private:
  //the smaller endpoint
  int u;
  //the larger endpoint
  int v;
};

class edge_hash {
public: 
  int operator()(const edge& e) const {
    return e.getU()*e.getV() + e.getV();
  }
};

template<typename T>
using edge_map = std::unordered_map<edge,T,edge_hash>;

typedef std::unordered_set<edge,edge_hash> edge_set;

#endif
