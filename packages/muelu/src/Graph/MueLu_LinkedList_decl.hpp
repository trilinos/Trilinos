#ifndef MUELU_LINKEDLIST_DECL_HPP
#define MUELU_LINKEDLIST_DECL_HPP

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION // Otherwise, class will be declared twice because _decl.hpp file also have the class definition (FIXME)

#include "MueLu_ConfigDefs.hpp"

/* ------------------------------------------------------------------------- */
/* linked list structures from ML for holding free node information          */
/* ------------------------------------------------------------------------- */

namespace MueLu {

  typedef struct MueLu_Node_Struct
  {
    int nodeId;
    struct MueLu_Node_Struct *next;
  } MueLu_Node;
  
  class LinkedList {

  public:    
    LinkedList() : nodeHead(NULL), nodeTail(NULL) ;

    ~LinkedList() ;

    bool IsEmpty() ;

    void Add(int iNode) ;

    int Pop() ;

  private:
    MueLu_Node *nodeHead;
    MueLu_Node *nodeTail;

    void DeleteHead() ;

  };

}

//TODO: nodeTail unused -> remove?

#endif // HAVE_MUELU_EXPLICIT_INSTANTIATION
#endif // MUELU_LINKEDLIST_DECL_HPP
