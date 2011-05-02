#ifndef MUELU_LINKEDLIST_HPP
#define MUELU_LINKEDLIST_HPP

/* ------------------------------------------------------------------------- */
/* linked list structures from ML for holding free node information          */
/* ------------------------------------------------------------------------- */

namespace MueLu {

// using Teuchos::ArrayView;
// using Teuchos::ArrayRCP;

  typedef struct MueLu_Node_Struct
  {
    int nodeId;
    struct MueLu_Node_Struct *next;
  } MueLu_Node;
  
  class LinkedList {

  public:    
    LinkedList() {
      MueLu_Node *newNode = new MueLu_Node;      
      newNode->nodeId = 0;
      nodeHead = newNode;
      nodeTail = newNode;
      newNode->next = NULL;
    }

    ~LinkedList() {
      MueLu_Node *newNode = NULL;
      while ( nodeHead != NULL )
        {
          newNode = nodeHead;
          nodeHead = newNode->next;
          delete newNode;
        }
    }

    bool IsEmpty() {
      return nodeHead == NULL;
    }

    void Add(int iNode) {
      MueLu_Node *newNode = new MueLu_Node;
      newNode->nodeId = iNode;
      newNode->next = NULL;
      if ( nodeHead == NULL ) {
          nodeHead = newNode;
          nodeTail = newNode;
        } else {
          nodeTail->next = newNode;
          nodeTail = newNode;
        }
    }

    int Pop() { // get head and remove first node
      if (IsEmpty()) throw(1);

      MueLu_Node *newNode = nodeHead;
      int iNode = newNode->nodeId;
      nodeHead = newNode->next;
      delete newNode;
      return iNode;
    }

  private:
    MueLu_Node *nodeHead;
    MueLu_Node *nodeTail;

  };

}
#endif
//TODO: nodeTail unused -> remove?
