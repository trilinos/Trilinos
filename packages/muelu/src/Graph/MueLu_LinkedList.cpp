#include "MueLu_LinkedList.hpp"

namespace MueLu {

  LinkedList::LinkedList() : nodeHead(NULL), nodeTail(NULL) { }

  LinkedList::~LinkedList() {
    while (nodeHead != NULL)
      DeleteHead();
  }

  bool LinkedList::IsEmpty() {
    return nodeHead == NULL;
  }

  void LinkedList::Add(int iNode) {
    MueLu_Node *newNode = new MueLu_Node;
    newNode->nodeId = iNode;
    newNode->next = NULL;
    if (nodeHead == NULL) {
      nodeHead = newNode;
      nodeTail = newNode;
    } else {
      nodeTail->next = newNode;
      nodeTail = newNode;
    }
  }

  int LinkedList::Pop() { // get head and remove first node
    if (IsEmpty()) throw(1);

    int iNode = nodeHead->nodeId;
    DeleteHead();
    return iNode;
  }

  void LinkedList::DeleteHead() {
    if (IsEmpty()) throw(1);
      
    MueLu_Node *newNode = nodeHead;
    nodeHead = newNode->next;
    delete newNode;
  }

}

//TODO: nodeTail unused -> remove?
