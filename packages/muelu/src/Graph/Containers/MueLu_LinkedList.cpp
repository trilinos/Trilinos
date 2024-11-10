// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <stddef.h>  // for NULL
#include "MueLu_LinkedList.hpp"

namespace MueLu {

LinkedList::LinkedList()
  : nodeHead(NULL)
  , nodeTail(NULL) {}

LinkedList::~LinkedList() {
  while (nodeHead != NULL)
    DeleteHead();
}

bool LinkedList::IsEmpty() {
  return nodeHead == NULL;
}

void LinkedList::Add(int iNode) {
  MueLu_Node *newNode = new MueLu_Node;
  newNode->nodeId     = iNode;
  newNode->next       = NULL;
  if (nodeHead == NULL) {
    nodeHead = newNode;
    nodeTail = newNode;
  } else {
    nodeTail->next = newNode;
    nodeTail       = newNode;
  }
}

int LinkedList::Pop() {  // get head and remove first node
  if (IsEmpty()) return -1;

  int iNode = nodeHead->nodeId;
  DeleteHead();
  return iNode;
}

void LinkedList::DeleteHead() {
  if (IsEmpty()) return;

  MueLu_Node *newNode = nodeHead;
  nodeHead            = newNode->next;
  delete newNode;
}

}  // namespace MueLu

// TODO: nodeTail unused -> remove?
