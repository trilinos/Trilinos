// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LINKEDLIST_HPP
#define MUELU_LINKEDLIST_HPP

/* ------------------------------------------------------------------------- */
/* linked list structures from ML for holding free node information          */
/* ------------------------------------------------------------------------- */

namespace MueLu {

typedef struct MueLu_Node_Struct {
  int nodeId;
  struct MueLu_Node_Struct *next;
} MueLu_Node;

class LinkedList {
 public:
  LinkedList();

  ~LinkedList();

  bool IsEmpty();

  void Add(int iNode);

  int Pop();

 private:
  MueLu_Node *nodeHead;
  MueLu_Node *nodeTail;

  void DeleteHead();
};

}  // namespace MueLu

#endif  // MUELU_LINKEDLIST_HPP
