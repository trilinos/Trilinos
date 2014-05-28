/**   ------------------------------------------------------------
 *    Copyright 1998-2007 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

// ---------------------------------------------------------------------
// Author:     H. Carter Edwards
//
// Purpose:    Associative container for allocated data objects
// ---------------------------------------------------------------------
// Acknowledgements:
//
//   Most all of the algorithms in this class were obtained from
// the Hewlett-Packard source for the Standard Template Library,
// thus the inclusion of Hewlett-Packard's copyright notice.
// Some minor modifications were obtained from Silicon Graphics'
// Standard Template Library source.
// ---------------------------------------------------------------------
/*
 * Copyright (c) 1996,1997
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */
/*
Red-black tree class, designed for use in implementing STL
associative containers (set, multiset, map, and multimap). The
insertion and deletion algorithms are based on those in Cormen,
Leiserson, and Rivest, Introduction to Algorithms (MIT Press, 1990),
except that

(1) the header cell is maintained with links not only to the root
but also to the leftmost node of the tree, to enable constant time
begin(), and to the rightmost node of the tree, to enable linear time
performance when used with the generic set algorithms (set_union,
etc.);

(2) when a node being deleted has two children its successor node is
relinked into its place, rather than copied, so that the only
iterators invalidated are those referring to the deleted node.
*/
// ---------------------------------------------------------------------

// The header

#include <stk_util/diag/Mapv.hpp>

#include <stdexcept>
#if defined(SIERRA_IA64_OPTIMIZER_WARN)
#include <stk_util/diag/Env.hpp>
#endif
#include <string>

namespace sierra {

// ---------------------------------------------------------------------

inline
MapvNodeBase * MapvBase::minimum( MapvNodeBase * x )
  { while ( x->left ) x = x->left ; return x ; }

inline
MapvNodeBase * MapvBase::maximum( MapvNodeBase * x )
  { while ( x->right ) x = x->right ; return x ; }

// ---------------------------------------------------------------------

inline void MapvBase::rotate_left( MapvNodeBase * x )
{
  MapvNodeBase * y = x->right ;

  x->right = y->left ;
  if ( y->left ) y->left->parent = x ;
  y->parent = x->parent ;

  if      ( x == nRoot() )         root(y);
  else if ( x == x->parent->left ) x->parent->left  = y ;
  else                             x->parent->right = y ;

  y->left   = x ;
  x->parent = y ;
}

inline void MapvBase::rotate_right( MapvNodeBase * x )
{
  MapvNodeBase * y = x->left ;

  x->left = y->right ;
  if ( y->right ) y->right->parent = x;
  y->parent = x->parent ;

  if      ( x == nRoot() )          root(y);
  else if ( x == x->parent->right ) x->parent->right = y ;
  else                              x->parent->left  = y ;

  y->right  = x;
  x->parent = y;
}

// ---------------------------------------------------------------------

void MapvBase::insert( MapvNodeBase * y , MapvNodeBase * z , bool z_lt_y )
{
  z->remove_from_container();

  if ( y == nEnd() ) { // First node inserted
    root(z);
    leftmost(z);
    rightmost(z);
    z->parent = header() ; // header is 'super-root'
  }
  else {
    if ( z_lt_y ) {
      y->left = z ;
      // maintain leftmost() pointing to minimum node
      if ( y == leftmost() ) leftmost(z);
    }
    else {
      y->right = z;
      // maintain rightmost() pointing to maximum node
      if ( y == rightmost() ) rightmost(z);
    }
    z->parent = y ;
  }
  z->left  = 0 ;
  z->right = 0 ;
  z->color = red ;
  ++Count ;

  // -------------------------------------------------------------------
  // Rebalance, 'y' and 'z' are reused as a local variable

  while ( z != nRoot() && z->parent->color == red ) {
    if ( z->parent == z->parent->parent->left ) {
      y = z->parent->parent->right ;
      if ( y && y->color == red ) {
	z->parent->color         = black;
	y->color                 = black;
	z->parent->parent->color = red;
	z = z->parent->parent ;
      }
      else {
	if ( z == z->parent->right ) {
	    z = z->parent ;
	    rotate_left(z);
	}
	z->parent->color         = black;
	z->parent->parent->color = red;
	rotate_right( z->parent->parent );
      }
    }
    else {
      y = z->parent->parent->left ;
      if ( y && y->color == red ) {
	z->parent->color         = black;
	y->color                 = black;
	z->parent->parent->color = red;
	z = z->parent->parent ;
      }
      else {
	if ( z == z->parent->left ) {
	    z = z->parent ;
	    rotate_right(z);
	}
	z->parent->color         = black;
	z->parent->parent->color = red;
	rotate_left(z->parent->parent);
      }
    }
  }
  nRoot()->color = black;
}

// ---------------------------------------------------------------------

void MapvBase::remove( MapvNodeBase * node )
{
  static const char method_name[] = "MapvBase::remove" ;

  if ( container(node) != this ) {
    std::string msg(method_name);
    msg.append(" given object not in this container");
    throw std::invalid_argument( msg );
  }

  if ( 1 == Count ) { // The last node ?

    if ( node != leftmost() || node != rightmost() || node != nRoot() ) {
      std::string msg(method_name);
      msg.append(" internal data structure corrupted" );
      throw std::runtime_error( msg );
    }

    leftmost( nREnd() );
    rightmost( nEnd() );
    root(0);
    Count = 0 ;
    header()->color = red ;
    node->left = node->right = node->parent = 0 ; node->color = 0 ;
    return ;
  }

  MapvNodeBase * z = node ;
  MapvNodeBase * y = node ;
  MapvNodeBase * x = 0 ;
  MapvNodeBase * x_parent = 0 ;

  // Ready to remove

  if ( y->left == 0 ) {       // z has at most one non-null child. y == z
    x = y->right ;            // x might be null
  }
  else if ( y->right == 0 ) { // z has exactly one non-null child. y == z
    x = y->left ;             // z is not null
  }
  else {                      // z has two non-null children.
     y = y->right ;           // Set y to z's successor.
     while ( y->left ) y = y->left ;
     x = y->right ;           // x might be null
  }

  if ( y != z ) { // relink y in place of z. y is z's successor
    z->left->parent = y ;
    y->left = z->left ;
    if ( y != z->right ) {
      x_parent = y->parent ;
      if ( x ) x->parent = x_parent ;
      y->parent->left = x;   // y must be a left child
      y->right = z->right;
      z->right->parent = y;
    } else {
      x_parent = y;  // needed in case x == 0
    }
    if ( nRoot() == z) {
      root(y);
    }
    else if ( z->parent->left == z) {
      z->parent->left = y;
    }
    else {
      z->parent->right = y;
    }
    y->parent = z->parent;
    { int c = y->color; y->color = z->color; z->color = c ; }
    y = z;
    // y points to node to be actually deleted
  }
  else {  // y == z
    x_parent = y->parent ;
    if ( x ) x->parent = x_parent ; // possibly x == 0
    if ( nRoot() == z) {
      root(x);
    }
    else if ( z->parent->left == z ) {
      z->parent->left = x;
    }
    else {
      z->parent->right = x;
    }
    if ( leftmost() == z )  {
      if ( z->right == 0 ) { // z->left must be null also
	// makes leftmost() == nEnd() if z == nRoot()
	leftmost( z->parent );
      }
      else {
	leftmost( minimum(x) );
      }
    }
    if ( rightmost() == z )  {
      if ( z->left == 0 ) { // z->right must be null also
	// makes rightmost() == nEnd() if z == nRoot()
	rightmost( z->parent );
      }
      else { // x == z->left
	rightmost( maximum(x) );
      }
    }
  }
  if ( y->color != red ) {
    while ( x != nRoot() && ( x == 0 || x->color == black ) ) {
      if ( x == x_parent->left ) {
	MapvNodeBase * w = x_parent->right ;
	if ( w->color == red ) {
	  w->color        = black;
	  x_parent->color = red;
	  rotate_left(x_parent);
	  w = x_parent->right ;
	}
	if ((w->left  == 0 || w->left->color  == black) &&
	    (w->right == 0 || w->right->color == black)) {
	  w->color = red ;
	  x = x_parent ;
	  x_parent = x_parent->parent ;
	}
	else {
	  if (w->right == 0 || w->right->color == black) {
	      if ( w->left ) w->left->color = black;
	      w->color = red;
	      rotate_right(w);
	      w = x_parent->right ;
	  }
	  w->color = x_parent->color ;
	  x_parent->color = black;
	  if ( w->right ) w->right->color = black;
	  rotate_left(x_parent);
	  break;
	}
      }
      else {  // same as then clause with "right" and "left" exchanged
	MapvNodeBase * w = x_parent->left ;
	if ( w->color == red ) {
	  w->color = black;
	  x_parent->color = red;
	  rotate_right(x_parent);
	  w = x_parent->left ;
	}
	if ((w->right == 0 || w->right->color == black) &&
	    (w->left  == 0 || w->left->color  == black)) {
	  w->color = red;
	  x = x_parent ;
	  x_parent = x_parent->parent ;
	}
	else {
	  if ( w->left == 0 || w->left->color == black ) {
	    if ( w->right ) w->right->color = black;
	    w->color = red;
	    rotate_left(w);
	    w = x_parent->left ;
	  }
	  w->color = x_parent->color ;
	  x_parent->color = black;
	  if ( w->left ) w->left->color = black;
	  rotate_right(x_parent);
	  break;
	}
      }
    }
    if ( x ) x->color = black;
  }

  y->left = y->right = y->parent = 0 ; y->color = 0 ;

  --Count ; // Decrement the tree's count
}

// ---------------------------------------------------------------------
// A reverse communicating method for deleting all entries

MapvNodeBase * MapvBase::unbalancing_removal( MapvNodeBase ** n )
{
  MapvNodeBase * t = *n ;

  while ( t != header() && t->parent ) {
    if      ( t->left  ) { t = t->left ; }
    else if ( t->right ) { t = t->right ; }
    else { // Move to parent and remove this leaf
      *n = t->parent ; t->parent = 0 ;
      if ( (*n)->left == t ) (*n)->left  = 0 ;
      else                   (*n)->right = 0 ;
    }
  }

  if ( t == header() ) {

    header()->parent = 0 ;
    header()->left   = 0 ;
    header()->right  = 0 ;
    header()->color  = red ;  /* Color the header node red */

    Count = 0 ;

    left_end.parent = 0 ;
    left_end.left   = 0 ;
    left_end.right  = 0 ;
    left_end.color  = black ;

    right_end.parent = 0 ;
    right_end.left   = 0 ;
    right_end.right  = 0 ;
    right_end.color  = black ;

    leftmost(  header()->left  = nREnd() ); // left end of the tree
    rightmost( header()->right = nEnd() );  // right end of the tree

    t = 0 ;
  }

  return t ;
}

// ---------------------------------------------------------------------
void MapvBase::WarnOptimize()
{
#if defined(SIERRA_IA64_OPTIMIZER_WARN) && defined(NDEBUG)
    static bool warn_once=true;
    int my_proc=Env::parallel_rank();

    if ( warn_once && my_proc==0){
        warn_once = false;
        std::cerr << "Optimizing previous versions of the intel compiler "
                  << "caused errors in Mapv.\n"
                  << "Results may be suspect.\n";
    }
#endif
}

// ---------------------------------------------------------------------
// Virtual destructor

MapvBase::~MapvBase()
{
  if ( Count || nRoot() != 0 ) {
    std::string msg("MapvBase destructor, container is not empty");
    throw std::logic_error( msg );
  }
}

}
