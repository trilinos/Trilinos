
#ifndef Epetra_HashTable_H_
#define Epetra_HashTable_H_

class Epetra_HashTable
{
  struct Node
  {
     int Key;
     int Value;
     Node * Ptr;

     Node( const int key = 0, const int value = 0, Node * ptr = 0 )
     : Key(key), Value(value), Ptr(ptr) {}
  };

  Node ** Container_;
  int Size_;
  unsigned int Seed_;

  int Func( const int key ) { return (Seed_ ^ key)%Size_; }
     
 public:

  Epetra_HashTable( const int size, const unsigned int seed = (2654435761U) )
  : Size_(size), Seed_(seed)
  {
    Container_ = new Node * [size];
    for( int i = 0; i < size; ++i ) Container_[i] = 0;
  }

  Epetra_HashTable( const Epetra_HashTable & obj )
  : Size_(obj.Size_),
    Seed_(obj.Seed_)
  {
    Container_ = new Node * [Size_];
    for( int i = 0; i < Size_; ++i ) Container_[i] = 0;
    {for( int i = 0; i < Size_; ++i )
    {
      Node * ptr = obj.Container_[i];
      while( ptr ) { Add( ptr->Key, ptr->Value ); ptr = ptr->Ptr; }
    }}
  }

  ~Epetra_HashTable()
  {
    Node * ptr1;
    Node * ptr2;
    for( int i = 0; i < Size_; ++i )
    {
      ptr1 = Container_[i];
      while( ptr1 ) { ptr2 = ptr1; ptr1 = ptr1->Ptr; delete ptr2; }
    }
  }

  void Add( const int key, const int value )
  {
    int v = Func(key);
    Node * n1 = Container_[v];
    Container_[v] = new Node(key,value,n1);
  }

  int Get( const int key )
  {
    Node * n = Container_[ Func(key) ];
    while( n && (n->Key != key) ) n = n->Ptr;
    if( n ) return n->Value;
    else    return -1;
  }

};

#endif
