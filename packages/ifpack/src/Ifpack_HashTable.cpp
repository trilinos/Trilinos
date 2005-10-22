#include "Ifpack_HashTable.h"

NodeArray* Ifpack_HashTable::Data_ = 0;
NodeArray* Ifpack_HashTable::FirstData_ = 0;
Node** Ifpack_HashTable::Container_ = 0;
int Ifpack_HashTable::ContainerSize_ = 0;
