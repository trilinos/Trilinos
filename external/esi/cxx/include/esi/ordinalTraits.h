#ifndef __ESI_ordinalTraits_h
#define __ESI_ordinalTraits_h

namespace esi {

template<class T>
struct ordinalTraits {
   static inline const char* name() {
     cout << "esi::ordinalTraits: unsupported ordinal type." << endl; abort(); 
     return(NULL);
   };
};

template<>
struct ordinalTraits<int> {
   typedef int ordinal_type;
   static inline const char* name() { return("int"); };
};

template<>
struct ordinalTraits<long> {
   typedef long ordinal_type;
   static inline const char* name() {return("int");};
};


#ifdef ESI_NO_TYPENAME_KEYWORD
#define TYPENAME
#else
#define TYPENAME typename
#endif

};     // esi namespace
#endif

