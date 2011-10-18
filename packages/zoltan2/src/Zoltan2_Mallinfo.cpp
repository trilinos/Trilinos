// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// A primitive utility for measuring memory usage at runtime.
//
// This utility is only built if -DZOLTAN2_ENABLE_MALLINFO:BOOL=ON
// is included in the cmake configuration command.
//

#include <malloc.h>  // for mallinfo - which isn't portable
#include <map> 
#include <string> 
#include <algorithm> 
#include <iostream> 

using namespace std;

namespace Zoltan2 {
// Map from a label to before/after memory use counts
static std::map<std::string, std::pair<int, int> > mem;

int getAllocatedMemory()
{
  struct mallinfo minfo = mallinfo();
  return minfo.arena;
}

// mallocCount looks at the current allocated heap memory and
// associates it with the argument label.  If this is the first
// occurence of the label, it's assumed to be a "before" count,
// and if it's the second or later, it's assumed to be the "after" 
// count.

int mallocCount(std::string label) {
  typedef std::map<std::string, std::pair<int, int> >::iterator it_t;
  struct mallinfo minfo = mallinfo();
  int memcount = minfo.arena;

  it_t curr = mem.find(label);
  if (curr == mem.end()){
    std::pair<int, int> results(memcount, 0);
    mem[label] = results;
  }
  else{
    curr->second.second = memcount;
  }
  return memcount;
}

// If a label is found, return the memory count
int getMallocCount(std::string label){
  int sum=0;
  typedef std::map<std::string, std::pair<int, int> >::iterator it_t;
  it_t curr = mem.find(label);
  if (curr != mem.end()){
    int before = curr->second.first;
    int after = curr->second.second;
    if (after > 0)
      sum = after - before; 
  }
  return sum;
}

// print out all labels and their memory counts, in
// no particular order

void printMallocCount()
{
  typedef std::map<std::string, std::pair<int, int> >::iterator it_t;
  for (it_t curr = mem.begin(); curr != mem.end(); ++curr)
  {
    int before = curr->second.first;
    int after = curr->second.second;
    if (after != 0)
      std::cout << curr->first << ": " << after-before << std::endl;
  }
  std::cout << std::endl;
}

void eraseMallocCount()
{
  mem.clear();
}
} // namespace
