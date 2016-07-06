#include <assert.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <Ioss_Sort.h>

namespace {
  template <typename INT> bool verify_sorted(std::vector<INT> v)
  {
    for (size_t i=1; i < v.size(); i++) {
      if (v[i-1] > v[i]) {
	std::cerr << "Unsorted at position " << i << "\n";
	return false;
      }
    }
    return true;
  }
}

int main() {
  std::random_device rd;
  std::mt19937_64 rng(rd());

  const int sawtooth = 1;
  const int do_rand = 2;
  const int stagger = 3;
  const int plateau = 4;
  const int shuffle = 5;

  for (size_t n : { 100, 1023, 1024, 1025, (2<<16)-1, 2<<16, (2<<16)+1}) {
    std::cerr << "\nSize: " << n << ": ";
    for (size_t m = 1; m < 2*n; m *= 2) {
      std::cerr << m;
      for (auto dist : { sawtooth, do_rand, stagger, plateau, shuffle }) {
	std::cerr << ".";
	std::vector<int64_t> x(n);
	
	size_t i=0,j=0,k=1;
	switch (dist) {
	case sawtooth:
	  for ( ; i < n; i++) {
	    x[i] =i%m;
	  }
	  break;
	case do_rand:
	  for ( ; i < n; i++) {
	    x[i] = rng() % m;
	  }
	  break;
	case stagger:
	  for ( ; i < n; i++) {
	    x[i] = (i*m + i) % n;
	  }
	  break;
	case plateau:
	  for ( ; i < n; i++) {
	    x[i] = std::min(i, m);
	  }
	  break;
	case shuffle:
	  for ( ; i < n; i++) {
	    x[i] = rng()%m ? (j+=2): (k+=2);
	  }
	  break;
	}

	Ioss::qsort(x); // Copy of x
	assert(verify_sorted(x));
	  
	std::reverse(x.begin(), x.end()); // Reversed
	Ioss::qsort(x);
	assert(verify_sorted(x));
	  
	std::reverse(&x[0], &x[n/2]); // Front half reversed
	Ioss::qsort(x);
	assert(verify_sorted(x));

	std::reverse(&x[n/2], &x[n]); // Back half reversed
	Ioss::qsort(x);
	assert(verify_sorted(x));

	Ioss::qsort(x); // Already sorted
	assert(verify_sorted(x));

	for (size_t p=0; p < n; p++) {
	  x[p] += p%5;
	}
	Ioss::qsort(x); // Dithered
	assert(verify_sorted(x));
      }
    }
  }
  std::cerr << "\nDone\n";
}
