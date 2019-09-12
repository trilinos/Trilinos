// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#define PR(a)  " " << #a << " = " << a << " "
#define PRINT(a) do { std::cout << #a << " = " << a ; } while(0)
#define PRINT2(a,b) do { std::cout << #a << " = " << a << " " << #b << " = " << b ; } while(0)
#define PRINTLN(a) do { std::cout << #a << " = " << a << std::endl; } while(0)
#define PRINTLN2(a,b) do { std::cout << #a << " = " << a << " " << #b << " = " << b << std::endl; } while(0)
#define PLN() do { std::cout <<  std::endl; } while(0)

#define DPRINT(a) do { if (debug_print) std::cout << #a << " = " << a ; } while(0)
#define DPRINT2(a,b) do { if (debug_print) std::cout << #a << " = " << a << " " << #b << " = " << b ; } while(0)
#define DPRINTLN(a) do { if (debug_print) std::cout << #a << " = " << a << std::endl; } while(0)
#define DPRINTLN2(a,b) do { if (debug_print) std::cout << #a << " = " << a << " " << #b << " = " << b << std::endl; } while(0)
#define DPLN() do { if (debug_print) std::cout <<  std::endl; } while(0)
