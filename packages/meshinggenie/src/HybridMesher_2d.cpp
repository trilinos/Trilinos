
//@HEADER
// ************************************************************************
// 
//               MeshingGenie: Fracture Meshing Services Package 
//                 Copyright 2011 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

// R5.0

#include "HybridMesher_2d.h"

using namespace std;

double cross_product(double x1, double y1, double x2, double y2)
{
  return x1*y2 - x2*y1;
}

void on_edge(std::vector<double> _x, std::vector<double> _y,
	     std::vector<size_t> &OnVorBound, size_t num_voronoi,
	     double ex1, double ey1, double ex2, double ey2)
{
  double ex_change = ex2-ex1, ey_change = ey2-ey1;
  size_t n;
  for (n=0; n<num_voronoi; ++n) {
    if (abs(cross_product(_x[n]-ex1, _y[n]-ey1, ex_change, ey_change)) <
	0.00000000000001 &&
	(ex1<=_x[n]&&_x[n]<=ex2 || ex2<=_x[n]&&_x[n]<=ex1) &&
	(ey1<=_y[n]&&_y[n]<=ey2 || ey2<=_y[n]&&_y[n]<=ey1)) {
      OnVorBound.push_back(n);
    }
  }
  
  size_t num_vor = OnVorBound.size()-1;
  for (n=num_vor; n>0; --n) {
    size_t swaps = 0;
    for (size_t m=0; m<n; ++m) {
      if ((_x[OnVorBound[m+1]]-_x[OnVorBound[m]])*ex_change<0 ||
	  (_y[OnVorBound[m+1]]-_y[OnVorBound[m]])*ey_change<0) {
	size_t tmp = OnVorBound[m];
	OnVorBound[m] = OnVorBound[m+1];
	OnVorBound[m+1] = tmp;
	++swaps;
      }
    }
    if (!swaps) {
      break;
    }
  }
}

void on_bound(std::vector<double> _x, std::vector<double> _y, size_t num_voronoi,
	      std::vector<double> _Holes, std::vector<size_t> &OnBound)
{
    std::vector<size_t> OnSegBound;
    size_t num_bound;
    size_t j;
    on_edge(_x,_y,OnSegBound,num_voronoi,_Holes[0],_Holes[1],_Holes[2],_Holes[3]);
    num_bound = OnSegBound.size();
    for (j=0; j<num_bound; ++j) {
      OnBound.push_back(OnSegBound[j]);
    }

    size_t num_segs = _Holes.size()-2;
    for (size_t i=2; i<num_segs; i+=2) {
      OnSegBound.clear();
      on_edge(_x, _y, OnSegBound, num_voronoi, _Holes[i], _Holes[i+1], _Holes[i+2],
	      _Holes[i+3]);
      num_bound = OnSegBound.size();
      for (j=1; j<num_bound; ++j) {
	OnBound.push_back(OnSegBound[j]);
      }
    }

    OnSegBound.clear();
    on_edge(_x, _y, OnSegBound, num_voronoi, _Holes[num_segs], _Holes[num_segs+1],
	    _Holes[0], _Holes[1]);
    num_bound = OnSegBound.size()-1;
    for (j=1; j<num_bound; ++j) {
      OnBound.push_back(OnSegBound[j]);
    }
}

void segment_intersection(double x1, double y1, double x2, double y2,
			 double x3, double y3, double x4, double y4,
			 double *x, double *y)
{
  double A1 = y2-y1, A2 = y4-y3;
  double B1 = x1-x2, B2 = x3-x4;
  double C1 = A1*x1 + B1*y1, C2 = A2*x3 + B2*y3;
  double det = A1*B2 - A2*B1;
  *x = (B2*C1-B1*C2) / det;
  *y = (A1*C2-A2*C1) / det;
}

void HybridMesher_2d::pave_region(size_t h, size_t i, size_t num_voronoi)
{
  double ex1 = _VorBound[h], hx1 = _StrBound[h];
  double ey1 = _VorBound[h+1], hy1 = _StrBound[h+1];
  double ex2 = _VorBound[i], hx2 = _StrBound[i];
  double ey2 = _VorBound[i+1], hy2 = _StrBound[i+1];
  double ex_change = ex2-ex1, ey_change = ey2-ey1;
  if (cross_product(hx2-ex1,hy2-ey1,ex_change,ey_change)<=0 ||
      cross_product(hx1-ex2,hy1-ey2,hx2-ex2,hy2-ey2)<=0 ||
      cross_product(ex1-hx2,ey1-hy2,hx1-hx2,hy1-hy2)<=0 ||
      cross_product(ex2-hx1,ey2-hy1,ex1-hx1,ey1-hy1)<=0) {
    printf("Structured mesh region not well formed!\n");
    exit(1);
  }
  
  size_t n;
  std::vector<size_t> OnVorBound;
  on_edge(_x, _y, OnVorBound, num_voronoi, ex1, ey1, ex2, ey2);
  size_t num_vor = OnVorBound.size()-1;
  
  if (!h) {
    if (_NumStr < 2) {
      printf("Too few structured layers requested.\n");
      exit(1);
    }
    double delta_h1x = (hx1-ex1)/_NumStr;
    double delta_h1y = (hy1-ey1)/_NumStr;
    for (n=1; n<=_NumStr; ++n) {
      _x.push_back(ex1 + n*delta_h1x);
      _y.push_back(ey1 + n*delta_h1y);
    }
  }

  double delta_h2x = (hx2-hx1)/num_vor;
  double delta_h2y = (hy2-hy1)/num_vor;
  double delta_h3x = (hx2-ex2)/_NumStr;
  double delta_h3y = (hy2-ey2)/_NumStr;
  std::vector<size_t> element;
  size_t num_nodes = _x.size() - _NumStr;
  for (n=1; n<num_vor; ++n) {
    double x, y;
    segment_intersection(_x[OnVorBound[n]], _y[OnVorBound[n]],
			 hx1 + n*delta_h2x, hy1 + n*delta_h2y,
			 _x[num_nodes], _y[num_nodes],
			 ex2+delta_h3x, ey2+delta_h3y,
			 &x, &y);
    _x.push_back(x);
    _y.push_back(y);
    element.clear();
    element.push_back(OnVorBound[n - 1]);
    element.push_back(num_nodes + (n-1) * _NumStr);
    element.push_back(num_nodes + n*_NumStr);
    element.push_back(OnVorBound[n]);
    _elements.push_back(element);

    for (size_t m=2; m<_NumStr; ++m) {
      segment_intersection(_x[OnVorBound[n]], _y[OnVorBound[n]],
			   hx1 + n*delta_h2x, hy1 + n*delta_h2y,
			   _x[num_nodes+m-1], _y[num_nodes+m-1],
			   ex2 + m*delta_h3x, ey2 + m*delta_h3y,
			   &x, &y);
      _x.push_back(x);
      _y.push_back(y);
      element.clear();
      element.push_back(num_nodes + (n-1)*_NumStr + m - 2);
      element.push_back(num_nodes + (n-1)*_NumStr + m - 1);
      element.push_back(num_nodes + n*_NumStr + m - 1);
      element.push_back(num_nodes + n*_NumStr + m - 2);
      _elements.push_back(element);
    }

    _x.push_back(hx1 + n*delta_h2x);
    _y.push_back(hy1 + n*delta_h2y);
    element.clear();
    element.push_back(num_nodes + (n-1)*_NumStr + _NumStr - 2);
    element.push_back(num_nodes + (n-1)*_NumStr + _NumStr - 1);
    element.push_back(num_nodes + n*_NumStr + _NumStr - 1);
    element.push_back(num_nodes + n*_NumStr + _NumStr - 2);
    _elements.push_back(element);
  }

  if (i) {
    _x.push_back(ex2 + delta_h3x);
    _y.push_back(ey2 + delta_h3y);
    element.clear();
    element.push_back(OnVorBound[num_vor - 1]);
    element.push_back(num_nodes + (num_vor-1) * _NumStr);
    element.push_back(num_nodes + num_vor*_NumStr);
    element.push_back(OnVorBound[num_vor]);
    _elements.push_back(element);

    for (n=2; n<=_NumStr; ++n) {
      _x.push_back(ex2 + n*delta_h3x);
      _y.push_back(ey2 + n*delta_h3y);
      element.clear();
      element.push_back(num_nodes + (num_vor-1)*_NumStr + n - 2);
      element.push_back(num_nodes + (num_vor-1)*_NumStr + n - 1);
      element.push_back(num_nodes + num_vor*_NumStr + n - 1);
      element.push_back(num_nodes + num_vor*_NumStr + n - 2);
      _elements.push_back(element);
    }
  } else {
    element.clear();
    element.push_back(OnVorBound[num_vor - 1]);
    element.push_back(num_nodes + (num_vor-1) * _NumStr);
    element.push_back(num_voronoi);
    element.push_back(OnVorBound[num_vor]);
    _elements.push_back(element);

    for (n=2; n<=_NumStr; ++n) {
      element.clear();
      element.push_back(num_voronoi + n - 2);
      element.push_back(num_nodes + (num_vor-1)*_NumStr + n - 2);
      element.push_back(num_nodes + (num_vor-1)*_NumStr + n - 1);
      element.push_back(num_voronoi + n - 1);
      _elements.push_back(element);
    }
  }
}

int HybridMesher_2d::execute()
{
  MeshingGenie_2d genie = MeshingGenie_2d(_dm,_VorBound,_Holes,_Cracks,1,false);

  genie.use_fixed_seed(_fixed_seed);

  genie.execute();

  size_t h;

  genie.get_Voronoi_Tessellation(_x, _y, _elements, 0.05 * _dm, false, true);

  size_t num_structured = _StrBound.size();
  if (num_structured > _VorBound.size()) {
    printf("Too many points on structured mesh boundary!\n");
    exit(1);
  }

  num_structured -= 2;
  size_t num_voronoi = _x.size();
  for (h=0; h<num_structured; h+=2) {
    pave_region(h, h+2, num_voronoi);
  }

  if (_StrBound.size() == _VorBound.size()) {
    pave_region(num_structured, 0, num_voronoi);
  }

  size_t num_holes = _Holes.size();
  for (h=0; h<num_holes; ++h) {
    std::vector<size_t> OnBound;
    on_bound(_x, _y, num_voronoi, _Holes[h], OnBound);

    size_t num_bound = OnBound.size();
    std::vector<double> Bound;
    int j;
    for (j=num_bound-1; j>=0; --j) {
      Bound.push_back(_x[OnBound[j]]);
      Bound.push_back(_y[OnBound[j]]);
    }

    std::vector< std::vector<double> > nholes, ncracks;
    MeshingGenie_2d subgenie = MeshingGenie_2d(_dm,Bound,nholes,ncracks,1,false);

    subgenie.use_fixed_seed(_fixed_seed);

    subgenie.execute();

    std::vector<double> x, y;
    std::vector< std::vector<size_t> > elements;
    subgenie.get_CDT_Tessellation(x, y, elements);

    std::vector<int> Remap;
    size_t num_2 = x.size();
    Remap.resize(num_2);
    for (j=0; j<num_2; ++j) {
     Remap[j] = -1;
    }

    std::vector<size_t> OnBound2;
    on_bound(x, y, num_2, _Holes[h], OnBound2);
    for (j=0; j<num_bound; ++j) {
      Remap[OnBound2[j]] = OnBound[j];
    }

    size_t num_change = _x.size();

    for (j=0; j<num_2; ++j) {
      if (Remap[j] < 0) {
	Remap[j] = num_change++;
	_x.push_back(x[j]);
	_y.push_back(y[j]);
      }
    }

    size_t num_elements = elements.size();
    for (j=0; j<num_elements; ++j) {
      size_t num = elements[j].size();
      std::vector<size_t> element;
      for (size_t k=0; k<num; ++k) {
	element.push_back(Remap[elements[j][k]]);
      }
      _elements.push_back(element);
    }
  }
  
  return 0;
}

void HybridMesher_2d::get_Tessellation(std::vector<double> &x,
				       std::vector<double> &y,
				       std::vector< std::vector<size_t> > &elements)
{
  x = _x;
  y = _y;
  elements = _elements;
}
