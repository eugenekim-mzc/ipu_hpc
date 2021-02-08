// Copyright (c) 2018 Graphcore Ltd. All rights reserved.

#include <poplar/Vertex.hpp>
using namespace poplar;

class plus4ways_Vertex : public Vertex
{
public:
  // Fields
  InOut<Vector<int>> a;
  Input<int> nx;
  Input<int> ny;
  Input<int> padding;
  Input<int> chunk;

  // Compute function
  bool compute()
  {
    for (int tid = 0; tid < chunk; tid++)
    {
		int tid1 = tid+padding;
        int i = tid1/nx;
        int j = tid1%nx;

		// Add '1' to 4 neighbor points (west,east,north,south)
		if((i>0 && i<ny-1) && (j>0 && j<nx-1))
		{
	    	a[tid-1] += 1;   //west 
	    	a[tid+1] += 1;   //east   
	    	a[tid+nx] += 1;   //north  
	    	a[tid-nx] += 1;   //south 
		}
	}
  	return true;
  }
};
