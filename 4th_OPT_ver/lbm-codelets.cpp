// Copyright (c) 2018 Graphcore Ltd. All rights reserved.
#include <poplar/Vertex.hpp>
#include "common.h"
#include <ipudef.h>

using namespace poplar;

class LBMcoreVertex : public Vertex
{
public:
	// Fields
	Input<Vector<float>> f;
	InOut<Vector<float>> ftmp;

	Input<Vector<float>> e;
	Input<Vector<float>> w;
	Input<float> tau;
	Input<int> chunk;

	// Compute function
	bool compute()
	{
		float tau1  = 1.0f/tau;
		float tau11 = (1.0f - tau1);

		for (unsigned tid = 0; tid < chunk; tid++)
		{

  			unsigned id = tid*sim_dir;

			// velocity components 
			const float speeds_0 = f[id+0];
			const float speeds_1 = f[id+1];
			const float speeds_2 = f[id+2];
			const float speeds_3 = f[id+3];
			const float speeds_4 = f[id+4];
			const float speeds_5 = f[id+5];
			const float speeds_6 = f[id+6];
			const float speeds_7 = f[id+7];
			const float speeds_8 = f[id+8];

			// local density total 
			const float rho =
				speeds_0 + speeds_1 + speeds_2 +
				speeds_3 + speeds_4 + speeds_5 +
				speeds_6 + speeds_7 + speeds_8;
    		
			// x velocity component 
			const float u_x = ((speeds_1 + speeds_5 + speeds_8)- 
							   (speeds_3 + speeds_6 + speeds_7)) /rho;

			// y velocity component 
			const float u_y = ((speeds_2 + speeds_5 + speeds_6)-
    		                   (speeds_4 + speeds_7 + speeds_8)) /rho;
				
			// square of velocity
  			float u_sq  = u_x*u_x + u_y*u_y;
  			float c_sq  = 1.0f - 1.5f*u_sq;
    		

			/*
  			for (unsigned l = 0; l < sim_dir; l++)
  			{
  				unsigned id1 = id+l;
    			float edotu 	 = e[l] * u_x + e[sim_dir + l] * u_y;
    			float tempfeq  = w[l] * rho * (edotu*(3.0f + 4.5f*edotu) + c_sq);
    			ftmp[id1] = tau11*f[id1] + tau1*tempfeq;
  			}
			*/

        	const float w2    = 1.0f/36.0f*rho*tau1;            // (1/36) *rho/tau
        	const float w1    = w2*4.0f;                        // (1/9) *rho/tau
        	const float w0    = w1*4.0f;                        // (4/9) *rho/tau

        	const auto relax2d = [&](const float speed, const float weight, const float edotu) -> float {
            	return tau11*speed + weight * ( edotu * (3.0f + 4.5f*edotu) + c_sq);
        	};
        	auto speeds_out_3 = relax2d(speeds_3, w1, -u_x);
        	auto speeds_out_4 = relax2d(speeds_4, w1, -u_y);
        	auto speeds_out_6 = relax2d(speeds_6, w2, -u_x + u_y);
        	auto speeds_out_7 = relax2d(speeds_7, w2, -u_x - u_y);
        	auto speeds_out_1 = relax2d(speeds_1, w1, u_x);
        	auto speeds_out_2 = relax2d(speeds_2, w1, u_y);
        	auto speeds_out_5 = relax2d(speeds_5, w2, u_x + u_y);
        	auto speeds_out_8 = relax2d(speeds_8, w2, u_x - u_y);
        	const float speeds_out_0 = tau11*speeds_0  + w0 * c_sq;

         	ftmp[id+0] = speeds_out_0;
         	ftmp[id+1] = speeds_out_1;
         	ftmp[id+2] = speeds_out_2;
     		ftmp[id+3] = speeds_out_3;
         	ftmp[id+4] = speeds_out_4;
        	ftmp[id+5] = speeds_out_5;
         	ftmp[id+6] = speeds_out_6;
         	ftmp[id+7] = speeds_out_7;
         	ftmp[id+8] = speeds_out_8;
		}
  		return true;
	}
};


class LBMstreamingVertex : public Vertex
{
public:
	// Fields
	InOut<Vector<float>> f;
	Input<Vector<float>> ftmp;
	Input<Vector<int>> eint;
	Input<unsigned> padding;
	Input<unsigned> dir;
	Input<unsigned> nx;
	Input<unsigned> ny;
	Input<int> chunk;

	// Compute function
	bool compute()
	{
  		const int imove = eint[sim_dir + dir];
  		const int jmove = eint[dir];

  		for (unsigned tid = 0; tid < chunk; tid++)
  		{
			const unsigned tid1 = tid + padding;

			//eugene 
      		//const int i = tid1 / nx + imove;
      		//const int j = tid1 % nx + jmove;

      		const int i = tid1 / WIDTH + imove;
      		const int j = tid1 % WIDTH + jmove;

      		if((i>-1 && i<ny) && (j>-1 && j<nx))
      		{
  				unsigned id = tid*sim_dir + dir;
          		f[id] = ftmp[id];
      		}
  		}
  		return true;
	}
};


class LBMboundaryVertex : public Vertex
{
public:
  	// Fields
  	InOut<Vector<float>> f;
  	Input<Vector<float>> ftmp;
  	Input<Vector<float>> e;
  	Input<Vector<float>> w;
  	Input<Vector<int>> eint;
  	Input<Vector<unsigned>> oppi;
  	Input<float> rho0;
  	Input<float> u0;
	Input<unsigned> nx;
	Input<unsigned> ny;
  	Input<unsigned> padding;
  	Input<int> chunk;

  	// Compute function
  	bool compute()
	{
		const float trho = rho0;
		const float tvel_x = u0;
		const float tvel_y = 0.0f;

		/*
		const float w0    = 6.0f*trho;
		const float w2    = w[2]*w0;    	// (1/36) *rho *6.0
   		const float w1    = w2*4.0f;     	// (1/9) *rho *6.0 
		*/

    	for (unsigned tid = 0; tid < chunk; tid++)
    	{
			unsigned tid1 = tid + padding;
        	unsigned i = tid1 / nx; // y축좌표
        	unsigned j = tid1 % nx; // x축좌표

        	// 맨 윗줄 벽 노드의 바로 밑줄노드
        	if (i == ny - 1)
        	{
            	if (j > 0 && j < nx - 1)  //  0 < x < nx-1 맨 왼쪽, 오른쪽 노드 제외 (윗쪽 양 모서리벽)
            	{

                	for (unsigned l = 0; l < sim_dir; ++l)
                	{
                    	int imove = eint[sim_dir + l];
                    	int jmove = eint[l];

                    	if (imove > 0)  // up
                    	{
                    		unsigned oppind = oppi[l];
                        	f[tid*sim_dir + oppind] = ftmp[tid*sim_dir + l] - 6.0f*w[l] * trho*(e[l] * tvel_x + e[sim_dir + l] * tvel_y);
                    	}
                	}

					/*
        			const auto bctop = [&](const float speed, const float weight, const float edotu) -> float {
            			return speed - weight*edotu;
					};

					// case of "imove > 0" : 2,5,6 directions  
					int id = tid*sim_dir;
        			f[id+4] = bctop(ftmp[id+2], w1, tvel_y);
        			f[id+7] = bctop(ftmp[id+5], w2, tvel_x + tvel_y);
        			f[id+8] = bctop(ftmp[id+6], w2, -tvel_x + tvel_y);
					*/
            	}
        	}

        	// 맨 밑줄 바닥 노드
        	else if (i == 0)
        	{
            	for (unsigned l = 0; l < sim_dir; ++l)
            	{
                	int imove = eint[sim_dir + l];
                	int jmove = eint[l];

                	if (j + jmove > nx-1 || j + jmove < 0 || imove < 0)  // right, left, down
                	{
                		unsigned oppind = oppi[l];
                    	f[tid*sim_dir+oppind] = ftmp[tid*sim_dir+l];
                	}
            	}
				
        	}

        	// 왼쪽 벽
        	else if (j == 0)
        	{
            	for (unsigned l = 0; l < sim_dir; ++l)
            	{
                	int imove = eint[sim_dir + l];
                	int jmove = eint[l];

                	if (i + imove > ny - 1 || jmove < 0 || i + imove < 0)  // up, left, right
                	{
                		unsigned oppind = oppi[l];
                    	f[tid*sim_dir + oppind] = ftmp[tid*sim_dir + l];
                	}
            	}
        	}

        	// 오른쪽 벽
        	else if (j == nx-1)
        	{
            	for (unsigned l = 0; l < sim_dir; ++l)
            	{
                	int imove = eint[sim_dir + l];
                	int jmove = eint[l];

                	if (jmove > 0 || i+imove > ny - 1 || i+imove < 0) // right, up, down
                	{
                		unsigned oppind = oppi[l];
                    	f[tid*sim_dir + oppind] = ftmp[tid*sim_dir + l];
                	}
            	}
        	}
    	}
    	return true;
  	}
};





/* 
class KepresaveVertex : public Vertex
{
public:
  // Fields
  InOut<Vector<float>> prho;
  InOut<Vector<float>> pvel;

  Input<Vector<float>> rho;
  Input<Vector<float>> vel;
  Input<int> chunk;

  // Compute function
  bool compute()
  {
    for (int tid = 0; tid < chunk; tid++)
    {
      prho[tid] 			= rho[tid];
      pvel[tid * dim_s] 	= vel[tid * dim_s];
      pvel[tid * dim_s + 1]	= vel[tid * dim_s + 1];
    }
    return true;
  }
};

// Print result on Host!!!!
class KeprintVertex : public Vertex
{
public:
  // Fields
  InOut<Vector<float>> rho;
  InOut<Vector<float>> vel;
  InOut<Vector<float>> prtrho;
  InOut<Vector<float>> prtvelx;
  InOut<Vector<float>> prtvely;
  Input<float> u0;
  Input<int> nx;
  Input<int> ny;

  bool compute()
  {
    for (int i = 0; i < ny; ++i)
    {
        int j = nx / 2;
        int id = j + i * nx;
        prtrho[i] = rho[id];
      	prtvelx[i] = vel[id * dim_s] / u0;
      	prtvely[i] = vel[id * dim_s+1] / u0;
    }
    return true;
  }
};


class KeErrorVertex : public Vertex
{
public:
  // Fields
  InOut<Vector<float>> error;
  Input<Vector<float>> rho;
  Input<Vector<float>> prho;
  Input<int> chunk;

  // Compute function
  bool compute()
  {
    float temp = 0.0;
    error[0] = 0.0;
    
    for (int tid = 0; tid < chunk; tid++)
    {
      temp = prho[tid] - rho[tid];
      if(temp < 0) temp = -1.0 * temp;
      error[0] += temp;
    }
    return true;
  }
};

*/
