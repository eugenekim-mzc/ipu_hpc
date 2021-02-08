// Copyright (c) 2018 Graphcore Ltd. All rights reserved.
#include <poplar/Vertex.hpp>
#include "common.h"

using namespace poplar;


class KemacroscopicVertex : public Vertex
{
public:
  // Fields
  InOut<Vector<float>> rho;
  InOut<Vector<float>> vel;
  Input<Vector<float>> f;
  Input<Vector<float>> e;
  Input<Vector<float>> w;
  Input<int> chunk;


  // Compute function
  bool compute()
  {
    for (int tid = 0; tid < chunk; tid++)
    {
        float temprho = 0;
        float tempvelx = 0;
        float tempvely = 0;

        for (int l = 0; l < sim_dir; ++l)
        {
            temprho += f[tid*sim_dir + l];
            tempvelx += f[tid*sim_dir + l] * e[l];
            tempvely += f[tid*sim_dir + l] * e[sim_dir + l];
        }
        tempvelx /= temprho;
        tempvely /= temprho;
        rho[tid] = temprho;
        vel[tid*dim_s] = tempvelx;
        vel[tid*dim_s+1] = tempvely;
    }
    return true;
  }


/*
  // Compute function
  bool compute()
  {
	int tid_sim_dir, tid_dim_s, tid_sim_dir1;
    float temprho, tempvelx, tempvely;

    for (int tid = 0; tid < chunk; tid++)
    {
      temprho = 0; tempvelx = 0; tempvely = 0;

	  tid_sim_dir = tid*sim_dir;
	  tid_dim_s   = tid*dim_s;
      for (int l = 0; l < sim_dir; ++l)
      {
        //temprho  += f[tid * sim_dir + l];
        //tempvelx += f[tid * sim_dir + l] * e[l];
        //tempvely += f[tid * sim_dir + l] * e[sim_dir + l];

	    tid_sim_dir1 = tid*sim_dir + l;
        temprho  += f[tid_sim_dir1];
        tempvelx += f[tid_sim_dir1] * e[l];
        tempvely += f[tid_sim_dir1] * e[sim_dir + l];
      }

      //tempvelx /= temprho;
      //tempvely /= temprho;
      rho[tid] = temprho;

      //vel[tid * dim_s] = tempvelx;
      //vel[tid * dim_s + 1] = tempvely;

      vel[tid_dim_s]     = tempvelx/temprho;
      vel[tid_dim_s + 1] = tempvely/temprho;
    }
    return true;
  }
*/

};

class KecollisionVertex : public Vertex
{
public:
  // Fields
  Input<Vector<float>> rho;
  Input<Vector<float>> vel;
  Input<Vector<float>> f;
  InOut<Vector<float>> ftmp;

  Input<Vector<float>> e;
  Input<Vector<float>> w;
  Input<float> tau;
  Input<int> chunk;

  // Compute function
  bool compute()
  {
    for (int tid = 0; tid < chunk; tid++)
    {
        float edotu = 0;
        float usquare = 0;
        float tempfeq = 0;
        float tempftmp = 0;
        float temprho = rho[tid];
        for (int l = 0; l < sim_dir; l++)
        {
            edotu = e[l] * vel[tid*dim_s] + e[sim_dir + l] * vel[tid*dim_s + 1];
            usquare = vel[tid*dim_s] * vel[tid*dim_s] + vel[tid*dim_s + 1] * vel[tid*dim_s + 1];
            tempfeq = w[l] * temprho * (1.0 + 3.0*edotu + 4.5*edotu * edotu - 1.5*usquare);
            tempftmp = f[tid*sim_dir + l] + 1.0 / tau * (tempfeq - f[tid*sim_dir + l]);
            ftmp[tid*sim_dir + l] = tempftmp;
        }
    }
    return true;
  }


/*
  // Compute function
  bool compute()
  {
	int   tid_sim_dir, tid_dim_s, tid_sim_dir1;
	float edotu, usquare, tempfeq, tempftmp, temprho;
	float tempvelx, tempvely;
	float tau1 = 1.0/tau;

    for (int tid = 0; tid < chunk; tid++)
    {
      edotu 	= 0;
      usquare 	= 0;
      tempfeq 	= 0;
      tempftmp 	= 0;
      temprho 	= rho[tid];
	  tempvelx 	= vel[tid_dim_s];
	  tempvely 	= vel[tid_dim_s+1];

	  tid_sim_dir = tid*sim_dir;
	  tid_dim_s   = tid*dim_s;

      usquare  = tempvelx*tempvelx + tempvely*tempvely;
      for (int l = 0; l < sim_dir; l++)
      {
        //edotu = e[l] * vel[tid * dim_s] + e[sim_dir + l] * vel[tid * dim_s + 1];
        //usquare = vel[tid * dim_s] * vel[tid * dim_s] + vel[tid * dim_s + 1] * vel[tid * dim_s + 1];
        //tempfeq = w[l] * temprho * (1.0 + 3.0 * edotu + 4.5 * edotu * edotu - 1.5 * usquare);
        //tempftmp = f[tid * sim_dir + l] + 1.0 / tau * (tempfeq - f[tid * sim_dir + l]);
        //ftmp[tid * sim_dir + l] = tempftmp;

	  	tid_sim_dir1 = tid_sim_dir+l;
        edotu 	 = e[l] * tempvelx + e[sim_dir + l] * tempvely;
        tempfeq  = w[l] * temprho * (1.0 + 3.0 * edotu + 4.5 * edotu * edotu - 1.5 * usquare);
        tempftmp = f[tid_sim_dir1] + tau1 * (tempfeq - f[tid_sim_dir1]);
        ftmp[tid_sim_dir1] = tempftmp;
      }
    }
    return true;
  }
*/
};

class KestreamingVertex : public Vertex
{
public:
  // Fields
  InOut<Vector<float>> f;
  Input<Vector<float>> ftmp;
  Input<Vector<int>> eint;
  Input<int> padding;
  Input<int> dir;
  Input<int> nx;
  Input<int> ny;
  Input<int> chunk;

  // Compute function
  bool compute()
  {
    int imove = eint[sim_dir + dir];
    int jmove = eint[dir];
    for (int tid = 0; tid < chunk; tid++)
    {
        int tid1 = tid + padding;
        int i = tid1/nx + imove;
        int j = tid1%nx + jmove;
        if((i>-1 && i<ny) && (j>-1 && j<nx))
        {
            f[tid*sim_dir + dir] = ftmp[tid*sim_dir + dir];
        }
    }
    return true;
  }
};


class KeboundaryVertex : public Vertex
{
public:
  // Fields
  InOut<Vector<float>> f;
  Input<Vector<float>> ftmp;
  Input<Vector<float>> e;
  Input<Vector<float>> w;
  Input<Vector<int>> eint;
  Input<Vector<int>> oppi;
  Input<float> rho0;
  Input<float> u0;
  Input<int> nx;
  Input<int> ny;
  Input<int> padding;
  Input<int> chunk;

  // Compute function
  bool compute()
  {
    for (int tid = 0; tid < chunk; tid++)
    {
		int tid1 = tid + padding;
        int i = tid1 / nx; // y축좌표
        int j = tid1 % nx; // x축좌표

        // 맨 윗줄 벽 노드의 바로 밑줄노드
        if (i == ny - 1)
        {
            if (j > 0 && j < nx - 1)  //  0 < x < nx-1 맨 왼쪽, 오른쪽 노드 제외 (윗쪽 양 모서리벽)
            {
                float trho = rho0;
                float tvel[dim_s];

                tvel[0] = u0;
                tvel[1] = 0.0;

                for (int l = 0; l < sim_dir; ++l)
                {
                    int imove = eint[sim_dir + l];
                    int jmove = eint[l];

                    if (imove > 0)  // up
                    {
                    	int oppind = oppi[l];
                        f[tid*sim_dir + oppind] = ftmp[tid*sim_dir + l] - 6.0*w[l] * trho*(e[l] * tvel[0] + e[sim_dir + l] * tvel[1]);
                    }
                }
            }
        }
        // 맨 밑줄 바닥 노드
        else if (i == 0)
        {
            for (int l = 0; l < sim_dir; ++l)
            {
                int imove = eint[sim_dir + l];
                int jmove = eint[l];

                if (j + jmove > nx-1 || j + jmove < 0 || imove < 0)  // right, left, down
                {
                	int oppind = oppi[l];
                    f[tid*sim_dir+oppind] = ftmp[tid*sim_dir+l];
                }
            }
        }
        // 왼쪽 벽
        else if (j == 0)
        {
            for (int l = 0; l < sim_dir; ++l)
            {
                int imove = eint[sim_dir + l];
                int jmove = eint[l];

                if (i + imove > ny - 1 || jmove < 0 || i + imove < 0)  // up, left, right
                {
                	int oppind = oppi[l];
                    f[tid*sim_dir + oppind] = ftmp[tid*sim_dir + l];
                }
            }
        }
        // 오른쪽 벽
        else if (j == nx-1)
        {
            for (int l = 0; l < sim_dir; ++l)
            {
                int imove = eint[sim_dir + l];
                int jmove = eint[l];

                if (jmove > 0 || i+imove > ny - 1 || i+imove < 0) // right, up, down
                {
                	int oppind = oppi[l];
                    f[tid*sim_dir + oppind] = ftmp[tid*sim_dir + l];
                }
            }
        }
    }
    return true;
  }



/*
  // Compute function
  bool compute()
  {

	int i, j, oppind, imove, jmove;
	float trho;
	float tvel[dim_s];

    for (int tid = 0; tid < chunk; tid++)
    {
      i = tid / nx; // y축좌표
      //j = tid % nx; // x축좌표
      j = tid - i; // x축좌표

      // 맨 윗줄 벽 노드의 바로 밑줄노드
      if (i == ny - 1)
      {
        if (j > 0 && j < nx - 1) //  0 < x < nx-1 맨 왼쪽, 오른쪽 노드 제외 (윗쪽 양 모서리벽)
        {
          trho = rho0;
          tvel[0] = u0;
          tvel[1] = 0.0;

          for (int l = 0; l < sim_dir; ++l)
          {
            oppind = oppi[l];
            imove = eint[sim_dir + l];
            jmove = eint[l];

            if (imove > 0) // up
            {
              f[tid * sim_dir + oppind] = ftmp[tid * sim_dir + l] - 6.0 * w[l] * trho * (e[l] * tvel[0] + e[sim_dir + l] * tvel[1]);
            }
          }
        }
      }
      // 맨 밑줄 바닥 노드
      else if (i == 0)
      {
        for (int l = 0; l < sim_dir; ++l)
        {
          oppind = oppi[l];
          imove = eint[sim_dir + l];
          jmove = eint[l];

          if (j + jmove > nx - 1 || j + jmove < 0 || imove < 0) // right, left, down
          {
            f[tid * sim_dir + oppind] = ftmp[tid * sim_dir + l];
          }
        }
      }
      // 왼쪽 벽
      else if (j == 0)
      {
        for (int l = 0; l < sim_dir; ++l)
        {
          oppind = oppi[l];
          imove = eint[sim_dir + l];
          jmove = eint[l];

          if (i + imove > ny - 1 || jmove < 0 || i + imove < 0) // up, left, right
          {
            f[tid * sim_dir + oppind] = ftmp[tid * sim_dir + l];
          }
        }
      }
      // 오른쪽 벽
      else if (j == nx - 1)
      {
        for (int l = 0; l < sim_dir; ++l)
        {
          oppind = oppi[l];
          imove = eint[sim_dir + l];
          jmove = eint[l];

          if (jmove > 0 || i + imove > ny - 1 || i + imove < 0) // right, up, down
          {
            f[tid * sim_dir + oppind] = ftmp[tid * sim_dir + l];
          }
        }
      }
    }
    return true;
  }
*/
};



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



/* 
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

*/

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
