#include <iostream>
#include <chrono>
#include <vector>
#include <memory>
#include <poplar/Graph.hpp>
#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poputil/TileMapping.hpp>
#include "common.h"

using namespace std;
using namespace poplar;
using namespace poplar::program;

#define  div1 8000

int slice_chunk(const unsigned idx, const int length, const unsigned chunk, int slice_in[2], int *chunksize)
{
    unsigned begin = idx*chunk;
    unsigned end   = (idx+1)*chunk;

    *chunksize = end - begin;
    if(*chunksize<1) return 1;

    if(end > length) end = length;
    begin *= sim_dir;
    end   *= sim_dir;
    slice_in[0]  = begin;
    slice_in[1]  = end;

    return 0;
}


int slice_streaming(const unsigned idx, const unsigned dir, const int length, const unsigned chunk, const unsigned nx, int slice_in[2], int slice_out[2], int *chunksize, unsigned *padding)
{
    int D2Q9eint[sim_dir*dim_s] = {0,1,0,-1,0,1,-1,-1,1,  0,0,1,0,-1,1,1,-1,-1};
    int margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];

    int begin = idx*chunk;
    int end   = (idx+1)*chunk;

    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;

    *chunksize = end - begin;
    if(*chunksize<1) return 1;

    *padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
    begin *= sim_dir;
    end   *= sim_dir;

    slice_in[0]  = begin;
    slice_in[1]  = end;
    slice_out[0] = begin_f;
    slice_out[1] = end_f;

    return 0;
}

Program buildMainProgram(Graph &graph, Tensor f, Tensor ftmp, Tensor e, Tensor w, Tensor eint, Tensor oppi, float rho0, float u0, float tau, int length, unsigned nx, unsigned ny)
{
	unsigned chunk = length/div1;
  	int slice_in[2];
  	int slice_out[2];
  	int chunksize;
	unsigned dir; 
	unsigned padding;


// *****************************************************************************************
// Macroscopic flow fields and Collision 
// *****************************************************************************************
	ComputeSet LBMcoreCS = graph.addComputeSet("LBMcoreCS");
	for (unsigned i = 0; i < div1 ; i++)
	{
		if(slice_chunk(i,length,chunk,slice_in,&chunksize)) continue;
  		auto v = graph.addVertex(LBMcoreCS,
                          		"LBMcoreVertex",
                          		{{"f", 		f.slice(slice_in[0], slice_in[1])		},
                          		 {"ftmp",	ftmp.slice(slice_in[0], slice_in[1])	},
                           		 {"e", 		e										},
                           		 {"w", 		w										},
                           		 {"tau", 	tau										},
                           		 {"chunk",	chunksize 								}});

  		unsigned tile = (i* tileNum)/div1;
  		graph.setTileMapping(v, 	tile);
  		graph.setTileMapping(f.slice(slice_in[0], slice_in[1]), tile);
  		graph.setTileMapping(ftmp.slice(slice_in[0], slice_in[1]), tile);
  	}



// *****************************************************************************************
// streaming in 9 directions : Separate ComputeSets to avoid Race Conditions 
// *****************************************************************************************

  	ComputeSet Stream_mid_CS = graph.addComputeSet("Stream_mid_CS");
  	for (unsigned i = 0; i < div1 ; i++)
  	{
  		dir = 0;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;
    	//printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,slice_in[0],slice_in[1],slice_out[0],slice_out[1],chunksize,padding);

    	auto v  = graph.addVertex(Stream_mid_CS,
                             	"LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    	unsigned tile = (i* tileNum)/div1;
    	graph.setTileMapping(v , 	 tile);
  	}


  	ComputeSet Stream_E_CS = graph.addComputeSet("Stream_E_CS");
  	for (unsigned i = 0; i < div1 ; i++)
  	{
  		dir = 1;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;
    	auto v  = graph.addVertex(Stream_E_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_N_CS = graph.addComputeSet("Stream_N_CS");
  for (unsigned i = 0; i < div1 ; i++)
  {
  		dir = 2;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;
    auto v  = graph.addVertex(Stream_N_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_W_CS = graph.addComputeSet("Stream_W_CS");
  for (unsigned i = 0; i < div1 ; i++)
  {
  		dir = 3;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;

    auto v  = graph.addVertex(Stream_W_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_S_CS = graph.addComputeSet("Stream_S_CS");
  for (unsigned i = 0; i < div1 ; i++)
  {
  		dir = 4;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;

    auto v  = graph.addVertex(Stream_S_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_NE_CS = graph.addComputeSet("Stream_NE_CS");
  for (unsigned i = 0; i < div1 ; i++)
  {
  		dir = 5;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;
    auto v  = graph.addVertex(Stream_NE_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_NW_CS = graph.addComputeSet("Stream_NW_CS");
  for (unsigned i = 0; i < div1 ; i++)
  {
  		dir = 6;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;

    auto v  = graph.addVertex(Stream_NW_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }


  ComputeSet Stream_SW_CS = graph.addComputeSet("Stream_SW_CS");
  for (unsigned i = 0; i < div1 ; i++)
  {
  		dir = 7;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;

    auto v  = graph.addVertex(Stream_SW_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_SE_CS = graph.addComputeSet("Stream_SE_CS");
  for (unsigned i = 0; i < div1 ; i++)
  {
  		dir = 8;
		if(slice_streaming(i,dir,length,chunk,nx,slice_in,slice_out,&chunksize,&padding)) continue;

    auto v  = graph.addVertex(Stream_SE_CS,
                             "LBMstreamingVertex",
                             	{{"f", 			f.slice(slice_out[0], slice_out[1])		},
                              	 {"ftmp", 		ftmp.slice(slice_in[0], slice_in[1])	},
                              	 {"eint", 		eint									},
                              	 {"padding", 	padding 								},
                              	 {"dir", 	 	dir										},
                              	 {"nx", 		nx										},
                              	 {"ny", 		ny										},
                              	 {"chunk",		chunksize 								}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }



// *****************************************************************************************
// *****************************************************************************************
// *****************************************************************************************

  ComputeSet LBMboundCS = graph.addComputeSet("LBMboundCS");
  for (int i = 0; i < div1; i++)
  {
    int begin  = i * chunk;
    int end    = (i + 1) * chunk;
    if(end > length) end = length;
    int chunksize = end - begin;
    if(chunksize<1) continue;
	
    int padding = begin;
	begin *= sim_dir;
	end   *= sim_dir;

//    printf("st=%4d,ed=%4d |size=%4d padding=%4d\n",begin,end,chunksize,padding);
    auto v = graph.addVertex(LBMboundCS,
                             "LBMboundaryVertex",
                             {{"f", 		f.slice(begin, end)		},
                              {"ftmp", 		ftmp.slice(begin, end)	},
                              {"e", 		e			},
                              {"w", 		w			},
                              {"eint", 		eint		},
                              {"oppi", 		oppi		},
                              {"rho0", 		rho0		},
                              {"u0", 		u0			},
                              {"nx", 		nx			},
                              {"ny", 		ny			},
                              {"padding", 	padding 	},
                              {"chunk",		chunksize 	}});
    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  //return Sequence(Execute(LBMcoreCS), Execute(LBMboundCS));
  return Sequence(Execute(LBMcoreCS), \
			  	Execute(Stream_mid_CS), Execute(Stream_E_CS), Execute(Stream_N_CS),Execute(Stream_W_CS), \
			  	Execute(Stream_S_CS), Execute(Stream_NE_CS), Execute(Stream_NW_CS), Execute(Stream_SW_CS), \
			  	Execute(Stream_SE_CS),  \
		        Execute(LBMboundCS));
}






int main()
{
  Inputval *pinput = new Inputval();
  HostNode *pnode = new HostNode();
  D2Q9model *pmodel = new D2Q9model();
  
  D2Q9.e[0][0] = 0.0; D2Q9.e[0][1] = 1.0; D2Q9.e[0][2] = 0.0;
  D2Q9.e[0][3] = -1.0; D2Q9.e[0][4] = 0.0; D2Q9.e[0][5] = 1.0;
  D2Q9.e[0][6] = -1.0; D2Q9.e[0][7] = -1.0; D2Q9.e[0][8] = 1.0;

  D2Q9.e[1][0] = 0.0; D2Q9.e[1][1] = 0.0; D2Q9.e[1][2] = 1.0;
  D2Q9.e[1][3] = 0.0; D2Q9.e[1][4] = -1.0; D2Q9.e[1][5] = 1.0;
  D2Q9.e[1][6] = 1.0; D2Q9.e[1][7] = -1.0; D2Q9.e[1][8] = -1.0;

  D2Q9.eint[0][0] = 0; D2Q9.eint[0][1] = 1; D2Q9.eint[0][2] = 0;
  D2Q9.eint[0][3] = -1; D2Q9.eint[0][4] = 0; D2Q9.eint[0][5] = 1;
  D2Q9.eint[0][6] = -1; D2Q9.eint[0][7] = -1; D2Q9.eint[0][8] = 1;

  D2Q9.eint[1][0] = 0; D2Q9.eint[1][1] = 0; D2Q9.eint[1][2] = 1;
  D2Q9.eint[1][3] = 0; D2Q9.eint[1][4] = -1; D2Q9.eint[1][5] = 1;
  D2Q9.eint[1][6] = 1; D2Q9.eint[1][7] = -1; D2Q9.eint[1][8] = -1;

  D2Q9.oppi[0] = 0; D2Q9.oppi[1] = 3; D2Q9.oppi[2] = 4;
  D2Q9.oppi[3] = 1; D2Q9.oppi[4] = 2; D2Q9.oppi[5] = 7;
  D2Q9.oppi[6] = 8; D2Q9.oppi[7] = 5; D2Q9.oppi[8] = 6;

  D2Q9.w[0] = 4.0 / 9.0; D2Q9.w[1] = 1.0 / 9.0; D2Q9.w[2] = 1.0 / 9.0;
  D2Q9.w[3] = 1.0 / 9.0; D2Q9.w[4] = 1.0 / 9.0; D2Q9.w[5] = 1.0 / 36.0;
  D2Q9.w[6] = 1.0 / 36.0; D2Q9.w[7] = 1.0 / 36.0; D2Q9.w[8] = 1.0 / 36.0;

  INPUT.h = 1000; INPUT.nx = INPUT.h; INPUT.ny = INPUT.h;
  INPUT.maxitr = 5000; INPUT.rho0 = 1.0; INPUT.u0 = 0.1;
  INPUT.Re = 100; INPUT.length = INPUT.nx * INPUT.ny;
  D2Q9.tau = 3.0 * INPUT.u0 * INPUT.h / INPUT.Re + 0.50;


  DeviceManager manager = DeviceManager::createDeviceManager();

  Device device;
  bool success = false;

  for (auto &hwDevice : manager.getDevices(TargetType::IPU, 1))
  {
    device = std::move(hwDevice);
    std::cerr << "Trying to attatch to IPU" << device.getId() << std::endl;

    if ((success = device.attach()))
      std::cerr << "Attached to IPU" << device.getId() << std::endl;
  }

  if (!success)
  {
    std::cerr << "Error attaching to device" << std::endl;
    return -1;
  }

  Target target = device.getTarget();
  Graph graph(target);

  // Add codelets to the graph
  graph.addCodelets("lbm-codelets.cpp");
  //graph.addCodelets("lbm-codelets.cpp", CodeletFileType::Auto,"-O3");

  //============================================================================================
  // Create and initialize variables on Host(CPU)
  std::cout << "Create and initialize variables on Host(CPU)\n";
  vector<float> hf(INPUT.length * sim_dir, 0);
  for (unsigned tid = 0; tid < INPUT.length; tid++)
  {
    float velx0 = 0.0;
    float vely0 = 0.0;
    // equilbrium 값으로 초기화
    for (int l = 0; l < sim_dir; l++)
    {
      float edotu = D2Q9.e[0][l] * velx0 + D2Q9.e[1][l] * vely0;
      float usquare = velx0*velx0 + vely0*vely0;
      hf[tid*sim_dir + l] = D2Q9.w[l]*INPUT.rho0 * (1.0 + 3.0*edotu + 4.5*edotu*edotu - 1.5*usquare);
    }
  }

  //============================================================================================    ​
  //============================================================================================
  // Add Variables to the graph
  std::cout<< "Add Variables to the graph\n";
  Tensor f 		= graph.addVariable(FLOAT, {unsigned(INPUT.length) * sim_dir}, "f");
  Tensor ftmp 	= graph.addVariable(FLOAT, {unsigned(INPUT.length) * sim_dir}, "ftmp");
  Tensor e 		= graph.addVariable(FLOAT, {dim_s * sim_dir}, "e");
  Tensor eint 	= graph.addVariable(INT, {dim_s * sim_dir}, "eint");
  Tensor w 		= graph.addVariable(FLOAT, {sim_dir}, "w");
  Tensor oppi 	= graph.addVariable(UNSIGNED_INT, {sim_dir}, "oppi");
  //============================================================================================

  graph.setTileMapping(e, tileNum);
  graph.setTileMapping(w, tileNum);
  graph.setTileMapping(eint, tileNum);
  graph.setTileMapping(oppi, tileNum);

  //============================================================================================
  // Create Host write handles for variables on the Device
  std::cout << "Create Host write handles for variables on the Device\n";
  graph.createHostWrite("H2D-f", f);
  graph.createHostWrite("H2D-ftmp", ftmp);

  graph.createHostWrite("H2D-e", e);
  graph.createHostWrite("H2D-eint", eint);
  graph.createHostWrite("H2D-w", w);
  graph.createHostWrite("H2D-oppi", oppi);

  graph.createHostRead("D2H-f", ftmp);
  //============================================================================================

  //============================================================================================
  // Copy Host data via the write handle to variables on the Device
  Sequence mainProg, iterProg;

  Program single_step = buildMainProgram(graph, f, ftmp, e, w, eint, oppi, INPUT.rho0, INPUT.u0, D2Q9.tau, INPUT.length, INPUT.nx, INPUT.ny);

  mainProg.add(single_step);
  iterProg.add(Repeat(INPUT.maxitr, mainProg));
  
  Engine engine(graph, iterProg);
  engine.load(device);

  engine.writeTensor("H2D-f", 		hf.data());
  engine.writeTensor("H2D-ftmp", 	hf.data());

  engine.writeTensor("H2D-e", 		D2Q9.e);
  engine.writeTensor("H2D-eint", 	D2Q9.eint);
  engine.writeTensor("H2D-w", 		D2Q9.w);
  engine.writeTensor("H2D-oppi", 	D2Q9.oppi);

  std::cout << "Copy Host data via the write handle to variables on the Device\n";
  
  //============================================================================================

    auto start_t = std::chrono::high_resolution_clock::now();
  std::cout << "Running program\n";
  engine.run(0);
  std::cout << "Program complete\n";

    auto end_t = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end_t - start_t;
    printf("Elapsed time(IPU) = %fs | h=%d\n",elapsed.count(), INPUT.h);


  std::cout << " >>> Copy velocity data from Device via read handle\n";
  engine.readTensor("D2H-f", hf.data());

	/*
 	printf("====================== VEL ====================\n");
	for (int i = INPUT.ny-1; i >= 0; i--)
  	{
 		for (int j = 0; j < INPUT.nx; j++) 
    	{
	  		int id = (j + i * INPUT.nx)*sim_dir;
			float rho  = hf[id+0]+hf[id+1]+hf[id+2]+hf[id+3]+hf[id+4]+hf[id+5]+hf[id+6]+hf[id+7]+hf[id+8];
           	float velx = hf[id+1] - hf[id+3] + hf[id+5] - hf[id+6] - hf[id+7] + hf[id+8];
           	float vely = hf[id+2] - hf[id+4] + hf[id+5] + hf[id+6] - hf[id+7] - hf[id+8];
	  		printf("%10.7f ",velx/(rho*INPUT.u0));
        }
	 	printf("\n");
  	}
  	printf("\n");
  	printf("\n");
	*/

	/*
  	printf("===================== Ftmp ====================\n");
  	for (int l=0; l<sim_dir; l++)
  	{
		for (int i = INPUT.ny-1; i >= 0; i--)
      	{
          	for (int j = 0; j < INPUT.nx; j++)
          	{
              	int id = j + i * INPUT.nx;
              	printf("%10.6f ",hf[id*sim_dir+l]);
          	}
          	printf("\n");
      	}
      	printf("\n");
  	}
  	printf("\n");
  	printf("\n");
	*/

  	//for (int i = 0 ; i< INPUT.ny; i++)
  	for (int i = INPUT.ny - 10 ; i< INPUT.ny; i++)
  	{
		int id = (INPUT.nx / 2 + i * INPUT.nx)*sim_dir;
		float rho  = hf[id+0]+hf[id+1]+hf[id+2]+hf[id+3]+hf[id+4]+hf[id+5]+hf[id+6]+hf[id+7]+hf[id+8];
       	float velx = hf[id+1] - hf[id+3] + hf[id+5] - hf[id+6] - hf[id+7] + hf[id+8];
		printf("%10.6f, %10.6f \n",(i + 0.5)/INPUT.ny, velx/(rho*INPUT.u0));
  	}

  return 0;
}
