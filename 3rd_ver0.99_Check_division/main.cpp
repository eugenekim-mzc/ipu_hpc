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


Program buildMainProgram(Graph &graph, Tensor rho, Tensor vel, Tensor prho, Tensor pvel, Tensor f, Tensor ftmp, Tensor e, Tensor w, Tensor eint, Tensor oppi, float rho0, float u0, float tau, int length, int nx, int ny)
{
  int D2Q9eint[sim_dir*dim_s] = {0,1,0,-1,0,1,-1,-1,1,  0,0,1,0,-1,1,1,-1,-1};
  int dir = 0;
  int margin = 0;
  
  /*
  int chunk = length/div;
  for (int i = 0; i < div ; i++)
  {
    int begin_rho  = i * chunk;
    int end_rho    = (i + 1) * chunk;
    if(end_rho > length) end_rho = length;
    int chunksize = end_rho - begin_rho;
    if(chunksize<1) continue;

    int begin_vel  = begin_rho * dim_s; 
    int end_vel    = end_rho * dim_s; 

    int begin_f    = begin_rho * sim_dir;
    int end_f      = end_rho * sim_dir;

    unsigned tile = (i* tileNum)/div;
    graph.setTileMapping(rho.slice(begin_rho, end_rho), tile);
    graph.setTileMapping(vel.slice(begin_vel, end_vel), tile);
    graph.setTileMapping(f.slice(begin_f, end_f), tile);
    graph.setTileMapping(ftmp.slice(begin_f, end_f), tile);
    //graph.setTileMapping(prho.slice(begin_rho, end_rho), tile);
    //graph.setTileMapping(pvel.slice(begin_vel, end_vel), tile);
  }
	*/

  int chunk = length/div1;
  ComputeSet kemacroCS = graph.addComputeSet("kemacroCS");
  for (int i = 0; i < div1 ; i++)
  {
    int begin_rho  = i * chunk;
    int end_rho    = (i + 1) * chunk;
    if(end_rho > length) end_rho = length;
    int chunksize = end_rho - begin_rho;
    if(chunksize<1) continue;

    int begin_vel  = begin_rho * dim_s; 
    int end_vel    = end_rho * dim_s; 

    int begin_f    = begin_rho * sim_dir;
    int end_f      = end_rho * sim_dir;


    auto v = graph.addVertex(kemacroCS,
                            "KemacroscopicVertex",
                            {{"rho", 	rho.slice(begin_rho, end_rho)  	},
                             {"vel", 	vel.slice(begin_vel, end_vel)  	},
                             {"f", 		f.slice(begin_f, end_f)		},
                             {"e", 		e		},
                             {"w", 		w		},
                             {"chunk",	chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v, 	tile);
    graph.setTileMapping(rho.slice(begin_rho, end_rho), tile);
    graph.setTileMapping(vel.slice(begin_vel, end_vel), tile);
    graph.setTileMapping(f.slice(begin_f, end_f), tile);
    graph.setTileMapping(ftmp.slice(begin_f, end_f), tile);
    graph.setTileMapping(prho.slice(begin_rho, end_rho), tile);
    graph.setTileMapping(pvel.slice(begin_vel, end_vel), tile);
  }


  ComputeSet kecollCS = graph.addComputeSet("kecollCS");
  for (int i = 0; i < div1; i++)
  {
    int begin_rho  = i * chunk;
    int end_rho    = (i + 1) * chunk;
    if(end_rho > length) end_rho = length;
    int chunksize = end_rho - begin_rho;
    if(chunksize<1) continue;

    int begin_vel  = begin_rho * dim_s; 
    int end_vel    = end_rho * dim_s; 

    int begin_f    = begin_rho * sim_dir;
    int end_f      = end_rho * sim_dir;

    auto v  = graph.addVertex(kecollCS,
                             "KecollisionVertex",
                             {{"rho", 		rho.slice(begin_rho, end_rho) }, 
                              {"vel", 		vel.slice(begin_vel, end_vel)  },
                              {"f", 		f.slice(begin_f, end_f)		},
                              {"ftmp", 		ftmp.slice(begin_f, end_f)		},
                              {"e", 		e		},
                              {"w", 		w		},
                              {"tau", 		tau		},
                              {"chunk",		chunksize  	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }


// *****************************************************************************************
// streaming in 9 directions : Separate ComputeSets to avoid Race Conditions 
// *****************************************************************************************
  ComputeSet Stream_mid_CS = graph.addComputeSet("Stream_mid_CS");
  dir = 0;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;
    if(chunksize<1) continue;

    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_mid_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_E_CS = graph.addComputeSet("Stream_E_CS");
  dir = 1;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;

    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_E_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_N_CS = graph.addComputeSet("Stream_N_CS");
  dir = 2;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;

    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_N_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_W_CS = graph.addComputeSet("Stream_W_CS");
  dir = 3;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;

    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_W_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_S_CS = graph.addComputeSet("Stream_S_CS");
  dir = 4;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;

    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_S_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_NE_CS = graph.addComputeSet("Stream_NE_CS");
  dir = 5;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;

    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_NE_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_NW_CS = graph.addComputeSet("Stream_NW_CS");
  dir = 6;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;

    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_NW_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }


  ComputeSet Stream_SW_CS = graph.addComputeSet("Stream_SW_CS");
  dir = 7;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;

    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_SW_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }

  ComputeSet Stream_SE_CS = graph.addComputeSet("Stream_SE_CS");
  dir = 8;
  margin =  nx*D2Q9eint[sim_dir + dir] + D2Q9eint[dir];
  //printf("margin=%d\n",margin);
  for (int i = 0; i < div1 ; i++)
  {
	int begin = i*chunk;
    int end   = (i+1)*chunk;
    if(begin < (0      - margin)) begin = 0      - margin;
    if(end   > (length - margin)) end   = length - margin;
    int chunksize = end - begin;
    if(chunksize<1) continue;
    int padding = begin;
    int begin_f = (begin + margin)*sim_dir;
    int end_f   = (end   + margin)*sim_dir;
	begin *= sim_dir;
	end   *= sim_dir;
//    printf("DIR=%1d |st=%4d,ed=%4d |st_f=%4d,ed_f=%4d |size=%4d padding=%4d\n",dir,begin,end,begin_f,end_f,chunksize,padding);

    auto v  = graph.addVertex(Stream_SE_CS,
                             "KestreamingVertex",
                             {{"f", 		f.slice(begin_f, end_f)			},
                              {"ftmp", 		ftmp.slice(begin, end)			},
                              {"eint", 		eint		},
                              {"padding", 	padding },
                              {"dir", 	 	dir		},
                              {"nx", 		nx		},
                              {"ny", 		ny		},
                              {"chunk",		chunksize 	}});

    unsigned tile = (i* tileNum)/div1;
    graph.setTileMapping(v , 	 tile);
  }



// *****************************************************************************************
// *****************************************************************************************
// *****************************************************************************************

  ComputeSet keboundCS = graph.addComputeSet("keboundCS");
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
    auto v = graph.addVertex(keboundCS,
                             "KeboundaryVertex",
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

/*
  ComputeSet kepresaveCS = graph.addComputeSet("kepresaveCS");
  //for (int i = 0; i < sublen; i++)
  {
    //auto rho_p  = rho.slice		(i*div, 		(i+1)*div);
    //auto vel_p  = vel.slice		(i*div*dim_s, 	(i+1)*div*dim_s);
    //auto prho_p = prho.slice	(i*div, 		(i+1)*div);
    //auto pvel_p = pvel.slice	(i*div*dim_s, 	(i+1)*div*dim_s);

    auto v = graph.addVertex(kepresaveCS,
                             "KepresaveVertex",
                            {{"rho", 	rho  	},
                             {"prho", 	prho  	},
                             {"vel", 	vel  	},
                             {"pvel", 	pvel  	},
                             {"sublen", sublen	}});

    graph.setTileMapping(v, 	 5);
   
    //graph.setTileMapping(v, 	 (i * tileNum) / sublen);
    //graph.setTileMapping(prho_p, (i * tileNum) / sublen);
    //graph.setTileMapping(pvel_p, (i * tileNum) / sublen);
  }
  //Need reduce of Error!!!!
  ComputeSet keErrorCS = graph.addComputeSet("keErrorCS");
  for (int i = 0; i < sublen; i++)
  {
    auto rho_p  = rho.slice		(i*div, 		i*div + div);
    auto prho_p = prho.slice	(i*div, 		i*div + div);

    auto v = graph.addVertex(keErrorCS,
                             "KeErrorVertex",
                            {{"rho", 	rho		},
                             {"prho", 	prho	},
                             {"sublen", sublen	}});

    graph.setTileMapping(v, 	 (i * tileNum) / sublen);
  }
  */	

//  return Sequence(Execute(kemacroCS), Execute(kecollCS), Execute(keboundCS));
  return Sequence(Execute(kemacroCS), Execute(kecollCS), \
			  	Execute(Stream_mid_CS), Execute(Stream_E_CS), Execute(Stream_N_CS),Execute(Stream_W_CS), \
			  	Execute(Stream_S_CS), Execute(Stream_NE_CS), Execute(Stream_NW_CS), Execute(Stream_SW_CS), \
			  	Execute(Stream_SE_CS),  \
		          Execute(keboundCS));
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
  INPUT.maxitr = 100000; INPUT.rho0 = 1.0; INPUT.u0 = 0.1;
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

  //============================================================================================
  // Create and initialize variables on Host(CPU)
  std::cout << "Create and initialize variables on Host(CPU)\n";
  vector<float> hrho(INPUT.length, 0);
  vector<float> hvel(INPUT.length * dim_s, 0);
  vector<float> hf(INPUT.length * sim_dir, 0);

  for (unsigned tid = 0; tid < INPUT.length; tid++)
  {
    hrho[tid] = INPUT.rho0;
    hvel[tid * dim_s] = 0.0;
    hvel[tid * dim_s + 1] = 0.0;
    // equilbrium 값으로 초기화
    for (int l = 0; l < sim_dir; l++)
    {
      float edotu = D2Q9.e[0][l] * hvel[tid] + D2Q9.e[1][l] * hvel[tid + 1];
      float usquare = hvel[tid] * hvel[tid] + hvel[tid + 1] * hvel[tid + 1];
      hf[tid * sim_dir + l] = D2Q9.w[l] * hrho[tid] * (1.0 + 3.0 * edotu + 4.5 * edotu * edotu - 1.5 * usquare);
    }
  }
  //============================================================================================    ​
  //============================================================================================
  // Add Variables to the graph
  std::cout<< "Add Variables to the graph\n";
  Tensor f 		= graph.addVariable(FLOAT, {unsigned(INPUT.length) * sim_dir}, "f");
  Tensor ftmp 	= graph.addVariable(FLOAT, {unsigned(INPUT.length) * sim_dir}, "ftmp");
  Tensor rho 	= graph.addVariable(FLOAT, {unsigned(INPUT.length)}, "rho");
  Tensor prho 	= graph.addVariable(FLOAT, {unsigned(INPUT.length)}, "prho");
  Tensor vel 	= graph.addVariable(FLOAT, {unsigned(INPUT.length) * dim_s}, "vel");
  Tensor pvel 	= graph.addVariable(FLOAT, {unsigned(INPUT.length) * dim_s}, "pvel");
  Tensor e 		= graph.addVariable(FLOAT, {dim_s * sim_dir}, "e");
  Tensor eint 	= graph.addVariable(INT, {dim_s * sim_dir}, "eint");
  Tensor w 		= graph.addVariable(FLOAT, {sim_dir}, "w");
  Tensor oppi 	= graph.addVariable(INT, {sim_dir}, "oppi");
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
  graph.createHostWrite("H2D-rho", rho);
  graph.createHostWrite("H2D-prho", prho);
  graph.createHostWrite("H2D-vel", vel);
  graph.createHostWrite("H2D-pvel", pvel);

  graph.createHostWrite("H2D-e", e);
  graph.createHostWrite("H2D-eint", eint);
  graph.createHostWrite("H2D-w", w);
  graph.createHostWrite("H2D-oppi", oppi);

  graph.createHostRead("D2H-vel", vel);
  graph.createHostRead("D2H-rho", rho);
  graph.createHostRead("D2H-ftmp", ftmp);
  //============================================================================================

  //============================================================================================
  // Copy Host data via the write handle to variables on the Device
  Sequence mainProg, iterProg;

  Program single_step = buildMainProgram(graph, rho, vel, prho, pvel, f, ftmp, e, w, eint, oppi, INPUT.rho0, INPUT.u0, D2Q9.tau, INPUT.length, INPUT.nx, INPUT.ny);

  mainProg.add(single_step);
  iterProg.add(Repeat(INPUT.maxitr, mainProg));
  
  Engine engine(graph, iterProg);
  engine.load(device);

  engine.writeTensor("H2D-f", 		hf.data());
  engine.writeTensor("H2D-ftmp", 	hf.data());
  engine.writeTensor("H2D-rho", 	hrho.data());
  engine.writeTensor("H2D-prho", 	hrho.data());
  engine.writeTensor("H2D-vel", 	hvel.data());
  engine.writeTensor("H2D-pvel", 	hvel.data());

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
  engine.readTensor("D2H-vel", hvel.data());
  engine.readTensor("D2H-rho", hrho.data());
  engine.readTensor("D2H-ftmp", hf.data());

	/*
  printf("====================== RHO ====================\n");
	for (int i = INPUT.ny-1; i >= 0; i--)
  //for (int i = 0; i < INPUT.ny; i++) 
  {
    for (int j = 0; j < INPUT.nx; j++) 
    {
	  int id = j + i * INPUT.nx;
	  printf("%10.7f ",hrho[id]);
    }
    printf("\n");
  }
  printf("\n");
  printf("\n");

  printf("====================== VEL ====================\n");
	for (int i = INPUT.ny-1; i >= 0; i--)
  //for (int i = 0; i < INPUT.ny; i++) 
  {
    for (int j = 0; j < INPUT.nx; j++) 
    {
	  int id = j + i * INPUT.nx;
	  printf("%10.7f ",hvel[id*dim_s]/INPUT.u0);
    }
    printf("\n");
  }

  printf("\n");
  printf("\n");
  printf("===================== Ftmp ====================\n");
  for (int l=0; l<sim_dir; l++)
  {
 		for (int i = INPUT.ny-1; i >= 0; i--)
      //for (int i = 0; i < INPUT.ny; i++)
      {
          for (int j = 0; j < INPUT.nx; j++)
          {
              int id = j + i * INPUT.nx;
              printf("%10.7f ",hf[id*sim_dir+l]);
          }
          printf("\n");
      }
      printf("\n");
  }

  for (int i = INPUT.ny-1; i >= 0; i--)
  {
	int id = INPUT.nx / 2 + i * INPUT.nx;
	printf("%10.7f, %10.7f \n",(i + 0.5)/INPUT.ny, hvel[id*dim_s]/INPUT.u0);
  }
	*/

  for (int i = 0 ; i< INPUT.ny; i++)
  {
	int id = INPUT.nx / 2 + i * INPUT.nx;
	printf("%10.7f, %10.7f \n",(i + 0.5)/INPUT.ny, hvel[id*dim_s]/INPUT.u0);
  }


  return 0;
}
