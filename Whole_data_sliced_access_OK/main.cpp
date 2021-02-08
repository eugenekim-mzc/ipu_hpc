#include <iostream>
#include <stdio.h>
#include <chrono>
#include <vector>
#include <memory>
#include <poplar/Graph.hpp>
#include <poplar/DeviceManager.hpp>
#include <poplar/Engine.hpp>
#include <poputil/TileMapping.hpp>
#define tileNum 1472
#define div 	20
//#define sim_dir			9
#define h              10
#define RED		   "\x1b[31m"
#define RESET	   "\x1b[0m"
using namespace std;
using namespace poplar;
using namespace poplar::program;


Program buildMainProgram(Graph &graph, Tensor oneset, Tensor sliced, int length, int nx, int ny)
{
	/*
  	for(int i=0; i<div; i++)
	{
		int chunk = length/div;
    	unsigned begin = i*chunk;
    	unsigned end   = (i+1)*chunk;
    	unsigned tile = (i*tileNum)/div;
    	graph.setTileMapping(sliced.slice(begin,end), tile);
    	//graph.setTileMapping(oneset.slice(begin,end), tile);
	}
	*/

   	graph.setTileMapping(sliced,0);

  	//============================================================================================
	// Handle whole array at once 
	ComputeSet plus4waysCS = graph.addComputeSet("plus4waysCS");
   	auto v = graph.addVertex(plus4waysCS,
                           	"plus4ways_Vertex",
                           {{"a", 		oneset	},
                           	{"nx",		nx 		},
                           	{"ny",		ny 		},
                           	{"padding",	0  		},
                           	{"chunk",	length 	}});
   	graph.setTileMapping(v,	0);
   	graph.setTileMapping(oneset, 0);

  	//============================================================================================
	// Handle sliced chunk array in loop
  	ComputeSet plus4ways_slicedCS = graph.addComputeSet("plus4ways_slicedCS");
	int div1 = 2;
	int chunk = length/div1;
  	//for(int i=0; i<div1; i++)
  	for(int i=0; i<1; i++)
  	{
    	unsigned begin = i*chunk;
    	unsigned end   = (i+1)*chunk;
		if(end>length) end = length;
		int size = end - begin;
		printf("begin=%d | end=%d | size=%d\n",begin,end,size);
    	for (int tid = 0; tid < size; tid++)
    	{
			int tid1 = tid+begin;
        	int ii = tid1/nx;
        	int jj = tid1%nx;
			if((ii>0 && ii<ny-1) && (jj>0 && jj<nx-1))
			{
				printf("%2d %2d|%2d|%2d %2d %2d %2d\n", ii, jj, tid1, tid-1,tid+1,tid+nx,tid-nx);
			}
		}
    	auto v1 = graph.addVertex(plus4ways_slicedCS,
                            	"plus4ways_Vertex",
                               {{"a", 		sliced.slice(begin,end)},
                             	{"nx",		nx 		},
                             	{"ny",		ny 		},
                             	{"padding",	begin	},
                             	{"chunk",	size  	}});
    	unsigned tile = (i*tileNum)/div1;
    	graph.setTileMapping(v1,	tile);
  	}

  	//============================================================================================
  	return Sequence(Execute(plus4waysCS), Execute(plus4ways_slicedCS));
}



int main()
{
  	int nx =h; int ny =h;
  	int length = nx*ny;

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

  	//============================================================================================
  	// Add codelets to the graph
  	graph.addCodelets("codelets.cpp");

  	//============================================================================================
  	// Create and initialize variables on Host(CPU)
  	std::cout << "Create and initialize variables on Host(CPU)\n";
  	vector<int> honeset(length, 0);
  	vector<int> hsliced(length, 0);

  	printf("================= Initial (whole data) ======================\n");
	for (int i = ny-1; i >= 0; i--)
  	{
    	for (int j = 0; j < nx; j++) 
    	{
	  		int id = j + i * nx;
	  		printf("%3d ",honeset[id]);
    	}
    	printf("\n");
  	}
  	printf("\n");

  	//============================================================================================
  	// Add Variables to the graph
  	std::cout<< "Add Variables to the graph\n";
  	Tensor oneset  	= graph.addVariable(INT, {unsigned(length)}, "oneset");
  	Tensor sliced 	= graph.addVariable(INT, {unsigned(length)}, "sliced");

  	//============================================================================================
  	// Create write $ read handles for variables HOST-Device
  	std::cout << "Create write $ read handles for variables HOST <--> Device\n";
  	graph.createHostWrite("H2D-oneset",	oneset);
  	graph.createHostWrite("H2D-sliced", sliced);
  	graph.createHostRead("D2H-oneset", 	oneset);
  	graph.createHostRead("D2H-sliced", 	sliced);

  	//============================================================================================
  	// Program and Engine 
  	Sequence mainProg;
  	Program one_step = buildMainProgram(graph, oneset, sliced, length, nx, ny);

  	mainProg.add(one_step);
  	Engine engine(graph, mainProg);
  	engine.load(device);

  	//============================================================================================
  	// Copy Host data via the write handle to variables on the Device
  	std::cout << "Copy Host data via the write handle to variables on the Device\n";
    // 초기값은 모두 '0'으로 세팅된 같은 값이 들어감!!!
  	engine.writeTensor("H2D-oneset", 	honeset.data());
  	engine.writeTensor("H2D-sliced", 	honeset.data());
  	//engine.writeTensor("H2D-sliced", 	hsliced.data());

  	//============================================================================================
  	// Start (Engine) Program
    auto start_t = std::chrono::high_resolution_clock::now();
  	std::cout << "Running program\n";
  	engine.run(0);
  	std::cout << "Program complete\n";

    auto end_t = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_t - start_t;
    printf("Elapsed time(IPU) = %fs | Width=%d\n",elapsed.count(), nx);

  	//============================================================================================
	// Copy Device data to Host, Print Result  
  	std::cout << " >>> Copy data from Device via read handle\n";
  	engine.readTensor("D2H-oneset", honeset.data());
  	engine.readTensor("D2H-sliced", hsliced.data());

  	printf("========================= Result ============================\n");
  	printf("======================== Whole set ==========================\n");
 	for (int i = ny-1; i >= 0; i--)
    {
       	for (int j = 0; j < nx; j++)
       	{
           	int id = j + i * nx;
			if(honeset[id]>0) 	printf(RED "%3d " RESET, honeset[id]);
			else       			printf("%3d ", honeset[id]);
      	}
      	printf("\n");
  	}
  	printf("\n");
  	printf("======================== Sliced set =========================\n");
 	for (int i = ny-1; i >= 0; i--)
    {
       	for (int j = 0; j < nx; j++)
       	{
           	int id = j + i * nx;
			if(hsliced[id]>0) 	printf(RED "%3d " RESET, hsliced[id]);
			else       			printf("%3d ", hsliced[id]);
      	}
      	printf("\n");
  	}
  	//============================================================================================

  	return 0;
}
