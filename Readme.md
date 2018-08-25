#SimpleMPL

> **SimpleMPL** stands for simple multiple patterning lithography framework.
### Pre-requisite

- [Limbo](https://github.com/limbo018/Limbo): require LIMBO_DIR environment variable to the path where Limbo is installed. SimpleMPL is based on Limbo library.

### How To Compile

```bash
$ git clone https://github.com/limbo018/SimpleMPL.git
$ cd simplempl/src/mpl/
$ make
```
Here are some optional parameters when make :

```bash
# default DBG (debug) is off
DBG = 0
# default GPROF is off, used to enable runtime profiling
GPROF = 0
# default GUROBI is off 
GUROBI = 0
# default LEMONCBC is off 
LEMONCBC = 0
# default CSDP is off 
CSDP = 0
```

### Features
 * Contact or metal layer decomposition 
 * No stitching
 * Support 3 or 4 coloring 
 * Density control
 * Multi-threading
 * Small memory usage
 * Multiple algorithms: 
     ILP (Gurobi or Lemon CBC), 
     SDP (Csdp API), 
     LP  (Gurobi API)
 * Dancing Links

The Csdp API used in SimpleMPL has been modified and built for threading safety at high level. 

### How To Execute

```bash
$ cd bin/
$ ./SimpleMPL
```

A table of options :

```bash
-help (false)                toggle printing help message
-in                          input gds file name
-out ()                      output gds file name
-coloring_distance (0)       a floating point number indicating number of coloring distance in nanometer
-color_num                   an integer indicating number of masks (colors) < 3|4 >
-simplify_level (3)          an integer indicating graph simplification level < 0|1|2|3 >
-thread_num (1)              an integer indicating maximum thread number
-path_layer                  an integer indicating layer for conflict edges
-precolor_layer              an integer indicating layer for pre-colored patterns
-uncolor_layer               an integer indicating layer for coloring
-algo (BACKTRACK)            algorithm type < ILP|BACKTRACK|LP|SDP >
-shape (RECTANGLE)           shape mode < RECTANGLE|POLYGON >
-verbose (false)             toggle controling screen messages
-dbg_comp_id (4294967295)    debug component id
```

Now Dancing Links is an independent component, and the source code is in /src/dlx. One example is in DL_main.cpp. 

### License

- BSD-3-clause License [[LINK](https://github.com/limbo018/SimpleMPL/blob/master/LICENSE)]

### Authors

| Name         | Affiliation         | email                                                     |
| ------------ | ------------------- | --------------------------------------------------------- |
| Yibo Lin     | ECE Dept, UT Austin | [yibolin@cerc.utexas.edu](mailto:yibolin@cerc.utexas.edu) |
| Bei Yu       | CSE Dept, CUHK      | [byu@cse.cuhk.edu.hk](mailto:byu@cse.cuhk.edu.hk)         |
| David Z. Pan | ECE Dept, UT Austin | [dpan@ece.utexas.edu](mailto:dpan@ece.utexas.edu)         |


