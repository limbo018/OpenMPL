# OpenMPL

> **OpenMPL** stands for open multiple patterning lithography framework.
### Pre-requisite

- [Limbo](https://github.com/limbo018/Limbo): require LIMBO_DIR environment variable to the path where Limbo is installed. OpenMPL is based on Limbo library.

### How To Compile

```bash
$ git clone https://github.com/limbo018/OpenMPL.git --recursive
$ mkdir build
$ cd build
$ cmake ..
$ make
$ make install
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

The Csdp API used in OpenMPL has been modified and built for threading safety at high level. 

### How To Execute

```bash
$ cd bin/
$ ./OpenMPL
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
-verbose (false)             toggle controlling screen messages
-stitch (false)              toggle controlling stitch insertion, default is false
-dbg_comp_id (4294967295)    debug component id
-use_stitch 				 use stitch to avoid conflict
gen_stitch					 generate stitch candidate
```

One exmaple : /bin/run.sh.

### Possible Compiler Problems

+ default CFLAGS of boost and gurobi are different in newest version
  + downgrade the boost version

+ ```
  SimpleMPL.cpp:461:5: error: ‘graph_simplification_type’ has no member named ‘set_isVDDGND’
  ```

  + checkout to ***stitch*** branch in your limbo directory

### License

- BSD-3-clause License [[LINK](https://github.com/limbo018/OpenMPL/blob/master/LICENSE)]

### Authors

| Name         | Affiliation         | email                                                     |
| ------------ | ------------------- | --------------------------------------------------------- |
| Yibo Lin     | School of EECS, PKU | [yibolin@pku.edu.cn](mailto:yibolin@pku.edu.cn) |
| Bei Yu       | CSE Dept, CUHK      | [byu@cse.cuhk.edu.hk](mailto:byu@cse.cuhk.edu.hk)         |
| Qi Sun       | CSE Dept, CUHK      | [qsun@cse.cuhk.edu.hk](mailto:qsun@cse.cuhk.edu.hk)       |
| David Z. Pan | ECE Dept, UT Austin | [dpan@ece.utexas.edu](mailto:dpan@ece.utexas.edu)         |


