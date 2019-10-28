# OpenMPL

> **OpenMPL** stands for open multiple patterning lithography framework.

| Stitch Insertion | Graph Simplification | Decomposition |
| ---------------- | -------------------- | ------------- | 
| <img src=/images/stitch-2.png width=150> | <img src=/images/biconnected.png width=250> | <img src=/images/total_c2.gif width=200> |

### Pre-requisite

- [GCC](https://gcc.gnu.org)
    - Recommend GCC 4.8 or later. 
    - Other compilers may also work, but not tested. 

- [CMake](https://cmake.org)
    - Require 3.8.2 or later. 

- [Boost](https://www.boost.org)
    - Require 1.55 or later. 
    - Need to install and visible for linking. 
    - Custom installation path may require to export ```BOOST_ROOT``` for [CMake](https://cmake.org/cmake/help/v3.8/module/FindBoost.html). 

- [Limbo](https://github.com/limbo018/Limbo)
    - Integrated as a git submodule.

- [Gurobi](https://www.gurobi.com) (Optional)
    - ILP solver. 

### Publications

* [Wei Li](https://wadmes.github.io/cv/), [Yuzhe Ma](http://yuzhe630.github.io/), [Qi Sun](http://qisunchn.top/), [Yibo Lin](http://yibolin.com), [Iris Hui-Ru Jiang](http://www.ee.ntu.edu.tw/profile1?id=1060726), [Bei Yu](http://www.cse.cuhk.edu.hk/~byu/index.html), [David Z. Pan](http://users.ece.utexas.edu/~dpan/), 
    "**OpenMPL: An Open Source Layout Decomposer**", 
    IEEE International Conference on ASIC (ASICON), Chongqing, China, Oct. 29–Nov. 1, 2019.
([preprint](https://arxiv.org/pdf/1809.07554v3.pdf))

### How To Compile

```bash
$ git clone https://github.com/limbo018/OpenMPL.git 
$ cd OpenMPL
$ git submodule update --init --recursive
$ mkdir build
$ cd build
$ cmake .. 
$ make
$ make install
```
The default installation path is the repo directory. It can bee changed via 
```
cmake .. -DCMAKE_INSTALL_PREFIX=your_installation_path
```

### Features
 * Contact or metal layer decomposition 
 * No stitching
 * Support 3 or 4 coloring 
 * Density control
 * Multi-threading
 * Small memory usage
 * Multiple algorithms: 
    * ILP (Gurobi API)
    * SDP (Csdp API)
    * LP  (Gurobi API)
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
| Yibo Lin     | School of EECS, PKU | [yibolin@pku.edu.cn](mailto:yibolin@pku.edu.cn)           |
| Bei Yu       | CSE Dept, CUHK      | [byu@cse.cuhk.edu.hk](mailto:byu@cse.cuhk.edu.hk)         |
| Wei Li       | CSE Dept, CUHK      | [werry715@gmail.com](mailto:werry715@gmail.com)           |
| Qi Sun       | CSE Dept, CUHK      | [qsun@cse.cuhk.edu.hk](mailto:qsun@cse.cuhk.edu.hk)       |
| David Z. Pan | ECE Dept, UT Austin | [dpan@ece.utexas.edu](mailto:dpan@ece.utexas.edu)         |


