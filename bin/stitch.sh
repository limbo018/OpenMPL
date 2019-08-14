#!

# ========================================================================
#                         SimpleMPL 1.X Usage                         
#  -help (false)                toggle printing help message
#  -in                          input gds file name
#  -out ()                      output gds file name
#  -coloring_distance (0)       a floating point number indicating number of coloring distance in nanometer
#  -color_num                   an integer indicating number of masks (colors) < 3|4 >
#  -simplify_level (3)          an integer indicating graph simplification level < 0|1|2|3 >
#  -thread_num (1)              an integer indicating maximum thread number
#  -path_layer                  an integer indicating layer for conflict edges
#  -precolor_layer              an integer indicating layer for pre-colored patterns
#  -uncolor_layer               an integer indicating layer for coloring
#  -algo (BACKTRACK)            algorithm type < ILP|BACKTRACK|LP|SDP|DL >
#  -shape (RECTANGLE)           shape mode < RECTANGLE|POLYGON >
#  -verbose (false)             toggle controlling screen messages
#  -stitch (false)              toggle controlling stitch insertion
#  -use_stitch
#  gen_stitch
#  -dbg_comp_id (4294967295)    debug component id
# ========================================================================

#-uncolor_layer  100  \ 
#-uncolor_layer 201  \
#-uncolor_layer 202  \
#-uncolor_layer 203 \
# if the benchmark contains polygon shapes, -shape must be set to POLYGON;
# otherwise, set -shape to RECTANGLE is more memory efficient
./OpenMPL  \
    -shape POLYGON \
    -in bench/sim_c5.gds \
    -out benchout/sim_c1_sti.gds  \
    -coloring_distance 120 \
	-uncolor_layer 1 \
	-uncolor_layer 101 \
    -color_num 3 \
    -algo DL\
    -thread_num 8 \
    -use_stitch \
    gen_stitch
