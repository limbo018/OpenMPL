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
#  -algo (BACKTRACK)            algorithm type < ILP|BACKTRACK|LP|SDP >
#  -shape (RECTANGLE)           shape mode < RECTANGLE|POLYGON >
#  -verbose (false)             toggle controling screen messages
#  -dbg_comp_id (4294967295)    debug component id
# ========================================================================


# if the benchmark contains polygon shapes, -shape must be set to POLYGON;
# otherwise, set -shape to RECTANGLE is more memory efficient
./OpenMPL  \
    -shape POLYGON \
	-in bench/total_c1.gds \
    -out benchout/total_c1_out.gds \
	-coloring_distance 160 \
    -uncolor_layer 101 \
    -uncolor_layer 1 \
	-color_num 4 \
    -algo LP \
    -thread_num 8 \
	-simplify_level 3
