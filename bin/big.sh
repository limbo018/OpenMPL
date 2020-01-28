#!
# ========================================================================
#                         SimpleMPL Usage                         
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
#  -algo (BACKTRACK)            algorithm type < ILP|BACKTRACK|LP|SDP|DL|ILP_UPDATED >
#  -shape (RECTANGLE)           shape mode < RECTANGLE|POLYGON >
#  -verbose (false)             toggle controlling screen messages
#  -use_stitch                  use stitch to avoid conflict
#  -gen_stitch                  generate stitch candidate
#  -dbg_comp_id                 debug component id
#  -weight_stitch               a floating point number indicating the weight of stitch
# ========================================================================

# if the benchmark contains polygon shapes, -shape must be set to POLYGON;
# otherwise, set -shape to RECTANGLE is more memory efficient
mkdir -p benchout
./OpenMPL  \
    -shape POLYGON \
    -in big/sparc_mul_top_85.gds \
    -out benchout/sim_c1_sti.gds  \
    -coloring_distance 100 \
	-uncolor_layer 15 \
	-uncolor_layer 16 \
    -color_num 3 \
    -algo DL \
    -thread_num 8 \
    -use_stitch \
    -gen_stitch \
    -record 0\
