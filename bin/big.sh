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

# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_90.gds -out benchout/out_sparc_exu_alu_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_90.gds -out benchout/out_sparc_exu_alu_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_90.gds -out benchout/out_sparc_exu_alu_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_90.gds -out benchout/out_sparc_exu_alu_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_702.gds -out benchout/out_sparc_exu_alu_702.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_702.gds -out benchout/out_sparc_exu_alu_702.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_702.gds -out benchout/out_sparc_exu_alu_702.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_702.gds -out benchout/out_sparc_exu_alu_702.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_70.gds -out benchout/out_sparc_exu_byp_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_70.gds -out benchout/out_sparc_exu_byp_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_70.gds -out benchout/out_sparc_exu_byp_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_70.gds -out benchout/out_sparc_exu_byp_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_80.gds -out benchout/out_sparc_exu_byp_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_80.gds -out benchout/out_sparc_exu_byp_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_80.gds -out benchout/out_sparc_exu_byp_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_80.gds -out benchout/out_sparc_exu_byp_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_85.gds -out benchout/out_sparc_exu_byp_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_85.gds -out benchout/out_sparc_exu_byp_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_85.gds -out benchout/out_sparc_exu_byp_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_85.gds -out benchout/out_sparc_exu_byp_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_90.gds -out benchout/out_sparc_exu_byp_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_90.gds -out benchout/out_sparc_exu_byp_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_90.gds -out benchout/out_sparc_exu_byp_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_byp_90.gds -out benchout/out_sparc_exu_byp_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_70.gds -out benchout/out_sparc_exu_div_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_70.gds -out benchout/out_sparc_exu_div_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_70.gds -out benchout/out_sparc_exu_div_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_70.gds -out benchout/out_sparc_exu_div_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_80.gds -out benchout/out_sparc_exu_div_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_80.gds -out benchout/out_sparc_exu_div_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_80.gds -out benchout/out_sparc_exu_div_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_80.gds -out benchout/out_sparc_exu_div_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_85.gds -out benchout/out_sparc_exu_div_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_85.gds -out benchout/out_sparc_exu_div_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_85.gds -out benchout/out_sparc_exu_div_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_85.gds -out benchout/out_sparc_exu_div_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_90.gds -out benchout/out_sparc_exu_div_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_90.gds -out benchout/out_sparc_exu_div_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_90.gds -out benchout/out_sparc_exu_div_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_div_90.gds -out benchout/out_sparc_exu_div_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_70.gds -out benchout/out_sparc_ffu_ctl_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_70.gds -out benchout/out_sparc_ffu_ctl_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_70.gds -out benchout/out_sparc_ffu_ctl_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_70.gds -out benchout/out_sparc_ffu_ctl_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_80.gds -out benchout/out_sparc_ffu_ctl_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_80.gds -out benchout/out_sparc_ffu_ctl_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_80.gds -out benchout/out_sparc_ffu_ctl_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_80.gds -out benchout/out_sparc_ffu_ctl_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_85.gds -out benchout/out_sparc_ffu_ctl_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_85.gds -out benchout/out_sparc_ffu_ctl_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_85.gds -out benchout/out_sparc_ffu_ctl_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_85.gds -out benchout/out_sparc_ffu_ctl_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_90.gds -out benchout/out_sparc_ffu_ctl_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_90.gds -out benchout/out_sparc_ffu_ctl_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_90.gds -out benchout/out_sparc_ffu_ctl_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_ffu_ctl_90.gds -out benchout/out_sparc_ffu_ctl_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_70.gds -out benchout/out_sparc_exu_alu_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_70.gds -out benchout/out_sparc_exu_alu_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_70.gds -out benchout/out_sparc_exu_alu_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_70.gds -out benchout/out_sparc_exu_alu_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_80.gds -out benchout/out_sparc_exu_alu_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_80.gds -out benchout/out_sparc_exu_alu_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_80.gds -out benchout/out_sparc_exu_alu_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_exu_alu_80.gds -out benchout/out_sparc_exu_alu_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_70.gds -out benchout/out_sparc_mul_top_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_70.gds -out benchout/out_sparc_mul_top_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_70.gds -out benchout/out_sparc_mul_top_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_70.gds -out benchout/out_sparc_mul_top_70.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_80.gds -out benchout/out_sparc_mul_top_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_80.gds -out benchout/out_sparc_mul_top_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_80.gds -out benchout/out_sparc_mul_top_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_80.gds -out benchout/out_sparc_mul_top_80.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_85.gds -out benchout/out_sparc_mul_top_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_85.gds -out benchout/out_sparc_mul_top_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_85.gds -out benchout/out_sparc_mul_top_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
# ./OpenMPL -shape POLYGON -in big/sparc_mul_top_85.gds -out benchout/out_sparc_mul_top_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

./OpenMPL -shape POLYGON -in big/sparc_mul_top_90.gds -out benchout/out_sparc_mul_top_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 1 -use_stitch -gen_stitch -record 1
./OpenMPL -shape POLYGON -in big/sparc_mul_top_90.gds -out benchout/out_sparc_mul_top_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 1 -use_stitch -gen_stitch -record 1
./OpenMPL -shape POLYGON -in big/sparc_mul_top_90.gds -out benchout/out_sparc_mul_top_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 1 -use_stitch -gen_stitch -record 1
./OpenMPL -shape POLYGON -in big/sparc_mul_top_90.gds -out benchout/out_sparc_mul_top_90.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 1 -use_stitch -gen_stitch -record 1

# ./OpenMPL \
#     -shape POLYGON \
#     -in big/sparc_exu_alu_85.gds \
#     -out benchout/cohen_sim_c1_sti.gds \
#     -coloring_distance 100 \
# 	-uncolor_layer 15 \
# 	-uncolor_layer 16 \
#     -color_num 3 \
#     -algo BACKTRACK \
#     -thread_num 1\
#     -use_stitch \
#     -gen_stitch \
#     -record 1\