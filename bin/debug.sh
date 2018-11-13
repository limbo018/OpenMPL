file OpenMPL
set args -shape POLYGON -in bench/sim_c1.gds -out benchout/sim_c1_SDP -simplify_level 3  -coloring_distance 120 -uncolor_layer 1 -uncolor_layer 101 -color_num 3 -algo SDP -thread_num 1 -stitch 
r
