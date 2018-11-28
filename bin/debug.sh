file OpenMPL
set args -shape POLYGON -in bench/total_c5.gds -out benchout/total_c5_SDP -simplify_level 3  -coloring_distance 120 -uncolor_layer 1 -uncolor_layer 101 -color_num 3 -algo SDP -thread_num 1 -stitch 
#break SDPColoringCsdp.h:574
r
