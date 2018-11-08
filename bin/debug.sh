file OpenMPL
set args -shape POLYGON -in bench/sim_c1.gds -out benchout/sim_c1_out.gds -simplify_level 2 -coloring_distance 120 -uncolor_layer 1 -uncolor_layer 101  -color_num 3 -algo ILP  -thread_num 1 -stitch -verbose
break LayoutDBPolygon.cpp:280
break LayoutDBPolygon.cpp:300
break LayoutDBPolygon.cpp:331
r
