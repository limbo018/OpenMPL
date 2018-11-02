file OpenMPL
set args -shape POLYGON -in bench/sim_c1.gds -out benchout/sim_c1.gds -coloring_distance 120 -uncolor_layer 1 -uncolor_layer 101 -color_num 3 -algo ILP -thread_num 8 -stitch
break LayoutDBPolygon.cpp:269
break LayoutDBPolygon.cpp:286
r
