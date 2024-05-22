mkdir -p benchout

./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo ILP -thread_num 8 -use_stitch -gen_stitch -record 1
./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo SDP -thread_num 8 -use_stitch -gen_stitch -record 1
./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo DL -thread_num 8 -use_stitch -gen_stitch -record 1
./OpenMPL -shape POLYGON -in big/sparc_exu_alu_85.gds -out benchout/out_sparc_exu_alu_85.gds -coloring_distance 100 -uncolor_layer 15 -uncolor_layer 16 -color_num 3 -algo BACKTRACK -thread_num 8 -use_stitch -gen_stitch -record 1