#!

# ========================================================================
#                         SimpleMPL 1.X Usage                         
#  "-in"                 : followed by input gds name. 
#  "-output"             : followed by output gds name (default: "output.gds"). 
#  "-coloring_distance"  : followed by floating point number indicating number of coloring distance in micron. 
#  "-color_num"          : followed by integer indicating number of masks (colors). 
#  "-thread_num"         : followed by integer, maximum thread number
#  "-precolor_layer"     : followed by an integer, pre-coloring layer
#  "-uncolor_layer"      : followed by an integer, layer for coloring
#  "-path_layer"         : followed by an integer, layer for conflict edges
# ========================================================================


./SimpleMPL  -in bench/via2_local_precolor.gds  -out out.gds  -coloring_distance 0.13 \
             -uncolor_layer  100 \
             -precolor_layer 201  -precolor_layer 202  -precolor_layer 203

