#########################################################################
# File Name: gdb.sh
# Author: Yibo Lin
# mail: yibolin@utexas.edu
# Created Time: Sat 23 May 2015 05:22:02 PM CDT
#########################################################################
#!/bin/bash

benchmark_dir="/home/usr1/shared_benchmarks/imec_7nm"
benchmark="Via2_local_precolor.gds"
color_num=3

output="${benchmark%.*}-out.gds"

# this parameter works for via2.gds 
#gdb \
#	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
#	--args ./bin/SimpleMPL \
#	-in "${benchmark_dir}/${benchmark}" \
#	-out "${output}" \
#	-uncolor_layer 208 \
#	-uncolor_layer 209 \
#	-uncolor_layer 210 \
#	-uncolor_layer 211 \
#	-uncolor_layer 216 \
#	-path_layer 207 \
#	-color_num ${color_num}

gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args ./bin/SimpleMPL \
	-in "${benchmark_dir}/${benchmark}" \
	-out "${output}" \
	-uncolor_layer 100 \
	-precolor_layer 201 \
	-precolor_layer 202 \
	-precolor_layer 203 \
	-coloring_distance 1300 \
	-color_num ${color_num}

#./bin/SimpleMPL \
#	-in "${benchmark_dir}/${benchmark}-flat.gds" \
#	-out "${benchmark}-out.gds" \
#	-uncolor_layer 208 \
#	-uncolor_layer 209 \
#	-uncolor_layer 210 \
#	-uncolor_layer 211 \
#	-uncolor_layer 216 \
#	-path_layer 207 \
#	-color_num ${color_num} > log

#$LIBRARIES_DIR/memusg \
#	./bin/SimpleMPL \
#	-in "${benchmark_dir}/${benchmark}-flat.gds" \
#	-out "${benchmark}-out.gds" \
#	-uncolor_layer 208 \
#	-uncolor_layer 209 \
#	-uncolor_layer 210 \
#	-uncolor_layer 211 \
#	-uncolor_layer 216 \
#	-path_layer 207 \
#	-color_num ${color_num}
