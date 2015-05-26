#########################################################################
# File Name: gdb.sh
# Author: Yibo Lin
# mail: yibolin@utexas.edu
# Created Time: Sat 23 May 2015 05:22:02 PM CDT
#########################################################################
#!/bin/bash

color_num=3
benchmark="output_4x4-flat.gds"
#benchmark="Via2_local_precolor.gds"

if [[ $benchmark == output_* ]]; then 
	benchmark_dir="/home/usr1/shared_benchmarks/imec_7nm/dpt_array"
elif [[ $benchmark == Via2_local_precolor* ]]; then
	benchmark_dir="/home/usr1/shared_benchmarks/imec_7nm"
fi

#output="${benchmark%.*}-out.gds"
output=""

if [[ $benchmark == output_* ]]; then 

# this parameter works for output_1x1.gds 
gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args ./bin/SimpleMPL \
	-in "${benchmark_dir}/${benchmark}" \
	-out "${output}" \
	-uncolor_layer 208 \
	-uncolor_layer 209 \
	-uncolor_layer 210 \
	-uncolor_layer 211 \
	-uncolor_layer 216 \
	-path_layer 207 \
	-color_num ${color_num}

elif [[ $benchmark == Via2_local_precolor* ]]; then

# this parameter works for Via2_local_precolor.gds
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

fi

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
