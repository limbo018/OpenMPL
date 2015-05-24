#########################################################################
# File Name: gdb.sh
# Author: Yibo Lin
# mail: yibolin@utexas.edu
# Created Time: Sat 23 May 2015 05:22:02 PM CDT
#########################################################################
#!/bin/bash

benchmark_dir="${BENCHMARKS_DIR}/imec_7nm/dpt_array"
benchmark="output_1x1"

gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args ./bin/SimpleMPL \
	-in "${benchmark_dir}/${benchmark}-flat.gds" \
	-out "${benchmark}-out.gds" \
	-uncolor_layer 208 \
	-uncolor_layer 209 \
	-uncolor_layer 210 \
	-uncolor_layer 211 \
	-uncolor_layer 216 \
	-path_layer 207 \
	-color_num 4
