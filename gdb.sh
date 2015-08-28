#########################################################################
# File Name: gdb.sh
# Author: Yibo Lin
# mail: yibolin@utexas.edu
# Created Time: Sat 23 May 2015 05:22:02 PM CDT
#########################################################################
#!/bin/bash

color_num=4
simplify_level=2
thread_num=8
algo=LP # BACKTRACK or ILP or LP
#benchmark="output_20x20-flat.gds"
#benchmark="Via2_local_precolor.gds"
#benchmark="via2_local_precolor.gds"
benchmark="sim_s3.gds"
#benchmark="mpl_sim_s5_c3_algo0.gds" # output from mpl 

if [[ $benchmark == output_* ]]; then 
	benchmark_dir="/home/local/eda03/shared_benchmarks/imec_7nm/dpt_array"
elif [[ $benchmark == Via2_local_precolor* ]]; then
	benchmark_dir="/home/local/eda03/shared_benchmarks/imec_7nm"
elif [[ $benchmark == via2_local_precolor* ]]; then
	benchmark_dir="./bin/bench"
elif [[ $benchmark == sim_* ]]; then
    benchmark_dir="${BENCHMARKS_DIR}/ISCAS_sim"
elif [[ $benchmark == total_* ]]; then
    benchmark_dir="${BENCHMARKS_DIR}/ISCAS_total"
elif [[ $benchmark == mpl_* ]]; then 
    benchmark_dir="${BENCHMARKS_DIR}/mpl_output/ISCAS_sim"
fi

#output="${benchmark%.*}-out.gds"
output=""

if [[ $benchmark == output_* ]]; then 

# this parameter works for output_1x1.gds 
gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args ./bin/SimpleMPL \
    -shape "RECTANGLE" \
	-in "${benchmark_dir}/${benchmark}" \
	-out "${output}" \
	-uncolor_layer 208 \
	-uncolor_layer 209 \
	-uncolor_layer 210 \
	-uncolor_layer 211 \
	-uncolor_layer 216 \
	-path_layer 207 \
	-color_num ${color_num} \
	-simplify_level ${simplify_level} \
	-thread_num ${thread_num} \
	-algo ${algo} #\
#	-verbose

elif [[ $benchmark == Via2_local_precolor* ]]; then

# this parameter works for Via2_local_precolor.gds
gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args ./bin/SimpleMPL \
    -shape "RECTANGLE" \
	-in "${benchmark_dir}/${benchmark}" \
	-out "${output}" \
	-uncolor_layer 100 \
	-precolor_layer 201 \
	-precolor_layer 202 \
	-precolor_layer 203 \
	-coloring_distance 130 \
	-color_num ${color_num} \
	-simplify_level ${simplify_level} \
	-thread_num ${thread_num} \
	-algo ${algo} \
	-verbose

elif [[ $benchmark == via2_local_precolor* ]]; then

# this parameter works for via2_local_precolor.gds
gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args ./bin/SimpleMPL \
    -shape "RECTANGLE" \
	-in "${benchmark_dir}/${benchmark}" \
	-out "${output}" \
	-uncolor_layer 100 \
	-precolor_layer 201 \
	-precolor_layer 202 \
	-precolor_layer 203 \
	-coloring_distance 300 \
	-color_num ${color_num} \
	-simplify_level ${simplify_level} \
	-thread_num ${thread_num} \
	-algo ${algo} \
	-verbose

elif [[ $benchmark == sim_* ]]; then

# this parameter works for sim_c1
    if [[ $color_num == 3 ]]; then
        coloring_distance=120
    else
        coloring_distance=160
    fi
gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args \
    ./bin/SimpleMPL \
    -shape "POLYGON" \
	-in "${benchmark_dir}/${benchmark}" \
	-out "${output}" \
	-uncolor_layer 1 \
    -uncolor_layer 101 \
	-coloring_distance ${coloring_distance} \
	-color_num ${color_num} \
	-simplify_level ${simplify_level} \
	-thread_num ${thread_num} \
	-algo ${algo} \
    -dbg_comp_id 583 \
    -verbose

elif [[ $benchmark == mpl_* ]]; then 

# this parameter works for mpl_sim_c9
gdb \
	-ex "source ${LIBRARIES_DIR}/gdb_container.sh" \
	--args \
    ./bin/SimpleMPL \
    -shape "POLYGON" \
	-in "${benchmark_dir}/${benchmark}" \
	-out "${output}" \
    -precolor_layer 3 \
    -precolor_layer 4 \
    -precolor_layer 5 \
    -uncolor_layer 0 \
	-coloring_distance 120 \
	-color_num ${color_num} \
	-simplify_level ${simplify_level} \
	-thread_num ${thread_num} \
	-algo ${algo} \
    -dbg_comp_id 1686000 \
    -verbose

fi

