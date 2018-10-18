#!/bin/bash

#if [ $# -lt 1 ]; then
#	echo "$0: usage: source gempl.sh < sim | total >"
#	return
#fi

benchmark_suit_hint=total
color_num_map=(3 4)
coloring_distance_map=(120 160)
algo_map=("LP" "ILP" "SDP")
simplify_level=3
thread_num=8
precolor_layer_map=()
uncolor_layer_map=(1 101)
shape="POLYGON"


current_date=$(date +'%m-%d-%Y')
gempl_dir="gempl_${current_date}"
mkdir -p ${gempl_dir}


BENCHMARKS_DIR="/research/byu2/qsun/project/ISCAS_bench"
insert1=""
insert2=""



if [ $benchmark_suit_hint == "sim" ]; then
	benchmark_dir="${BENCHMARKS_DIR}${insert1}/ISCAS_sim"
	benchmark_map=()
elif [ $benchmark_suit_hint == "total" ]; then
	benchmark_dir="${BENCHMARKS_DIR}${insert1}/ISCAS_total"
	benchmark_map=(\
		"total_c1" \
		"total_c2" \
		"total_c3" \
		"total_c4" \
		"total_c5" \
		)
fi

for i in "${!color_num_map[@]}"; do
	color_num="${color_num_map[$i]}"
	coloring_distance="${coloring_distance_map[$i]}"

	for algo in "${algo_map[@]}"; do
		for benchmark in "${benchmark_map[@]}"; do
			cmd="./OpenMPL -shape ${shape} -simplify_level ${simplify_level} -thread_num ${thread_num}"
			for precolor_layer in "${precolor_layer_map[@]}"; do
				cmd="${cmd} -precolor_layer ${precolor_layer}"
			done
			for uncolor_layer in "${uncolor_layer_map[@]}"; do
				cmd="${cmd} -uncolor_layer ${uncolor_layer}" 
			done
			cmd="${cmd} -algo ${algo}"
			cmd="${cmd} -color_num ${color_num} -coloring_distance ${coloring_distance}"
			cmd="${cmd} -in ${benchmark_dir}/${benchmark}.gds"

			# run command
			echo "time (${cmd}) > ${gempl_dir}/${benchmark}_c${color_num}_${algo}.rpt 2>&1"	
			time (${cmd}) > ${gempl_dir}/${benchmark}_c${color_num}_${algo}.rpt 2>&1
		done
	done
done

