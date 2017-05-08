#########################################################################
# File Name: run.sh
# Author: Yibo Lin
# mail: yibolin@utexas.edu
# Created Time: Thu 27 Aug 2015 09:48:41 AM CDT
#########################################################################
#!/bin/bash

if [ $# -lt 1 ]; then
    echo "$0: usage: source run.sh < ISCAS_sim|ISCAS_total > < mpl|nompl > <suffix (optional)>"
    return
fi

benchmark_suit_hint=$1
color_num_map=(3 4)
coloring_distance_map=(120 160) # coloring_distance is paired with color_num
algo_map=("LP") # BACKTRACK, ILP, LP, MIS
simplify_level=3
thread_num=8
use_mpl_output="false"
suffix=""

if [ $# -ge 2 ]; then
    if [ $2 == "mpl" ]; then
        use_mpl_output="true"
    fi
fi 

if [ $# -ge 3 ]; then
    suffix=".$3"
fi 

# report directory is determined with date 
current_date=$(date +'%m-%d-%Y')
rpt_dir="rpt_${current_date}"
mkdir -p ${rpt_dir}

shape="POLYGON"
if [ $use_mpl_output == "true" ]; then
    simplify_level=0 # forbid simplification 
    thread_num=1 # use single thread 
    insert1="/mpl_output" # for benchmark_dir
    insert2="mpl_" # for benchmark 
    precolor_layer_map=(3 4 5 6)
    uncolor_layer_map=() # dummy 
else 
    insert1="" # for benchmark_dir
    insert2="" # for benchmark 
    precolor_layer_map=()
    uncolor_layer_map=(1 101)
fi
# run ISCAS_sim
if [ $benchmark_suit_hint == "ISCAS_sim" ]; then 
    benchmark_dir="${BENCHMARKS_DIR}${insert1}/ISCAS_sim"
    benchmark_map=(\
        "${insert2}sim_c1" \
        "${insert2}sim_c2" \
        "${insert2}sim_c3" \
        "${insert2}sim_c4" \
        "${insert2}sim_c5" \
        "${insert2}sim_c6" \
        "${insert2}sim_c7" \
        "${insert2}sim_c8" \
        "${insert2}sim_c9" \
        "${insert2}sim_c10" \
        "${insert2}sim_s1" \
        "${insert2}sim_s2" \
        "${insert2}sim_s3" \
        "${insert2}sim_s4" \
        "${insert2}sim_s5" \
        )
#    benchmark_map=(\
#        "${insert2}sim_c10" \
#        )
# run ISCAS_total
elif [ $benchmark_suit_hint == "ISCAS_total" ]; then 
    benchmark_dir="${BENCHMARKS_DIR}${insert1}/ISCAS_total"
    benchmark_map=(\
        "${insert2}total_c1" \
        "${insert2}total_c2" \
        "${insert2}total_c3" \
        "${insert2}total_c4" \
        "${insert2}total_c5" \
        "${insert2}total_c6" \
        "${insert2}total_c7" \
        "${insert2}total_c8" \
        "${insert2}total_c9" \
        "${insert2}total_c10" \
        "${insert2}total_s1" \
        "${insert2}total_s2" \
        "${insert2}total_s3" \
        "${insert2}total_s4" \
        "${insert2}total_s5" \
        )
fi

for i in "${!color_num_map[@]}"; do
    color_num="${color_num_map[$i]}"
    coloring_distance="${coloring_distance_map[$i]}"
#    if [[ $color_num == 4 ]]; then
#        continue
#    fi

    for algo in "${algo_map[@]}"; do

        for benchmark in "${benchmark_map[@]}"; do
            cmd="./bin/SimpleMPL -shape ${shape} -simplify_level ${simplify_level} -thread_num ${thread_num}"
            for precolor_layer in "${precolor_layer_map[@]}"; do
                cmd="${cmd} -precolor_layer ${precolor_layer}"
            done
            for uncolor_layer in "${uncolor_layer_map[@]}"; do
                cmd="${cmd} -uncolor_layer ${uncolor_layer}"
            done
            cmd="${cmd} -algo ${algo}"
            cmd="${cmd} -color_num ${color_num} -coloring_distance ${coloring_distance}"

            if [ $use_mpl_output == "false" ]; then
                cmd="${cmd} -in ${benchmark_dir}/${benchmark}.gds"
                # no output 
                #output="${benchmark}-out.gds"
                #cmd="${cmd} -out ${output}"

                # run command 
                echo "time (${cmd}) > ${rpt_dir}/${benchmark}_c${color_num}_${algo}${suffix}.rpt 2>&1"
                time (${cmd}) > ${rpt_dir}/${benchmark}_c${color_num}_${algo}${suffix}.rpt 2>&1
            else 
                # 0 for mpl ILP, 1 for mpl SDP
                for mpl_algo in `seq 0 1`; do 
                    cmd="${cmd} -in ${benchmark_dir}/${benchmark}_c${color_num}_algo${mpl_algo}.gds"
                    echo "time (${cmd}) > ${rpt_dir}/${benchmark}_c${color_num}_${algo}_mpl${mpl_algo}${suffix}.rpt 2>&1"
                    time (${cmd}) > ${rpt_dir}/${benchmark}_c${color_num}_${algo}_mpl${mpl_algo}${suffix}.rpt 2>&1
                done
            fi 
        done
    done

done
