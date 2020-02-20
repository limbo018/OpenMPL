import os

algorithms = ["ILP",
              "SDP",
              "DL",
              "BACKTRACK"]
shape = "POLYGON"
file_path = ""
out = ""
coloring_distance = "100"
layer1 = "15"
layer2 = "16"
color_num = "3"
algo = ""
thread_num = "1"
record = "1"


input_folder = "big/"
output_folder = "benchout/"
file_names = os.listdir(input_folder)


fp = open('experiments.sh', 'w')

for file_name in file_names:
    file_path = input_folder + file_name
    out = output_folder + "out_" + file_name
    for algo in algorithms:
        sh = "./OpenMPL " + "-shape " + shape + " -in " + file_path + " -out " + out + " -coloring_distance " + coloring_distance + " -uncolor_layer " + layer1 + " -uncolor_layer " + layer2 + " -color_num " + color_num + " -algo " + algo + " -thread_num " + thread_num + " -use_stitch -gen_stitch" + " -record " + record
        fp.write(sh)
        fp.write("\n")
    fp.write("\n")

fp.close()
