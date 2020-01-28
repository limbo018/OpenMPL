import subprocess
import os,glob
files = [f for f in glob.glob('./cell/*.gds') if os.path.isfile(f)]
for i in files:
    for j in files:
        for k in range(2):
            print(i+" "+j+ " "+str(k))
            if k == 0:
                subprocess.run(["/uac/gds/wli/big_home/repository/OpenMPL/bin/OpenMPL","-shape","POLYGON","-in", i,"-in2",j ,"-out",
                "benchout/sim_c1_sti.gds","-coloring_distance","100","-uncolor_layer","15","-uncolor_layer","16",
                "-color_num","3","-algo","ILP","-thread_num","8","-use_stitch","-gen_stitch","-record","2"])
            else:
                subprocess.run(["/uac/gds/wli/big_home/repository/OpenMPL/bin/OpenMPL","-shape","POLYGON","-in", i,"-in2",j ,"-in2_flip","-out",
                "benchout/sim_c1_sti.gds","-coloring_distance","100","-uncolor_layer","15","-uncolor_layer","16",
                "-color_num","3","-algo","ILP","-thread_num","8","-use_stitch","-gen_stitch","-record","2"])
