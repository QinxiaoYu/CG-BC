import commands
import os

# --------------- Instancias DE --------------- #
fo = open("DEresults-CGBC.txt", "w");
fo.write("Inst LB UB %Gap #cols %EF t(s) #nodes t(s) z^* T(s)\n");
fo.close();
for nv in [100, 200, 300, 400, 500, 600]:
    for num_inst in [1, 2, 3]:
        inst = 'dcmst' + str(nv) + "_" + str(num_inst)
        command =  './CG-BC ~/Dropbox/DCMST/allinstances/DE/' + inst + ".dat" + " DEresults-CGBC.txt" + " 2"
        print command
        print "\n"
        os.system(command)
# --------------------------------------------- #

# --------------- Instancias DR --------------- #
fo = open("DRresults-CGBC.txt", "w");
fo.write("Inst LB UB %Gap #cols %EF t(s) #nodes t(s) z^* T(s)\n");
fo.close();
for nv in [100, 200, 300, 400, 500, 600]:
    for num_inst in [1, 2, 3]:
        inst = 'R123_dcmst' + str(nv) + "_" + str(num_inst)
        command =  './CG-BC ~/Dropbox/DCMST/allinstances/R123/' + inst + ".dat" + " DRresults-CGBC.txt" + " 2"
        print command
        print "\n"
        os.system(command)
# --------------------------------------------- #

# --------------- Instancias LH-EUC --------------- #
fo = open("LHEucresults-CGBC.txt", "w");
fo.write("Inst LB UB %Gap #cols %EF t(s) #nodes t(s) z^* T(s)\n");
fo.close();

for nv in ["100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"]:
    for d in ["3", "4", "5"]:
	inst = 'inst_euc_' + str(nv) + "_" + str(d)
	command =  './CG-BC ~/Dropbox/DCMST/allinstances/LH/Euc/' + inst + " LHEucresults-CGBC.txt" + " 2"
	print command
	print "\n"
	os.system(command)
# --------------------------------------------- #

# --------------- Instancias LH-Random --------------- #
fo = open("LHRandomresults-CGBC.txt", "w");
fo.write("Inst LB UB %Gap #cols %EF t(s) #nodes t(s) z^* T(s)\n");
fo.close();

for nv in ["100", "200", "300", "400", "500", "600", "700", "800", "900", "1000"]:
    for d in ["3", "4", "5"]:
	inst = 'inst_noneuc_' + str(nv) + "_" + str(d)
	command =  './CG-BC ~/Dropbox/DCMST/allinstances/LH/Random/' + inst + " LHRandomresults-CGBC.txt" + " 2"
	print command
	print "\n"
	os.system(command)
# --------------------------------------------- #

# --------------- Instancias ANDINST --------------- #
fo = open("ANDINSTresults-CGBC.txt", "w");
fo.write("Inst LB UB %Gap #cols %EF t(s) #nodes t(s) z^* T(s)\n");
fo.close();

for nv in [100, 200, 300, 400]:
    for num_inst in [1, 2, 3]:
        inst = 'tb1ct' + str(nv) + "_" + str(num_inst)
        command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/' + inst + ".txt" + " ANDINSTresults-CGBC.txt" + " 1"
        print command
        print "\n"
        os.system(command)

for num_inst in [1, 2, 3]:
    inst = 'tb2ct500_' + str(num_inst)
    command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/' + inst + ".txt" + " ANDINSTresults-CGBC.txt" + " 1"
    print command
    print "\n"
    os.system(command)

command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb2ct600_1.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb2ct700_1.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb2ct800_1.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb2ct900_1.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb3ct1000_1.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb3ct2000_1.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb3ct2000_2.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb3ct2000_3.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb3ct2000_4.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
command =  './CG-BC ~/Dropbox/DCMST/allinstances/ANDINST/tb3ct2000_5.txt' + " ANDINSTresults-CGBC.txt" + " 1"
os.system(command)
# --------------------------------------------- #
