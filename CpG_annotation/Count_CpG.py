#!/usr/bin/env python3
import os

f_1 = open("Candidate_probes_annotated.txt", "r", encoding = "utf-8")
f_2 = open("isoform_candidate.txt", "r", encoding = "utf-8")
o_1 = open("tmp.txt", "w", encoding = "utf-8")

i = 0
meth_list = []
meth_index_dict = {}

for file_1 in f_1:
    line_1 = file_1.strip().split('\t')

    if line_1[14] != "NotAnnotated":
        meth_list.append([line_1[1], line_1[3], line_1[4], line_1[5], line_1[9], line_1[10], line_1[11], line_1[12], line_1[13], line_1[14]])

        for item in line_1[14].split(";"):
            if item not in meth_index_dict:
                meth_index_dict[item] = str(i)
            elif item in meth_index_dict:
                if meth_index_dict[item] != str(i):
                    meth_index_dict[item] = meth_index_dict[item] + ";" + str(i)

        i += 1

for file_2 in f_2:
    line_2 = file_2.strip().split('\t')

    tmp = [line_2[0],line_2[1]]

    if line_2[1] in meth_index_dict:
        for item in meth_index_dict[line_2[1]].split(";"):
            tmp.append("*".join(meth_list[int(item)]))
    else:
        tmp = [line_2[0],line_2[1],"None_meth_candidate"]

    list_w = '\t'.join(tmp)
    list_w2 = list_w + "\n"
    o_1.writelines(list_w2)


#---------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------Count CpGs
#---------------------------------------------------------------------------------------------------------------------
f_1 = open("tmp.txt", "r", encoding = "utf-8")
o_1 = open("Annotated_methylationsites_count_result.txt", "w", encoding = "utf-8")

tmp = ["No", "gene", "No.annotatedCpG", "N_Shelf", "N_Shore", "Island", "S_Shore", "S_Shelf", "TSS1500", "TSS200", "UTR5", "Exon1at","ExonBnd",\
"UTR3", "DNase_N", "OpenChromatin_N", "TFBS", "Enhancer"]

list_w = '\t'.join(tmp)
list_w2 = list_w + "\n"
o_1.writelines(list_w2)

i = 0
for file_1 in f_1:
    line_1 = file_1.strip().split('\t')

    if line_1[2] == "None_meth_candidate":
        annotated_CpG = "0"
        tmp = [line_1[0], line_1[1], str(annotated_CpG)]
    else:
        annotated_CpG = len(line_1) - 2
        CpG = []
        Group = []
        DNase = []
        OpenChromatin = []
        TFBS = []
        Enhancer = []
        Gene = []

        tmp = [line_1[0], line_1[1], str(annotated_CpG)]


        for item in range(2,len(line_1)):
            CpG.append(line_1[item].split("*")[3])

            if line_1[item].split("*")[5] != "":
                DNase.append(line_1[item].split("*")[5])
            if line_1[item].split("*")[6] != "":
                OpenChromatin.append(line_1[item].split("*")[6])
            if line_1[item].split("*")[7] != "":
                TFBS.append(line_1[item].split("*")[7])
            if line_1[item].split("*")[8] != "":
                Enhancer.append(line_1[item].split("*")[8])
            Gene = line_1[item].split("*")[9].split(";")

            Group1 = line_1[item].split("*")[4].split(";")

            j = 0
            tmp = []
            for item1 in Gene:
                if item1 == line_1[1]:
                    tmp.append(j)
                j += 1

            Group2 = []
            for item2 in tmp:
                Group2.append(Group1[item2])

            Group3 = list(set(Group2))

            for item_g2 in Group3:
                Group.append(item_g2)


        #Count
        #CpG
        N_Shelf = CpG.count('N_Shelf')
        N_Shore = CpG.count('N_Shore')
        Island = CpG.count('Island')
        S_Shore = CpG.count('S_Shore')
        S_Shelf = CpG.count('S_Shelf')

        #Group
        TSS1500 = Group.count('TSS1500')
        TSS200 = Group.count('TSS200')
        UTR5 = Group.count("5'UTR")
        Exon1st = Group.count('1stExon')
        ExonBnd = Group.count('ExonBnd')
        UTR3 = Group.count("3'UTR")

        #Others
        DNase_N = len(DNase)
        OpenChromatin_N = len(OpenChromatin)
        TFBS_N = len(TFBS)
        Enhancer_N = len(Enhancer)

        tmp = [line_1[0], line_1[1], str(annotated_CpG), str(N_Shelf), str(N_Shore), str(Island), str(S_Shore), str(S_Shelf), str(TSS1500), str(TSS200),\
        str(UTR5), str(Exon1st), str(ExonBnd), str(UTR3), str(DNase_N), str(OpenChromatin_N), str(TFBS_N), str(Enhancer_N)]

    if i != 0:
        list_w = '\t'.join(tmp)
        list_w2 = list_w + "\n"
        o_1.writelines(list_w2)
    i += 1

os.remove("tmp.txt")
