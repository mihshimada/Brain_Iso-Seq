#!/usr/bin/env python3
import os
import sys

f_1 = open(sys.argv[1], "r", encoding = "utf-8")
f_2 = open("candidate_probes.txt", "r", encoding = "utf-8")
o_1 = open("tmp.txt", "w", encoding = "utf-8")

i = 0
probe_dict = {}
for file_1 in f_1:
    line_1 = file_1.strip().split(',')

    if len(line_1) >= 52:
        if line_1[49] != "NA" and i != 0:
            probe_dict[line_1[0]] = line_1[48] + "*" + str(int(line_1[49])+1) + "*" + line_1[14] + "*" + line_1[19] + "*" + line_1[26] + "*" + line_1[27]\
             + "*" + line_1[28] + "*" + line_1[29] + "*" + line_1[33] + "*" + line_1[35] + "*" + line_1[37] + "*" + line_1[23]
        elif line_1[49] == "NA" and i != 0:
            probe_dict[line_1[0]] = line_1[48] + "*" + "NA" + "*" + line_1[14] + "*" + line_1[19] + "*" + line_1[26] + "*" + line_1[27]\
             + "*" + line_1[28] + "*" + line_1[29] + "*" + line_1[33] + "*" + line_1[35] + "*" + line_1[37] + "*" + line_1[23]

        i += 1
        print(i)

i = 0
for file_2 in f_2:
    line_2 = file_2.strip().split('\t')

    if i == 0:
        tmp = [line_2[0], line_2[1], "chr", "mapinfo", "strand", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Group", "GencodeBasicV12_NAME", "GencodeBasicV12_Accession",\
         "GencodeBasicV12_Group","DNase_Hypersensitivity_NAME", "OpenChromatin_NAME","TFBS_NAME","450k_Enhancer"]
    else:
        tmp = [line_2[0], line_2[1]]

        for item in probe_dict[line_2[1]].split("*"):
            tmp.append(item)

    list_w = ','.join(tmp)
    list_w2 = list_w + "\n"
    o_1.writelines(list_w2)

    i += 1


#-----------------------------------------------------------------------------------------------------

f_1 = open("gene_transcript.txt", "r", encoding = "utf-8")
f_2 = open("tmp.txt", "r", encoding = "utf-8")
o_1 = open("Candidate_probes_annotated.txt", "w", encoding = "utf-8")


ENSG_dict = {}
for file_1 in f_1:
    line_1 = file_1.strip().split('\t')

    if line_1[1].split(".")[0] not in ENSG_dict:
        ENSG_dict[line_1[1].split(".")[0]] = line_1[0].split(".")[0]
    else:
        ENSG_dict[line_1[1].split(".")[0]] = ENSG_dict[line_1[1].split(".")[0]] + ";" + line_1[0].split(".")[0]

i = 0
for file_2 in f_2:
    line_2 = file_2.strip().split(',')

    if len(line_2[8]) != "":
        tmp = []
        for item in line_2[8].split(";"):
            if item.split(".")[0] in ENSG_dict:
                tmp.append(ENSG_dict[item.split(".")[0]])
            else:
                tmp.append("NA")

            tmp3 = ';'.join(tmp)
        tmp4 = line_2
        if i == 0:
            tmp4.append("ENSG_ID")
        else:
            tmp4.append(tmp3)
    else:
        tmp4 = line_2
        tmp4.append("NotAnnotated")
        print(line_2)

    list_w = '\t'.join(tmp4)
    list_w2 = list_w + "\n"
    o_1.writelines(list_w2)

    i += 1

os.remove("tmp.txt")
