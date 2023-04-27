#!/usr/bin/env python3
import sys

f_1 = open(sys.argv[1], "r", encoding = "utf-8")
o_1 = open("gene_transcript.txt", "w", encoding = "utf-8")


tmp_list = []

for file_1 in f_1:
    line_1 = file_1.strip().split('\t')
    if len(line_1) >= 3 and line_1[2] == "transcript":
        gene_name = line_1[8].split(";")[0].split(" ")[1].replace('"', '')
        transcript_name = line_1[8].split(";")[1].split(" ")[2].replace('"', '')

        tmp_list.append(gene_name + ";" + transcript_name)

tmp_list2 = list(set(tmp_list))
print(tmp_list2)

for item in tmp_list2:

    tmp = item.split(";")

    list_w = '\t'.join(tmp)
    list_w2 = list_w + "\n"
    o_1.writelines(list_w2)

f_1.close
