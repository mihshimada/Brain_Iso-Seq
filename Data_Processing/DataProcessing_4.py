#!/usr/bin/env python3
import pandas as pd
import collections
import os


f_1 = open("files.txt", "r", encoding = "utf-8")

filename = []
i = 1
for file_1 in f_1:
    line_1 = file_1.strip().split('\t') #改行コートをはがしてタブ区切
    filename.append(line_1[0])

f_1.close

j = 0
list = []
for file in filename:
    samplename = file.split('.')[0]

    filename2 = "./talon_result_95/" + samplename + ".talon.classification_isoform.combined_count.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")

    filename3 = "./talon_result_95/" + samplename + ".talon.classification_filtered_F.tsv"
    o_1 = open(filename3, "w", encoding = "utf-8")

    for file_2 in f_2:
        line_2 = file_2.strip().split('\t') #改行コートをはがしてタブ区切

        if line_2[1] != "1":
            list_w = '\t'.join(line_2)
            list_w2 = list_w + "\n"
            o_1.writelines(list_w2)

    print(j)
    j += 1
