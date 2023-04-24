#!/usr/bin/env python3

import pandas as pd
import collections
import os


f_1 = open("files.txt", "r", encoding = "utf-8")

filename = []
All = []
i = 1
for file_1 in f_1:
    line_1 = file_1.strip().split('\t')
    filename.append(line_1[0])

f_1.close

j = 0
list = []
for file in filename:
    samplename = file.split('.')[0]
    All.append([samplename])

    filename2 = "./talon_result_95/" + samplename + ".talon.classification.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")

    i = 1
    for file_2 in f_2:
        line_2 = file_2.strip().split('\t')
        #print(len(line_2))
        list.append([line_2[71],samplename,1])
        All[j].append(line_2[71])

        i+=1

    j += 1
    print(j)

list_pd = pd.DataFrame(list, columns=['category','sample','n'])



result = list_pd.groupby(['category']).count()
result.to_csv("tmp.csv")

#List IDs that are expressed in two or more samples
f_3 = open("tmp.csv", "r", encoding = "utf-8")

more2sample = []
sample_list = []
sample_list_index = []

k = 0
for line_3 in f_3:
    item_3 = line_3.strip().split(',')

    if item_3[2] != "1":
        more2sample.append(item_3[0])
        sample_list_index.append(item_3[0])

        tmp = [item_3[0]]
        for i in All:
            for j in i:
                if item_3[0] == j:
                    tmp.append(i[0])
        sample_list.append(tmp)
    k += 1


#Reflected in the classicication file
for file in filename:
    samplename = file.split('.')[0]

    filename2 = "./talon_result_95/" + samplename + ".talon.classification.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")
    o_1 = open("./talon_result_95/" + samplename + ".talon.classification_more2.tsv", "w", encoding = "utf-8")

    i = 1
    for file_2 in f_2:
        line_2 = file_2.strip().split('\t') #改行コートをはがしてタブ区切

        if i == 1:
            tmp = ["SuperID","More.2.sample","Samples"]
            for j in line_2:
                tmp.append(j)
            tmp.append("\n")
            list_w = '\t'.join(tmp)
            o_1.writelines(list_w)

        else:
            if line_2[71] in more2sample:
                flg = "TRUE"
                index_no = more2sample.index(line_2[71])

                l = 0
                tmp2 = []
                for m in sample_list[index_no]:
                    if l != 0:
                        tmp2.append(m)
                    l += 1

                tmp3 = ','.join(tmp2)

            else:
                flg = "FALSE"
                tmp3 = "NA"

            tmp = [line_2[71],flg,tmp3]
            for j in range(len(line_2)-1):
                tmp.append(line_2[j])
            tmp.append("\n")
            list_w = '\t'.join(tmp)
            o_1.writelines(list_w)
        i += 1
    o_1.close
