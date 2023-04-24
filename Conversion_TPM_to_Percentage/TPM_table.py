#!/usr/bin/env python3
import os
import pandas as pd

f_1 = open("files.txt", "r", encoding = "utf-8")


filename = []
i = 1
for file_1 in f_1:
    line_1 = file_1.strip().split('\t')
    filename.append(line_1[0])

f_1.close

j = 0

list_all = []
list_title = []

for file in filename:
    samplename = file.split('.')[0]
    list_title.append(samplename)

    filename2 = samplename + ".talon.classification_filtered_F.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")

    k = 0
    dict_tmp = {}
    for file_2 in f_2:
        if k != 0:
            line_2 = file_2.strip().split('\t')
            gene_name = line_2[0]
            if gene_name not in dict_tmp:
                dict_tmp[gene_name] = line_2[3]
            elif gene_name in dict_tmp:
                num = float(dict_tmp[gene_name])
                num_re = num + float(line_2[3])
                dict_tmp[gene_name] = float(num_re)
        k += 1
    list_tmp = list(dict_tmp.items())

    list_all.append(list_tmp)


for i in range(len(list_all)):
    list_pd = pd.DataFrame(list_all[i], columns=['gene_name',list_title[i]])
    if i == 0:
        list_pd_combined = list_pd
    elif i != 0:
        list_pd_combined = pd.merge(list_pd_combined, list_pd,how='outer')
    print(i)

print(list_pd_combined)
list_pd_combined.to_csv("tmp.csv")


#annotation情報の付加
f_3 = open("tmp.csv", "r", encoding = "utf-8")
TPM_all = {}
j = 0
tmp_sample = []
for file_3 in f_3:
    line_3 = file_3.strip().split(',')
    if j == 0:
        for item in range(len(line_3)):
            if item >= 2:
                tmp_sample.append(line_3[item])
    else:
        tmp2 = []
        for k in range(len(line_3)):
            if k >= 2:
                if line_3[k] != "":
                    tmp2.append(line_3[1] + "/" + line_3[k])
                else:
                    tmp2.append(line_3[1] + "/" + "0")
        tmp_line = '*'.join(tmp2)

        if line_3[1].split("_")[0] not in TPM_all:
            TPM_all[line_3[1].split("_")[0]] = tmp_line

        elif line_3[1].split("_")[0] in TPM_all:
            tmp_line2 = TPM_all[line_3[1].split("_")[0]] + ";" + tmp_line
            TPM_all[line_3[1].split("_")[0]] = tmp_line2
    j += 1
f_3.close()


f_4 = open("tmp.csv", "r", encoding = "utf-8")
o_1 = open("isoform_TPM_withAnnotation.txt", "w", encoding = "utf-8")
j = 0
for file_4 in f_4:
    line_4 = file_4.strip().split(',')

    if j == 0:
        tmp_line = '*'.join(tmp_sample)
        isoform = "isoform_count"
    else:
        if line_4[1].split("_")[0] in TPM_all:
            tmp = []
            if len(TPM_all[line_4[1].split("_")[0]].split(";")) == 1:
                tmp_line = "NA"
                isoform = "1"
            elif len(TPM_all[line_4[1].split("_")[0]].split(";")) != 1:
                for item in TPM_all[line_4[1].split("_")[0]].split(";"):
                    if item.split("/")[0] != line_4[1]:
                        tmp.append(item)
                tmp_line = ';'.join(tmp)
                isoform = str(len(TPM_all[line_4[1].split("_")[0]].split(";")))

    line_4.append(isoform)
    line_4.append(tmp_line)
    list_w = '\t'.join(line_4)
    list_w2 = list_w + "\n"
    o_1.writelines(list_w2)

    j += 1

os.remove("tmp.csv")
