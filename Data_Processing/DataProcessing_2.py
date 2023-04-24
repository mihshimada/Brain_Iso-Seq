#!/usr/bin/env python3

f_1 = open("files.txt", "r", encoding = "utf-8")

filename = []
All = []
i = 1
for file_1 in f_1:
    line_1 = file_1.strip().split('\t') #改行コートをはがしてタブ区切
    filename.append(line_1[0])

f_1.close

j = 0
for file in filename:
    samplename = file.split('.')[0]
    All.append([samplename])

    filename2 = "./abundance/" + samplename + ".hifi.collapsed.abundance.txt"
    f_2 = open(filename2, "r", encoding = "utf-8")

    norm_dict = {}
    for file_2 in f_2:
        line_2 = file_2.strip().split('\t') #改行コートをはがしてタブ区切

        if len(line_2) >= 3:
            count_norm = line_2[1] + ":" + line_2[2]
            norm_dict[line_2[0]] = count_norm

    f_2.close

    #print(norm_dict)



    filename3 = "./talon_result_95/" + samplename + ".talon.classification_more2.tsv"
    o_1 = open("./talon_result_95/" + samplename + ".talon.classification_more2_normcount.tsv", "w", encoding = "utf-8")


    f_3 = open(filename3, "r", encoding = "utf-8")

    i = 1
    for file_3 in f_3:
        line_3 = file_3.strip().split('\t') #改行コートをはがしてタブ区切

        if i == 1:
            tmp = ["SuperID", "count", "norm_count","More.2.sample","More.2.sample_NONredundant"]
            for j in range(2, len(line_3)):
                tmp.append(line_3[j])
        else:
            tmp_s = line_3[2].split(",")
            tmp_s2 = set(tmp_s)
            tmp_s3 = ','.join(tmp_s2)

            tmp = [line_3[0], norm_dict[line_3[3]].split(":")[0], str(float(norm_dict[line_3[3]].split(":")[1])*1000000), str(len(tmp_s2)),tmp_s3]
            for j in range(2, len(line_3)):
                tmp.append(line_3[j])

        tmp.append("\n")
        list_w = '\t'.join(tmp)
        o_1.writelines(list_w)

        i += 1

    f_3.close
    o_1.close
