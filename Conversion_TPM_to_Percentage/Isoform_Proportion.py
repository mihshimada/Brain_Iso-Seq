#!/usr/bin/env python3

f_1 = open("isoform_TPM_withAnnotation.txt", "r", encoding = "utf-8")
o_1 = open("isoform_proportion_individual.txt", "w", encoding = "utf-8")

i = 0
for file_1 in f_1:
    line_1 = file_1.strip().split('\t')

    if i != 0:
        percent_C = []

        if line_1[len(line_1) - 2] == "1":
            for j in range(2,len(line_1) - 2):
                if line_1[j] == "":
                    percent_C.append(0)
                else:
                    percent_C.append(100)

        else:
            list_C = []
            k = 0
            for j in range(2,len(line_1) - 2):
                if line_1[j] == "":
                    list_C.append([0])
                else:
                    list_C.append([float(line_1[j])])

            for k in line_1[len(line_1) - 1].split(";"):
                for l in range(len(k.split("*"))):
                    if l <= 3:
                        list_C[l].append(float(k.split("*")[l].split("/")[1]))

            percent_C = []
            for m in list_C:
                if sum(m) != 0:
                    percent_C.append(m[0]/sum(m) * 100)
                else:
                    percent_C.append(0)


        tmp = [line_1[1], line_1[len(line_1) - 2]]
        for item in percent_C:
            tmp.append(str(item))


        list_w = '\t'.join(tmp)
        list_w2 = list_w + "\n"
        o_1.writelines(list_w2)

    if i == 0:
        tmp = ["ID", "Number_of_isoform"]
        for j in range(1, len(line_1) - 2):
            tmp.append(line_1[j])
            print(j)
        list_w = '\t'.join(tmp)
        list_w2 = list_w + "\n"
        o_1.writelines(list_w2)

    i += 1

f_1.close
