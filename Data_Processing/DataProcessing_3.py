#!/usr/bin/env python3
import pandas as pd

samplename = input('Enter your sample name : ')

filename3 = "./talon_result_95/" + samplename + ".talon.classification_isoform.combined_count.tsv"
o_1 = open(filename3, "w", encoding = "utf-8")
#-----------------------------------------------------------------------------function1 (addition)
def syukei_sum(list_no,colname):
    filename2 = "./talon_result_95/" + samplename + ".talon.classification_more2_normcount.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")

    list = []
    list_index = []

    for file_2 in f_2:
        line_2 = file_2.strip().split('\t')

        if line_2[0] not in list_index:
            list_index.append(line_2[0])
            list.append([line_2[0],line_2[list_no]])

        elif line_2[0] in list_index:
            index_no = list_index.index(line_2[0])
            list[index_no].append(line_2[list_no])
    f_2.close

    list_countN = []
    for i in list:
        if len(i) == 2:
            list_countN.append(i)
        else:
            countN = 0
            for j in range(len(i)):
                if j == 0:
                    tmp = [i[j]]
                else:
                    countN = countN + float(i[j])
            tmp.append(str(countN))
            list_countN.append(tmp)
    #print(list_countN)
    list_countN_pd = pd.DataFrame(list_countN,columns=["SuperID",colname+"_sum"])
    return list_countN_pd


#-----------------------------------------------------------------------------function2 (regard as a character string and summarize)
def syukei_str(list_no,collapse,colname):
    filename2 = "./talon_result_95/" + samplename + ".talon.classification_more2_normcount.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")

    list = []
    list_index = []

    for file_2 in f_2:
        line_2 = file_2.strip().split('\t')

        if line_2[0] not in list_index:
            list_index.append(line_2[0])
            list.append([line_2[0],line_2[list_no]])

        elif line_2[0] in list_index:
            index_no = list_index.index(line_2[0])
            list[index_no].append(line_2[list_no])
    f_2.close

    list2 = []
    if collapse == "Y":
        for i in list:
            for j in range(len(i)):
                if j == 0:
                    tmp = [i[j]]
                    tmp2 = []
                else:
                    tmp2.append(i[j])
            #print(tmp2)
            tmp3 = set(tmp2)
            tmp4 = ';'.join(tmp3)
            tmp.append(tmp4)

            list2.append(tmp)
    elif collapse == "N":
        for i in list:
            for j in range(len(i)):
                if j == 0:
                    tmp = [i[j]]
                    tmp2 = []
                else:
                    tmp2.append(i[j])
            #print(tmp2)
            tmp4 = ';'.join(tmp2)
            tmp.append(tmp4)

            list2.append(tmp)
    #print(list_countN)

    list2_pd = pd.DataFrame(list2,columns=["SuperID",colname])
    return list2_pd

#-----------------------------------------------------------------------------function3 (Output only the average (summary=N) || Output the average and original values)
def syukei_average(list_no,summary,colname):
    filename2 = "./talon_result_95/" + samplename + ".talon.classification_more2_normcount.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")

    list = []
    list_index = []

    for file_2 in f_2:
        line_2 = file_2.strip().split('\t')

        if line_2[0] not in list_index:
            list_index.append(line_2[0])
            list.append([line_2[0],line_2[list_no]])

        elif line_2[0] in list_index:
            index_no = list_index.index(line_2[0])
            list[index_no].append(line_2[list_no])
    f_2.close

    list2 = []
    list3 = []
    if summary == "Y":
        k = 0
        for i in list:
            if k == 0:
                list2.append(i)
                list3.append(i)
            else:
                for j in range(len(i)):
                    if j == 0:
                        tmp = [i[j]]
                        tmp5 = [i[j]]
                        tmp2 = []
                        tmp2_2 = []
                    else:
                        if i[j] != "NA":
                            tmp2.append(float(i[j]))
                            tmp2_2.append(i[j])
                        elif i[j] == "NA":
                            tmp2.append(i[j])
                            tmp2_2.append(i[j])
                if "NA" not in tmp2:
                    average = sum(tmp2)/len(tmp2)
                else:
                    average = "NA"
                tmp3 = set(tmp2_2)
                tmp4 = ';'.join(tmp3)
                tmp5.append(tmp4)
                tmp.append(str(average))


                list2.append(tmp5)
                list3.append(tmp)
            k += 1

        list2_pd = pd.DataFrame(list2,columns=["SuperID",colname])
        #print(list2_pd)
        colname2 = str(colname) + "_average"
        list3_pd = pd.DataFrame(list3,columns=["SuperID",colname2])
        #print(list3_pd)
        list4_pd = pd.merge(list2_pd, list3_pd, on='SuperID')
        #print(list4_pd)

        return list4_pd
        #return list3

    elif summary == "N":
        k = 0
        for i in list:
            if k == 0:
                list2.append(i)
                list3.append(i)
            else:
                for j in range(len(i)):
                    if j == 0:
                        tmp = [i[j]]
                        tmp2 = []
                    else:
                        tmp2.append(float(i[j]))
                #print(tmp2)
                average = sum(tmp2)/len(tmp2)
                tmp.append(str(average))
                list3.append(tmp)
        colname2 = str(colname) + "_average"
        list3_pd = pd.DataFrame(list3,columns=["SuperID",colname2])
        return list3_pd
    #print(list_countN)


#-----------------------------------------------------------------------------function4 (Tabulate duplicate samples）
def syukei_sample(list_no):
    filename2 = "./talon_result_95/" + samplename + ".talon.classification_more2_normcount.tsv"
    filename3 = "./talon_result_95/" + samplename + ".talon.classification_isoform.combined_count.tsv"
    f_2 = open(filename2, "r", encoding = "utf-8")
    o_1 = open(filename3, "w", encoding = "utf-8")

    list = []
    list_index = []

    for file_2 in f_2:
        line_2 = file_2.strip().split('\t')
        temp = line_2[list_no].split(',')

        if line_2[0] not in list_index:
            list_index.append(line_2[0])
            index_no = list_index.index(line_2[0])
            list.append([line_2[0]])
            for item in temp:
                if item == "NA":
                    list[index_no].append(samplename)
                else:
                    list[index_no].append(item)
        elif line_2[0] in list_index:
            index_no = list_index.index(line_2[0])
            for item in temp:
                if item == "NA":
                    list[index_no].append(samplename)
                else:
                    list[index_no].append(item)
    f_2.close

    list2 = []

    for i in list:
        for j in range(len(i)):
            if j == 0:
                tmp = [i[j]]
                tmp2 = []
            else:
                tmp2.append(i[j])
        #print(tmp2)
        tmp3 = set(tmp2)
        tmp4 = ';'.join(tmp3)
        tmp.append(tmp4)
        list2.append(tmp)


    list3 = []
    for item in list2:
        count = len(item[1].split(";"))

        for item2 in range(len(item)):
            if item2 == 0:
                tmp = [item[item2]]
                tmp.append(str(count))
            else:
                tmp.append(item[item2])
        list3.append(tmp)



    list3_pd = pd.DataFrame(list3,columns=["SuperID","sample_count","sample"])
    return list3_pd



#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
data_pd_0 = syukei_sample(4)
data_pd_1 = syukei_sum(2,"norm_count")
data_pd_2 = syukei_str(19,"Y","annot_gene_name")
data_pd_3 = syukei_str(20,"Y","annot_transcript_name")
data_pd_4 = syukei_str(21,"Y","gene_novelty")
data_pd_5 = syukei_str(22,"N","transcript_novelty")
data_pd_6 = syukei_str(23,"N","ISM_subtype")
data_pd_7 = syukei_str(29,"N","isoform")
data_pd_8 = syukei_str(30,"Y","chrom.y")
data_pd_9 = syukei_str(31,"Y","strand.y")
data_pd_10 = syukei_average(32,"Y","length")
data_pd_11 = syukei_str(33,"Y","exons")
data_pd_12 = syukei_str(34,"Y","structural_category")
data_pd_13 = syukei_str(35,"Y","associated_gene")
data_pd_14 = syukei_str(36,"Y","associated_transcript")
data_pd_15 = syukei_str(37,"Y","ref_length")
data_pd_16 = syukei_str(38,"Y","ref_exons")
data_pd_17 = syukei_average(39,"Y","diff_to_TSS")
data_pd_18 = syukei_average(40,"Y","diff_to_TTS")
data_pd_19 = syukei_average(41,"Y","diff_to_gene_TSS")
data_pd_20 = syukei_average(42,"Y","diff_to_gene_TTS")
data_pd_21 = syukei_str(43,"Y","subcategory")
data_pd_22 = syukei_str(44,"N","RTS_stage")
data_pd_23 = syukei_str(45,"N","all_canonical")
data_pd_24 = syukei_str(57,"N","FSM_class")
data_pd_25 = syukei_str(58,"Y","coding")
data_pd_26 = syukei_average(59,"Y","ORF_length")
data_pd_27 = syukei_average(60,"Y","CDS_length")
data_pd_28 = syukei_str(61,"N","CDS_start")
data_pd_29 = syukei_str(62,"N","CDS_end")
data_pd_30 = syukei_average(63,"Y","CDS_genomic_start")
data_pd_31 = syukei_average(64,"Y","CDS_genomic_end")
data_pd_32 = syukei_str(65,"Y","predicted_NMD")
data_pd_33 = syukei_average(66,"Y","perc_A_downstream_TTS")
data_pd_34 = syukei_str(67,"Y","seq_A_downstream_TTS")
data_pd_35 = syukei_average(68,"Y","dist_to_CAGE_peak")
data_pd_36 = syukei_str(69,"Y","within_CAGE_peak")
data_pd_37 = syukei_str(70,"N","dist_to_polyA_site")
data_pd_38 = syukei_str(71,"Y","within_polyA_site")
data_pd_39 = syukei_str(72,"Y","polyA_motif")
data_pd_40 = syukei_average(73,"Y","polyA_dist")
data_pd_41 = syukei_str(74,"Y","polyA_motif_found")
data_pd_42 = syukei_str(75,"Y","ORF_seq")
#data_pd_43 = syukei_str(77,"N","POS_MLprob")


last_0 = pd.merge(data_pd_0, data_pd_1, on='SuperID')
last_1 = pd.merge(last_0, data_pd_2, on='SuperID')
last_2 = pd.merge(last_1, data_pd_3, on='SuperID')
last_3= pd.merge(last_2 , data_pd_4, on='SuperID')
last_4= pd.merge(last_3 , data_pd_5, on='SuperID')
last_5= pd.merge(last_4 , data_pd_6, on='SuperID')
last_6= pd.merge(last_5 , data_pd_7, on='SuperID')
last_7= pd.merge(last_6 , data_pd_8, on='SuperID')
last_8= pd.merge(last_7 , data_pd_9, on='SuperID')
last_9= pd.merge(last_8 , data_pd_10, on='SuperID')
last_10= pd.merge(last_9 , data_pd_11, on='SuperID')
last_11= pd.merge(last_10 , data_pd_12, on='SuperID')
last_12= pd.merge(last_11 , data_pd_13, on='SuperID')
last_13= pd.merge(last_12 , data_pd_14, on='SuperID')
last_14= pd.merge(last_13 , data_pd_15, on='SuperID')
last_15= pd.merge(last_14 , data_pd_16, on='SuperID')
last_16= pd.merge(last_15 , data_pd_17, on='SuperID')
last_17= pd.merge(last_16 , data_pd_18, on='SuperID')
last_18= pd.merge(last_17 , data_pd_19, on='SuperID')
last_19= pd.merge(last_18 , data_pd_20, on='SuperID')
last_20= pd.merge(last_19 , data_pd_21, on='SuperID')
last_21= pd.merge(last_20 , data_pd_22, on='SuperID')
last_22= pd.merge(last_21 , data_pd_23, on='SuperID')
last_23= pd.merge(last_22 , data_pd_24, on='SuperID')
last_24= pd.merge(last_23 , data_pd_25, on='SuperID')
last_25= pd.merge(last_24 , data_pd_26, on='SuperID')
last_26= pd.merge(last_25 , data_pd_27, on='SuperID')
last_27= pd.merge(last_26 , data_pd_28, on='SuperID')
last_28= pd.merge(last_27 , data_pd_29, on='SuperID')
last_29= pd.merge(last_28 , data_pd_30, on='SuperID')
last_30= pd.merge(last_29 , data_pd_31, on='SuperID')
last_31= pd.merge(last_30 , data_pd_32, on='SuperID')
last_32= pd.merge(last_31 , data_pd_33, on='SuperID')
last_33= pd.merge(last_32 , data_pd_34, on='SuperID')
last_34= pd.merge(last_33 , data_pd_35, on='SuperID')
last_35= pd.merge(last_34 , data_pd_36, on='SuperID')
last_36= pd.merge(last_35 , data_pd_37, on='SuperID')
last_37= pd.merge(last_36 , data_pd_38, on='SuperID')
last_38= pd.merge(last_37 , data_pd_39, on='SuperID')
last_39= pd.merge(last_38 , data_pd_40, on='SuperID')
last_40= pd.merge(last_39 , data_pd_41, on='SuperID')
last_41= pd.merge(last_40 , data_pd_42, on='SuperID')
#last_42= pd.merge(last_41 , data_pd_43, on='SuperID')

print(last_41)


last_list = last_41.values.tolist()

last_list_name = last_41.columns.tolist()
list_w = '\t'.join(last_list_name)
list_w2 = list_w + "\n"
o_1.writelines(list_w2)

i = 0
for item in last_list:
    if i != 0:
        list_w = '\t'.join(item)
        list_w2 = list_w + "\n"
        o_1.writelines(list_w2)
    i += 1
