#!/usr/bin/env python3

filename_1  = "exon_position.txt"
filename_2  = "combination.txt"

#Define a function to detect AS

def file_writing(out_name, tmp, count):
    if count >= 1:
        o_1 = open(out_name, 'a')
        for item in tmp:
            tmp2 = [str(i) for i in item]
            tmp_w_2 = '\t'.join(tmp2)
            tmp_w_3 = tmp_w_2 + "\n"
            o_1.writelines(tmp_w_3)
        o_1.close()


def AS_call(list_1,list_2,keys_1,keys_2,chr):
    if len(list_1) >= len(list_2):
        list1 = list_1
        list2 = list_2
    else:
        list1 = list_2
        list2 = list_1

    overlap_exon_no = []
    for i in range(len(list2)):
        list2_tmp = list2[i]
        len_list2 = list2_tmp[1] - list2_tmp[0] + 1


        for k in range(len(list1)):
            len_list1 = list1[k][1] - list1[k][0] + 1
            if list1[k][0] >= list2_tmp[0] and list1[k][0] <= list2_tmp[1]:
                over_s = list1[k][0]
                over_e = min(list2_tmp[1], list1[k][1])
                len_over = over_e - over_s + 1
            elif list2_tmp[0] >= list1[k][0] and list2_tmp[0] <= list1[k][1]:
                over_s = list2_tmp[0]
                over_e = min(list2_tmp[1], list1[k][1])
                len_over = over_e - over_s + 1
            else:
                len_over = 0


            over_per_for2 = (len_over / min(len_list1,len_list2)) * 100
            if over_per_for2 >= 50:
                tmp = [k,i]
                overlap_exon_no.append(tmp)


    #---------------------------------------------------------------------Variables
    completely_different = 0
    diff_transcription_start = 0
    diff_polyadenyration_site = 0
    diff_5p_firstexon = 0
    diff_3p_lastexon = 0
    diff_5p_ss = 0
    diff_3p_ss = 0
    intron_expression = 0
    exon_skip = 0
    cassette_exon = 0
    complecated_exon = 0
    #-----------------------------------------------------------Completely different
    if len(overlap_exon_no) == 0:
        completely_different = 1
    #-------------------------------------------------------diff_firstexon
    else:
        if [0,0] not in overlap_exon_no:
            diff_transcription_start = 1
            if list1[0][0] >= 0:
                diff_transcription_start_500bp_up = [[chr, int(list1[0][0]) - 500, int(list1[0][0]), keys_1, keys_2]]
                diff_transcription_start_500bp_up.append([chr, int(list2[0][0]) - 500, int(list2[0][0]), keys_1, keys_2])
                diff_transcription_start_500bp_down = [[chr, int(list1[0][0]), int(list1[0][0]) + 500 ,  keys_1, keys_2]]
                diff_transcription_start_500bp_down.append([chr, int(list2[0][0]), int(list2[0][0]) + 500, keys_1, keys_2])
                diff_transcription_start_1kbp_up = [[chr, int(list1[0][0]) - 1000, int(list1[0][0]), keys_1, keys_2]]
                diff_transcription_start_1kbp_up.append([chr, int(list2[0][0]) - 1000, int(list2[0][0]), keys_1, keys_2])
                diff_transcription_start_1kbp_down = [[chr, int(list1[0][0]), int(list1[0][0]) + 1000 ,  keys_1, keys_2]]
                diff_transcription_start_1kbp_down.append([chr, int(list2[0][0]), int(list2[0][0]) + 1000, keys_1, keys_2])
            elif list1[0][0] < 0:
                diff_transcription_start_500bp_up = [[chr, abs(int(list1[0][0])), abs(int(list1[0][0])) + 500 , keys_1, keys_2]]
                diff_transcription_start_500bp_up.append([chr, abs(int(list2[0][0])), abs(int(list2[0][0])) + 500, keys_1, keys_2])
                diff_transcription_start_500bp_down = [[chr, abs(int(list1[0][0])) - 500 , abs(int(list1[0][0])), keys_1, keys_2]]
                diff_transcription_start_500bp_down.append([chr, abs(int(list2[0][0])) - 500, abs(int(list2[0][0])), keys_1, keys_2])
                diff_transcription_start_1kbp_up = [[chr,  abs(int(list1[0][0])), abs(int(list1[0][0])) + 1000 , keys_1, keys_2]]
                diff_transcription_start_1kbp_up.append([chr, abs(int(list2[0][0])), abs(int(list2[0][0])) + 1000, keys_1, keys_2])
                diff_transcription_start_1kbp_down = [[chr, abs(int(list1[0][0])) - 1000 , abs(int(list1[0][0])), keys_1, keys_2]]
                diff_transcription_start_1kbp_down.append([chr, abs(int(list2[0][0])) - 1000, abs(int(list2[0][0])), keys_1, keys_2])


    #------------------------------------------------------diff_lastexon
        if [len(list1)-1,len(list2)-1] not in overlap_exon_no:
            diff_polyadenyration_site = 1
            if list1[0][0] >= 0:
                diff_polyadenyration_site_500bp_up = [[chr, int(list1[len(list1)-1][1]) - 500 , int(list1[len(list1)-1][1]), keys_1, keys_2]]
                diff_polyadenyration_site_500bp_up.append([chr, int(list2[len(list2)-1][1]) - 500, int(list2[len(list2)-1][1]), keys_1, keys_2])
                diff_polyadenyration_site_500bp_down = [[chr, int(list1[len(list1)-1][1]), int(list1[len(list1)-1][1]) + 500 , keys_1, keys_2]]
                diff_polyadenyration_site_500bp_down.append([chr, int(list2[len(list2)-1][1]), int(list2[len(list2)-1][1]) + 500, keys_1, keys_2])
                diff_polyadenyration_site_1kbp_up = [[chr, int(list1[len(list1)-1][1]) - 1000 , int(list1[len(list1)-1][1]), keys_1, keys_2]]
                diff_polyadenyration_site_1kbp_up.append([chr, int(list2[len(list2)-1][1]) - 1000, int(list2[len(list2)-1][1]), keys_1, keys_2])
                diff_polyadenyration_site_1kbp_down = [[chr, int(list1[len(list1)-1][1]), int(list1[len(list1)-1][1]) + 1000, keys_1, keys_2]]
                diff_polyadenyration_site_1kbp_down.append([chr, int(list2[len(list2)-1][1]), int(list2[len(list2)-1][1]) + 1000, keys_1, keys_2])
            elif list1[0][0] < 0:
                diff_polyadenyration_site_500bp_up = [[chr, abs(int(list1[len(list1)-1][1])), abs(int(list1[len(list1)-1][1])) + 500 , keys_1, keys_2]]
                diff_polyadenyration_site_500bp_up.append([chr, abs(int(list2[len(list2)-1][1])), abs(int(list2[len(list2)-1][1])) + 500, keys_1, keys_2])
                diff_polyadenyration_site_500bp_down = [[chr, abs(int(list1[len(list1)-1][1]))  - 500 , abs(int(list1[len(list1)-1][1])), keys_1, keys_2]]
                diff_polyadenyration_site_500bp_down.append([chr, abs(int(list2[len(list2)-1][1])) - 500,abs(int(list2[len(list2)-1][1])), keys_1, keys_2])
                diff_polyadenyration_site_1kbp_up = [[chr, abs(int(list1[len(list1)-1][1])), abs(int(list1[len(list1)-1][1]))  + 1000 , keys_1, keys_2]]
                diff_polyadenyration_site_1kbp_up.append([chr, abs(int(list2[len(list2)-1][1])), abs(int(list2[len(list2)-1][1])) + 1000, keys_1, keys_2])
                diff_polyadenyration_site_1kbp_down = [[chr, abs(int(list1[len(list1)-1][1])) - 1000 , abs(int(list1[len(list1)-1][1])), keys_1, keys_2]]
                diff_polyadenyration_site_1kbp_down.append([chr, abs(int(list2[len(list2)-1][1])) - 1000, abs(int(list2[len(list2)-1][1])), keys_1, keys_2])
    #--------------------------------------------------------------diff_5p_firstexon
        if [0,0] in overlap_exon_no:
            if list1[0][0] != list2[0][0]:
                diff_5p_firstexon = 1
    #---------------------------------------------------------------diff_3p_lastexon
        if [len(list1)-1,len(list2)-1] in overlap_exon_no:
            if list1[-1][1] != list2[-1][1]:
                diff_3p_lastexon = 1
    #----------------------------------------------------------diff_5p_ss/diff_3p_ss
        diff_5p_ss = 0
        diff_3p_ss = 0
        for i in range(len(overlap_exon_no)):
            if (overlap_exon_no[i][0] != 0 and overlap_exon_no[i][1] != 0) and (list1[overlap_exon_no[i][0]][0] != list2[overlap_exon_no[i][1]][0]) and \
            diff_5p_ss == 0:
                if list1[overlap_exon_no[i][0]][0] >= 0: #+鎖
                    min_1 = min(int(list1[overlap_exon_no[i][0]][0]),int(list2[overlap_exon_no[i][1]][0]))
                    max_1 = max(int(list1[overlap_exon_no[i][0]][0]),int(list2[overlap_exon_no[i][1]][0]))
                    diff_5p_ss_list = [[chr, min_1 - 500, max_1 + 500, keys_1, keys_2]]
                elif list1[overlap_exon_no[i][0]][0] < 0: #-鎖
                    min_1 = min(abs(int(list1[overlap_exon_no[i][0]][0])),abs(int(list2[overlap_exon_no[i][1]][0])))
                    max_1 = max(abs(int(list1[overlap_exon_no[i][0]][0])),abs(int(list2[overlap_exon_no[i][1]][0])))
                    diff_5p_ss_list = [[chr, min_1 - 500, max_1 + 500, keys_1, keys_2]]
                diff_5p_ss += 1

            elif (overlap_exon_no[i][0] != 0 and overlap_exon_no[i][1] != 0) and (list1[overlap_exon_no[i][0]][0] != list2[overlap_exon_no[i][1]][0]) and \
            diff_5p_ss != 0:
                if list1[overlap_exon_no[i][0]][0] >= 0: #+鎖
                    min_1 = min(int(list1[overlap_exon_no[i][0]][0]),int(list2[overlap_exon_no[i][1]][0]))
                    max_1 = max(int(list1[overlap_exon_no[i][0]][0]),int(list2[overlap_exon_no[i][1]][0]))
                    diff_5p_ss_list.append([chr, min_1 - 500, max_1 + 500, keys_1, keys_2])
                elif list1[overlap_exon_no[i][0]][0] < 0: #-鎖
                    min_1 = min(abs(int(list1[overlap_exon_no[i][0]][0])),abs(int(list2[overlap_exon_no[i][1]][0])))
                    max_1 = max(abs(int(list1[overlap_exon_no[i][0]][0])),abs(int(list2[overlap_exon_no[i][1]][0])))
                    diff_5p_ss_list.append([chr, min_1 - 500, max_1 + 500, keys_1, keys_2])

                diff_5p_ss += 1

            if (overlap_exon_no[i][0] != len(list1)-1 and overlap_exon_no[i][1] != len(list2)-1) and (list1[overlap_exon_no[i][0]][1] != list2[overlap_exon_no[i][1]][1]) and \
            diff_3p_ss == 0:
                if list1[overlap_exon_no[i][0]][0] >= 0: #+鎖
                    min_1 = min(int(list1[overlap_exon_no[i][0]][1]),int(list2[overlap_exon_no[i][1]][1]))
                    max_1 = max(int(list1[overlap_exon_no[i][0]][1]),int(list2[overlap_exon_no[i][1]][1]))
                    diff_3p_ss_list = [[chr, min_1 - 500, max_1 + 500, keys_1, keys_2]]
                elif list1[overlap_exon_no[i][0]][0] < 0: #-鎖
                    min_1 = min(abs(int(list1[overlap_exon_no[i][0]][1])),abs(int(list2[overlap_exon_no[i][1]][1])))
                    max_1 = max(abs(int(list1[overlap_exon_no[i][0]][1])),abs(int(list2[overlap_exon_no[i][1]][1])))
                    diff_3p_ss_list = [[chr, min_1 - 500, max_1 + 500, keys_1, keys_2]]
                diff_3p_ss += 1
            elif (overlap_exon_no[i][0] != len(list1)-1 and overlap_exon_no[i][1] != len(list2)-1) and (list1[overlap_exon_no[i][0]][1] != list2[overlap_exon_no[i][1]][1]) and \
            diff_3p_ss != 0:
                if list1[overlap_exon_no[i][0]][0] >= 0: #+鎖
                    min_1 = min(int(list1[overlap_exon_no[i][0]][1]),int(list2[overlap_exon_no[i][1]][1]))
                    max_1 = max(int(list1[overlap_exon_no[i][0]][1]),int(list2[overlap_exon_no[i][1]][1]))
                    diff_3p_ss_list.append([chr, min_1 - 500, max_1 + 500, keys_1, keys_2])
                elif list1[overlap_exon_no[i][0]][0] < 0: #-鎖
                    min_1 = min(abs(int(list1[overlap_exon_no[i][0]][1])),abs(int(list2[overlap_exon_no[i][1]][1])))
                    max_1 = max(abs(int(list1[overlap_exon_no[i][0]][1])),abs(int(list2[overlap_exon_no[i][1]][1])))
                    diff_3p_ss_list.append([chr, min_1 - 500, max_1 + 500, keys_1, keys_2])
                diff_3p_ss += 1
    #--------------------------------------------------------------intron_expression
        #--------------------------intron_sequence↓
        intron_1 = []
        intron_2 = []
        z = 0
        for i in range(len(list1)):
            if z == 0:
                start = list1[i][1] + 1
            else:
                end = list1[i][0] - 1
                intron_1.append([start,end])
                start = list1[i][1] + 1
            z+=1
        z = 0
        for i in range(len(list2)):
            if z == 0:
                start = list2[i][1] + 1
            else:
                end = list2[i][0] - 1
                intron_2.append([start,end])
                start = list2[i][1] + 1
            z+=1
        #-------------------------intron_sequence↑
        for item in list1:
            for item2 in intron_2:
                if item[0] <= item2[0] and item[1]>= item2[1] and intron_expression == 0:
                    if item[0] >= 0: #+鎖
                        intron_retention_up_5 = [[chr, int(item2[0]) - 500 , int(item2[0]) + 500, keys_1, keys_2]]
                        intron_retention_down_5 = [[chr, int(item2[1]) - 500 , int(item2[1]) + 500, keys_1, keys_2]]
                        intron_retention_expressing = [[chr, int(item2[0]), int(item2[1]), keys_1, keys_2]]
                    if item[0] < 0: #-鎖
                        intron_retention_up_5 = [[chr, abs(int(item2[0])) - 500 , abs(int(item2[0])) + 500, keys_1, keys_2]]
                        intron_retention_down_5 = [[chr, abs(int(item2[1])) - 500 , abs(int(item2[1])) + 500, keys_1, keys_2]]
                        intron_retention_expressing = [[chr, abs(int(item2[1])), abs(item2[0]), keys_1, keys_2]]
                    intron_expression += 1
                elif item[0] <= item2[0] and item[1]>= item2[1] and intron_expression != 0:
                    if item[0] >= 0: #+鎖
                        intron_retention_up_5.append([chr, int(item2[0]) - 500 , int(item2[0]) + 500, keys_1, keys_2])
                        intron_retention_down_5.append([chr, int(item2[1]) - 500 , int(item2[1]) + 500, keys_1, keys_2])
                        intron_retention_expressing.append([chr, int(item2[0]), int(item2[1]), keys_1, keys_2])
                    if item[0] < 0: #-鎖
                        intron_retention_up_5.append([chr, abs(int(item2[0])) - 500 , abs(int(item2[0])) + 500, keys_1, keys_2])
                        intron_retention_down_5.append([chr, abs(int(item2[1])) - 500 , abs(int(item2[1])) + 500, keys_1, keys_2])
                        intron_retention_expressing.append([chr, abs(int(item2[1])), abs(item2[0]), keys_1, keys_2])
                    intron_expression += 1
        for item in list2:
            for item2 in intron_1:
                if item[0] <= item2[0] and item[1]>= item2[1] and intron_expression == 0:
                    if item[0] >= 0: #+鎖
                        intron_retention_up_5 = [[chr, int(item2[0]) - 500 , int(item2[0]) + 500, keys_1, keys_2]]
                        intron_retention_down_5 = [[chr, int(item2[1]) - 500 , int(item2[1]) + 500, keys_1, keys_2]]
                        intron_retention_expressing = [[chr, int(item2[0]), int(item2[1]), keys_1, keys_2]]
                    if item[0] < 0: #-鎖
                        intron_retention_up_5 = [[chr, abs(int(item2[0])) - 500 , abs(int(item2[0])) + 500, keys_1, keys_2]]
                        intron_retention_down_5 = [[chr, abs(int(item2[1])) - 500 , abs(int(item2[1])) + 500, keys_1, keys_2]]
                        intron_retention_expressing = [[chr, abs(int(item2[1])), abs(item2[0]), keys_1, keys_2]]
                    intron_expression += 1
                elif item[0] <= item2[0] and item[1]>= item2[1] and intron_expression != 0:
                    if item[0] >= 0: #+鎖
                        intron_retention_up_5.append([chr, int(item2[0]) - 500 , int(item2[0]) + 500, keys_1, keys_2])
                        intron_retention_down_5.append([chr, int(item2[1]) - 500 , int(item2[1]) + 500, keys_1, keys_2])
                        intron_retention_expressing.append([chr, int(item2[0]), int(item2[1]), keys_1, keys_2])
                    if item[0] < 0: #-鎖
                        intron_retention_up_5.append([chr, abs(int(item2[0])) - 500 , abs(int(item2[0])) + 500, keys_1, keys_2])
                        intron_retention_down_5.append([chr, abs(int(item2[1])) - 500 , abs(int(item2[1])) + 500, keys_1, keys_2])
                        intron_retention_expressing.append([chr, abs(int(item2[1])), abs(item2[0]), keys_1, keys_2])
                    intron_expression += 1

    #---------------------------------------exon_skip/cassette_exon/complecated_exon
        z = 0
        for i in range(len(overlap_exon_no)):
            if z == 0:
                tmp1_1 = overlap_exon_no[i][0]
                tmp2_1 = overlap_exon_no[i][1]
            else:
                tmp1_2 = overlap_exon_no[i][0]
                tmp2_2 = overlap_exon_no[i][1]

                if tmp1_2 - tmp1_1 == 1 and tmp2_2 - tmp2_1 == 1:
                    structure_normal = 1
                elif tmp1_2 - tmp1_1 == 0 or tmp2_2 - tmp2_1 == 0:
                    structure_intron_expressing = 1
                elif (tmp1_2 - tmp1_1 == 2 and tmp2_2 - tmp2_1 == 1)\
                     or (tmp1_2 - tmp1_1 == 1 and tmp2_2 - tmp2_1 == 2):#exon_skip


                    pos1 = min(int(list1[tmp1_1][1]), int(list2[tmp2_1][1]))
                    pos2 = max(int(list1[tmp1_2][0]), int(list2[tmp2_2][0]))

                    if exon_skip == 0:
                        if pos1 >= 0: #+strand
                            exon_skip_pos_5p = [[chr, pos1 - 500, pos1 + 500, keys_1, keys_2]]
                            exon_skip_pos_3p = [[chr, pos2 - 500, pos2 + 500, keys_1, keys_2]]
                            exon_skip_pos_all = [[chr, pos1, pos2, keys_1, keys_2]]
                        elif pos1 < 0: #-strand
                            exon_skip_pos_5p = [[chr, abs(pos2) - 500, abs(pos2) + 500, keys_1, keys_2]]
                            exon_skip_pos_3p = [[chr, abs(pos1) - 500, abs(pos1) + 500, keys_1, keys_2]]
                            exon_skip_pos_all = [[chr, abs(pos2), abs(pos1), keys_1, keys_2]]
                        exon_skip += 1
                    elif exon_skip != 0:
                        if pos1 >= 0: #+strand
                            exon_skip_pos_5p.append([chr, pos1 - 500, pos1 + 500, keys_1, keys_2])
                            exon_skip_pos_3p.append([chr, pos2 - 500, pos2 + 500, keys_1, keys_2])
                            exon_skip_pos_all.append([chr, pos1, pos2, keys_1, keys_2])
                        elif pos1 < 0: #-strand
                            exon_skip_pos_5p.append([chr, abs(pos2) - 500, abs(pos2) + 500, keys_1, keys_2])
                            exon_skip_pos_3p.append([chr, abs(pos1) - 500, abs(pos1) + 500, keys_1, keys_2])
                            exon_skip_pos_all.append([chr, abs(pos2), abs(pos1), keys_1, keys_2])
                        exon_skip += 1

                elif (tmp1_2 - tmp1_1 == 2) and (tmp2_2 - tmp2_1 == 2):
                    pos1 = min(int(list1[tmp1_1][1]), int(list2[tmp2_1][1]))
                    pos2 = max(int(list1[tmp1_2][0]), int(list2[tmp2_2][0]))

                    if cassette_exon == 0:
                        if pos1 >= 0: #+strand
                            cassette_pos_5p = [[chr, pos1 - 500, pos1 + 500, keys_1, keys_2]]
                            cassette_pos_3p = [[chr, pos2 - 500, pos2 + 500, keys_1, keys_2]]
                            cassette_pos_all = [[chr, pos1, pos2, keys_1, keys_2]]
                        elif pos1 < 0: #-strand
                            cassette_pos_5p = [[chr, abs(pos2) - 500, abs(pos2) + 500, keys_1, keys_2]]
                            cassette_pos_3p = [[chr, abs(pos1) - 500, abs(pos1) + 500, keys_1, keys_2]]
                            cassette_pos_all = [[chr, abs(pos2), abs(pos1), keys_1, keys_2]]
                        cassette_exon += 1
                    elif cassette_exon != 0:
                        if pos1 >= 0: #+strand
                            cassette_pos_5p.append([chr, pos1 - 500, pos1 + 500, keys_1, keys_2])
                            cassette_pos_3p.append([chr, pos2 - 500, pos2 + 500, keys_1, keys_2])
                            cassette_pos_all.append([chr, pos1, pos2, keys_1, keys_2])
                        elif pos1 < 0: #-strand
                            cassette_pos_5p.append([chr, abs(pos2) - 500, abs(pos2) + 500, keys_1, keys_2])
                            cassette_pos_3p.append([chr, abs(pos1) - 500, abs(pos1) + 500, keys_1, keys_2])
                            cassette_pos_all.append([chr, abs(pos2), abs(pos1), keys_1, keys_2])
                        cassette_exon += 1
                else:
                    complecated_exon += 1

                tmp1_1 = overlap_exon_no[i][0]
                tmp2_1 = overlap_exon_no[i][1]
            z += 1


    tmp = [[str(keys_1),str(keys_2),str(completely_different),str(diff_transcription_start),str(diff_polyadenyration_site),str(diff_5p_firstexon)\
    ,str(diff_3p_lastexon),str(diff_5p_ss),str(diff_3p_ss),str(intron_expression),str(exon_skip),str(cassette_exon)\
    ,str(complecated_exon)]]


    file_writing("Count_result.txt", tmp, 1)
    if diff_transcription_start != 0:
        file_writing("diff_firstexon_500bp_up.txt", diff_transcription_start_500bp_up, diff_transcription_start)
        file_writing("diff_firstexon_500bp_down.txt", diff_transcription_start_500bp_down, diff_transcription_start)
        file_writing("diff_firstexon_1kbp_up.txt", diff_transcription_start_1kbp_up, diff_transcription_start)
        file_writing("diff_firstexon_1kbp_down.txt", diff_transcription_start_1kbp_down, diff_transcription_start)
    if diff_polyadenyration_site != 0:
        file_writing("diff_lastexon_500bp_up.txt", diff_polyadenyration_site_500bp_up, diff_polyadenyration_site)
        file_writing("diff_lastexon_500bp_down.txt", diff_polyadenyration_site_500bp_down, diff_polyadenyration_site)
        file_writing("diff_lastexon_1kbp_up.txt", diff_polyadenyration_site_1kbp_up, diff_polyadenyration_site)
        file_writing("diff_lastexon_1kbp_down.txt", diff_polyadenyration_site_1kbp_down, diff_polyadenyration_site)
    if diff_5p_ss != 0:
        file_writing("diff_5p_ss_list.txt", diff_5p_ss_list, diff_5p_ss)
    if diff_3p_ss != 0:
        file_writing("diff_3p_ss_list.txt", diff_3p_ss_list, diff_3p_ss)
    if intron_expression != 0:
        file_writing("intron_retention_up_5.txt", intron_retention_up_5, intron_expression)
        file_writing("intron_retention_down_5.txt", intron_retention_down_5, intron_expression)
        file_writing("intron_retention_expressing.txt", intron_retention_expressing, intron_expression)
    if exon_skip != 0:
        file_writing("exon_skip_pos_5p.txt", exon_skip_pos_5p, exon_skip)
        file_writing("exon_skip_pos_3p.txt", exon_skip_pos_3p, exon_skip)
        file_writing("exon_skip_pos_all.txt", exon_skip_pos_all, exon_skip)
    if cassette_exon != 0:
        file_writing("cassette_pos_5p.txt", cassette_pos_5p, cassette_exon)
        file_writing("cassette_pos_3p.txt", cassette_pos_3p, cassette_exon)
        file_writing("cassette_pos_all.txt", cassette_pos_all, cassette_exon)


tmp_1st = [["isoform_1","isoform_2","completely_different","diff_first_exon","diff_last_exon","diff_5p_firstexon","diff_3p_lastexon",\
"diff_5p_ss", "diff_3p_ss","intron_retention","exon_skip","cassette_exon","complecated_exon"]]
file_writing("Count_result.txt", tmp_1st, 1)


def get_keys_from_value(d, val):
    return [k for k, v in d.items() if v == val]



pos_dict = {}
f_1 = open(filename_1, "r", encoding = "utf-8")
for line_1 in f_1:
    item_1 = line_1.strip().split('\t') #改行コートをはがしてタブ区切り

    pos_dict[item_1[0]] = item_1[3] + "*" + item_1[1] + "*" + item_1[2]

f_1.close()


f_2 = open(filename_2, "r", encoding = "utf-8")
for line_2 in f_2:
    item_2 = line_2.strip().split('\t')

    keys_1 = item_2[0]
    keys_2 = item_2[1]

    exon_list_SpV1 = []
    exon_list_SpV2 = []
    exon_list_SpV1_pre = []
    exon_list_SpV2_pre = []

    for item in pos_dict[item_2[0]].split("*")[2].split(","):
        if pos_dict[item_2[0]].split("*")[0] == "+":
            exon_list_SpV1.append([int(item.split(":")[0]), int(item.split(":")[1])])
        elif pos_dict[item_2[0]].split("*")[0] == "-":
            exon_list_SpV1_pre.append([int(item.split(":")[1])*(-1), int(item.split(":")[0])*(-1)])

    if pos_dict[item_2[0]].split("*")[0] == "-":
            exon_list_SpV1 = exon_list_SpV1_pre[::-1]

    for item in pos_dict[item_2[1]].split("*")[2].split(","):
        if pos_dict[item_2[1]].split("*")[0] == "+":
            exon_list_SpV2.append([int(item.split(":")[0]), int(item.split(":")[1])])
        elif pos_dict[item_2[1]].split("*")[0] == "-":
            exon_list_SpV2_pre.append([int(item.split(":")[1])*(-1), int(item.split(":")[0])*(-1)])

    if pos_dict[item_2[1]].split("*")[0] == "-":
            exon_list_SpV2 = exon_list_SpV2_pre[::-1]

    AS_call(exon_list_SpV1,exon_list_SpV2,keys_1,keys_2,pos_dict[item_2[0]].split("*")[1])

f_2.close()
