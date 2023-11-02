#!/usr/bin/python

import numpy as np
import sklearn as skl
import re
import sys, math


def read_fimo_file(file1):
    """
    HM05227 HUM     23      38      -       15.3496 2.79e-06                GGCGGAGCGCCGGGGA
    HM00132 HUM     24      38      +       14.3669 1.79e-06                CCCCGGCGCTCCGCC
    """
    f = open(file1, 'r')
    i = 0
    shape_list = []
    sdict = {}
    hdict, mdict = {}, {}
    lh, lm = 0, 0
    for line in f:
        line = line.strip('\n')
        sent = line.split('\t')
        sdict[sent[1]] = 1
        if len(sdict) == 1:
            hdict[lh] = sent
            hdict[lh][2], hdict[lh][3] = int(hdict[lh][2]), int(hdict[lh][3])
            lh = lh + 1
        else:
            mdict[lm] = sent
            mdict[lm][2], mdict[lm][3] = int(mdict[lm][2]), int(mdict[lm][3])
            lm = lm + 1
    f.close()
    return(hdict, mdict)


def in_block(st1, en1, st2, en2, cutoff = 0.01):
    flag = 0
    if((st2 < en1) and 1.0*(en1 - st2 + 1)/(en1 - st1 + 1) >= cutoff and 1.0*(en1 - st2 + 1)/(en2 - st2 + 1) >= cutoff):
        flag = 1
    return(flag)


def motif_to_block(mdict1, cutoff = 0.01):
    num1 = len(mdict1)
    bdict1, bdict2 = {}, {}
    index1, st1, en1 = -1, -1, -1
    blist = []
    for i in range(num1):
        sent = mdict1[i]
        st2, en2 = int(sent[2]), int(sent[3])
        if index1 == -1:
            index1 = index1 + 1
            bdict1[index1] = [i]
            st1, en1 = st2, en2
        else:
            flag = 1
            for k in bdict1[index1]:
                sent2 = mdict1[k]
                st3, en3 = int(sent2[2]), int(sent2[3])
                if in_block(st3, en3, st2, en2, cutoff) == 0:
                    flag = 0
                    break
            if flag == 1:
                bdict1[index1].append(i)
                en1 = max(en1, en2)
            else:
                bdict2[index1] = [st1, en1]
                index1 = index1 + 1
                bdict1[index1] = [i]
                st1, en1 = st2, en2
    if en1 >= 0:
        bdict2[index1] = [st1, en1]
    return(bdict1, bdict2)


def motif_to_block0(mdict1, cutoff = 0.01):
    num1 = len(mdict1)
    bdict1, bdict2 = {}, {}
    index1, st1, en1 = -1, -1, -1
    blist = []
    for i in range(num1):
        sent = mdict1[i]
        st2, en2 = int(sent[2]), int(sent[3])
        if in_block(st1, en1, st2, en2, cutoff) == 0:
            if en1 >= 0:
                bdict2[index1] = [st1, en1]
            index1 = index1 + 1
            if index1 >= 0:
                bdict1[index1] = [i]
                st1, en1 = st2, en2
        else:
            bdict1[index1].append(i)
            en1 = max(en1, en2)
    if en1 >= 0:
        bdict2[index1] = [st1, en1]
    return(bdict1, bdict2)


def cal_motif_score(p1, p2):
    score1 = -(math.log10(p1) + math.log10(p2))/2.0
    return(score1)


def score_block(blist1, blist2, mdict1, mdict2):
    msdict1, msdict2 = {}, {}
    for ke1 in blist1:
        sent1 = mdict1[ke1]
        sid = sent1[0] + "|" + sent1[4]
        if sid in msdict1:
            if - math.log10(float(sent1[6])) > msdict1[sid]:
                msdict1[sid] = - math.log10(float(sent1[6]))
        else:
            msdict1[sid] = - math.log10(float(sent1[6]))
    for ke2 in blist2:
        sent2 = mdict2[ke2]
        sid = sent2[0] + "|" + sent2[4]
        if sid in msdict2:
            if - math.log10(float(sent2[6])) > msdict2[sid]:
                msdict2[sid] = - math.log10(float(sent2[6]))
        else:
            msdict2[sid] = - math.log10(float(sent2[6]))
    score = 0.0
    flag, num = 0, 0
    for ke1 in msdict1:
        if ke1 in msdict2:
            flag = 1
            num = num + 1
            if (msdict1[ke1] + msdict2[ke1])/2 > score:
                score = (msdict1[ke1] + msdict2[ke1])/2
                #score = score + (msdict1[ke1] + msdict2[ke1])/2
    if num == 0:
        score = 0.0
        for ke1 in msdict1:
            if msdict1[ke1] > score:
                score = msdict1[ke1]
        for ke2 in msdict2:
            if msdict2[ke2] > score:
                score = msdict2[ke2]
        score = - score
    else:
        score = score + num
    return(score)



def cal_sigmoid(x):
    max1, min1 = 10.0, -10.0
    M1, m1 = 6.0, -6.0
    Slop = 1
    n_score = 2*1/(1 + math.exp(- Slop * ((x - 4.0 - min1) *(M1 - m1)/(max1 - min1) + m1))) - 1.0;
    return(n_score)


def cal_matrix(bdict1, bdict2, mdict1, mdict2):
    len1, len2 = len(bdict1), len(bdict2)
    mat1 = [[0 for j in range(len2)] for i in range(len1)]
    mat2 = [[0 for j in range(len2+1)] for i in range(len1+1)]
    for i in range(len1):
        for j in range(len2):
            mat1[i][j] = cal_sigmoid(score_block(bdict1[i], bdict2[j], mdict1, mdict2))
    for i in range(len1):
        for j in range(len2):
            mat2[i+1][j+1] = max(mat2[i][j+1], mat2[i+1][j], mat1[i][j] + mat2[i][j])
    return(mat1, mat2)


def trace_back(mat1, mat2):
    len1, len2 = len(mat1), len(mat1[0])
    i, j = len1, len2
    V, W = mat1, mat2
    p1, p2, r = 0, 0, 0
    pa, pb, la, lb = [], [], [], []
    if W[i][j] == W[i-1][j]:
        while W[i][j] == W[i-1][j]:
            i = i - 1 
    if W[i][j] == W[i][j-1]:
        while W[i][j] == W[i][j-1]:
            j = j - 1
    while (i > 0 and j > 0):
        if( W[i][j] == W[i-1][j-1] + V[i-1][j-1] ):
            pa.append(i - 1); pb.append(j - 1); p1 = p1 + 1; i = i - 1; j = j - 1
        elif( W[i][j] == W[i-1][j] ):
            pa.append(i - 1); pb.append(-1); p1 = p1 + 1; i = i - 1
        else:
            pa.append(-1); pb.append(j - 1); p1 = p1 + 1; j = j - 1
    while( pa[p1-1] == -1 or pb[p1-1] == -1):
        p1 = p1 - 1
    while( pa[r] == -1 or pb[r] == -1 ):
        r = r + 1
    for p2 in range(r, p1):
        la.append(pa[p2]); lb.append(pb[p2])
    la.reverse()
    lb.reverse()
    return(la, lb)


def common_block(blist1, blist2, mdict1, mdict2):
    """
    HM05227 HUM     23      38      -       15.3496 2.79e-06                GGCGGAGCGCCGGGGA
    HM00132 HUM     24      38      +       14.3669 1.79e-06                CCCCGGCGCTCCGCC
    """
    msdict1, msdict2 = {}, {}
    nblist1, nblist2 = [], []
    st1, en1, st2, en2 = -1, -1, -1, -1
    for ke1 in blist1:
        sent1 = mdict1[ke1]
        sid = sent1[0] + "|" + sent1[4]
        if sid in msdict1:
            if - math.log10(float(sent1[6])) > msdict1[sid][0]:
                msdict1[sid] = [- math.log10(float(sent1[6])), ke1, int(sent1[2]), int(sent1[3])]
        else:
            msdict1[sid] = [- math.log10(float(sent1[6])), ke1, int(sent1[2]), int(sent1[3])]
    for ke2 in blist2:
        sent2 = mdict2[ke2]
        sid = sent2[0] + "|" + sent2[4]
        if sid in msdict2:
            if - math.log10(float(sent2[6])) > msdict2[sid][0]:
                msdict2[sid] = [- math.log10(float(sent2[6])), ke2, int(sent2[2]), int(sent2[3])]
        else:
            msdict2[sid] = [ - math.log10(float(sent2[6])), ke2, int(sent2[2]), int(sent2[3])]
    score = 0.0
    flag, num = 0, 0
    for ke1 in msdict1:
        if ke1 in msdict2:
            flag = 1
            num = num + 1
            sent1 = msdict1[ke1]
            sent2 = msdict2[ke1]
            nblist1.append(sent1[1])
            nblist2.append(sent2[1])
            if st1 < 0:
                st1, en1, st2, en2 = sent1[2], sent1[3], sent2[2], sent2[3]
            else:
                st1, en1, st2, en2 = min(sent1[2], st1), max(sent1[3], en1), min(sent2[2], st2), max(sent2[3], en2)
    return(nblist1, nblist2, st1, en1, st2, en2)


def find_com(la, lb, bdict1, bdict2, mdict1, mdict2):
    bdict1_list, bdict2_list = {}, {}
    bdict1_corr, bdict2_corr = {}, {}
    nla, nlb = [], []
    len1 = len(la)
    for i in range(len1):
        if( la[i] > -1 and lb[i] > -1 ):
            nla.append(la[i])
            nlb.append(lb[i])
            nblist1, nblist2, st1, en1, st2, en2 = common_block(bdict1[la[i]], bdict2[lb[i]], mdict1, mdict2)
            bdict1_list[la[i]] = nblist1
            bdict2_list[lb[i]] = nblist2
            bdict1_corr[la[i]] = [st1, en1]
            bdict2_corr[lb[i]] = [st2, en2]
    return(nla, nlb, bdict1_list, bdict2_list, bdict1_corr, bdict2_corr)    


def cal_oud(listx, listy):
    dis = 0.0
    len1 = len(listx)
    for i in range(1, len1):
        dis = dis + math.pow((listx[i] - listx[i - 1]) - (listy[i] - listy[i - 1]),2)
    return(1.0 * dis)

def cal_oud2(listx, listy):
    dis = 0.0
    len1 = len(listx)
    if(len1 >= 2):
        for i in range(1, len1):
            dis = dis + math.pow((listx[i] - listx[i - 1]) - (listy[i] - listy[i - 1]),2)
        dis = math.pow(dis / (len1 - 1), 0.5)
    return(1.0 * dis)

def cal_ham(listx, listy):
    dis = 0.0
    len1 = len(listx)
    for i in range(1, len1):
        dis = dis + abs((listx[i] - listx[i - 1]) - (listy[i] - listy[i - 1]))
    return(1.0 * dis)


def Sscore(la, lb, mat1):
    score = 0.0
    len1 = len(la)
    for i in range(len1):
        if( la[i] > -1 and lb[i] > -1 ):
            score = score + mat1[la[i]][lb[i]]
    return(1.0 * score)


def trimmer(la, lb, bdict1, bdict2, mdict1, mdict2, mat1, trimper):
    num1 = 0
    len1 = len(la)
    ta, tb = [], []
    for i in range(len1):
        #print(str(la[i])+":"+str(lb[i])+"|")
        if( la[i] > -1 and lb[i] > -1 ):
            num1 = num1 + 1
    #print("\n")
    if num1 <= 3:
        ta = [k for k in la]
        tb = [k for k in lb]
        nla, nlb, bdict1_list, bdict2_list, bdict1_corr, bdict2_corr = find_com(la, lb, bdict1, bdict2, mdict1, mdict2)
        return(ta, tb, bdict1_list, bdict2_list, bdict1_corr, bdict2_corr)
    else:
        trinum = 1.0 * num1 * trimper;
        nla, nlb, bdict1_list, bdict2_list, bdict1_corr, bdict2_corr = find_com(la, lb, bdict1, bdict2, mdict1, mdict2)
        if trinum <= 1:
            trinum = 1
        sla = [bdict1_corr[nla[i]][0] for i in range(num1)]
        slb = [bdict2_corr[nlb[i]][0] for i in range(num1)]
        sumdis = cal_ham(sla, slb); sumsco = Sscore(nla, nlb, mat1)
        #print(str(sumdis) + "|" + str(sumsco) + "|" + str(len(sla)))
        outpri = ""
        for k in range(1,len(sla)):
            outpri = outpri + "|" + str(sla[k]-sla[k-1]) + "|" + str(slb[k]-slb[k-1])
        #print(outpri)
        outpri = ""
        for k in range(len(nla)):
            outpri = outpri + "|" + str(mat1[nla[k]][nlb[k]])
        #print(outpri)
        i = 0; j = num1 - 1
        for k in range(int(trinum)):
            if ( cal_ham(sla[0:(i+2)], slb[0:(i+2)])/Sscore(nla[0:(i+1)], nlb[0:(i+1)], mat1) > cal_ham(sla[(j-1):num1], slb[(j-1):num1])/Sscore(nla[j:num1], nlb[j:num1], mat1)):
                i = i + 1
            else:
                j = j - 1
        #print(str(i) + "|" + str(j) + "|" + str(num1))
        i = i - 1
        while i >= 0:
            s01 = cal_ham(sla[0:(i+2)], slb[0:(i+2)]) / sumdis
            s02 = Sscore(nla[0:(i+1)], nlb[0:(i+1)], mat1) / sumsco
            s03 = cal_ham(sla[i:(i+2)], slb[i:(i+2)]) / sumdis
            s04 = Sscore(nla[i:(i+1)], nlb[i:(i+1)], mat1) / sumsco
            #print(str(s01) + "|" + str(s02) + "|" + str(s03) + "|" + str(s04))
            if (cal_ham(sla[0:(i+2)], slb[0:(i+2)]) / sumdis  > Sscore(nla[0:(i+1)], nlb[0:(i+1)], mat1) / sumsco) and (cal_ham(sla[i:(i+2)], slb[i:(i+2)]) / sumdis  > Sscore(nla[i:(i+1)], nlb[i:(i+1)], mat1) / sumsco):
                break
            else:
                i = i - 1
        #print(str(i) + "|" + str(j) + "|" + str(num1))
        j = j + 1
        while j <= num1 - 1:
            s01 = cal_ham(sla[(j-1):num1], slb[(j-1):num1]) / sumdis
            s02 = Sscore(nla[j:num1], nlb[j:num1], mat1) / sumsco
            s03 = cal_ham(sla[(j-1):(j+1)], slb[(j-1):(j+1)]) / sumdis
            s04 = Sscore(nla[j:(j+1)], nlb[j:(j+1)], mat1) / sumsco
            #print(str(s01) + "|" + str(s02) + "|" + str(s03) + "|" + str(s04))
            if (cal_ham(sla[(j-1):num1], slb[(j-1):num1]) / sumdis  > Sscore(nla[j:num1], nlb[j:num1], mat1) / sumsco) and (cal_ham(sla[(j-1):(j+1)], slb[(j-1):(j+1)]) / sumdis  > Sscore(nla[j:(j+1)], nlb[j:(j+1)], mat1) / sumsco):
                break
            else:
                j = j + 1
        #print(str(i) + "|" + str(j) + "|" + str(num1))
        tris = i + 1
        trie = num1 - j
        #print(str(tris) + "--" + str(trie))
        i = 0; j = len1 - 1; k = 0
        while k < tris:
            if( la[i] > -1) and (lb[i] > -1 ):
                k = k + 1
            i = i + 1
        while (la[i] == -1) or (lb[i] == -1):
            i = i + 1
        k = 0
        while k < trie:
            if( la[j] > -1) and (lb[j] > -1 ):
                k = k + 1
            j = j - 1
        while (la[j] == -1) or (lb[j] == -1):
            j = j - 1
        for k in range(i,j+1):
            ta.append(la[k])
            tb.append(lb[k])
        return(ta, tb, bdict1_list, bdict2_list, bdict1_corr, bdict2_corr)


def block_name(name1, index1, num1):
    return(name1 + "0"*(num1-len(str(index1))) + str(index1))


def output_file(ta, tb, bdict1, bdict2, bdict1_list, bdict2_list, bdict1_corr, bdict2_corr, mdict1, mdict2, mat1, out_pre):
    firname, secname = mdict1[0][1], mdict2[0][1]
    outfile1 = out_pre + "." + firname + ".fimo.s"
    outfile2 = out_pre + "." + secname + ".fimo.s"
    outfile3 = out_pre + ".align"
    outfile4 = out_pre + ".dyn"

    sla, slb = [], []
    for i in range(len(ta)):
        if( ta[i] > -1 and tb[i] > -1 ):
            sla.append(bdict1_corr[ta[i]][0])
            slb.append(bdict2_corr[tb[i]][0])
    #sumdis = cal_ham(sla, slb)
    sumdis = cal_oud2(sla, slb)

    Score1 = Sscore(ta, tb, mat1)
    out1 = open(outfile1, 'w')
    out2 = open(outfile2, 'w')
    out3 = open(outfile3, 'w')
    out4 = open(outfile4, 'w')

    out4.write(">>>Score%f\n" % (round(Score1,3)))
    outstr4 = ""
    for i in range(len(ta)):
        if ta[i] == -1:
            outstr4 = outstr4 + "------"
        else:
            outstr4 = outstr4 + block_name("HUM", ta[i], 3)
    out4.write(outstr4 + "\n")
    outstr4 = ""
    for i in range(len(tb)):
        if tb[i] == -1:
            outstr4 = outstr4 + "------"
        else:
            outstr4 = outstr4 + block_name("MOU", tb[i], 3)
    out4.write(outstr4 + "\n\n")
    
    for i in range(len(ta)):
        outstr4 = ""
        if ta[i] == -1:
            outstr4 = outstr4 + "------:------\t"
        else:
            outstr4 = outstr4 + block_name("HUM", ta[i], 3) + ":"
            sent1 = bdict1[ta[i]]
            for k in sent1:
                sent2 = mdict1[k]
                outstr4 = outstr4 + sent2[0] + "_" + sent2[4] + ","
            outstr4 = outstr4 + "\t"
        if tb[i] == -1:
            outstr4 = outstr4 + "------:------\n"
        else:
            outstr4 = outstr4 + block_name("MOU", tb[i], 3) + ":"
            sent1 = bdict2[tb[i]]
            for k in sent1:
                sent2 = mdict2[k]
                outstr4 = outstr4 + sent2[0] + "_" + sent2[4] + ","
            outstr4 = outstr4 + "\n"
        out4.write(outstr4)
    out4.close()
    for i in range(len(ta)):
        if( ta[i] > -1 and tb[i] > -1 ):
            sent1 = bdict1_list[ta[i]]
            block_id = block_name("HUM", ta[i], 3)
            for k in sent1:
                sent2 = mdict1[k]
                outstr1 = sent2[0] + "\t" + sent2[1] + "\t" + str(sent2[2]) + "\t" + str(sent2[3]) + "\t" + sent2[4] + "\t" + sent2[5] + "\t" + sent2[6] + "\t" + block_id + "\t" + sent2[7] + "\n"
                out1.write(outstr1)
            sent1 = bdict2_list[tb[i]]
            block_id = block_name("MOU", tb[i], 3)
            for k in sent1:
                sent2 = mdict2[k]
                outstr2 = sent2[0] + "\t" + sent2[1] + "\t" + str(sent2[2]) + "\t" + str(sent2[3]) + "\t" + sent2[4] + "\t" + sent2[5] + "\t" + sent2[6] + "\t" + block_id + "\t" + sent2[7] + "\n"
                out2.write(outstr2)
    out1.close()
    out2.close()
    outstr3 = out_pre + "\t" + "0001" + "\t" + str(round(Score1,3)) + "\t" + str(round(sumdis,3))
    st1, en1, st2, en2, num1, k1, max_num = -1, -1, -1, -1, 0, 0, 0
    for i in range(len(ta)):
        if( ta[i] > -1 and tb[i] > -1 ):
            k1 = k1 + 1
            sent1 = bdict1_list[ta[i]]
            sent2 = bdict2_list[tb[i]]
            num1 = num1 + len(sent1)
            max_num = max(max_num, len(sent1))
            if st1 == -1:
                st1, en1, st2, en2 = bdict1_corr[ta[i]][0], bdict1_corr[ta[i]][1], bdict2_corr[tb[i]][0], bdict2_corr[tb[i]][1]
            else:
                st1, en1, st2, en2 = min(bdict1_corr[ta[i]][0], st1), max(bdict1_corr[ta[i]][1], en1), min(bdict2_corr[tb[i]][0], st2), max(bdict2_corr[tb[i]][1], en2)
            sla.append(bdict1_corr[ta[i]][0])
            slb.append(bdict2_corr[tb[i]][0])
    outstr3 = outstr3 + "\t" + str(st1) + "_" + str(en1) + "\t" + str(st2) + "_" + str(en2) + "\t" + str(num1) + "\t" + str(k1) + "\t" + str(max_num) + "\n"
    out3.write(outstr3)
    outstr1, outstr2 = "", ""
    for i in range(len(ta)):
        if( ta[i] > -1 and tb[i] > -1 ):
            sent1 = bdict1_list[ta[i]]
            block_id = block_name("HUM", ta[i], 3)
            for k in sent1:
                sent2 = mdict1[k]
                outstr1 = outstr1 + block_id + "_" + sent2[0] + "_" + sent2[4] + "_" + str(round(-math.log10(float(sent2[6])),3)) + " "    
            sent1 = bdict2_list[tb[i]]
            block_id = block_name("MOU", tb[i], 3)
            for k in sent1:
                sent2 = mdict2[k]
                outstr2 = outstr2 + block_id + "_" + sent2[0] + "_" + sent2[4] + "_" + str(round(-math.log10(float(sent2[6])),3)) + " "
    out3.write(outstr1 + "\n")
    out3.write(outstr2 + "\n")
    outstr1, outstr2 = "", ""
    for i in range(len(ta)):
        if( ta[i] > -1 and tb[i] > -1 ):
            sent1 = bdict1_corr[ta[i]]
            block_id = block_name("HUM", ta[i], 3)
            #for k in sent1:
            #sent2 = mdict1[k]
            outstr1 = outstr1 + block_id + "_" + str(sent1[0]) + " "    
            sent2 = bdict2_corr[tb[i]]
            block_id = block_name("MOU", tb[i], 3)
            #for k in sent1:
            #sent2 = mdict2[k]
            outstr2 = outstr2 + block_id + "_" + str(sent2[0]) + " "
    out3.write(outstr1 + "\n")
    out3.write(outstr2 + "\n")
    out3.close()

if __name__ == '__main__':
    input_file = sys.argv[1]
    block_rate = float(sys.argv[2])
    trim_rate = float(sys.argv[3])
    #output_pre = sys.argv[3]
    hdict1, mdict1 = read_fimo_file(input_file)
    hblock1, hblock2 = motif_to_block(hdict1, block_rate)
    mblock1, mblock2 = motif_to_block(mdict1, block_rate)
    #out1 = open("test1.txt", 'w')
    """
    HM05227 HUM     23      38      -       15.3496 2.79e-06                GGCGGAGCGCCGGGGA
    HM00132 HUM     24      38      +       14.3669 1.79e-06                CCCCGGCGCTCCGCC
    """
    #for i in hblock1:
        #sent1 = hblock1[i]
        #outstr = str(i) + "|"
        #outstr2 = str(i) + "|"
        #for k in sent1:
        #    sent2 = hdict1[k]
        #    outstr = outstr + sent2[0] + "_" + str(sent2[2]) + "_" + str(sent2[3]) + "_" + sent2[4] + "|"
        #sent1 = hblock2[i]
        #outstr2 = outstr2 + str(sent1[0]) + "|" + str(sent1[1])
        #out1.write(outstr + "\n" + outstr2 + "\n")
    #out1.close()
    score_mat1, score_mat2 = cal_matrix(hblock1, mblock1, hdict1, mdict1)
    laba, labb = trace_back(score_mat1, score_mat2)
    ta, tb, hdict_list, mdict_list, hdict_corr, mdict_corr = trimmer(laba, labb, hblock1, mblock1, hdict1, mdict1, score_mat1, trim_rate)
    output_file(ta, tb, hblock1, mblock1, hdict_list, mdict_list, hdict_corr, mdict_corr, hdict1, mdict1, score_mat1, input_file)



