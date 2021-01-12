from STree4CS import STree4CS as STree
from datavisual import *
from alignment import *

data = read_from_file_with_enter(data_path)
keys,samples= statistics_length(data)

#['MW411947.1', 'MW411948.1', 'MW411949.1', 'MT232662.1', 'MT232664.1', 'MW403692.1', 'MW403693.1', 'MW403694.1', 'MW403695.1', 'MW403696.1', 'MW406485.1', 'MW406487.1', 'MW406488.1', 'MW406490.1', 'MW406491.1', 'MW320729.1', 'MW320730.1', 'MW320731.1', 'MW320733.1', 'MW320735.1', 'MW064259.1', 'MW064260.1', 'MW064263.1', 'MW064264.1', 'MW342706.1', 'MW404672.1', 'MW404673.1', 'MW365356.1', 'MW365357.1']
#print(samples[0])
#print(list(samples[0]))
# st = STree.STree4CS([list(samples[0]),list(samples[1])])

# print(st.lcsEndStart())
#[list(samples[0]),list("ACCTTCCC")]

#找出俩个序列最大公共字符串 起始序列对
def find_commom_substrings(seq1,seq2):
    seqindex=0
    index=0
    resultAlignList=[]
    while(index<len(seq2) and seqindex<len(seq1)):
        st=STree.STree4CS([seq1[seqindex:-1],seq2[index:-1]])
        start,end=st.lcsEndStart()
        length=end-start
        startindex=seqindex+start
        seqindex=seqindex+end

        st1 = STree.STree4CS([seq2[index:-1]])
        start1=st1.find(seq1[startindex:seqindex])+index
        if start1>=index and length>0:
            print("start",startindex,"length",length,"seqindex",seqindex)
            print("start1",start1,"length",length,"index",start1+length)
            resultAlignList.append([startindex,start1,length])
            index=start1+length
        else:
            index=index+1
    return resultAlignList






def listToString(seq):
    return "".join(seq)

#对齐后 对应相同字符 分数为10， score=len(seq)*10
def calculateSameSequence(seq):
    return len(seq)*10

def calculate_score_align(seq1,seq2):
    resultList=find_commom_substrings(list(samples[0]),list(samples[1]))
    score=0
    for i in resultList[1:-1]:
        start1=i[0]
        start2=i[1]
        length=i[2]
        str1="".join(seq1[start1:start1+length])
        str2="".join(seq2[start2:start2+length])
        a,b,c=align(str1,str2)
        print("seq1",b,"\n seq2",c,"\n score",a)
        score=score+a
    return score

print(keys)
#print(find_commom_substrings(list(samples[0]),list(samples[1])))
#print(calculate_score_align(list(samples[0]),list(samples[1])))