from alignment import *
from STree4CS import STree4CS as STree
from datavisual import *
import psutil
import os
from time import clock


class PairsequenceItem:
    """
    the base pair genome sequence ,offer function as follows:
    function: 
        getPairScoreAlignSequence:
            return score of two aligned sequence, and the aligned two sequence
    """
    def __init__(self,seq1="",seq2="",start1=0,start2=0,common=False):
        self.start1=start1
        self.start2=start2
        if len(seq1)==0:
            self.seq1=""
        else:
            self.seq1=seq1
        if len(seq2)==0:
            self.seq2=""
        else:
            self.seq2=seq2
        self.common=common

    # return score, s1,s2 aligned
    def getPairScoreAlignSequence(self):
        if self.common:
            return len(self.seq1)*10,self.seq1,self.seq1
        else:
            return align("".join(self.seq1),"".join(self.seq2))


class Node:
    """
        the basic node of trees, and the data is a PairsequenceItem class.
    """
    def __init__(self, data):
        self.left = None
        self.right = None
        self.data = data

def createTree(seq1,seq2):
    """
        create a binary tree using digui algorithm, and first find longest comom sequence 
        and using the LCS to split original sequence.
        exit requirement: then the sequence length less than 200;
    """
    if len(seq1)==0 and len(seq2)==0:
        return None
    elif len(seq1)>200 and len(seq2)>200:  
        try:
            #print("seq1",len(seq1),"seq2",len(seq2))
            st=STree.STree4CS([seq1,seq2])
            #print("FindLSC")
            start,end=st.lcsEndStart()
            st1 = STree.STree4CS([seq2])
            start1=st1.find(seq1[start:end])
            #print("construct PairsequenceItem")
            pair=PairsequenceItem(seq1[start:end],seq2[start1:start1+end-start],start,start1,common=True)
            root=Node(data=pair)
            #print("construct left node")
            root.left=createTree(seq1[0:start],seq2[0:start1])
            #print("construct right node")
            root.right=createTree(seq1[end:-1],seq2[start1+end-start:-1])
            #print("return node")
            return root
        except :
            print("error seq1",len(seq1),"seq2",len(seq2))
    else:
        #print("uncommon sequence:\tseq1",seq1,"seq2",seq2)
        pa =PairsequenceItem(seq1,seq2,0,0,common=False)
        return Node(data=pa)


def inorderTraversal(root):   #先把根节点入栈，如果左子树一直不为空，就一直入栈，直到 
                                ##把所有左节点入栈，然后pop栈顶元素，指针指向栈顶元素的右子树
    """
    :type root: TreeNode
    :rtype: List[int]
    :return: ordered sequence aligned segment list
    """
    stack=[]
    res=[]
    if not root:
        return []
    while root or stack:
        while root:
            stack.append(root)
            root=root.left
        if stack:
            a=stack.pop()
            root=a.right
            res.append(a.data)
    return res

def getScoreAlignSequence(datalist):
    """
    :type datalist: List[]
    :rtype: Node
    :return: 
        aligned sequenceA; aligned sequenceB, score
    """
    strA=""
    strB=""
    score=0
    for pair in datalist:
        s,a,b=pair.getPairScoreAlignSequence()
        strA=strA+"".join(a)
        strB=strB+"".join(b)
        score=score+s
    return strA,strB,score
def _time_analyze_():
    """
    :calculate time cost
    """
    start = clock()
    pass
    finish = clock()
    
def memoryMonitor():
    """
    :calculate memory cost
    """
    #print(u'当前进程的内存使用：%.4f GB' % (psutil.Process(os.getpid()).memory_info().rss / 1024 ) )
    #当前进程的内存使用：0.6119 KB
    memomybefore=psutil.Process(os.getpid()).memory_info().rss / 1024
    pass
    memomyAfter=psutil.Process(os.getpid()).memory_info().rss / 1024
    used=memomyAfter-memomybefore
def saveToFile(des,score,seqa,seqb):
    """
    :saveResult to file: filename- ./data/alignResult.txt, mode=a
    MW411949.1VSMT232662.1	 score:295819	 memory:92564.0	 time:24.631434100000092
    seqa: 
    seqb:
    """
    file = open('./data/alignResult.txt', mode='a')
    # 方法1  write 写入
    file.write(des+score+'\n')
    file.write("seqa:\t"+seqa+'\n')
    file.write("seqb:\t"+seqb+'\n')
    file.close()

def handleSingle(seq1,seq2):
    """
        hande a pair of sequence
    """
    root=createTree(seq1,seq2)
    pairs=inorderTraversal(root)
    return getScoreAlignSequence(pairs)

def handleAll():
    """
        hande all a pair of sequence
    """
    data = read_from_file_with_enter(data_path)
    keys,samples= statistics_length(data)
    size=len(samples)
    for i in range(0,size):
        for j in range(i+1,size):
            print("handle sequence i=",i,"\t j=",j,"\t",keys[i],"\t",keys[j])
            memomybefore=psutil.Process(os.getpid()).memory_info().rss / 1024
            start = clock()
            strA,strB,score=handleSingle(list(samples[i]),list(samples[j]))
            memomyAfter=psutil.Process(os.getpid()).memory_info().rss / 1024
            finish = clock()
            analyse="\t score:"+str(score)+"\t memory:"+str(memomyAfter-memomybefore)+"\t time:"+str(finish-start)
            saveToFile(keys[i]+"VS"+keys[j],analyse,strA,strB)
            print("save to file:",keys[i],"VS",keys[j],"result:",analyse)


def testSinglePair():
    """
        test a pair of sequence
    """
    data = read_from_file_with_enter(data_path)
    keys,samples= statistics_length(data)

    memomybefore=psutil.Process(os.getpid()).memory_info().rss / 1024
    start = clock()

    strA,strB,score=handleSingle(list(samples[0]),list(samples[1]))
    memomyAfter=psutil.Process(os.getpid()).memory_info().rss / 1024
    finish = clock()
    analyse="\t score:"+str(score)+"\t memory:"+str(memomyAfter-memomybefore)+"\t time:"+str(finish-start)

    print("strA:",len(strA),"\n",strA)
    print("strB:",len(strB),"\n",strB)
    print("analyse",analyse)
    saveToFile(keys[0]+"VS"+keys[1],analyse,strA,strB)
if __name__ == "__main__":
   handleAll()






