
import pandas as pd
import matplotlib.pyplot as plt
 
import xlwt
# import seaborn as sns

 
data_path = "./data/sars-cov2.fasta"

 
def read_from_file_with_enter(filename):
    fr = open(filename,'r')
    key=""
    value=""
    samples = {}
    for line in fr:
        if line.startswith('>'):
            if key !="":
                samples[key]=value
            key=line
            value=""
            continue
        value += line[:-1]
    return samples

 
def read_from_file(filename):
    fr = open(filename, 'r')
    sample = ""
    samples = []
    for line in fr:
        if line.startswith('>'):
            if sample != "":
                samples.append(sample)
                sample=""
            continue
        sample+=line[:-1]
    return samples
 
def statistics_lengths(samples):
    keys=[]
    length=[]
    for key in samples:
        k=key.split(" ")[0].split(">")[1]
        l=len(samples[key])
        keys.append(k)
        length.append(l)
    return keys,length
def statistics_length(samples):
    keys=[]
    sample=[]
    for key in samples:
        k=key.split(" ")[0].split(">")[1]
        keys.append(k)
        sample.append(samples[key])
    return keys,sample

def getResultData():
    fr = open("./data/result.txt",'r')
    score=[]
    memory=[]
    time=[]
    for line in fr:
        #save to file: MW411947.1 VS MW411949.1 result:   score:298219    memory:31964.0  time:11.2022803
        if line.find("save to file")!=-1:
            s=float(line.split("score:")[-1].split(" ")[0])
            m=float(line.split("memory:")[-1].split(" ")[0])
            t=float(line.split("time:")[-1].split(" ")[0])
            score.append(abs(s))
            memory.append(abs(m))
            time.append(abs(t))
    # print("score:",score)
    # print("memory:",memory)
    # print("time:",time)
    return score,memory,time
#
def writeToExcel():
    itemname=['MW411947.1', 'MW411948.1', 'MW411949.1', 'MT232662.1', 'MT232664.1', 'MW403692.1', 'MW403693.1', 'MW403694.1', 'MW403695.1', 'MW403696.1', 'MW406485.1', 'MW406487.1', 'MW406488.1', 'MW406490.1', 'MW406491.1', 'MW320729.1', 'MW320730.1', 'MW320731.1', 'MW320733.1', 'MW320735.1', 'MW064259.1', 'MW064260.1', 'MW064263.1', 'MW064264.1', 'MW342706.1', 'MW404672.1', 'MW404673.1', 'MW365356.1', 'MW365357.1']
    count=len(itemname)
    book = xlwt.Workbook() 
    score, memory,time=getResultData()
    datal=["score", "memory","time"]
    #创建表单
    for index in range(0,3):
        sheet = book.add_sheet(datal[index],cell_overwrite_ok=True)
        for i in range(0,len(itemname)):
            sheet.write(0,i+1,itemname[i])
            sheet.write(i+1,0,itemname[i])
        data=[]
        if index==0:
            data=score
        elif index==1:
            data=memory
        else:
            data=time
        k=0
        for i in range(0,count):
            for j in range(i+1,count):
                sheet.write(i+1,j+1,data[k])
                sheet.write(j+1,i+1,data[k])
                k=k+1
        for i in range(0,count):
            sheet.write(i+1,i+1,0)
    book.save('./data/resultMatrix.xls')


if __name__ == "__main__":
    # data = read_from_file_with_enter(data_path)
    # keys,length= statistics_lengths(data)
    # print(keys)
    # print(length)
    writeToExcel()
    
    # plt.figure(figsize=(6,4))

    # x = range(len(keys))
    # y = length
    # plt.plot(x, y, 'ro-')
    # plt.xticks(x, keys, rotation=45)
    # plt.margins(0.08)
    # plt.subplots_adjust(bottom=0.15)
    # plt.show()


# plt.figure(figsize=(6,4))

# g = sns.distplot(a=length,label="GenData",kde=False,bins=20)

# plt.ylabel("Number")
# plt.title("Length of sample")
# plt.legend(loc='best')
 
# plt.figure(figsize=(6,4))
# g = sns.kdeplot(data=length,shade=True)
# plt.title("Kernel density estimation")
 
# plt.show(g)