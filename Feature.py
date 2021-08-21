#BKF
import pandas as pd
import numpy as np
import itertools
import os
import sys


gene_type="DNA"
type_value="T"
if gene_type=="RNA":
    type_value="U"
elif gene_type=="DNA":
    type_value="T"

def read_fasta_file():
    '''
    used for load fasta data and transformd into numpy.array format
    '''
    fh = open('S3.txt', 'r')
    seq = []
    for line in fh:
        if line.startswith('>'):
            continue
        else:
            seq.append(line.replace('\n', '').replace('\r', ''))
    fh.close()
    matrix_data = np.array([list(e) for e in seq])
    print(matrix_data)
    return matrix_data

def AthMethPre_extract_one_line(data_line):
   
    #A=[0,0,0,1]
    #T=[0,0,1,0]
    #C=[0,1,0,0]
    #G=[1,0,0,0]
    #N=[0,0,0,0]
    feature_representation={"A":A,"C":C,"G":G,"N":N}
    feature_representation[type_value]=U
    beginning=0
    end=len(data_line)-1
    one_line_feature=[]
    alphabet='ACNG'
    alphabet+=type_value
    matrix_two=["".join(e) for e in itertools.product(alphabet, repeat=2)] 
    matrix_three=["".join(e) for e in itertools.product(alphabet, repeat=3)
    matrix_four=["".join(e) for e in itertools.product(alphabet, repeat=4)]
    feature_two=np.zeros(25)
    feature_three=np.zeros(125)
    feature_four=np.zeros(625)
    for index,data in enumerate(data_line):
        if index==beginning or index==end:
            one_line_feature.extend(feature_representation["N"])
        elif data in feature_representation.keys():
            one_line_feature.extend(feature_representation["N"])
            one_line_feature.extend(feature_representation[data])
        if "".join(data_line[index:(index+2)]) in matrix_two and index <= end-1:
            feature_two[matrix_two.index("".join(data_line[index:(index+2)]))]+=1
        if "".join(data_line[index:(index+3)]) in matrix_three and index <= end-2:
            feature_three[matrix_three.index("".join(data_line[index:(index+3)]))]+=1
        if "".join(data_line[index:(index+4)]) in matrix_four and index <=end-3:
            feature_four[matrix_four.index("".join(data_line[index:(index+4)]))]+=1
    sum_two=np.sum(feature_two)
    sum_three=np.sum(feature_three)
    sum_four=np.sum(feature_four)
    one_line_feature.extend(feature_two/sum_two)
    one_line_feature.extend(feature_three/sum_three)
    one_line_feature.extend(feature_four/sum_four)
    return one_line_feature



def AthMethPre_extract_one_line_without(data_line):
   
    A=[1,0,0,0]
    C=[0,1,0,0]
    G=[0,0,1,0]
    T=[0,0,0,1]
    feature_representation={"A":A,"C":C,"G":G,"T":T}
    beginning=0
    end=len(data_line)-1
    one_line_feature=[]
    alphabet='ACGT'
    matrix_two=["".join(e) for e in itertools.product(alphabet, repeat=2)] 
    matrix_three=["".join(e) for e in itertools.product(alphabet, repeat=3)]
    matrix_four=["".join(e) for e in itertools.product(alphabet, repeat=4)]
    feature_two=np.zeros(16)
    feature_three=np.zeros(64)
    feature_four=np.zeros(256)
    for index,data in enumerate(data_line):
        one_line_feature.extend(feature_representation[data])
        if "".join(data_line[index:(index+2)]) in matrix_two and index <= end-1:
            feature_two[matrix_two.index("".join(data_line[index:(index+2)]))]+=1
        if "".join(data_line[index:(index+3)]) in matrix_three and index <= end-2:
            feature_three[matrix_three.index("".join(data_line[index:(index+3)]))]+=1
        if "".join(data_line[index:(index+4)]) in matrix_four and index <=end-3:
            feature_four[matrix_four.index("".join(data_line[index:(index+4)]))]+=1
    sum_two=np.sum(feature_two)
    sum_three=np.sum(feature_three)
    sum_four=np.sum(feature_four)
    one_line_feature.extend(feature_two/sum_two)
    one_line_feature.extend(feature_three/sum_three)
    one_line_feature.extend(feature_four/sum_four)
    return one_line_feature

def AthMethPre_feature_extraction(matrix_data,fill_NA):
    if fill_NA=="1":
        final_feature_matrix=[AthMethPre_extract_one_line(e) for e in matrix_data]
    elif fill_NA=="0":
        final_feature_matrix=[AthMethPre_extract_one_line_without(e) for e in matrix_data]
    return final_feature_matrix


fill_NA="0"
matrix_data=read_fasta_file()
final_feature_matrix=AthMethPre_feature_extraction(matrix_data,fill_NA)
print(np.array(final_feature_matrix).shape)
pd.DataFrame(final_feature_matrix).to_excel('D:/DNA起始位点/S3/BKF.xlsx',header=None,index=False)

#RFHCP
import pandas as pd
import numpy as np
import os
import sys
import itertools




gene_type='DNA'
fill_NA='0'


if gene_type=="RNA":
    gene_value="U"
elif gene_type=="DNA":
    gene_value="T"


def convert_with(dataPath,outputPath):
    """RFH feature"""
    lines=open(dataPath).readlines()
    finally_text = open(outputPath, 'w')
    finnaly_lines=""
    for line in lines:
        if line.strip()=="":continue
        if line.strip()[0] in ['A','G','C',gene_value,'N']:
            position_mark=0
            count_AGCT=[0,0,0,0,0]
            temp = ""
            for x in list(line.strip()):
                position_mark+=1
                if x=="A" or x=="G":temp+="1,"
                else:temp+="0,"
                if x=="A" or x==gene_value:temp+="1,"
                else:temp+="0,"
                if x == "A" or x == "C":temp+="1,"
                else:temp+="0,"
                if x == "A":
                    count_AGCT[0] += 1
                    temp +=str(round(count_AGCT[0] / position_mark*1.0,2))
                    temp+=','
                elif x == "G":
                    count_AGCT[1] += 1
                    temp +=str(round(count_AGCT[1] / position_mark*1.0,2))
                    temp += ','
                elif x == "C":
                    count_AGCT[2] += 1
                    temp +=str(round(count_AGCT[2] / position_mark*1.0,2))
                    temp += ','
                elif x == gene_value:
                    count_AGCT[3] += 1
                    temp +=str(round(count_AGCT[3] / position_mark*1.0,2))
                    temp += ','
                elif x == "N":
                    count_AGCT[4] += 1
                    temp +=str(round(count_AGCT[4] / position_mark*1.0,2))
                    temp += ','

            finnaly_lines+=((temp[:len(temp)-1])+'\n')
            #finally_text.write(temp+'\n')
    finally_text.writelines(finnaly_lines)
    finally_text.close()
    
def convert_without(dataPath,outputPath):
    """RFH feature"""
    lines=open('S1.txt').readlines()
    finally_text = open('D:/DNA起始位点/S1/RFHCP1.csv', 'w')
    finnaly_lines=""
    for line in lines:
        if line.strip()=="":continue
        if line.strip()[0] in ['A','G','C',gene_value]:
            position_mark=0
            count_AGCT=[0,0,0,0]
            temp = ""
            for x in list(line.strip()):
                position_mark+=1
                if x=="A" or x=="G":temp+="1,"
                else:temp+="0,"
                if x=="A" or x==gene_value:temp+="1,"
                else:temp+="0,"
                if x == "A" or x == "C":temp+="1,"
                else:temp+="0,"
                if x == "A":
                    count_AGCT[0] += 1
                    temp +=str(round(count_AGCT[0] / position_mark*1.0,2))
                    temp+=','
                elif x == "G":
                    count_AGCT[1] += 1
                    temp +=str(round(count_AGCT[1] / position_mark*1.0,2))
                    temp += ','
                elif x == "C":
                    count_AGCT[2] += 1
                    temp +=str(round(count_AGCT[2] / position_mark*1.0,2))
                    temp += ','
                elif x == gene_value:
                    count_AGCT[3] += 1
                    temp +=str(round(count_AGCT[3] / position_mark*1.0,2))
                    temp += ','

            finnaly_lines+=((temp[:len(temp)-1])+'\n')
    finally_text.writelines(finnaly_lines)
    finally_text.close()
    

if fill_NA=="1":
    convert_with(' .txt',path+outputname)
    data=pd.read_csv(path+outputname,header=None,index_col=False)
    print(data.values.shape)
elif fill_NA=="0":
convert_without('S1.txt','D:/DNA起始位点/S1/RFHCP.csv')
