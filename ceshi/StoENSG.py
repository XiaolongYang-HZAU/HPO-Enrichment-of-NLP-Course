#生信1802杨晓龙 原创
#2021/5/12
#将symbol编号转化为ENSG编号
import mygene
mg = mygene.MyGeneInfo()

symb = []
with open("symbol_Parkinson.txt") as f: #读入symbol文件
    for l in f:
        symb.append(l.strip())
#print(symb)
out = mg.querymany(symb, scopes='symbol', fields='ensembl.gene', species='human') 
#print(out)

result = []
for i in out:                      #获取对应的ENSG编号列表
    if 'ensembl' in i.keys():
        genelist = i['ensembl']
        if type(genelist) == list:
            for j in genelist:
                result.append(j['gene'])
        else :
            result.append(genelist['gene'])

with open("ENSGresult.txt", "w") as output: #保存为txt文件
    output.write(str(result))