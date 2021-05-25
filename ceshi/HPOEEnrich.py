#生信1802杨晓龙 原创
#2021/5/12
#这个函数能生成main函数需要的前置文件
import pandas as pd
import matplotlib.pyplot as plt
from bokeh.plotting import show
from bokeh.io import output_notebook
output_notebook()
from hpoea.enrich import GSEA
from hpoea.plot import LineagePlot, dot_plot

input_txt = "D:/linuxvs/ziranyuyan/xiangmu/ceshi/ceshi/ENSGresult.txt"
ENSG=[]
with open(input_txt) as f: #读入symbol文件
    for l in f:
        ENSG.append(l.strip())
print(ENSG[0:3])

for i in range(len(ENSG)//5 + 1):
    print(" ".join(ENSG[i*5:(i+1)*5]))

from hpoea.utils.idconvert import EntrezEnsemblConvert
cvt = EntrezEnsemblConvert()
ENTRZ = cvt.ensembl2entrez(ENSG)
print(ENTRZ[0:3])
print(len(ENTRZ))

gsea = GSEA()
gsea.enrich(ENTRZ)
gsea.multiple_test_corretion(method='fdr_bh')
print(gsea.enrichment_table.head(1))
gsea.enrichment_table.shape[0]
gsea.filter(by='padj', threshold=0.01)#筛选padj小于0.01
#print(type(gsea.enrichment_table))
t=gsea.enrichment_table
t.to_csv("enrichiment_table.csv")#保存分析结果表格为enrichiment_table.csv