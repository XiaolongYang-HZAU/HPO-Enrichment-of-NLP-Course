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

input_csv = "D:/linuxvs/ziranyuyan/xiangmu/ceshi/ceshi/SAGD_00055.csv"
exp = pd.read_csv(input_csv)
exp.columns = ["ensembl_id"] + list(exp.columns)[1:]
print(exp.head(3))

sig = exp[exp.padj <= 0.05]
test_genes = list(sig.ensembl_id)
print(len(test_genes))#选择具有显著差异表达的条件基因:padj<= 0.05

for i in range(len(test_genes)//5 + 1):
    print(" ".join(test_genes[i*5:(i+1)*5]))

from hpoea.utils.idconvert import EntrezEnsemblConvert
cvt = EntrezEnsemblConvert()
test_entrez = cvt.ensembl2entrez(test_genes)
#print(test_entrez)
print(len(test_entrez))

gsea = GSEA()
gsea.enrich(test_entrez)
gsea.multiple_test_corretion(method='fdr_bh')
print(gsea.enrichment_table.head(1))
gsea.enrichment_table.shape[0]
gsea.filter(by='padj', threshold=0.01)#筛选padj小于0.01
#print(type(gsea.enrichment_table))
t=gsea.enrichment_table
t.to_csv("enrichiment_table0.01.csv")#保存分析结果表格为enrichiment_table0.01.csv
terms = list(gsea.enrichment_table.HPO_term_ID)
lin = LineagePlot()
fig, ax = plt.subplots(figsize=(20, 10))
lin.plot(terms, ax=ax)
