 

### HPO富集分析的一些参考网站和包

> 往届学长作业：https://github.com/tongliu-liu/HPO-enrichment-analysis
>
> R包：HPOSim：https://pubmed.ncbi.nlm.nih.gov/25664462/#:~:text=The%20Human%20Phenotype%20Ontology%20%28HPO%29%20provides%20a%20standardized,used%20offline%20and%20provide%20only%20few%20similarity%20measures.
>
> **python HPO基因集合富集分析示例：**https://nanguage.github.io/examples/hpo_enrich/example_sagd_00055.html
>
> ![image-20210422184443263](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422184443263.png)
>
> 
>
> 官网，下载HPO数据，查HPO编号和注释用：https://hpo.jax.org/app/
>
> SAGD示例数据库（配合Python包食用）：http://bioinfo.life.hust.edu.cn/SAGD#!/download
>
> 人类疾病数据库（配合R包）：https://omim.org/

老师推荐：HPOterms在文本中得到，与水稻的难度类似？

挑战性：实验设计，获取基因集合，要把研究的问题立起来，要有基因的具体来源。

## 使用python包进行HPO富集分析

> 环境：Windows10下的vscode
>
> python 3.7.4 64-bit based on conda

#### 依赖包的安装

首先要安装git（windows10系统下）

> 参考网址：https://git-scm.com/download/win
>
> https://blog.csdn.net/qq_32786873/article/details/80570783

按照以上网址安装完毕后，将git加入环境变量：![image-20210421142420742](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210421142420742.png)

> 依赖包网址：https://github.com/Nanguage/BioTMCourse/tree/master/HPO%20enrich

``` powershell
git clone https://github.com/Nanguage/BioTMCourse.git
cd  D:\linuxvs\ziranyuyan\xiangmu\ceshi\ceshi\BioTMCourse\HPO enrich
python setup.py install #报错了，找不到request包
pip install request -i http://pypi.douban.com/simple/ --trusted-host pypi.douban.com#终端中安装request
pip install bokeh
#安装完后再运行python setup.py install仍然报相同错误，先不管，尝试运行后续代码
```

#### 依赖包准备

``` python
import pandas as pd
import matplotlib.pyplot as plt
from bokeh.plotting import show
from bokeh.io import output_notebook
output_notebook()

from hpoea.enrich import GSEA
from hpoea.plot import LineagePlot, dot_plot#导入包
```

#### 数据准备

> 本笔记仅展示代码的可行性，不再对输入数据和数据背景，研究价值和具体的结果分析做出介绍，需要结合具体情况进行分析。
>
> 在这个例子中，我们使用性别相关的基因作为输入例子，它从SAGD数据库下载。选择数据组SAGD_00055(人类下丘脑组织)作为输入。

SAGD_00055.csv格式如图，这是基因表达差异分析的结果

![image-20210421153653767](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210421153653767.png)

其中log2FoldChange是对差异倍数取log2的值，这个值为负则对应基因为下调基因（相对表达量低），反之则为上调基因。

> http://www.pinlue.com/article/2019/07/1303/489298644023.html

padj是矫正过后的p值，越小则差异表达越显著。

``` python
input_csv = "D:/linuxvs/ziranyuyan/xiangmu/ceshi/ceshi/SAGD_00055.csv"#读入数据
exp = pd.read_csv(input_csv)
exp.columns = ["ensembl_id"] + list(exp.columns)[1:]
print(exp.head(3))
'''        ensembl_id     baseMean  ...     FPKM_M    FPKM_F
0  ENSG00000012817   850.482453  ...   6.015774  0.020867
1  ENSG00000067048   908.726851  ...  10.234487  0.026688
2  ENSG00000129824  1493.628628  ...  39.948216  0.155030

[3 rows x 9 columns]'''
```

#### 数据筛选

``` python
#选择具有显著差异表达的条件基因:padj<= 0.05
sig = exp[exp.padj <= 0.05]
test_genes = list(sig.ensembl_id)
len(test_genes)#54
```

共选择了54个基因，以下是基因列表：

``` python
for i in range(len(test_genes)//5 + 1):
    print(" ".join(test_genes[i*5:(i+1)*5]))
'''ENSG00000012817 ENSG00000067048 ENSG00000129824 ENSG00000114374 ENSG00000131002
ENSG00000198692 ENSG00000165246 ENSG00000183878 ENSG00000229807 ENSG00000233864
ENSG00000099725 ENSG00000067646 ENSG00000099715 ENSG00000176728 ENSG00000154620
ENSG00000260197 ENSG00000206159 ENSG00000241859 ENSG00000278847 ENSG00000273906
ENSG00000227289 ENSG00000228764 ENSG00000232226 ENSG00000229236 ENSG00000169953
ENSG00000225117 ENSG00000270641 ENSG00000092377 ENSG00000267793 ENSG00000232348
ENSG00000224060 ENSG00000229163 ENSG00000215580 ENSG00000002586 ENSG00000231535
ENSG00000227494 ENSG00000133048 ENSG00000258484 ENSG00000184895 ENSG00000233070
ENSG00000228787 ENSG00000005889 ENSG00000214717 ENSG00000135245 ENSG00000064886
ENSG00000215301 ENSG00000229238 ENSG00000147050 ENSG00000198535 ENSG00000130600
ENSG00000217896 ENSG00000126012 ENSG00000261600 ENSG00000169248'''
```

#### 格式转化

基因名称是ENS格式，但HPO需要Entrez格式。

把基因转换为正确的格式

``` python
from hpoea.utils.idconvert import EntrezEnsemblConvert
cvt = EntrezEnsemblConvert() #报错，缺了个包，在终端中： pip3 install biothings_client解决
test_entrez = cvt.ensembl2entrez(test_genes)
print(len(test_entrez))#36
#转换格式后缺失了一些基因（原本有54个）
```

#### 富集分析

这里使用GSEA类中的一部分函数功能完成富集分析

![image-20210422173529632](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422173529632.png)

在导入“类”GSEA时，出现以上报错，浏览enrich.py后发现是由于该函数自动下载的HPO数据格式出错导致，手动修改函数导入自己下载并修改为正确格式后的HPO文件，修改位置如图C:\Users\Yangxiaolong\AppData\Local\Programs\Python\Python37\Lib\site-packages\hpoea-0.0.0-py3.7.egg\hpoea

![image-20210422184836786](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422184836786.png)

完成以上错误处理后开始进行富集分析：

``` python
gsea = GSEA()
gsea.enrich(test_entrez) #富集分析
gsea.multiple_test_corretion(method='fdr_bh')
print(gsea.enrichment_table.head(1)) #查看结果的第一行
```

结果第一行如下

![image-20210422191743734](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422191743734.png)

HPOID：即对应概念基因集在HPO数据库的ID

HPONAME：简单的注释（注释的题目），在官网搜索对应的ID可以看到更详细的注释内容

relatedgene：当前HPO基因集下与输入的基因集有关联的基因

pvalue：反映了输入基因集与对应HPO基因集及其概念的关联程度，关联越密切，p值越小

padj：p的矫正值

``` python
#gsea.enrichment_table.shape[0]
gsea.filter(by='padj', threshold=0.05)#筛选padj小于0.05
#print(type(gsea.enrichment_table))
t=gsea.enrichment_table
t.to_csv("enrichiment_table.csv")#保存分析结果表格为enrichiment_table.csv
```

结果中的前五行如图所示（共有54行）：

![image-20210422195338531](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422195338531.png)

由于padj<=0.05的结果太多，为了方便展示，筛选padj<0.01的结果

``` python
gsea.filter(by='padj', threshold=0.01)#筛选padj小于0.01
t.to_csv("enrichiment_table0.01.csv")#保存分析结果表格为enrichiment_table0.01.csv
```



得到如下图18行结果：

![image-20210422195643166](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422195643166.png)



#### 结果可视化

绘图函数dotplot.py：

``` python
import pandas as pd
import matplotlib.pyplot as plt
from bokeh.plotting import show
from bokeh.io import output_notebook
output_notebook()
from hpoea.enrich import GSEA
from hpoea.plot import LineagePlot, dot_plot
def dotp(path):
    f = open(path, encoding='utf-8')
    data = pd.read_csv(f)
    p = dot_plot(data, size=20, x='pvalue')
    return p
```

对应的主函数main.py：

``` python
import dotplot
from PIL import Image
import matplotlib.pyplot as plt
from bokeh.plotting import show
p=dotplot.dotp("D:/linuxvs/enrichiment_table0.01.csv")
#p.savefig("D:/linuxvs/dot.png")
show(p)
#print(type(p))
```

结果如下：





![textplot](D:\linuxvs\textplot.png)



绘制网络图：

绘图函数mapplot.py

``` python
import pandas as pd
import matplotlib.pyplot as plt
from bokeh.plotting import show
from bokeh.io import output_notebook
output_notebook()
from hpoea.enrich import GSEA
from hpoea.plot import LineagePlot, dot_plot
def mapp(path):
    f = open(path, encoding='utf-8')
    data = pd.read_csv(f)
    terms = list(data.HPO_term_ID)
    lin = LineagePlot()
    fig, ax = plt.subplots(figsize=(20, 10))
    lin.plot(terms, ax=ax)
```

对应主函数部分：

``` python
mapplot.mapp("D:/linuxvs/enrichiment_table0.01.csv")
```

运行时以上函数时报错，经过排查，是由于撞墙而不能下载绘制网络图所需文件导致

> 手动下载：https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.obo
>
> 保存为hpo.obo

修改plot.py中的如图部分，导入自己下载的文件

![image-20210422210525616](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422210525616.png)



再执行提示缺少部分依赖包：

``` powershell
pip install graphviz
pip install pygraphviz
```

解决以上问题再运行主函数，获得图片如下：

![testmap](D:\linuxvs\testmap.png)

> 代码到此为之运行完毕。
>
> 进一步分析贼需要在官网查询更加详细的注释并进行推理分析，本笔记仅作举例：
>
> ![image-20210422213321979](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422213321979.png)
>
> HP0001450注释信息：一种与Y染色体编码的基因相关的遗传模式。
>
> 由此可以推理输入的基因集与Y染色体可能有关，而实际上我们输入的基因集就是下丘脑组织中与性别相关的基因，因此分析较为合理。

### R包HPOSim的使用测试（可能是一种更深度，更有价值的HPO分析）：

> HPOSim is an R package for analyzing phenotypic similarity for genes and diseases based on HPO data. Seven commonly used semantic similarity measures are implemented in HPOSim. Enrichment analysis of gene sets and disease sets are also implemented, including hypergeometric enrichment analysis and network ontology analysis (NOA).
>
> HPOSim是一个基于HPO数据分析基因和疾病表型相似性的R包。在HPOSim中实现了七种常用的语义相似性度量。实现了基因集和疾病集的富集分析，包括超几何富集分析和网络本体分析。
>
> HPOSim consists of two parts: (i) the similarity measures between phenotypes (HPO terms), between human genes (Entrez IDs) and between diseases (OMIM IDs), and (ii) HPO-based enrichment analysis (NOA and the hypergeometric method) for gene set and disease set.
>
> HPOSim包括两个部分:(i)表型之间(HPO)、人类基因之间(Entrez id)和疾病之间(OMIM id)的相似性度量，以及(ii)基于HPO的基因集和疾病集的富集分析(NOA和超几何方法)。
>
> 结果示例：
>
> ![image-20210422214744071](C:\Users\Yangxiaolong\AppData\Roaming\Typora\typora-user-images\image-20210422214744071.png)

``` R
install.packages("D:/linuxvs/ziranyuyan/xiangmu/ceshi/HPO.db_1.9.tar.gz", repos = NULL, type = "source")

```

我没有找到这个包的代码示例，只能依靠论文和help（）文件获取帮助，进度缓慢，不再赘述。