#生信1802杨晓龙 原创
#2021/5/12
#主函数
import dotplot
from PIL import Image
import matplotlib.pyplot as plt
from bokeh.plotting import show
import mapplot
p=dotplot.dotp("D:/linuxvs/ziranyuyan/xiangmu/ceshi/ceshi/head10-19.csv")
#p.savefig("D:/linuxvs/dot.png")
show(p)
#print(type(p))
mapplot.mapp("D:/linuxvs/ziranyuyan/xiangmu/ceshi/ceshi/enrichment_head15.csv")