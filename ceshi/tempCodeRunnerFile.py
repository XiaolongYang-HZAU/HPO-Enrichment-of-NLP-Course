import dotplot
from PIL import Image
import matplotlib.pyplot as plt
from bokeh.plotting import show
p=dotplot.dotp("D:/linuxvs/enrichiment_table0.01.csv")
#p.savefig("D:/linuxvs/dot.png")
plt.savefig(p,"D:/linuxvs/dot.png")