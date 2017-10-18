import numpy as np

n = np.load("intz.npz")
nd = n['lt']
o = np.load("gold_standard_00050.npz")
od = o['lt']

coco = np.corrcoef(nd,od)
print("correlation coefficient:\n {}".format(coco))