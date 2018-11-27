# -*- coding: utf-8 -*-
import numpy as np
import math
from matplotlib import mlab
import matplotlib.pyplot as plt
from array import array

n = np.loadtxt("data/t7.txt", delimiter='\n', dtype=np.float)
n1 = np.loadtxt("data/t17.txt", delimiter='\n', dtype=np.float)

#m = np.loadtxt("data/x.txt", delimiter='\n', dtype=np.float)
#m1 = np.loadtxt("data/y.txt", delimiter='\n', dtype=np.float)

z = np.loadtxt("data/x17.txt", delimiter='\n', dtype=np.float)
z1 = np.loadtxt("data/y17.txt", delimiter='\n', dtype=np.float)

#l = np.loadtxt("data/xlong.txt", delimiter='\n', dtype=np.float)
#l1 = np.loadtxt("data/ylong.txt", delimiter='\n', dtype=np.float)

plt.plot(n, n1, color='blue')
#plt.plot(m, m1, color='green')
plt.plot(z, z1, color='orange')
#plt.plot(l, l1, color='m')
plt.title("Неявная схема")
plt.show()