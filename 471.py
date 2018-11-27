# -*- coding: utf-8 -*-
import numpy as np
import math
from matplotlib import mlab
import matplotlib.pyplot as plt
from array import array

n = np.loadtxt("data/t7.txt", delimiter='\n', dtype=np.float)
n1 = np.loadtxt("data/t17.txt", delimiter='\n', dtype=np.float)


z = np.loadtxt("data/x17.txt", delimiter='\n', dtype=np.float)
z1 = np.loadtxt("data/y17.txt", delimiter='\n', dtype=np.float)

#построение графика точного и аналитического решения задачи в заданный момент времени
plt.plot(n, n1, color='blue')
plt.plot(z, z1, color='orange')
plt.title("Неявная схема")
plt.show()
