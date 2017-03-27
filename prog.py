import numpy as np
import matplotlib.pyplot as plt
import copy
from Pgraphs import Pgraphs


control = Pgraphs.Barabasi(1000)

u = copy.deepcopy(control)
d = copy.deepcopy(control)
c = copy.deepcopy(control)
m = copy.deepcopy(control)

power = 20
average = 1
floor = 0.5
dr = 0.01 # 1%
r = 1.
n = int((1.-floor)/dr)

pear_data = []
reci_data = []

upear_var = [[] for i in range(n)]
dpear_var = [[] for i in range(n)]
cpear_var = [[] for i in range(n)]
mpear_var = [[] for i in range(n)]

upear = np.float32([0]*n)
ureci = np.float32([0]*n)

dpear = np.float32([0]*n)
dreci = np.float32([0]*n)

cpear = np.float32([0]*n)
creci = np.float32([0]*n)

mpear = np.float32([0]*n)
mreci = np.float32([0]*n)

for j in range(average):
	pear_data = []
	reci_data = []
	for i in range(n):
		#u = copy.deepcopy(control)
		reciprocity = r-i*dr
		u.decrease_uniform_reci(reciprocity)
		degree, activation = u.matricial_walk(power , returnData=True)
		pear = np.corrcoef(degree, activation)
		pear = pear[1][0]
		pear_data.append(pear)
		reci_data.append(reciprocity)
		upear_var[i].append(pear)

	upear = upear + pear_data
	ureci = ureci + reci_data

upear = upear/average
ureci = ureci/average
upear_var = [np.array(i) for i in upear_var]
upear_var = [i.var() for i in upear_var]

u = plt.scatter(ureci, upear, color='black')
plt.plot(ureci, upear, color='black')

for j in range(average):
	pear_data = []
	reci_data = []
	for i in range(n):
		#d = copy.deepcopy(control)
		reciprocity = r-i*dr
		d.decrease_degree_reci(reciprocity)
		degree, activation = d.matricial_walk(power , returnData=True)
		pear = np.corrcoef(degree, activation)
		pear = pear[1][0]
		pear_data.append(pear)
		reci_data.append(reciprocity)
		dpear_var[i].append(pear)

	dpear = dpear + pear_data
	dreci = dreci + reci_data

dpear = dpear/average
dreci = dreci/average
dpear_var = [np.array(i) for i in dpear_var]
dpear_var = [i.var() for i in dpear_var]

d = plt.scatter(dreci, dpear,  color='red')
plt.plot(dreci, dpear, color='red')

for j in range(average):
	pear_data = []
	reci_data = []
	for i in range(n):
		#c = copy.deepcopy(control)
		reciprocity = r-i*dr
		c.decrease_clus_reci(reciprocity)
		degree, activation = c.matricial_walk(power , returnData=True)
		pear = np.corrcoef(degree, activation)
		pear = pear[1][0]
		pear_data.append(pear)
		reci_data.append(reciprocity)
		cpear_var[i].append(pear)

	cpear = cpear + pear_data
	creci = creci + reci_data

cpear = cpear/average
creci = creci/average
cpear_var = [np.array(i) for i in cpear_var]
cpear_var = [i.var() for i in cpear_var]

c = plt.scatter(creci, cpear, color='blue')
plt.plot(creci, cpear, color='blue')


for j in range(average):
	pear_data = []
	reci_data = []
	for i in range(n):
		#m = copy.deepcopy(control)
		reciprocity = r-i*dr
		m.decrease_MI_reci(reciprocity)
		degree, activation = m.matricial_walk(power , returnData=True)
		pear = np.corrcoef(degree, activation)
		pear = pear[1][0]
		pear_data.append(pear)
		reci_data.append(reciprocity)
		mpear_var[i].append(pear)

	mpear = mpear + pear_data
	mreci = mreci + reci_data

mpear = mpear/average
mreci = mreci/average
mpear_var = [np.array(i) for i in mpear_var]
mpear_var = [i.var() for i in mpear_var]

m = plt.scatter(mreci, mpear, color='green')
plt.plot(mreci, mpear, color='green')

a = [min(upear), min(cpear), min(dpear), min(mpear)]
a = min(a)
b = [min(ureci), min(creci), min(dreci), min(mreci)]
b = min(b) 

plt.grid(True)
plt.xlim(b*0.9, 1.1)
plt.ylim(a*0.9, 1.1)
plt.title('%s.png'%control.model, fontsize=16)
plt.xlabel('Reciprocidade', fontsize=16)
plt.ylabel('Coef Pearson', fontsize=16)

plt.legend((c,d, u, m), ('clustering', 'degree', 'uniforme', 'Matchin Index'), loc=2, scatterpoints=1)
plt.show()


