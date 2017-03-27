import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
import matplotlib.axes as ax
import random
import networkx as nx
import copy

from decimal import Decimal

class Pgraphs:
	def __init__(self, model, graph, n, dire):
		self.model = model
		self.graph = graph
		self.n = n
		self.dire = dire
		self.edcome_go = set()
		self.edcome = set()
		self.probabily_matrix = None

	@staticmethod
	def Read_GML_Graph(filename):
		graph = ig.Graph.Read_GML(filename)
		n = graph.vcount()
		filenamec = []
		for i in filename:	
			if i == '.':
				break
			else:
				filenamec.append(i)
		filename = ''.join(filenamec)
		if graph.is_directed() == True:
			clus = graph.clusters()
			graph = clus.giant()
			directed = False
		else:
			redl = []
			edl	= graph.get_edgelist()
			for i in edl:
				redl.append(tuple(reversed(i)))
			graph = ig.Graph(directed=True)
			graph.add_vertices(n)
			graph.add_edges(edl + redl)
			clus = graph.clusters()
			graph = clus.giant()	
			directed = True
		n = graph.vcount()
		return Pgraphs(filename, graph, n, directed)
			
		
	@staticmethod
	def Barabasi(n=1000, m=5, directed=False):
		if directed == False:
			g = ig.Graph.Barabasi(n, m)
			redl = []
			edl	= g.get_edgelist()
			for i in edl:
				redl.append(tuple(reversed(i)))
			graph = ig.Graph(directed=True)
			graph.add_vertices(n)
			graph.add_edges(edl + redl)
			clus = graph.clusters()
			graph = clus.giant()	
		elif directed == True:
			g = ig.Graph.Barabasi(n, m, directed=True)
			edge_list = g.get_edgelist()
			edge = random.choice(edge_list)
			#down here i incrase the reciprocity, because the barabasi model ins really low reciprocal
			while g.reciprocity() < 0.3:
					edge = random.choice(edge_list)
					if edge in edge_list and tuple(reversed(edge)) in edge_list:
						edge_list.pop(edge_list.index(edge))	
					elif edge in edge_list or tuple(reversed(edge)) not in edge_list:
						new_edge = tuple(reversed(edge))
						edge_list.append(new_edge)
						g.add_edge(edge[1], edge[0])
			clus = g.clusters(mode="STRONG")
			graph = clus.giant()
		model = 'Barabasi-Albert'
		n = graph.vcount()
		return Pgraphs(model , graph, n, directed)
	
	@staticmethod
	def Waxman(n=1000, a=20, bbox=1, directed=False):
		model = 'Waxman'
		matrix = np.zeros((n, n))
		layout = [((np.random.uniform(0, np.sqrt(bbox))), (np.random.uniform(0, np.sqrt(bbox)))) for i in np.arange(n)]
		distance_matrix = np.zeros((n,n))
		for i in range(n):
			for j in range(n):
				Dij = np.sqrt((layout[i][0] - layout[j][0]) ** 2 + (layout[i][1] - layout[j][1]) ** 2)
				distance_matrix[i][j] = Dij

			# above we raffle the positions
		for i in range(n):
			for j in range(n):
				if i == j:
					continue
				d = distance_matrix[i][j]
				p = np.exp(-(a * d))
				if np.random.sample() < p:
					matrix[i][j] = 1
		if directed == False:
			g = ig.Graph.Adjacency(list(matrix), mode='UPPER')
			redl = []
			edl	= g.get_edgelist()
			for i in edl:
				redl.append(tuple(reversed(i)))
			graph = ig.Graph(directed=True)
			graph.add_vertices(n)
			graph.add_edges(edl + redl)
			clus = graph.clusters()
			graph = clus.giant()	
		elif directed == True:
			g = ig.Graph.Adjacency(list(matrix))
			clus = g.clusters(mode="STRONG")
			graph = clus.giant()
		else:
			raise ValueError('directed need to bem a boolean')
		n = graph.vcount()
			
		graph.vs['coordinates'] = layout
		g['a'] = a
		return Pgraphs(model , graph, n, directed)
	
	@staticmethod
	def Erdos_Renyi(n=1000, p=0.02, directed=False):
			model = 'Erdos-Renyi'
			if directed == False:
				g = ig.Graph.Erdos_Renyi(n, p)
				redl = []
				edl	= g.get_edgelist()
				for i in edl:
					redl.append(tuple(reversed(i)))
				graph = ig.Graph(directed=True)
				graph.add_vertices(n)
				graph.add_edges(edl + redl)
				clus = graph.clusters()
				graph = clus.giant()		
			elif directed == True:
				g = ig.Graph.Erdos_Renyi(n, p, directed=True)
				clus = g.clusters(mode="STRONG")
				graph = clus.giant()
			n = graph.vcount()
			return Pgraphs(model , graph, n, directed)
		

	@staticmethod
	def Watts_Strogatz(n=1000, k=10, p=0.03, directed=False):
		model = 'Watts-Strogatz'
		if directed == False:
			g = ig.Graph.Watts_Strogatz(1, n, k, p)	
			redl = []
			edl = g.get_edgelist()
			for i in edl:
				redl.append(tuple(reversed(i)))
			g = ig.Graph(directed=True) 
			g.add_vertices(n)
			g.add_edges(edl + redl)
			clus = g.clusters()
			graph = clus.giant()
		elif directed == True:
			g = ig.Graph.Watts_Strogatz(1, n, k, 0)	
			redl = []
			edl = g.get_edgelist()
			for i in edl:
				redl.append(tuple(reversed(i)))
			g = ig.Graph(directed=True) 
			g.add_vertices(n)
			g.add_edges(edl + redl)
			adj_list = g.get_adjlist()
			for i in range(len(adj_list)):
				for j in range(len(adj_list[i])):
					if np.random.uniform() < p:
						g.delete_edges((i,j))
						newnode = np.random.randint(0,n)
						while newnode != i:
							newnode = np.random.randint(0,n)
						g.add_edges([(i, newnode)])
			clus = g.clusters()
			graph = clus.giant()
		n = graph.vcount()
		return Pgraphs(model , graph, n, directed)

	def adj_listf(self):
		self.adj_list = self.graph.get_adjlist()
		return self.adj_list

	def plot(self):
		ig.plot(self.graph)
		
	def average_short_path(self):
		a = self.graph.shortest_paths_dijkstra() 
		s = 0 
		control = 0
		for i in a:
			for j in i:
				if j > 10**100:
						continue
				else:	
					s = j + s
					control += 1
		s = s/control
		return s

	def shortest_path(self):
		a = self.graph.shortest_paths_dijkstra() 
		s_min = 1000 
		for i in a:
			for j in i:
				if j > 10**100:
						continue
				else:	
					if j < s_min and j != 0:
						s_min = j	
		return s_min

	def reciprocityf(self):
		return self.graph.reciprocity()

	def clustering(self):
		return self.graph.transitivity_avglocal_undirected()

	def edge_listf(self):
		self.edge_list = self.graph.get_edgelist()
		return self.edge_list
		
	def average_degreef(self):
		self.edge_listf()
		self.average_degree = len(self.edge_list)/self.n
		return self.average_degree

	def att_StrongComponent(self):
		clus = self.graph.clusters()
		g = clus.giant()
		self.graph = g
		self.n = len(self.graph.get_adjlist())	
		self.disjoint_edges()

	def disjoint_edges(self):
		self.edge_listf()
		aux = set(self.edge_list)
		l = len(aux)
		for i in self.edge_list:
			if tuple(reversed(i)) in aux:
				self.edcome_go.add(i)
			else:
				self.edcome.add(i)

	def get_prob_distance_matrix(self):
		if self.probabily_matrix == None:
			self.probabily_matrix = np.zeros((self.n, self.n))
			layout = self.graph.vs['coordinates']
			i = 0
			j = 0
			while i < self.n:
				while j < self.n:
					if i == j:
						j += 1
						continue
					Dij = np.sqrt((layout[i][0] - layout[j][0]) ** 2 + (layout[i][1] - layout[j][1]) ** 2)
					a = g['a']
					p = np.exp(-(a * Dij))
					self.probabily_matrix[i][j] = p
					j += 1
				i += 1
				j = i
			self.probabily_matrix = self.probabily_matrix + self.probabily_matrix.T
	
	def MatchingIndex_list(self):
		self.disjoint_edges()
		self.edge_listf()
		self.adj_listf()
		set_edges = set(self.edge_list)
		matilist = []
		for i in self.edge_list: 
			if i in self.edcome:
				matilist.append(0)
				continue
			else: 
				n = 0
				d = 0
				aux = set(self.adj_list[i[0]]).intersection(self.adj_list[i[1]])
				n = len(aux)
				d = len(self.adj_list[i[0]]) + len(self.adj_list[i[1]]) - 2
				mi = n/d
				matilist.append(mi)
		return matilist
	
	def ave_MatchingIndex(self):
		aux = self.MatchingIndex_list()
		return np.sum(aux)/len(aux)	
	
	### Dynamics
	def random_walk(self, steps=10**6):
		self.adj_listf()
		if len(self.adj_list) == 1:
			return 0
		self.steps = steps
		self.state_list = [0 for i in range(self.n)]
		self.start = int(np.random.uniform(0, self.n))
		self.position = self.start
		for i in range(self.steps):
			while True:
				if not self.adj_list[self.position]:
					self.position = random.randint(0,self.n)
				else:
					break
			self.position = random.choice(self.adj_list[self.position])
			self.state_list[self.position] += 1
		return self.state_list

	def transmission_matrix(self):
		self.adj_listf()
		self.tmatrix = np.zeros((self.n, self.n))
		for j in np.arange(self.n):
			for k in self.adj_list[j]:
				p = 1/len(self.adj_list[j])
				self.tmatrix[j][k] = p

	def matricial_walk(self, steps, plot=False, lg=None, returnData=False):
		self.adj_listf()
		self.steps = steps
		self.transmission_matrix()
		tn = np.linalg.matrix_power(self.tmatrix, steps)   
		self.p0 = self.tmatrix[0]
		self.pN = np.dot(tn.transpose() , self.p0)
		self.pN = self.pN/np.sum(self.pN)
		d = [len(self.adj_list[i]) for i in range(len(self.adj_list))]
		if plot == True:
			self.plot_AtivationDegree(list(d), list(self.pN), steps, lg)
		elif returnData== True:
			return list(d), list(self.pN)
		else: pass
	
	###plots
	def plot_AtivationDegree(self, degree, ativation, steps, lg):
		fig, ax = plt.subplots(1)
		pearson = np.corrcoef(degree , ativation)
		self.average_degreef()
		textstr = ('pearson coefficient = %f\nreciprocity = %f\nnodes = %i\nsteps on the network=%.1E\naverage degree = %f'%(pearson[1][0],
				self.reciprocityf(), self.n, Decimal(steps), self.average_degree))
		props = dict(boxstyle='round', facecolor='white', alpha=1)
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
		pl = plt.scatter(degree, ativation)
		plt.xlim((0, 1.1*max(degree)))
		plt.ylim((0,1.1*max(ativation)))
		plt.xlabel('degree')
		plt.ylabel('total activation')
		#plt.title('correlation - activation degree, %s'%lg)
		plt.savefig('%s - activation-degree %f,%s.jpeg'%(self.model, self.reciprocityf(), lg))
		#return pl
		plt.show()
		plt.cla()
		plt.clf()
		plt.close()

	def plot(self):
		ig.plot(self.graph)

	def ReciprocityPearson(self, reciprocity_min, reciprocity_max, type_decai, interval, power, subplot=False, returnData=False, plot=False):
		if self.reciprocityf() > reciprocity_max:
			self.decrease_uniform_reci(reciprocity_max)
		pearson_data = []
		reciprocity_data = []
		if type_decai=='clustering':
			reciprocity = self.reciprocityf()
			degree, ativation = self.matricial_walk( power, '', returnData=True)
			pearson = np.corrcoef(degree , ativation)
			pearson_data.append(pearson[1][0])
			reciprocity_data.append(reciprocity)
			for i in range(int((reciprocity_max-reciprocity_min)/interval)):
				self.decrease_clus_reci(reciprocity - interval)
				reciprocity = self.reciprocityf()
				degree, ativation = self.matricial_walk( power, '', returnData=True)
				pearson = np.corrcoef(degree , ativation)
				pearson_data.append(pearson[1][0])
				reciprocity_data.append(reciprocity)

		elif type_decai=='degree':
			reciprocity = self.reciprocityf()
			degree, ativation = self.matricial_walk( power, '', returnData=True)
			pearson = np.corrcoef(degree , ativation)
			pearson_data.append(pearson[1][0])
			reciprocity_data.append(reciprocity)
			for i in range(int((reciprocity_max-reciprocity_min)/interval)):
				self.decrease_degree_reci(reciprocity - interval)
				reciprocity = self.reciprocityf()
				degree, ativation = self.matricial_walk( power, '', returnData=True)
				pearson = np.corrcoef(degree , ativation)
				pearson_data.append(pearson[1][0])
				reciprocity_data.append(reciprocity)

		elif type_decai=='uniform':
			reciprocity = self.reciprocityf()
			degree, ativation = self.matricial_walk( power, '', returnData=True)
			pearson = np.corrcoef(degree , ativation)
			pearson_data.append(pearson[1][0])
			reciprocity_data.append(reciprocity)
			for i in range(int((reciprocity_max-reciprocity_min)/interval)):
				self.decrease_uniform_reci(reciprocity - interval)
				reciprocity = self.reciprocityf()
				degree, ativation = self.matricial_walk( power, '', returnData=True)
				pearson = np.corrcoef(degree , ativation)
				pearson_data.append(pearson[1][0])
				reciprocity_data.append(reciprocity)

		elif type_decai == 'MI':
			reciprocity = self.reciprocityf()
			degree, ativation = self.matricial_walk( power, '', returnData=True)
			pearson = np.corrcoef(degree , ativation)
			pearson_data.append(pearson[1][0])
			reciprocity_data.append(reciprocity)
			for i in range(int((reciprocity_max-reciprocity_min)/interval)):
				self.decrease_MI_reci(reciprocity - interval)
				reciprocity = self.reciprocityf()
				degree, ativation = self.matricial_walk( power, '', returnData=True)
				pearson = np.corrcoef(degree , ativation)
				pearson_data.append(pearson[1][0])
				reciprocity_data.append(reciprocity)

		else: raise TypeError('is need to specify the prefence of the edges reorganization')
		if plot == True:
			plt.plot(pearson_data, reciprocity_data, label='%s'%type_decai)
			plt.xlim(0,1)
			plt.ylim(0,1)
			plt.title('%s - Pearson X reciprocity'%self.model)
			plt.show()	
		elif returnData == True:
			return pearson_data, reciprocity_data
		else:
			subplot.plot(pearson_data, reciprocity_data, label='%s'%type_decai)
		
	def interative_dynamics(self, steps, lg, returnData=False, plot=False):
		self.adj_listf()
		pn = self.random_walk(steps)
		d = [len(self.adj_list[i]) for i in range(len(self.adj_list))]
		if returnData==True:
			return d, pn	
		elif plot==True:
			self.plot_AtivationDegree(d, pn, steps, lg)
		else:
			pass

	@staticmethod
	def Plot_Reciprocity_Pearson(power, average_n, interval, dt, Pgraph_obj, grid=True, save=False, plot=True):
		g = int(abs(interval[1]-interval[0])/dt)
		upear = np.float64([0]*g)
		ureci = np.float64([0]*g) 
		cpear = np.float64([0]*g)
		creci = np.float64([0]*g)

		dpear = np.float64([0]*g)
		dreci = np.float64([0]*g)

		mpear = np.float64([0]*g)
		mreci = np.float64([0]*g)
		
		control = Pgraph_obj

		u = copy.deepcopy(control)
		c = copy.deepcopy(control)
		d = copy.deepcopy(control)
		m = copy.deepcopy(control)

		for i in range(average_n):
			u = copy.deepcopy(control)
			c = copy.deepcopy(control)
			d = copy.deepcopy(control)
			m = copy.deepcopy(control)

			p, r = u.ReciprocityPearson( interval[0], interval[1], 'uniform', dt, power=power, returnData=True)
			upear = upear + p
			ureci = ureci + r
			p, r = c.ReciprocityPearson( interval[0], interval[1], 'clustering', dt, power=power, returnData=True)
			cpear = cpear + p
			creci = creci + r
			p, r = d.ReciprocityPearson( interval[0], interval[1], 'degree', dt, power=power, returnData=True)
			dpear = dpear + p
			dreci = dreci + r
			p, r = m.ReciprocityPearson( interval[0], interval[1], 'MI', dt, power=power, returnData=True)
			mpear = mpear + p
			mreci = mreci + r

		ureci = ureci/average_n
		upear = upear/average_n

		creci = creci/average_n
		cpear = cpear/average_n

		dreci = dreci/average_n
		dpear = dpear/average_n

		mreci = mreci/average_n
		mpear = mpear/average_n

		a = [min(upear), min(cpear), min(dpear), min(mpear)]
		a = min(a)
		b = [min(ureci),min(creci),min(dreci), min(mreci)]
		b = min(b)

		u = plt.scatter(ureci, upear, color='black', s=20)
		plt.plot(ureci, upear, color='black')
		c = plt.scatter(creci, cpear, color='blue', s=20)
		plt.plot(creci, cpear, color='blue')
		d = plt.scatter(dreci, dpear, color='red', s=20)
		plt.plot(dreci, dpear, color='red')
		m = plt.scatter(mreci, mpear, color='green', s=20)
		plt.plot(mreci, mpear, color='green')

		plt.grid(grid)
		plt.legend((u,c,d,m), ('uniform', 'clustering cof.', 'degree', 'matchin In'), fontsize=14, loc=2)
		plt.ylabel('Pearson Cof', fontsize=16)
		plt.xlabel('Reciprocidade', fontsize=16)
		plt.title('%s'%control.model)
		plt.xlim(b*0.9, 1.1)
		plt.ylim(a*0.9, 1.1)

		if save==True:
			plt.savefig('%s.jpeg'%control.model)
		if plot==True:
			plt.show()

		print(control.n)
		print(control.average_degreef())
		print(control.model)

	
	### Reorganization
	def model_connect_probability(self, node):
		if self.model == 'Waxman':
			self.get_prob_distance_matrix()
			pi = self.get_prob_distance_matrix[node]
			return pi, pi/np.sum(pi)

		elif self.model == 'Barabasi':
			degree_list = self.adj_list
			degree_weight = degree_list/np.sum(degree_list)
			return degree_list, degree_weight	

		elif self.model == 'Erdos-Renyi':
			p = np.array([1 for i in range(self.n)])
			pw = p/self.n
			return p, pw
	
		elif self.model == 'Watts-Strogatz':
			pass		

	def decrease_uniform_reci(self, rec, keepModel=False): 
		#the network should be a sparse one to increase the algorithm speed
		self.dire = True
		self.edge_listf()
		self.adj_listf()
		self.disjoint_edges()
		aux_adj_list = copy.copy(self.adj_list)
		aux_edge_list = copy.copy(self.edge_list)
		set_edges = set(aux_edge_list)
		to_add_edge = []
		to_del_edge = []
		aux_range = range(self.n)
		reciprocity = self.reciprocityf()
		while reciprocity > rec:
			node = np.random.choice(aux_range)
			if len(self.adj_list[node]) == 0: 
				continue
			target = np.random.choice(self.adj_list[node])
			while (node, target) not in set_edges:
				print('oila')
				target = np.random.choice(self.adj_list[node])

			if (node, target) in self.edcome:
				self.edcome.remove((node, target))
			elif (node, target) in self.edcome_go:
				self.edcome_go.remove((node, target))
				self.edcome_go.remove((target, node))
				self.edcome.add((target, node))

			to_del_edge.append((node, target))
			aux_adj_list[node].pop(aux_adj_list[node].index(target))

			new_target = random.randrange(len(self.adj_list))
			while new_target == node:
				new_target = random.randrange(len(self.adj_list))
				print('oika')
			to_add_edge.append((node, new_target))
			if (new_target, node) in self.edcome:
				self.edcome.remove((new_target, node))
				self.edcome_go.add((new_target, node))
				self.edcome_go.add((node, new_target))
			else:
				self.edcome.add((node, new_target))
			set_edges.add((node, new_target))
				
			reciprocity = len(self.edcome_go)/(len(self.edcome_go) + len(self.edcome))

		self.graph.delete_edges(to_del_edge)
		self.graph.add_edges(to_add_edge)
		self.edge_listf()
		self.adj_listf()
		self.att_StrongComponent()

	
	def decrease_clus_reci(self, rec):
		self.dire = True
		aux_range = None
		self.edge_listf()
		self.adj_listf()
		self.disjoint_edges()
		aux_adj_list = copy.copy(self.adj_list)
		aux_edge_list = copy.copy(self.edge_list)
		set_edges = set(aux_edge_list)
		to_add_edge = []
		to_del_edge = []
		aux_range = range(self.n)
		reciprocity = self.reciprocityf()
		while reciprocity > rec:
			clus_cof_list = np.float64(self.graph.transitivity_local_undirected(mode='zero'))
			clus_cof_weight = clus_cof_list/np.sum(clus_cof_list)
			
			node = np.random.choice(aux_range, p=clus_cof_weight, replace=False)
			while len(self.adj_list[node]) == 0: 
				node = np.random.choice(aux_range, p=clus_cof_weight, replace=False)

			plist = [clus_cof_list[i] for i in self.adj_list[node]]
			if np.sum(plist) == 0:
				target = np.random.choice(self.adj_list[node], replace=False)
				while (node, target) not in set_edges:
					target = np.random.choice(self.adj_list[node], replace=False)
			else:
				plist = plist/np.sum(plist)
				target = np.random.choice(self.adj_list[node], p=plist, replace=False)
				while (node, target) not in set_edges:
					target = np.random.choice(self.adj_list[node], p=plist, replace=False)

			if (node, target) in self.edcome:
				self.edcome.remove((node, target))
			elif (node, target) in self.edcome_go:
				self.edcome_go.remove((node, target))
				self.edcome_go.remove((target, node))
				self.edcome.add((target, node))

			to_del_edge.append((node, target))
			aux_adj_list[node].pop(aux_adj_list[node].index(target))

			new_target = random.randrange(len(self.adj_list))
			while new_target == node:
				new_target = random.randrange(len(self.adj_list))
			to_add_edge.append((node, new_target))
			if (new_target, node) in self.edcome:
				self.edcome.remove((new_target, node))
				self.edcome_go.add((new_target, node))
				self.edcome_go.add((node, new_target))
			else:
				self.edcome.add((node, new_target))
			set_edges.add((node, new_target))
				
			reciprocity = len(self.edcome_go)/(len(self.edcome_go) + len(self.edcome))
			

		self.graph.delete_edges(to_del_edge)
		self.graph.add_edges(to_add_edge)
		self.edge_listf()
		self.adj_listf()
		self.att_StrongComponent()

	def decrease_degree_reci(self, rec):
		self.dire = True
		aux_range = None
		self.adj_listf()
		self.disjoint_edges()

		aux_adj_list = copy.copy(self.adj_list)

		set_edges = set(self.edge_listf())
		to_add_edge = []
		to_del_edge = []
		aux_range = range(self.n)

		reciprocity = self.reciprocityf()
		degree_list = [len(i) for i in self.adj_list]
		sum_degree = np.sum(degree_list)
		degree_weight = degree_list/sum_degree
		len_zero_node = []

		while reciprocity > rec:
			node = np.random.choice(aux_range, p=degree_weight, replace=False)
			plist = [len(self.adj_list[i]) for i in self.adj_list[node]]
			while np.sum(plist) == 0: 
				node = np.random.choice(aux_range, p=degree_weight, replace=False)
				degree_weight[node] = 0
				if len(aux_adj_list[node]) == 0: 
					continue
				plist = [len(aux_adj_list[i]) for i in aux_adj_list[node]]
			plist = plist/np.sum(plist)
			target = np.random.choice(aux_adj_list[node], p=plist, replace=False)
			while (node, target) not in set_edges:
				target = np.random.choice(aux_adj_list[node], p=plist, replace=False)

			if (node, target) in self.edcome:
				self.edcome.remove((node, target))
			elif (node, target) in self.edcome_go:
				self.edcome_go.remove((node, target))
				self.edcome_go.remove((target, node))
				self.edcome.add((target, node))
					
			
			to_del_edge.append((node, target))
			set_edges.remove((node, target))
			degree_weight[target] = len(aux_adj_list[target])/sum_degree
			aux_adj_list[node].pop(aux_adj_list[node].index(target))

			new_target = random.randrange(len(self.adj_list))
			while new_target == node:
				new_target = random.randrange(len(self.adj_list))
			if (new_target, node) in self.edcome:
				self.edcome.remove((new_target, node))
				self.edcome_go.add((new_target, node))
				self.edcome_go.add((node, new_target))
			else:
				self.edcome.add((node, new_target))
			aux_adj_list[node].append(new_target)
			degree_weight[new_target] = len(aux_adj_list[new_target])/sum_degree
			set_edges.add((node, new_target))
			to_add_edge.append((node, new_target))
				
			reciprocity = len(self.edcome_go)/(len(self.edcome_go) + len(self.edcome))
		self.graph.add_edges(list(to_add_edge))
		self.graph.delete_edges(list(to_del_edge))
		self.edge_listf()
		self.adj_listf()
		self.att_StrongComponent()

	def decrease_MI_reci(self, reci):
		matilist = self.MatchingIndex_list()
		matilist = matilist/np.sum(matilist)
		aux_range = len(self.edge_listf())
		while self.reciprocityf() > reci:
			matilist = self.MatchingIndex_list()
			if np.sum(matilist) == 0:
				break
			matilist = matilist/np.sum(matilist)
			edge = np.random.choice(list(range(len(self.edge_listf()))), p=matilist, replace=False)	
			edge = self.edge_list[edge]
			self.graph.delete_edges([edge])
			n1 = np.random.randint(0, self.n)
			n2 = np.random.randint(0, self.n)
			while n2 == n1:
				n2 = np.random.randint(0, self.n)
			if (n1,n2) in self.edcome:
				self.edcome.remove((n1,n2))	
				self.edcome_go.add((n1,n2))	
				self.edcome_go.add((n2,n1))	
			self.graph.add_edges([(n1, n2)])
		self.att_StrongComponent()
			

	def chose_weighted_element(self, flist):
		r = random.uniform(0, sum(flist))
		s = 0.0
		for i in range(len(flist)):
			s += flist[i]
			if r < s:
				return	i	





