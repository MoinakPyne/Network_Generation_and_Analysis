# VENKATRAMAN RENGANATHAN - vxr131730
# MOINAK PYNE - mxp132030
# MECH 6317/EECS 6302 Dynamics of Complex Networks & Systems
# PROJECT, Due 28th April 2017
######################################################################################################

import zen
import numpy
import random
import urllib
import time
import unicodedata
from io import open
import matplotlib.pyplot as plt
plt.ioff()
from numpy import *
import numpy.linalg as la
from numpy.linalg import eig,norm
import sys
sys.path.append('../zend3js/')
import d3js
from time import sleep
import colorsys
import scipy.sparse as sparse
from scipy.sparse import csgraph
from scipy.integrate import odeint
from scipy import sparse
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import norm

######################################################################################
############################ FUNCTION DEFINITIONS ####################################
######################################################################################


def propagate(G,A,x,steps,slp=0.5,keep_highlights=False,update_at_end=False):
	# interactive = d3.interactive
	# d3.set_interactive(False)
	A = A.T  # adjacency matrix of the network G
	# d3.highlight_nodes_(list(where(x>0)[0]))
	# d3.update()
	# sleep(slp)
	# cum_highlighted = sign(x)
	for i in range(steps): # the brains
		x = sign(numpy.dot(A,x)) # the brains
		# cum_highlighted = sign(cum_highlighted+x)
		# if not update_at_end:
		# 	if not keep_highlights:
		# 		d3.clear_highlights()
		# 	d3.highlight_nodes_(list(where(x>0)[0]))
		# 	d3.update()
		# 	sleep(slp)
	# if update_at_end:
	# 	if not keep_highlights:
	# 		d3.clear_highlights()
	# 		d3.highlight_nodes_(list(where(x>0)[0]))
	# 	else:
	# 		d3.highlight_nodes_(list(where(cum_highlighted>0)[0]))
	# 	d3.update()
	# d3.set_interactive(interactive)
	if keep_highlights:
		return cum_highlighted
	else:
		return x

# prints the top five (num) nodes according to the centrality vector v
# v takes the form: v[nidx] is the centrality of node with index nidx
def print_top(G,v, num=5):
	idx_list = [(i,v[i]) for i in range(len(v))]
	idx_list = sorted(idx_list, key = lambda x: x[1], reverse=True)
	for i in range(min(num,len(idx_list))):
		nidx, score = idx_list[i]
		print '  %i. %s (%1.4f)' % (i+1,G.node_object(nidx),score)
		#print '  %i. %s' % (i+1,G.node_object(idx))

# returns the index of the maximum of the array
# if two or more indices have the same max value, the first index is returned
def index_of_max(v):
	return numpy.where(v == max(v))[0]

########### POWER LAW ALPHA CALCULATION ######################################
def calculate_alpha(G, kmin):
	summation = 0
	alpha = 0
	G_nodes = G.num_nodes
	# print '\n Number of Nodes in Graph = %i' % G_nodes
	deg = numpy.zeros((G_nodes))
	for i in range(0, G_nodes):
		deg[i] = G.degree_(i)
	deg_seq = numpy.sort(deg)
	N = sum(deg >= kmin)
	denominator = kmin - 0.5	
	for i in range(0, G_nodes):
		if(deg_seq[i] >= kmin):
			numerator = deg_seq[i]
			summation = summation + numpy.log(numerator/denominator)
	alpha = 1 + (N / summation)
	return alpha

## Plots the degree distribution and calculates the power law coefficent
def calc_powerlaw(G,kmin):
	ddist = zen.degree.ddist(G,normalize=False)
	cdist = zen.degree.cddist(G,inverse=True)
	k = numpy.arange(len(ddist))

	# Degree Distribution - Bar Plot
	plt.figure(figsize=(8,12))
	plt.bar(k,ddist, width=0.8, bottom=0, color='b')
	plt.xlabel('k')
	plt.ylabel('p(k)')
	plt.title('Degree Distribution Bar Plot')

	# Cumulative Degree Distribution - Log-Log Plot
	plt.figure(figsize=(8,12))
	plt.loglog(k, cdist)
	plt.xlabel('log(k)')
	plt.ylabel('log(cdf)')
	plt.title('Cumulative Degree Distribution Log-Log Plot')
	
	# Need to visualize and determine kmin and suppy it to below function.
	alpha = calculate_alpha(G, kmin)
	print 'Kmin = %i, Alpha = %1.4f' % (kmin, alpha)
	
	plt.show()

################ MODULARITY CALCULATION ################################################################
def modularity(G,c):
	d = dict()
	for k,v in c.iteritems():
		for n in v:
			d[n] = k
	Q, Qmax = 0,1
	for u in G.nodes_iter():
		if (("Pasadena" in G.node_data(u)['zenData']) or ("Massachusetts" in G.node_data(u)['zenData']) or ("Berkeley" in G.node_data(u)['zenData']) or ("Stanford" in G.node_data(u)['zenData']) or ("Ann Arbor" in G.node_data(u)['zenData']) or ("Independent" in G.node_data(u)['zenData'])):
			for v in G.nodes_iter():
				if (("Pasadena" in G.node_data(v)['zenData']) or ("Massachusetts" in G.node_data(v)['zenData']) or ("Berkeley" in G.node_data(v)['zenData']) or ("Stanford" in G.node_data(v)['zenData']) or ("Ann Arbor" in G.node_data(v)['zenData']) or ("Independent" in G.node_data(v)['zenData'])):
					if d[u] == d[v]:
						Q += ( int(G.has_edge(v,u)) - G.degree(u)*G.degree(v)/float(G.num_edges) )/float(G.num_edges)
						Qmax -= ( G.degree(u)*G.degree(v)/float(G.num_edges) )/float(G.num_edges)
	return Q, Qmax

######################### FRIENDSHIP PARADOX CHECKING ########################################################
def check_friendship_paradox(G):
	G_nodes = G.num_nodes
	G_edges = G.num_edges
	k_sum = 0
	k_sqr_sum = 0
	for i in range(0, G_nodes):
		i_deg = G.degree_(i)
		k_sum = k_sum + i_deg
		k_sqr_sum = k_sqr_sum + i_deg * i_deg
	k = k_sum / G_nodes
	k_sqr = k_sqr_sum / G_nodes
	print '\n <k^2> = %i' %k_sqr
	print '\n <k> = %i' %k
	neighbor_average = float(k_sqr)/float(k)
	node_average = float(k)
	if(float(k_sqr)/float(k) > float(k)):
		print 'Neighbor of a node has more neighbors (%1.3f) than the node itself (%1.3f) - Friendship Paradox is observed' %(neighbor_average, node_average)
	else:
		print 'Friendship Paradox is not observed'

## HELPER FUNCTIONS =======================================
def RGBToHTMLColor(rgb_tuple):
	""" convert an (R, G, B) tuple to #RRGGBB """
	hexcolor = '#%02x%02x%02x' % rgb_tuple
	# that's it! '%02x' means zero-padded, 2-digit hex values
	return hexcolor

def HTMLColorToRGB(colorstring):
	""" convert #RRGGBB to an (R, G, B) tuple """
	colorstring = colorstring.strip()
	if colorstring[0] == '#': colorstring = colorstring[1:]
	if len(colorstring) != 6:
		raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
	r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
	r, g, b = [int(n, 16) for n in (r, g, b)]
	return (r, g, b)

def color_interp(color1,color2,v,m=0,M=1):
	c1 = array(HTMLColorToRGB(color1))
	c2 = array(HTMLColorToRGB(color2))
	if v > M:
		c = tuple(c2)
	elif v < m:
		c = tuple(c1)
	else:
		#c = tuple( c1 + (c2-c1)/(M-m)*(v-m) ) # linear interpolation of color
		c = tuple( c1 + (c2-c1)*(1 - exp(-2*(v-m)/(M-m))) ) # logistic interpolation of color
	return RGBToHTMLColor(c)

def color_by_value(d3,G,x,color1='#77BEF5',color2='#F57878'):
	d3.set_interactive(False)
	m = x.min()
	M = x.max()
	for i in G.nodes_iter_():
		d3.stylize_node_(i, d3js.node_style(fill=color_interp(color1,color2,x[i])))
	d3.update()
	d3.set_interactive(True)

# G: Graph
# p: uniform probability to activate across an edge
# x: initial active seed set (as a list/array)
def influence_cascade(G,p,x):
	G = G.copy()
	x = x.copy()
	activated_nodes = set([])
	for i,xi in enumerate(x):
		if xi > 0:
			activated_nodes.add(G.node_object(i))

	while len(activated_nodes) > 0:			
		newly_activated = set([])
		for u in activated_nodes:
			x[G.node_idx(u)] = 1
			for v in G.neighbors(u):
				if random.random() <= p:
					newly_activated.add(v)
				G.rm_edge(u,v)
		activated_nodes = newly_activated	
	return x

# G: Graph
# edge_prob: uniform probability to activate across an edge
# x: initial active seed set (as a list/array)
def independent_cascade(G, p, x):
	G = G.copy()
	x = x.copy()
	activated_nodes = set([])
	for i,xi in enumerate(x):
		if xi > 0:
			activated_nodes.add(G.node_object(i))

	while len(activated_nodes) > 0:			
		newly_activated = set([])
		for u in activated_nodes:
			x[G.node_idx(u)] = 1
			for v in G.neighbors(u):
				if random.random() <= p:
					newly_activated.add(v)		
				if (random.randint(0,2) == 1):				
					G.rm_edge(u,v)
		activated_nodes = newly_activated	
	return x
## ========================================================

def get_sparse_adjacency_matrix(G):
	nodelist = G.nodes()
	nlen = len(nodelist)
	undirected = not G.is_directed()
	index=dict(zip(nodelist,range(nlen)))
	M = sparse.lil_matrix((nlen,nlen))
	nodeset = set(nodelist)
	for u,v,w in G.edges_iter(weight=True):
		if (u in nodeset) and (v in nodeset):
			i,j = index[u],index[v]
			M[i,j] = M[i,j] + w
			if undirected:
				M[j,i] = M[i,j]
	return M

def spectral_modularity_maximization(G, A):
	nodelist = G.nodes()
	nlen = len(nodelist)
	undirected = not G.is_directed()
	index=dict(zip(nodelist,range(nlen)))
	nodeset = set(nodelist)
	B = sparse.lil_matrix((nlen,nlen))
	for u,v,w in G.edges_iter(weight=True):
		if (u in nodeset) and (v in nodeset):
			i,j = index[u],index[v]
			B[i,j] = A[i,j] - (G.degree_(i) * G.degree_(j))/(2*float(G.num_edges))
			if undirected:
				B[j,i] = B[i,j]
	return B

######################################################################################
############################ NETWORK ANALYSIS ####################################
######################################################################################

G = zen.Graph()
################## AFTER NETWORK IS BUILT, IT SHOULD BE SUPPLIED HERE FOR ANALYSIS################
G = zen.io.gml.read('author_nw.gml',weight_fxn = lambda x: x['weight'])
print 'Loaded the network'
#G = zen.edgelist.read('author_nw.edgelist',weighted=True)
#d3 = d3js.D3jsRenderer(G, interactive=False, autolaunch=False)

G_Nodes = G.num_nodes
G_Edges = G.num_edges

#Printing the number of nodes and edges in the network (G) given in Fig 6.1a
print '#Nodes: %i, #Edges: %i' % (G_Nodes,G_Edges)

c = 2*float(G_Edges)/float(G_Nodes)
print 'Average Degree of the network = %1.4f' %float(c)

#M is the adjacency matrix for the network G
A = get_sparse_adjacency_matrix(G)
print 'Got Sparse Adjacency Matrix successfully'
(A_row, A_col) = A.shape
print '\nA matrix - #rows: %i, #cols: %i' % (A_row,A_col)

#d = zen.diameter(G)
d = 36 # Obtained from running code from Dr. Ruths code
# print '\nDiameter of the network = %i' %d

# # Find Global Clustering Coefficient
global_clustering = zen.algorithms.clustering.gcc(G) 
print 'Global Clustering Coefficient = %1.4f' %global_clustering


# Degree Centrality
print '\nDegree Centrality:'
v1 = numpy.zeros((G_Nodes,1))

for i in range(0,G_Nodes): #loops takes the total number of nodes in the entire network
	a = G.neighbors_(i) #takes the neighbors of each node and puts them in an array
	N1 = len(a)
	wt = 0
	for j in range(0,N1): #loops through the neighbors of a particular node
		w2 = a[j]
		wt = wt + G.weight(G.node_object(a[j]),G.node_object(i))
	v1[i] = wt

print '\nThe top five active authors/professors based on degree centrality are:'
print_top(G,v1, num=5)


print '\n============================================='
# Eigenvector Centrality
print '\nEigenvector Centrality (by Zen):'
G1 = zen.algorithms.centrality.eigenvector_centrality_(G,weighted=True, max_iter = 1000)
print_top(G,G1, num=5)


print '\n============================================='
# Katz Centrality
print '\nThe top five authors based on Katz Centrality :'
num_steps = 10 # should be 10
alpha = 10
x = sparse.lil_matrix((G_Nodes, 1)) # the state vector - initial centrality vector
node_index = G.node_idx('Manfred Morari')
x[node_index] = 1
one_vec = numpy.ones((G_Nodes,1))
for i in range(num_steps):
	if i > 0:
		x = x/numpy.linalg.norm(x) # at each step we need to normalize the centrality vector
	elif i == 0:
		x = x/norm(x)
	x = alpha * A*x + one_vec # "pass" the centrality one step forward
print_top(G,x)

print '\n============================================='
# PageRank
print '\nThe top five authors based on PageRank'
x = sparse.lil_matrix((G_Nodes, 1)) # the state vector - initial centrality vector
D_inv = sparse.lil_matrix((G_Nodes,G_Nodes))
beta = 1
one_vec = numpy.ones((G_Nodes,1))
num_steps = 10 # should be 10
for i in range(0, G_Nodes):
	node_i_out_degree = G.degree_(i)
	if node_i_out_degree < 1:
		node_i_out_degree = 1
	D_inv[i,i] = 1/node_i_out_degree
AD_inv = A * D_inv
alpha = 0.85 # alpha between 0 and lambda_max, beta = 1
for i in range(num_steps):
	if i > 0:
		x = x/numpy.linalg.norm(x) # at each step we need to normalize the centrality vector
	elif i == 0:
		x = x/norm(x)
	x = alpha * AD_inv *x + one_vec # "pass" the centrality one step forward
print '\n Page Rank Centrality with alpha = %1.2f, beta = %i :' % (alpha, beta)
print_top(G,x)


print '\n============================================='
# Betweenness Centrality
print '\nBetweenness Centrality'
G4 = zen.algorithms.centrality.betweenness_centrality_(G)
print_top(G,G4)


# POWER LAW ==============================================
calc_powerlaw(G,10)  # need to change kmin appropriately

MODULARITY ========================================================================================
mit_list = []
ucb_list = []
stanford_list = []
caltech_list = []
umich_list = []
independent_list = []

for i in range(0, G_Nodes):
	if("Pasadena" in G.node_data_(i)['zenData']):
		caltech_list.append(G.node_object(i))
	elif("Massachusetts" in G.node_data_(i)['zenData']):
		mit_list.append(G.node_object(i))
	elif("Berkeley" in G.node_data_(i)['zenData']):
		ucb_list.append(G.node_object(i))
	elif("Stanford" in G.node_data_(i)['zenData']):
		stanford_list.append(G.node_object(i))
	elif("Ann Arbor" in G.node_data_(i)['zenData']):
		umich_list.append(G.node_object(i))
	elif("Independent" in G.node_data_(i)['zenData']):
		independent_list.append(G.node_object(i))

caltech = len(caltech_list)
mit = len(mit_list)
stanford = len(stanford_list)
berkeley = len(ucb_list)
michigan = len(umich_list)
independent = len(independent_list)

print '\nCaltech accounts for %i authors out of %i total authors in the network' %(caltech, G_Nodes)
print '\nMIT accounts for %i authors out of %i total authors in the network' %(mit, G_Nodes)
print '\nStanford accounts for %i authors out of %i total authors in the network' %(stanford, G_Nodes)
print '\nBerkeley accounts for %i authors out of %i total authors in the network' %(berkeley, G_Nodes)
print '\nMichigan accounts for %i authors out of %i total authors in the network' %(michigan, G_Nodes)
print '\nIndependent Consultants account for %i authors out of %i total authors in the network' %(independent, G_Nodes)

########### spectral_modularity_maximization - Didn't finish #########################################################
# B = spectral_modularity_maximization(G, A)
# print '\nSuccessfully Obtained B Matrix for Spectral Modularity Maximization'
# B_max_eig, B_max_eig_vec = sparse.linalg.eigs(B, k=1, which='LM', return_eigenvectors=True)
# Q = B_max_eig.real * G_Nodes / (4 * G_Edges)
# print 'Spectral Modularity Maximization Q_max: %1.4f' %Q

c = {
	'Caltech': caltech_list,
	'MIT': mit_list,
	'UCB': ucb_list,
	'Stanford': stanford_list,
	'Michigan': umich_list,
	'Independent': independent_list
	}

Q, Qmax = modularity(G,c)
print 'Modularity: %1.4f / %1.4f' % (Q,Qmax)

# ########################################################################################
# # Configuraton Model of Given Sample Networks
# ########################################################################################
check_friendship_paradox(G)

########################################################################################
# DIFFUSION DYNAMICS, INFECTION AND INFLUENCE MODELS
########################################################################################

dt = 0.001 # the "infintesimal" size steps we take to integrate
T = 5 # the end of the simulation time
time = linspace(0,T,int(T/dt)) # the array of time points spaced by dt

# # DIFFUSION ==============================================
# print '============================\nDIFFUSION\n'
# x = numpy.zeros((G_Nodes, 1))
# node_index = G.node_idx('Manfred Morari')
# x[node_index] = 10
# I = sparse.eye(G_Nodes)
# D = I*A.sum(1)
# L = sparse.csgraph.laplacian(A)
# print 'Successfully obtained Laplacian Matrix'
# min_eig, zero_eig_vec = sparse.linalg.eigs(L, k=1, which='SM', return_eigenvectors=True)
# norm_zero_eig_vec = numpy.linalg.norm(zero_eig_vec)
# normal_zero_eig_vec = zero_eig_vec / norm_zero_eig_vec
# print 'Simulating Diffusion of Idea From Professor Manfred Morari...'
# diff_const = 1*dt;
# error = numpy.zeros(len(time))
# mean_x = numpy.zeros(len(time))
# for i,t in enumerate(time):
# 	# at each time point update the value of x
# 	x = x - diff_const * L * x
# 	x = x / numpy.linalg.norm(x)
# 	mean_x[i] = mean(x)
# 	error_value = x-normal_zero_eig_vec
# 	error[i] = numpy.linalg.norm(error_value)
# 	#error[i] = numpy.sqrt(numpy.sum((x-normal_zero_eig_vec)**2))
# 	x = x * numpy.linalg.norm(x)

# print '\nSimulating steady state behaviour of diffusion process using laplacian...'
# plt.figure()
# plt.plot(time,error)
# plt.xlabel('time')
# plt.ylabel('error') 
# plt.title('Euclidean Distance between Actual State and Equilibrium State')
# plt.savefig('diffusion_steady_state.png')

## SI MODEL ===============================================
print '============================\nSI MODEL\n'
#x = sparse.lil_matrix((G_Nodes, 1)) # the state vector
x = numpy.zeros((G_Nodes, 1))
dx = numpy.zeros(len(time))
node_index = G.node_idx('U Shaked')
x[node_index] = 10 
s = 1 - x
beta = 1.0
constant = beta*dt
for i,t in enumerate(time):
	# at each time point update the value of x
	term = constant*s*A.dot(x)
	x = x + term
	s = s - term
	dx[i] = sum(x)
plt.figure()
plt.plot(time, dx.transpose(), 'r') 
plt.xlabel('t')
plt.ylabel('Infected Population') 
plt.title('SI Model')
plt.savefig('si_model.png')


## SIR MODEL WITH EIGENVALUE CALCULATION==========================================
max_eig, max_eig_vec = sparse.linalg.eigs(A, k=1, which='LM', return_eigenvectors=True)
print 'Maximum Eigen Value = %1.4f' %max_eig.real
x = numpy.zeros(G_Nodes) # the state vector
r = numpy.zeros(G_Nodes) # the state vector
node_index = G.node_idx('U Shaked')
x[node_index] = 1
dx = numpy.zeros(len(time))
ds = numpy.zeros(len(time))
dr = numpy.zeros(len(time))
s = numpy.ones((G_Nodes,1))
s = 1 - x
beta = 1.0
gamma = beta * max_eig
gamma = gamma + 1 # just to make sure -> gamma > beta * eig_max so that network is not epidemic
for i,t in enumerate(time):
	# at each time point update the value of x
	term = beta * s*A.dot(x) * dt
	x = x + term - gamma*x*dt
	s = s - term
	#r = gamma*x*dt
	r = 1 - x - s
	dx[i] = numpy.sum(x)
	ds[i] = numpy.sum(s)
	dr[i] = numpy.sum(r)

plt.figure()
plt.plot(time, dx.transpose(), 'r', label='Infected State') 
plt.plot(time, ds.transpose(), 'b', label='Susceptible State')  
plt.plot(time, dr.transpose(), 'g', label='Recovered State') 
plt.xlabel('t')
plt.ylabel('Infected, Susceptible & Recovered Population') 
plt.title('SIR Model')
plt.legend()
plt.savefig('sir_model.png')

# INDEPENDENT CASCADE ====================================
print '============================\nINDEPENDENT CASCADE MODEL\n'
max_eig, max_eig_vec = sparse.linalg.eigs(A, k=1, which='LM', return_eigenvectors=True)
beta = 1.0
gamma = beta * max_eig.real
gamma = gamma + 1 # just to make sure -> gamma > beta * eig_max so that network is not epidemic
repeats = 1000
x = zeros(G_Nodes) # the state vector
p = 1 - exp(-beta/gamma) # Uniform Probability 
print '\n Probability = %1.4f' %p
node_index = G.node_idx('U Shaked')
x[node_index] = 1
x_prob=numpy.zeros((repeats,G_Nodes))
for i in range(0,repeats):
    x_prob[i,:]= influence_cascade(G,p,x)    
x_i_mean_1 = numpy.zeros( (G_Nodes,) )
for i in range(G_Nodes):
    for j in range(repeats):
        x_i_mean_1[i] = x_i_mean_1[i] + x_prob[j,i]/repeats

print '\n Influence Models using Percolation Process'
x = zeros(G_Nodes) # the state vector
node_index = G.node_idx('U Shaked')
x[node_index] = 1
x_prob=numpy.zeros((repeats,G_Nodes))
for i in range(0,repeats):
    x_prob[i,:]= independent_cascade(G,p,x)    
x_i_mean_2 = numpy.zeros( (G_Nodes,) )
for i in range(G_Nodes):
    for j in range(repeats):
        x_i_mean_2[i] = x_i_mean_2[i] + x_prob[j,i]/repeats

node_index_array = numpy.arange(0,G_Nodes,1)
plt.figure()
plt.plot(node_index_array, x_i_mean_1, 'r', label='Online Cascade') 
plt.plot(node_index_array, x_i_mean_2, 'b', label='Coin Toss Based Percolation') 
plt.plot(node_index_array, r, 'g', label='Late Time SIR') 
plt.xlabel('Node Indices')
plt.ylabel('Probability')
plt.title('Expected Node Activation Probability in Influence Models') 
plt.legend()
plt.savefig('influence_models.png')

#d3.stop_server()