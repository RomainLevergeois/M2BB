#!/usr/bin/env python3 

#----------------------------------------
# gitLab : https://gitlab.univ-nantes.fr/E18F227W/m2bb/-/tree/main/Language_de_script
#
# This script allows you to align 
# two by two of all the sequences contained 
# in a fasta file '<test2>.fasta'.
#
# The output is a network graph representing 
# sequences as nodes connected by
# edges with a color representing the alignment score.
#
# Author : Romain Levergeois
# School : University of Nantes
# Curriculum: Master 2 Bioinformatics for Biologists
#
# Script set up in the framework of the Master and free of rights.
#----------------------------------------

# We import the necessary modules
import functions as fct
import sys 
import argparse
import re 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import networkx as nx 
from netgraph import Graph
from Bio import pairwise2


# MAIN
## Add parser for arguments.
parser = argparse.ArgumentParser(
	description='This script allows you to perform a two by two alignment of all sequences contained in a fasta file.\nThe output is a network graph representing the sequences as nodes connected by edges with a color representing the alignment score.')

parser.add_argument("fasta_file",
					metavar='Path_to_file',
					type=str,
					action='store',
					help='Fasta file for alignement : .fasta or .fa')

parser.add_argument('--edge_label',
					metavar='bool',
					type=str,action='store',
					help='Draw edge label or not', 
					required=False, 
					default='False')

parser.add_argument('--node_label',
					metavar='bool',
					type=str, 
					action='store', 
					help='Draw node label or not',
					required=False,  
					default='False')

parser.add_argument('--colorbar',
					metavar='bool',
					type=str, 
					action='store', 
					help='Draw node label or not',
					required=False,  
					default='True')

parser.add_argument('--threshold',
					metavar='float',
					type=float ,
					action='store', 
					help='Set a threshold for drawing edges that have a label >= threshold',
					required=False,  
					default=0.0)

args = parser.parse_args()

## Definition of variables for the alignment
dict_sequences = fct.getfasta(args.fasta_file)
list_sequences = list(dict_sequences.values())
result_align = {}
header = list(dict_sequences.keys())

## Double loop to align all sequences 2 by 2
## contained in the fasta file. The result of the alignments (score of the
## best alignment) is stored in a dictionary whose key is the 
## combination in tuple of the 2 names of the 2 aligned sequences.


for i in range(0,len(list_sequences)-1) :
	for j in range(i+1,len(list_sequences)) :
		result_align[str(header[i]),str(header[j])] = pairwise2.align.globalxx(list_sequences[i],list_sequences[j],score_only = True)

if (args.threshold>max(result_align.values())):
	print("\n---------------------------TP1_languages_script.py------------------------------\nError:[4] <threshold> parameter is higher than the value of the highest alignment score, nothing to display.\n", file=sys.stderr)
	exit(4)


## Construction of the network graph
network = nx.DiGraph()

### The nodes to be displayed are selected according to the threshold entered in parameter.
dict_nodes_displayed = {}
node_color = {}
if args.node_label == 'False' :
	node_label2 = False
else :
	node_label2= True
	print("\n---------------------------TP1_languages_script.py------------------------------\nWarning: --node_label parameter set to true, pay attention to the length of sequences header (ID should be <= 6 chars)\n", file=sys.stderr)

color_name = ["black","gray","darkgray","lightgray","rosybrown","lightcoral","indianred","brown","firebrick","maroon","red","tomato","lightsalmon","orangered","sienna","bisque","darkorange","tan","orange","khaki","olive","yellow","yellowgreen","darkalivegreen","lawngreen","forestgreen","limegreen","darkgreen","lime","springgreen","aquamarine","lightseagreen","darkslategray","darkcyan","cyan", "steelblue","slategray","royalblue","midnightblue","madiumblue","blueviolet","darkviolet","plum","purple","magenta","deeppink","crimson"]
dic_color_to_rgb= {'aliceblue': '#F0F8FF', 'antiquewhite': '#FAEBD7', 'aqua': '#00FFFF', 'aquamarine': '#7FFFD4', 'azure': '#F0FFFF', 'beige': '#F5F5DC', 'bisque': '#FFE4C4', 'black': '#000000', 'blanchedalmond': '#FFEBCD', 'blue': '#0000FF', 'blueviolet': '#8A2BE2', 'brown': '#A52A2A', 'burlywood': '#DEB887', 'cadetblue': '#5F9EA0', 'chartreuse': '#7FFF00', 'chocolate': '#D2691E', 'coral': '#FF7F50', 'cornflowerblue': '#6495ED', 'cornsilk': '#FFF8DC', 'crimson': '#DC143C', 'cyan': '#00FFFF', 'darkblue': '#00008B', 'darkcyan': '#008B8B', 'darkgoldenrod': '#B8860B', 'darkgray': '#A9A9A9', 'darkgreen': '#006400', 'darkgrey': '#A9A9A9', 'darkkhaki': '#BDB76B', 'darkmagenta': '#8B008B', 'darkolivegreen': '#556B2F', 'darkorange': '#FF8C00', 'darkorchid': '#9932CC', 'darkred': '#8B0000', 'darksalmon': '#E9967A', 'darkseagreen': '#8FBC8F', 'darkslateblue': '#483D8B', 'darkslategray': '#2F4F4F', 'darkslategrey': '#2F4F4F', 'darkturquoise': '#00CED1', 'darkviolet': '#9400D3', 'deeppink': '#FF1493', 'deepskyblue': '#00BFFF', 'dimgray': '#696969', 'dimgrey': '#696969', 'dodgerblue': '#1E90FF', 'firebrick': '#B22222', 'floralwhite': '#FFFAF0', 'forestgreen': '#228B22', 'fuchsia': '#FF00FF', 'gainsboro': '#DCDCDC', 'ghostwhite': '#F8F8FF', 'gold': '#FFD700', 'goldenrod': '#DAA520', 'gray': '#808080', 'green': '#008000', 'greenyellow': '#ADFF2F', 'grey': '#808080', 'honeydew': '#F0FFF0', 'hotpink': '#FF69B4', 'indianred': '#CD5C5C', 'indigo': '#4B0082', 'ivory': '#FFFFF0', 'khaki': '#F0E68C', 'lavender': '#E6E6FA', 'lavenderblush': '#FFF0F5', 'lawngreen': '#7CFC00', 'lemonchiffon': '#FFFACD', 'lightblue': '#ADD8E6', 'lightcoral': '#F08080', 'lightcyan': '#E0FFFF', 'lightgoldenrodyellow': '#FAFAD2', 'lightgray': '#D3D3D3', 'lightgreen': '#90EE90', 'lightgrey': '#D3D3D3', 'lightpink': '#FFB6C1', 'lightsalmon': '#FFA07A', 'lightseagreen': '#20B2AA', 'lightskyblue': '#87CEFA', 'lightslategray': '#778899', 'lightslategrey': '#778899', 'lightsteelblue': '#B0C4DE', 'lightyellow': '#FFFFE0', 'lime': '#00FF00', 'limegreen': '#32CD32', 'linen': '#FAF0E6', 'magenta': '#FF00FF', 'maroon': '#800000', 'mediumaquamarine': '#66CDAA', 'mediumblue': '#0000CD', 'mediumorchid': '#BA55D3', 'mediumpurple': '#9370DB', 'mediumseagreen': '#3CB371', 'mediumslateblue': '#7B68EE', 'mediumspringgreen': '#00FA9A', 'mediumturquoise': '#48D1CC', 'mediumvioletred': '#C71585', 'midnightblue': '#191970', 'mintcream': '#F5FFFA', 'mistyrose': '#FFE4E1', 'moccasin': '#FFE4B5', 'navajowhite': '#FFDEAD', 'navy': '#000080', 'oldlace': '#FDF5E6', 'olive': '#808000', 'olivedrab': '#6B8E23', 'orange': '#FFA500', 'orangered': '#FF4500', 'orchid': '#DA70D6', 'palegoldenrod': '#EEE8AA', 'palegreen': '#98FB98', 'paleturquoise': '#AFEEEE', 'palevioletred': '#DB7093', 'papayawhip': '#FFEFD5', 'peachpuff': '#FFDAB9', 'peru': '#CD853F', 'pink': '#FFC0CB', 'plum': '#DDA0DD', 'powderblue': '#B0E0E6', 'purple': '#800080', 'rebeccapurple': '#663399', 'red': '#FF0000', 'rosybrown': '#BC8F8F', 'royalblue': '#4169E1', 'saddlebrown': '#8B4513', 'salmon': '#FA8072', 'sandybrown': '#F4A460', 'seagreen': '#2E8B57', 'seashell': '#FFF5EE', 'sienna': '#A0522D', 'silver': '#C0C0C0', 'skyblue': '#87CEEB', 'slateblue': '#6A5ACD', 'slategray': '#708090', 'slategrey': '#708090', 'snow': '#FFFAFA', 'springgreen': '#00FF7F', 'steelblue': '#4682B4', 'tan': '#D2B48C', 'teal': '#008080', 'thistle': '#D8BFD8', 'tomato': '#FF6347', 'turquoise': '#40E0D0', 'violet': '#EE82EE', 'wheat': '#F5DEB3', 'white': '#FFFFFF', 'whitesmoke': '#F5F5F5', 'yellow': '#FFFF00', 'yellowgreen': '#9ACD32'}

for seq in header :
	k = 0
	for key, value in result_align.items() :
		if (key[0]==seq or key[1]==seq) :
			try :
				dict_nodes_displayed[str(seq)].append(value)
			except KeyError :
				dict_nodes_displayed[str(seq)]=[value]
		k=k+1

if node_label2==False :
	if len(color_name)<len(dict_nodes_displayed.keys()) :
		print("\n---------------------------TP1_languages_script.py------------------------------\nWarning: --node_label parameter set to true. (to many nodes for the colorset)\n", file=sys.stderr)
		print("\n---------------------------TP1_languages_script.py------------------------------\nWarning: --node_label parameter set to true, pay attention to the length of sequences header (ID should be <= 6 chars)\n", file=sys.stderr)
		node_label2=True
else :
	node_color='green'
	
for key, value in dict_nodes_displayed.items() :
	cpt = 0
	for item in value :
		if item < args.threshold :
			cpt = cpt+1

	if cpt!=len(dict_nodes_displayed.keys())-1 :
		if node_label2==True :	
			network.add_node(key[0:7])
		else :
			network.add_node(key)
			node_color[key] = dic_color_to_rgb[color_name[0]]
	color_name.pop(0)

### We add the edges to the graph
dict_edges_displayed = {}  

#### We save the keys of edges that we will display in relation
#### to the limit imposed by the user and save the score.
for key in list(tuple(result_align.keys())) :
	if result_align[key]>=args.threshold :
		if node_label2==True :
			dict_edges_displayed[(str(key[0])[0:7],str(key[1])[0:7])]=result_align[key]
		else :
			dict_edges_displayed[key]=result_align[key]

for key in dict_edges_displayed.keys() :
		network.add_edge(*key, label=dict_edges_displayed[key])


### Definitions of graph parameters
pos = nx.circular_layout(network)
if args.edge_label=='True' :	
	edge_labels = nx.get_edge_attributes(network, 'label')
else :
	edge_labels = False

### We associate a color to the edges corresponding to the alignment score
valeurs_uniques = set(result_align.values())
valeurs_uniques = list(valeurs_uniques)
valeurs_uniques_sorted = sorted(valeurs_uniques, key = lambda x:float(x))

colors = mpl.colormaps['jet'].resampled(len(valeurs_uniques_sorted))
colors = list(colors(range(len(valeurs_uniques_sorted))))
edge_color={}

for h in range(0,len(dict_edges_displayed.keys())) :
	edge_color[list(dict_edges_displayed.keys())[h]]=colors[valeurs_uniques_sorted.index(list(dict_edges_displayed.values())[h])]

### We create the figure with matplotlib.
plt.subplots(3,2 ,sharey='col', gridspec_kw={'width_ratios': [40, 1],'height_ratios': [1,5,1]})

#### MatplotlibDeprecationWarning: Auto-removal of overlapping axes is deprecated since 3.6 
#### and will be removed two minor releases later; explicitly call ax.remove() as needed :
plt.subplot(3,2,1).remove()
plt.subplot(3,2,3).remove()
plt.subplot(3,2,5).remove()

subplot_map = plt.subplot(3,2,(1,5)) 

### We draw the grace network with netgraph.Graph() 
### github : "https://github.com/paulbrodersen/netgraph"
graph = Graph(network, edge_layout='curved',
	origin=(-1, -1),
	scale=(2, 2),
	node_color=node_color,
	node_layout=pos,
	node_size=8,
    node_labels=node_label2,
	node_label_fontdict=dict(size=6),
    edge_labels=edge_labels,
	edge_label_fontdict=dict(size=6),
	edge_alpha=1,
	edge_label_position=0.7,
	edge_color=edge_color,
	edge_width=1
)

### We adjust the borders of the figure
plt.subplots_adjust(
	top=1.0,
	bottom=0.0,
	left=0.0,
	right=0.910,
	hspace=0.0,
	wspace=0.0
)

### Construction of the edges legend : color gradient presented in colormap.
gradient1 = np.linspace(int(min(valeurs_uniques_sorted)),int(max(valeurs_uniques_sorted)),256)
gradient2 = np.linspace(int(min(valeurs_uniques_sorted)),int(max(valeurs_uniques_sorted)),256)
gradient = np.vstack((gradient1,gradient2))

### We create the subplot for the edges legend which is a colormap corresponding to the score of the best alignment.
subplot_legend = plt.subplot(324)

if args.colorbar=='True' :
	norm = mpl.colors.Normalize(vmin=min(result_align.values()), vmax=max(result_align.values()))
	mpl.colorbar.ColorbarBase(subplot_legend, cmap='jet',
									norm=norm,
									orientation='vertical',                                
									ticks=[min(result_align.values()),(((max(result_align.values())-min(result_align.values()))/2)+min(result_align.values())),max(result_align.values())])

	subplot_legend.text(
		1,
		1.02,
		'Score',
		va='center',
		ha='right',
		fontsize=10,
		transform=subplot_legend.transAxes
	)
else : 
	subplot_legend.set_axis_off()
	plt.subplots_adjust(
		top=1.0,
		bottom=0.0,
		left=0.0,
		right=1.0,
		hspace=0.0,
		wspace=0.0)

### We remove the axes, the borders of the box from the legended subplot
subplot_legend.spines.right.set_visible(False)
subplot_legend.spines.top.set_visible(False)
subplot_legend.get_xaxis().set_visible(False)

plt.subplot(322).set_axis_off()
plt.subplot(326).set_axis_off()

### We save the network figure in svg format
plt.savefig('network_netgraph.svg')

### Construction of the nodes legend if the parameter --node_label is equal to False 
if args.node_label=='False' :
	fig_legend = plt.figure()
	subplot_legend = plt.subplot()
	subplot_legend.set_axis_off()

	text_kwargs = dict(ha='left', va='center', fontsize=10)
	x_coord = []
	y_coord = []
	color_scatter = list(node_color.values())
	for i in range (1, len(node_color.keys())+1) :
		x_coord.append(1)
		y_coord.append(i)
		subplot_legend.text(1.2, i, str(list(node_color.keys())[i-1]), **text_kwargs)

	subplot_legend.scatter(x_coord, y_coord, s = 100, color = color_scatter )
	subplot_legend.set_aspect(0.4)
	plt.subplots_adjust(
		top=1.0,
		bottom=0.0,
		left=0.0,
		right=0.9,
		hspace=0.2,
		wspace=0.2
	)
	plt.xlim(0.9,6)

	### We save the nodes legend
	plt.savefig('legend.svg')

## We show figures
plt.show()

## Everything went well so we return the error code 0
exit(0)
# END