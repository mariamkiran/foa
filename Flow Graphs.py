
# coding: utf-8

# In[1]:



import csv
import networkx as nx
import matplotlib.pyplot as plt
import pygraphviz as pgvz

import pandas as pd
import numpy  as np

G  = nx.DiGraph()

pos = {}
node_label_pos = {}

with open('esnet6_proposed_footprint.txt') as f:
    reader = csv.DictReader(f , delimiter='\t')
            
    for row in reader :
        
        pos[row['SITE']]=(int(row['X']) , int(row['Y']))
        node_label_pos[row['SITE']] = (pos[row['SITE']][0] + 3 , pos[row['SITE']][1] + 3)

        G.add_node(row['SITE'] , SITE=row['SITE'],PE=float(row['PE']))

        total_bw = 0;
        if (row['LINK1'] is not None): 
            G.add_edges_from([ (row['SITE'] , row['LINK1']) , (row['LINK1'] , row['SITE']) ],
                       capacity=float(row['LINK1BW']) , 
                       impedance=1/(float(row['LINK1BW']))
                      )
            total_bw += float(row['LINK1BW'])

        if (row['LINK2'] is not None): 
            G.add_edges_from([ (row['SITE'] , row['LINK2']) , (row['LINK2'] , row['SITE']) ],
                       capacity=float(row['LINK2BW']) ,
                       impedance=1/(float(row['LINK2BW']))                       
                      )
            total_bw += float(row['LINK2BW'])
            
        if (row['LINK3'] is not None): 
            G.add_edges_from([ (row['SITE'] , row['LINK3']) , (row['LINK3'] , row['SITE']) ],
                       capacity=float(row['LINK3BW']) ,
                       impedance=1/(float(row['LINK3BW']))                       
                      ) 
            total_bw += float(row['LINK3BW'])
            
        if (row['LINK4'] is not None): 
            G.add_edges_from([ (row['SITE'] , row['LINK4']) , (row['LINK4'] , row['SITE']) ],
                       capacity=float(row['LINK4BW']) ,
                       impedance=1/(float(row['LINK4BW']))                       
                      )  
            total_bw += float(row['LINK4BW'])    
        
        G.node[row['SITE']]['bw'] = total_bw
        G.node[row['SITE']]['type'] = row['TYPE']
#       nx.set_node_attributes(G,'bw',{row['SITE']:total_bw})
#       nx.set_node_attributes(G,'type',{row['SITE']:row['TYPE']})        


types = nx.get_node_attributes(G,'type')
end_sites  = [ k for k,v in types.items() if v == 'S']
core_sites = [ k for k,v in types.items() if v == 'C']        

CG = G.copy()
CG.remove_nodes_from(end_sites)

nx.write_graphml(G,'network.graphml')


# In[2]:


#--------------------------------------------------------------
# calculate some global variables that can be used at any time
# in order to determine node and edge sizes based on capacity
#--------------------------------------------------------------

edges = G.edges()
link_bw = nx.get_edge_attributes(G,'capacity')


# In[3]:


NodeSizeFactor  = 100    #  NodeSize = BW in Tbits per second x NodeSizeFactor 
EdgeWidthFactor = 4      #  EdgeWidth = BW in Tbits per second x EdgeWidthFactor


#----------------------------------------------------------------------
# Shortcut functions to draw maps.  Sometimes they use global variables
# for expediancy.  
#----------------------------------------------------------------------
    
def drawBaseMap(G , node_color='grey' , edges=True , alpha=0.5 , arrows=False , campus=True) :

    types = nx.get_node_attributes(G,'type')
    end_sites  = [ k for k,v in types.items() if v == 'S']
    core_sites = [ k for k,v in types.items() if v == 'C']
    
    bw = nx.get_node_attributes(G,'bw')
    pe = nx.get_node_attributes(G,'PE')
       
    bw_size = [ NodeSizeFactor * bw[site] for site in core_sites ]
    nx.draw_networkx_nodes(G , nodelist=core_sites, pos=pos , node_size=bw_size , 
                           node_color='white', edgecolors='black',alpha=alpha)

    if campus :
        pe_size = [ NodeSizeFactor * pe[site] for site in end_sites ]
        nx.draw_networkx_nodes(G , nodelist=end_sites,  pos=pos , node_size=pe_size , 
                               node_color='lightblue', edgecolors='black', alpha=alpha)

    pe_size = [ NodeSizeFactor * pe[site] for site in core_sites ]    
    nx.draw_networkx_nodes(G , nodelist=core_sites, pos=pos , node_size=pe_size , 
                           node_color=node_color, edgecolors='black' , alpha=alpha)

    nx.draw_networkx_labels(G , pos=node_label_pos , font_size=10)
    
    if (edges) : 
        if (arrows) : 
            width = [EdgeWidthFactor*link_bw[e] for e in G.edges()]            
            nx.draw_networkx_edges (G , pos=pos, width=width, edge_color='lightgrey', alpha=alpha)
        else :
            width = [5*link_bw[e] for e in G.edges()]
            nx.draw_networkx_edges (G , pos=pos, width=width, edge_color='lightgrey', alpha=alpha , arrows=False)  
        
        edge_bw = nx.get_edge_attributes(G,'capacity')
        nx.draw_networkx_edge_labels(G , pos=pos, edge_labels=edge_bw, label_pos=0.75, 
                                     edge_color='lightgrey', font_size=9, alpha=alpha)


def highlightNodes(G , nodelist=None , node_color='red' , alpha=0.5 , node_shape='o') :
    if (not nodelist) : nodelist = G.nodes()
    nx.draw_networkx_nodes(G , pos=pos , nodelist=nodelist , 
                           node_size=[1200 , 1200] , 
                           node_color=node_color, edgecolors='black',
                           alpha=alpha , node_shape=node_shape
                          )
    
def hightlightEdges( G , edgelist=None , edge_color='red' , rate=0.1 , alpha=0.5, labels=False) :
    if (not edgelist) : edgelist = G.edges()
    nx.draw_networkx_edges (G , pos=pos, 
                                   edgelist=edgelist, 
                                   edge_color=edge_color,
                                   width=EdgeWidthFactor*rate , 
                                   alpha = alpha
                                 )
    if labels :
        edge_labels = []
        edge_labels = edge_labels.append(rate)
        nx.draw_networkx_edge_labels(G , edgelist=edgelist, pos=pos, edge_labels=rate, label_pos=0.75, 
                                     edge_color='lightgrey', font_size=9, alpha=alpha)






# First visualize the topology and capacity based on the model that was read in from the external CSV file.

# In[26]:


fig = plt.figure(figsize=(18,15))
plt.axis('off')

#drawBaseMap(CG , node_color='yellow' , arrows=False , alpha=0.5 )
drawBaseMap(G  , node_color='yellow' , arrows=False , alpha=0.5 )

fig.suptitle('Link and Site Bandwidth Summary for ESnet 6 : L2/L3 only', fontsize=20)

plt.show()


# In[27]:


from networkx.algorithms.flow import shortest_augmenting_path
from networkx.algorithms.flow import edmonds_karp

flow_source = 'BNL-2'
all_flows = {}
all_flow_rates = {}
for site in end_sites:
    flow_dest = site
    flow_rate , flows = nx.maximum_flow(G, flow_source, flow_dest , flow_func=edmonds_karp)
    all_flows[(flow_source,flow_dest)] = flows
    all_flow_rates[(flow_source,flow_dest)] = flow_rate

flow_dest = 'UCSD'
fig = plt.figure(figsize=(18,15))
plt.subplot(2,1,1)
plt.axis('off')

drawBaseMap(G,alpha=0.1,edges=False)   # Draw the background image of the network
#highlightNodes(G , nodelist=[flow_source, flow_dest], node_color=['yellowgreen','purple'])

# Plot the FLOW paths.  See nx.maximum_flows() to understand the "flows" datastructure.


flow_dest = 'UCSD'
flows = all_flows[(flow_source,flow_dest)]
flow_rate = all_flow_rates[(flow_source,flow_dest)]

for site , links in flows.items() :
    for link , rate in links.items() :
        if (rate > 0) :
            hightlightEdges(G , edgelist=[(site,link)] , rate=rate ,  alpha=0.5)
            
plt.title('Max FLOW possible from : ' + flow_source + " to " + flow_dest + " is %2.2f Tbps" % flow_rate, fontsize=20)



plt.subplot(2,1,2)
plt.axis('off')
drawBaseMap(G,alpha=0.1,edges=False)   # Draw the background image of the network

flow_dest = 'LBNL'            
flows = all_flows[(flow_source,flow_dest)]
flow_rate = all_flow_rates[(flow_source,flow_dest)]

for site , links in flows.items() :
    for link , rate in links.items() :
        if (rate > 0) :
            hightlightEdges(G , edgelist=[(site,link)] , rate=rate ,  
                            alpha=0.5 , edge_color='green' , labels=False)
                    
plt.title('Max FLOW possible from : ' + flow_source + " to " + flow_dest + " is %2.2f Tbps" % flow_rate, fontsize=20)

                
plt.show()


# In[30]:



flow_dest_list = end_sites

FG = nx.DiGraph()

for flow_dest in flow_dest_list:
    flows = all_flows[(flow_source, flow_dest)] 

#    print ("Processing target site " , flow_dest , "\n")
     
    for site , links in flows.items():
        for link_dest , rate in links.items() :
            if (rate > 0) : 
#                print (site , ":" , link_dest , " " , rate)
                FG.add_edge(site,link_dest, rate=rate)
                FG.node[site]['depth'] = 0
                FG.node[link_dest]['depth'] = 0
                FG.node[site]['height'] = 0
                FG.node[link_dest]['height'] = 0
                        
dfs_edges = [ (u,v) for u,v in nx.dfs_edges(FG,flow_source)]


height = 0
depth  = 0
for edge in dfs_edges:
    sn = edge[0]
    tn = edge[1]
    
    if (FG.node[sn]['depth'] < depth) : height += 1
    depth   = FG.node[sn]['depth'] + 1
    if FG.node[tn]['depth'] < depth : 
        FG.node[tn]['depth']  = depth
        FG.node[tn]['height'] = height 
        
#    print (tn , "=" , FG.node[tn] , "\t", sn , "=", FG.node[sn])

dfs_pos = {}
for node in FG.nodes() :
    dfs_pos[node] = (FG.node[node]['depth'],FG.node[node]['height'])

#print ("\n\n")
#print(dfs_pos)


# In[31]:


fig = plt.figure(figsize=(18,18))
plt.axis('off')

plt.title('All flow paths from %s to all campuses' % flow_source)

nx.draw_networkx_nodes(FG , pos=dfs_pos , node_size=1500 , node_shape='s',
                           node_color='white', edgecolors='black',alpha=1,labels=True)
nx.draw_networkx_edges (FG , pos=dfs_pos)
nx.draw_networkx_labels(FG , pos=dfs_pos, font_size=10)

plt.show()

