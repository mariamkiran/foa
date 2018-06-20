#  Set of routines useful for first order analysis of networks
#
#  Yatish Kumar
#  March 2018


import networkx as nx
import matplotlib.pyplot as plt

import pandas as pd
import numpy  as np

NodeSizeFactor  = 100    #  NodeSize = BW in Tbits per second x NodeSizeFactor 
EdgeWidthFactor = 4      #  EdgeWidth = BW in Tbits per second x EdgeWidthFactor

#----------------------------------------------------------------------

def importCSV(filename) :

    import csv

    G  = nx.DiGraph()
    pos = {}

    with open(filename) as f:
        reader = csv.DictReader(f , delimiter='\t')

        for row in reader :

            pos[row['SITE']]=(int(row['X']) , int(row['Y']))

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
            
    return (G , pos)


#----------------------------------------------------------------------
# Shortcut functions to draw maps.  Sometimes they use global variables
# for expediancy.  
#----------------------------------------------------------------------
    
def drawBaseMap(G , pos , node_color='grey' , edges=True , alpha=0.5 , arrows=False , campus=True) :

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

    node_label_pos = {}
    for k,v in pos.items():
        node_label_pos[k] = ( v[0] + 3, v[1] + 3 )
        
    nx.draw_networkx_labels(G , pos=node_label_pos , font_size=10)
    
    if (edges) : 
        if (arrows) : 
            width = [EdgeWidthFactor*G.edges[e]['capacity'] for e in G.edges()]            
            nx.draw_networkx_edges (G , pos=pos, width=width, edge_color='lightgrey', alpha=alpha)
        else :
            width = [EdgeWidthFactor*G.edges[e]['capacity'] for e in G.edges()]
            nx.draw_networkx_edges (G , pos=pos, width=width, edge_color='lightgrey', alpha=alpha , arrows=False)  
        
        edge_bw = nx.get_edge_attributes(G,'capacity')
        nx.draw_networkx_edge_labels(G , pos=pos, edge_labels=edge_bw, label_pos=0.75, 
                                     edge_color='lightgrey', font_size=9, alpha=alpha)


#----------------------------------------------------------------------
        
def highlightNodes(G , pos , nodelist=None , node_color='red' , alpha=0.5 , node_shape='o' , node_size=1200) :
    if (not nodelist) : nodelist = G.nodes()
    nx.draw_networkx_nodes(G , pos=pos , nodelist=nodelist , 
                           node_size=node_size , 
                           node_color=node_color, edgecolors='black',
                           alpha=alpha , node_shape=node_shape
                          )

#----------------------------------------------------------------------
    
def highlightEdges( G , pos , edgelist=None , edge_color='red' , font_weight='normal', rate=0.1 , alpha=0.5, labels=False) :

    if (not edgelist) : edgelist = G.edges()

    if isinstance(rate,dict):
        width = [EdgeWidthFactor*rate[e] for e in edgelist]
    elif isinstance(rate,list) :
        width = rate
    else :
        width = EdgeWidthFactor*rate
        
    nx.draw_networkx_edges (G , pos=pos, 
                                   edgelist=edgelist, 
                                   edge_color=edge_color,
                                   width=width , 
                                   alpha = alpha ,
                                   arrowstyle='simple'
                            )
    edge_labels = {}
    if labels :
        if isinstance(rate, dict):
            edge_labels = rate
        else :
            for edge in edgelist : edge_labels[edge] = rate

        nx.draw_networkx_edge_labels(G , edgelist=edgelist, pos=pos, edge_labels=edge_labels, label_pos=0.5, 
                                     font_color=edge_color, font_weight=font_weight, font_size=9, alpha=alpha)

#----------------------------------------------------------------------

def highlightEdgeVector(G,pos,edgelist=False,distance=[-0.5 , -0.25 , 0.25 , 0.5 ] , 
                        tail=0.2 , head = 0.2 , edge_color='black' , alpha = 1 ,
                        rate = 1 ,
                        style='tapered'
                       ):
    
    import math
    from math import cos
    from math import sin
    from math import pi
    
    def rotate(p , theta) :
        x = p[0]
        y = p[1]
        x1 =  x*cos(theta) + y*sin(theta)
        y1 = -x*sin(theta) + y*cos(theta)
        return (x1,y1)
    
    for delta in distance :
        if (edgelist == False) : edgelist = G.edges()
        for edge in edgelist :
            s = edge[0]
            t = edge[1]
            
            spos = pos[s]
            tpos = pos[t]

            spost = (0,0)
            tpost = (tpos[0] - spos[0],tpos[1]-spos[1])
            
            theta = math.acos( abs(tpost[0] - spost[0]) / math.sqrt((tpost[0]-spost[0])**2 + (tpost[1]-spost[1])**2) )


            if ((tpost[0] > spost[0]) & (tpost[1] < spost[1]))   : sgn = -1 
            elif ((tpost[0] < spost[0]) & (tpost[1] > spost[1])) : sgn = -1 
            else : sgn = 1

            theta = theta * sgn
            
            srpos = rotate(spost,theta)
            trpos = rotate(tpost,theta)

            if (srpos[0] > trpos[0]) : sgn = -1
            else :                     sgn = 1

            srpos = (srpos[0] + sgn*abs(tail) , srpos[1] + delta)
            trpos = (trpos[0] - sgn*abs(head) , trpos[1] + delta)

            spos1 = rotate(srpos,-theta)
            tpos1 = rotate(trpos,-theta)
 
            new_pos={}
            if (style == 'tapered') :
                new_pos[s] = spos
            else :
                new_pos[s] = (spos1[0]+spos[0],spos1[1]+spos[1])
                
            new_pos[t] = (tpos1[0]+spos[0],tpos1[1]+spos[1])
        
            nx.draw_networkx_edges(G,pos=new_pos , edgelist=[(s,t)] , width=EdgeWidthFactor*rate , edge_color=edge_color , alpha=alpha)   

#---------------------------------------------------------------------

def highlightFlows(G , pos , flow_sources, flow_dests, all_flows , edge_color='red') :

    flow_source = flow_sources[0]
    flow_dest = flow_dests[0]

    highlightNodes(G , pos=pos, nodelist=[flow_source, flow_dest], node_color=['yellowgreen','purple'])

    # Plot the FLOW paths.  See nx.maximum_flows() to understand the "flows" datastructure.

    flows = all_flows[(flow_source,flow_dest)]

    for site , links in flows.items() :
        for link , rate in links.items() :
            if (rate > 0) :
                highlightEdges(G , pos=pos, edgelist=[(site,link)] , rate=rate ,  
                                    alpha=0.5 , edge_color=edge_color , labels=True)

#---------------------------------------------------------------------
                
def computeMaxFlows (G , flow_sources , flow_dests) :

    from networkx.algorithms.flow import shortest_augmenting_path
    from networkx.algorithms.flow import edmonds_karp

    all_flows = {}
    all_flow_rates = {}

    for flow_source in flow_sources :
        for site in flow_dests:
            if ( flow_source != site ) :
                flow_dest = site
                flow_rate , flows = nx.maximum_flow(G, flow_source, flow_dest , flow_func=edmonds_karp)
                all_flows[(flow_source,flow_dest)] = flows
                all_flow_rates[(flow_source,flow_dest)] = flow_rate


    all_flow_rates_pd = pd.DataFrame(0.0 , columns = flow_sources , index = flow_dests)

    for link , rate in all_flow_rates.items() :
        all_flow_rates_pd[link[0]][link[1]] = rate

    total_flows_pd = pd.DataFrame(0.0 , columns = flow_sources , index = G.edges())
    for flow , sites in all_flows.items():
        for site , links in sites.items() :
            for link , rate in links.items() :
                if (rate > 0.0) :
                    total_flows_pd[flow[0]][(site,link)] += rate

    total_flows_pd['total'] = total_flows_pd.sum(axis=1) # add up total load on a link from all sources
    total_flows_pd = total_flows_pd[(total_flows_pd['total'] > 0)] # drop links that have absolutely no traffic

    return (total_flows_pd , all_flows , all_flow_rates , all_flow_rates_pd )

