#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_networkx.py
#  
#  Copyright 2021 123 <123@DESKTOP-GGMH6CF>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import networkx as nx
oo = float('inf')
def createGraph(filename) :
    G = nx.DiGraph()
    for line in open(filename) :
        strlist = line.split()
        n1 = int(strlist[0])
        n2 = int(strlist[1])
        G.add_edges_from([(n1, n2)]) #G.add_edges_from([(n1, n2)])
    return G
# create directed graph from file.
G = createGraph("out.txt")
print(G.edges)



d_gen = nx.dfs_edges(G,0)  #  按边深度搜索, 1为起点
b_gen = nx.bfs_edges(G,0)
print(list(d_gen), list(b_gen))
import matplotlib.pyplot as plt

pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos)
nx.draw_networkx_edges(G, pos)
nx.draw_networkx_labels(G, pos)
plt.show()
#print(nx.dfs_tree(G,1).nodes())  # 按点深搜
