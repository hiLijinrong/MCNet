import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from random import shuffle
import random
from math import ceil
from collections import defaultdict


# 设置点分布函数
def assign_positions_by_sector(G,rs,sector_dict):
    pos = {}
    sector_angle = 2 * np.pi / len(sector_dict)
    # 根据node的sector属性和node类型来设置位置
    for node,attr in G.nodes(data=True):
        node_sector = attr['sector']
        radius = np.random.uniform(*rs[attr['type']])
        angle = sector_angle * node_sector + sector_angle / 2
        angle_offset = np.random.uniform(-0.4*sector_angle , 0.4*sector_angle)
        pos[node] = (radius * np.cos(angle+angle_offset), radius * np.sin(angle+angle_offset))
    return pos

def extract_elements(lst): # 提取网络中的信息
    result = []
    for item in lst:
        if isinstance(item, (list, tuple)):
            result.extend(extract_elements(item))
        else:
            result.append(item)
    return result

def sort_breakdown(dmg_list,net,pop=True): 
    # 对损坏节点根据不同模式进行修复排序
    # 其中：pop=Ture表示默认以服务分口为排序依据
    dmg_weight = []
    for dmg in dmg_list:
        if pop:
            try:
                dmg_weight.append(net.nodes[dmg]['srv_popu'])
            except:
                dmg_weight.append(net.edges[dmg]['srv_popu'])
        else:
            degree = nx.degree_centrality(net)
            if isinstance(dmg,tuple):
                dmg_weight.append(np.mean([degree[item] for item in dmg]))
            else:
                dmg_weight.append(degree[dmg])
    return dmg_weight

def dmg2groups(lst, n): # 对顺坏节点进行重新分组
    group_size = ceil(len(lst) / n) 
    groups = []
    start = 0
    for i in range(n):
        end = start + group_size
        groups.append(lst[start:end])
        start = end
    return groups

def dmg_occur_time(t_seq,Nt,dt,dmg_list,groups,at):
    # 根据不同时间点进行采样
    dmg_occur = []
    dmg_groups = []
    grouped = dmg2groups(dmg_list,groups)
    for group in grouped:
        if group != []:
            dmg_groups.append(group)
    # 生成每个Groups所对应的时间点
    time_occur = random.sample(list(range(0,int((Nt-at)/dt))),Nt)
    if 0 not in time_occur:
        time_occur.append(0)
    time_occur = sorted([round(item,1) for item in time_occur])
    for idx in range(len(dmg_groups)):
        group_t = time_occur[idx]
        for dmg in dmg_groups[idx]:
            dmg_occur.append(t_seq[group_t])
    return dmg_occur

def sort_by_weights(elements, weights): # 根据权重对节点元素进行排序
    zipped = list(zip(elements, weights))
    sorted_zipped = sorted(zipped, key=lambda x: x[1], reverse=True)
    sorted_elements = [item[0] for item in sorted_zipped]
    return sorted_elements

def sort_by_occur_and_weights(elements, occur, weights):
    # 根据节点损坏时间和权重进行重新排序
    zipped = list(zip(elements, weights, occur))
    groups = defaultdict(list)
    for elem, weight, time in zipped:
        groups[time].append((elem, weight))
    sorted_groups = sorted(groups.items(), key=lambda x: x[0])
    sorted_elements = []
    for _, group in sorted_groups:
        sorted_group = sorted(group, key=lambda x: x[1], reverse=True)
        sorted_elements.extend([item[0] for item in sorted_group])
    return sorted_elements

def seq_dmgs_by_weight_occur(G,t_seq,Nt,dt,at,groups,dmgs,at_sametime=True,pop=True):
    # 根据节点损坏时间和权重得出修复顺序
    dmg_genes,dmg_edges,dmg_loads = dmgs # 电源节点，电源线和荷载节点三种类型损坏
    if at_sametime:
        dmgs = [dmg_genes, dmg_edges, dmg_loads]
        dmgs_w = [sort_breakdown(dmg,G,pop) for dmg in dmgs]
        dmgs_s = [sort_by_weights(dmg,dmg_w) for dmg,dmg_w in zip(dmgs,dmgs_w)]
        dmgs_seq = []
        for dmgs in dmgs_s:
            for dmg in dmgs:
                dmgs_seq.append(dmg)
        dmgs_occur = [0 for _ in range(len(dmgs_seq))]
    else:
        dmgs = dmg_genes+dmg_edges+dmg_loads
        dmgs_occur = dmg_occur_time(t_seq,Nt,dt,dmgs,groups,at)
        dmgs_w = sort_breakdown(dmgs,G,pop)
        dmgs_seq = sort_by_occur_and_weights(dmgs,dmgs_occur,dmgs_w)
    return dmgs_seq,dmgs_occur,dmgs_w

def whether_overloads(net, ab_loads, gnode, genes):
    # 超负荷判定
    back_up = []
    for ab_load in ab_loads:
        for gene in genes:
            if gene != gnode and nx.has_path(net,ab_load,gene):
                st_path = nx.shortest_path(net,ab_load,gene)
                count = sum(1 for node in st_path if ('R' in node or 'G' in node))
                if count<2 and gnode not in st_path:
                    back_up.append(st_path)
    return back_up

def get_nbs_roads(G,dmg_roads): 
    # 基于近邻关系的损失扩散模型
    nbs_roads = {}
    roadcount =[]
    roads = [road for roadlist in dmg_roads.values() for road in roadlist]
    for gene in dmg_roads.keys():
        road_edges = dmg_roads[gene]
        nbs_road = []
        for edge in road_edges:
            _,road_node = edge
            nbs = set(nx.neighbors(G,road_node))
            for nb in nbs:
                if 'L' not in nb and 'G' not in nb:
                    if all([(road_node,nb) not in roads,(nb,road_node) not in roads]):
                        nbs_road.append((road_node,nb))
        nbs_roads[gene] = nbs_road
        roadcount.append(len(nbs_road))
    return nbs_roads,roadcount

def refresh_color(G,eles,color,wd):
    # 对网络节点和连边根据属性来重新作色
    try:
        for edge in eles:
            G.edges[edge]['color'] = color
            G.edges[edge]['wd'] = wd
            G.edges[edge]['ls'] = 'solid'
    except:
        for node in eles:
            G.nodes[node]['color'] = color
            G.nodes[node]['ns'] = G.nodes[node]['ns']*1.5
    return True


def get_damges_by_zone(G,probs): 
    # 根据灾害类型进行绘图区域的分区
    # 初始化故障事件
    dmg_genes = []
    dmg_loads = []
    dmg_edges = []
    # 电源节点损坏
    for node, attr in G.nodes(data=True):
        if attr['type'] == 'generator':
            zone = attr['zone']
            if random.random() < probs[zone]['generator']:
                dmg_genes.append(node)

    # 荷载节点损坏
    for node, attr in G.nodes(data=True):
        if attr['type'] == 'load':
            zone = attr['zone']
            if random.random() < probs[zone]['load']:
                dmg_loads.append(node)
    # 线路损坏
    for u, v, attr in G.edges(data=True):
        if 'dmg_prob' in attr.keys():
            if random.random() < attr['dmg_prob']:
                dmg_edges.append((u, v))

    return dmg_genes,dmg_edges,dmg_loads



def damage_evluation(G,dmgs,gene_loads):
    # 根据异常类型对电网的故障节点进行移除
    dmg_nodes = extract_elements(dmgs)
    dmg_genes,dmg_edges,dmg_loads = dmgs
    dmg_nodes = dmg_genes + dmg_loads
    for load in dmg_nodes:
        try:
            G.remove_node(load)
        except:
            pass
    for edge in dmg_edges:
        try:
            G.remove_edge(edge)
        except:
            pass
    # 计算电网中所有的异常节点
    bugloads = {}
    for gene in gene_loads.keys():
        bugnodes = []
        if gene not in dmg_genes:
            for load in gene_loads[gene]:
                try:
                    all_paths = list(nx.shortest_path(G,load,gene))
                except:
                    bugnodes.append(load)
            if bugnodes:
                bugloads[gene] = bugnodes
    for gene in dmg_genes:
        bugloads[gene] = gene_loads[gene]
    return bugloads



def fixed_evluation(G,dmgs,gene_loads,dmg_genes):
    # 根据异常类型对电网的故障节点进行移除
    dmg_nodes = extract_elements(dmgs)
    for load in dmg_nodes:
        try:
            G.remove_node(load)
        except:
            pass
    for edge in dmg_nodes:
        try:
            G.remove_edge(edge)
        except:
            pass
    # 计算电网中所有的异常节点
    bugloads = {}
    for gene in gene_loads.keys():
        bugnodes = []
        if gene not in dmg_genes:
            for load in gene_loads[gene]:
                try:
                    all_paths = list(nx.shortest_path(G,load,gene))
                except:
                    bugnodes.append(load)
            if bugnodes:
                bugloads[gene] = bugnodes
    for gene in dmg_genes:
        bugloads[gene] = gene_loads[gene]
    
    return bugloads


# 计算节点的跨区域供电特性
def network_overload(G,bugloads,gene_loads):
    # 计算出当前电网所有潜在的超负荷运行情况
    backups = {}
    for gene1 in gene_loads.keys():
        backup = []
        for load in gene_loads[gene1]:
            for gene2 in gene_loads.keys():
                if gene2 != gene1:
                    try:
                        all_paths = list(nx.all_shortest_paths(G,load,gene2))
                        for st_path in all_paths:
                            count = sum(1 for item in st_path if 'G' in item)
                            if gene1 not in st_path and count<2: # 
                                backup.append([gene1,*st_path])
                    except:
                        pass
        if backup:
            backups[gene1] = backup

    # 计算当前边移除情况下的超负荷运行能力
    bug_nodes = {}
    inbackups = {}
    all_backs = []
    for gene in gene_loads.keys():
        withbackup = []
        if gene in bugloads.keys():
            for bugload in bugloads[gene]:
                inbackup = []
                if gene in backups.keys():
                    for back_up in backups[gene]:
                        if bugload in back_up:
                            inbackup.append(True)
                            all_backs.append(back_up)
                    if any(inbackup):
                        withbackup.append(bugload)
            if withbackup:
                inbackups[gene] = list(set(withbackup))
                # 考虑节点损坏的相互抵消效应
                bug_nodes[gene] = list(set(bugloads[gene])-set(withbackup))
                # 不考虑节点损坏的相互抵消效应
                #bug_nodes[gene] = list(set(bugloads[gene])) 
    return inbackups, backups, bug_nodes


def overload_evluation(rho,inbackups,backups,bug_nodes,bugloads,gene_loads):
    # 超负荷运行的多维度评估
    rho = rho  if rho ==2 else  1
    back_genes = {}
    for gene in inbackups.keys():
        back_paths =[]
        for load in inbackups[gene]:
            for st_path in backups[gene]:
                if load in st_path:
                    back_paths.append(st_path[1:])
        back_genes[gene] = back_paths
    overload_genes = {}
    for gene1 in gene_loads.keys():
        backs = []
        for gene2 in back_genes.keys():
            for back in back_genes[gene2]:
                if back[-1] == gene1:
                    backs.append(back[:-1])
        if backs:
            overload_genes[gene1] = list(set(extract_elements(backs)))
    # 更新最终的损坏节点
    init_dmg_loads = bugloads.copy()
    for gene in bug_nodes.keys():
        init_dmg_loads[gene] = bug_nodes[gene]
    # 计算超负荷时间,抵消超负荷节点本身的故障点
    orgs = {}
    overload_times = {}
    for gene in overload_genes.keys():
        de_nodes = 0  # 超负荷节点本身的故障点
        if gene in init_dmg_loads.keys():
            de_nodes = len(init_dmg_loads[gene])
        bs_nodes = len(gene_loads[gene])
        ov_nodes = len(overload_genes[gene])
        orgs[gene] = (ov_nodes-de_nodes)/bs_nodes
        overload_times[gene] = round(rho*np.exp(np.cos(orgs[gene])),2)
    return overload_genes, overload_times

def toDiGraph(G): # 转化为双向图
    # 创建一个空的有向图
    di_G = nx.DiGraph()
    for s, t, attr in G.edges(data=True):
        di_G.add_edge(s, t, **attr)
        di_G.add_edge(t, s, **attr)
    return di_G


def node_fix_cost(G,fix,dmg,trans_phi,fix_costs,type):
    # 电源节点修复所需的资源消耗和最大流动态评估
    zone = G.nodes[dmg]['zone']
    trans = trans_phi[zone]
    fixcost = (1+trans)*fix_costs[type]
    di_G = toDiGraph(G)
    maxflow = nx.maximum_flow_value(di_G, fix, dmg, capacity='afc')
    return fixcost,maxflow


def node_fix_cost1(G,fix,dmg,trans_phi,fix_costs,type):
    # 电源节点修复所需的资源消耗和最大流动态评估
    zone = G.nodes[dmg]['zone']
    trans = trans_phi[zone]
    fixcost = fix_costs[type]
    di_G = toDiGraph(G)
    maxflow = nx.maximum_flow_value(di_G, fix, dmg, capacity='afc')
    return fixcost,maxflow


def edge_fix_cost(G,fix,dmg,trans_phi,fix_costs):
    # 电线修复所需的资源消耗和最大流动态评估
    source,target  = dmg
    s_zone = G.nodes[source]['zone']
    t_zone = G.nodes[target]['zone']
    s_trans = trans_phi[s_zone]
    t_trans = trans_phi[t_zone]
    trans = 0.5*(s_trans + t_trans)
    di_G = toDiGraph(G)
    s_max_flow = nx.maximum_flow_value(di_G, fix, source, capacity='afc')
    t_max_flow = nx.maximum_flow_value(di_G, fix, target, capacity='afc')
    maxflow = 0.5*(s_max_flow+t_max_flow)
    fixcost = (1+trans)*fix_costs['edge']
    return fixcost,maxflow

def edge_fix_cost1(G,fix,dmg,trans_phi,fix_costs):
    # 电线修复所需的资源消耗和最大流动态评估
    source,target  = dmg
    di_G = toDiGraph(G)
    s_max_flow = nx.maximum_flow_value(di_G, fix, source, capacity='afc')
    t_max_flow = nx.maximum_flow_value(di_G, fix, target, capacity='afc')
    maxflow = 0.5*(s_max_flow+t_max_flow)
    fixcost = fix_costs['edge']
    return fixcost,maxflow


def flow_loss_path(G,bugloads,road_links,order):
    # 并根据顺坏的荷载节点计算影响的道路节点
    dmg_roads = {}
    init_dmg_loads = bugloads.copy()
    for gene in init_dmg_loads.keys():
        init_dmg_load  = init_dmg_loads[gene]
        dmg_road = []
        for load in init_dmg_load:
            sel_road = road_links[road_links['target']==load].values
            for road in sel_road:
                if 'S' not in road :
                    dmg_road.append((load,road[0]))
        if dmg_road:
            dmg_roads[gene] = dmg_road
    # 并根据顺坏的荷载节点计算影响的道路节点1-6阶
    nbs_roads = [dmg_roads]
    counts = []
    roads_list = []
    for i in range(order):
        if i > 0:
            nbs_roads.append({})
            nbs_roads[i], count = get_nbs_roads(G, nbs_roads[i-1])
            counts.append(count)
        roads_list.append([road for road_list in nbs_roads[i].values() for road in road_list])
    nbs_road_counts = [sum(count) for count in counts]
    return dmg_roads,nbs_roads,roads_list


def overload_computed(G,dmgs,gene_loads,rho):
    # 1.评估当前损坏状态下的实际存在故障的荷载点(考虑中介节点损坏造成的关联损坏效应)
    bugloads = damage_evluation(G,dmgs,gene_loads)
    # 2.评估当前损坏状态下的荷载节点超负荷运行情况
    inbackups, backups, bug_nodes = network_overload(G,bugloads,gene_loads)
    # 3.评估当前的超负荷运行状态和超负荷运行时间
    overload_genes, overload_times = overload_evluation(rho,inbackups,backups,bug_nodes,bugloads,gene_loads)

    return bugloads,overload_genes, overload_times


def max_cost_cal(G,dmg,fix,trans_phi,fix_costs):
    # 最大损失和最大流计算
    if 'G' in dmg: 
        fixcost,maxflow = node_fix_cost1(G,fix,dmg,trans_phi,fix_costs,'gene')
    elif isinstance(dmg, tuple):
        fixcost,maxflow = edge_fix_cost1(G,fix,dmg,trans_phi,fix_costs)
    else:
        fixcost,maxflow = node_fix_cost1(G,fix,dmg,trans_phi,fix_costs,'load')
    return fixcost,maxflow



def road_afc_loss(G,freq,loss_paths,flow_rate): # 动态路网的流量损失
    # 设置网络的损失速率
    loss_roads = loss_paths[freq] #获取当前时刻的受到影响的节点
    # 初始化损坏节点的损失效率
    loss_nums = {}
    loss_rate = {}
    if freq == 0:
        for source,target in loss_roads:
            loss_nums[source] = 1
            loss = flow_rate[G.nodes[source]['zone']]
            loss_rate[source] = loss
            G.nodes[source]['loss'] = loss
            # 计算当前异常道路的通行降低值
            ifc = G.edges[(source,target)]['ifc']
            afc = G.edges[(source,target)]['afc']
            if afc > ifc*(1-loss):
                G.edges[(source,target)]['afc'] = ifc*(1-loss)
    else:
        for source,target in loss_roads:
            if 'S' not in source:
                nbs_s = list(set(nx.neighbors(G,source)))
                nbs_n= sum([1 for item in nbs_s if 'G' not in item]) # 计算除去电源节点
                loss_nums[source] = nbs_n
                loss = flow_rate[G.nodes[source]['zone']]
                loss_rate[source] = loss/nbs_n 
                G.nodes[source]['loss'] = loss/nbs_n 
                # 计算当前异常道路的通行降低值
                ifc = G.edges[(source,target)]['ifc']
                afc = G.edges[(source,target)]['afc']
                if afc > ifc*(1-loss):
                    G.edges[(source,target)]['afc'] = ifc*(1-loss)


def flow_gain_path(G,fixloads,road_links,order):
    # 并根据损坏的荷载节点计算影响的道路节点
    fix_roads = {}
    for gene in fixloads.keys():
        init_fix_load  = fixloads[gene]
        fix_road = []
        for load in init_fix_load:
            sel_road = road_links[road_links['target']==load].values
            for road in sel_road:
                if 'S' not in road:
                    fix_road.append((load,road[0]))
        if fix_road:
            fix_roads[gene] = fix_road
    # 并根据顺坏的荷载节点计算影响的道路节点1-6阶
    nbs_roads = [fix_roads]
    counts = []
    roads_list = []
    for i in range(order):
        if i > 0:
            nbs_roads.append({})
            nbs_roads[i], count = get_nbs_roads(G, nbs_roads[i-1])
            counts.append(count)
        roads_list.append([road for road_list in nbs_roads[i].values() for road in road_list])
    return fix_roads,nbs_roads,roads_list


def road_afc_gain(G,dt,rho,fix_dmgs,aff,dmg_genes): # 动态的路网流量恢复关联dt
    # 判断存在损失的节点属于哪项修复过的dmg
    fixnodes = extract_elements(fix_dmgs)
    for node in fixnodes:
        if node in dmg_genes:
            # 修复电源节点下所有的节点
            sel_aff = aff[aff['generator node']==node]
            sel_road = sel_aff['road node']
            sel_load = sel_aff['load node']
            sel_nodes = list(set(sel_load).union(set(sel_road)))
            # 修复节点管辖的流量损失
            for source,target in G.edges():
                selected = [source in sel_nodes, target in sel_nodes]
                if any(selected):
                    ifc = G.edges[(source,target)]['ifc']
                    afc = G.edges[(source,target)]['afc']
                    if afc < ifc:
                        new_afc = afc*(1+dt*rho)
                        if new_afc < ifc:
                            G.edges[(source,target)]['afc'] = new_afc
                        else:
                            G.edges[(source,target)]['afc'] = ifc
        elif 'L' in node:
            # 修复荷载节点下所有的节点
            sel_aff = aff[aff['load node']==node]
            sel_load = sel_aff['road node']
            sel_nodes = list(sel_load)
            # 修复节点管辖的流量损失
            for source,target in G.edges():
                selected = [source in sel_nodes, target in sel_nodes]
                if any(selected):
                    ifc = G.edges[(source,target)]['ifc']
                    afc = G.edges[(source,target)]['afc']
                    if afc < ifc:
                        new_afc = afc*(1+dt*rho)
                        if new_afc < ifc:
                            G.edges[(source,target)]['afc'] = new_afc
                        else:
                            G.edges[(source,target)]['afc'] = ifc
    

def total_fix_gain(G,dt,rho):
    # 修复存在流量损失的边进行动态修复
    for source, target, attr in G.edges(data=True):
        ifc = G.edges[(source,target)]['ifc']
        afc = G.edges[(source,target)]['afc']
        if ifc-afc > 0: # 流量损失随时间衰减
            new_afc = afc*(1+dt*rho)
            if new_afc <ifc:
                G.edges[(source,target)]['afc'] = new_afc
            else:
                G.edges[(source,target)]['afc'] = ifc


def get_dmg_zone(G,dmg):
    # 获取当前损失节点所处分区
    if isinstance(dmg,tuple):
        zone = G.nodes[dmg[0]]['zone']
    else:
        zone = G.nodes[dmg]['zone']
    return zone


def dynamic_vis_net(G,dmgs,dmg_roads,freq,save_path,roads_path,save,bugloads,overload_genes,sector_dict):
    init_dmg_loads = bugloads.copy()
    dmg_genes,dmg_edges,_, = dmgs
    # 将G中所有的节点和道路都设置会灰色
    for edge in G.edges():
        G.edges[edge]['color'] = 'lightgray'
    for node in G.nodes():
        G.nodes[node]['color'] = 'gray'
    # 可视化电网顺坏与路网流量损失传递
    for gene in dmg_genes:
        G.nodes[gene]['color'] = 'blue'
    for gene in init_dmg_loads.keys():
        for load in init_dmg_loads[gene]:
            G.nodes[load]['color'] = 'blue'
    for gene in overload_genes.keys():
        overloads = overload_genes[gene]
        G.nodes[gene]['color'] = 'cyan' 
        for load in overloads:
            G.nodes[load]['color'] = 'green'
    for edge in dmg_edges:
        G.edges[edge]['color'] = 'blue'
        G.edges[edge]['wd'] = 3
        G.edges[edge]['ls'] = 'solid'
    for gene in dmg_roads.keys():
        edges = dmg_roads[gene]
        for edge in edges:
            G.edges[edge]['color'] = 'red'
            G.edges[edge]['wd'] = 5
            G.edges[edge]['ls'] = 'solid'
    # 绘制超负荷
    # 根据路网流量扩散对路网着色
    colors = ['#ff0100','#fe6603','#fd9802','#ffcc03','#ffff00','#ace601']
    for i in range(freq):
        wd = 5-i if i<5 else 1
        refresh_color(G,roads_path[i],colors[i],wd)
    # 设置可视化参数
    rs = {'mainten':(0,0.8),'generator':(1.2,1.8),'load':(2.2,3.8),'road':(4.2,7.8)}
    pos_sectors = assign_positions_by_sector(G,rs,sector_dict)

    # 提取节点和边的特性
    node_colors = [attr['color'] for node, attr in G.nodes(data=True)]
    node_sizes = [attr['ns'] for node, attr in G.nodes(data=True)]
    edge_colors = [attr['color'] for source,target,attr in G.edges(data=True)]
    styles = [attr['ls'] for source,target,attr in G.edges(data=True)]
    widths = [attr['wd'] for source,target,attr in G.edges(data=True)]
    # 绘制网络
    plt.figure(figsize=(10, 8))
    nx.draw(G, pos_sectors, node_color=node_colors, edge_color=edge_colors, with_labels=False,\
            node_size=node_sizes,alpha=0.6,style=styles,width=widths,font_size=6)
    plt.title('Multi-layer coupled networks visualization',fontsize=20)
    #plt.show()
    if save:
        plt.savefig(save_path, dpi=600)
    else:
        plt.show()
