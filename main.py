import pandas as pd
import numpy as np
import networkx as nx
from api import *
from info import *
import warnings
warnings.filterwarnings('ignore')
import multiprocessing as mp
from functools import partial

# 解析数据
generator_nodes = pd.read_csv('./static/data/mcn/generator nodes.csv')
power_links = pd.read_csv('./static/data/mcn/power links.csv')
prob_file = pd.read_csv('./static/data/mcn/probability of electrical grid .csv')
load_nodes = pd.read_csv('./static/data/mcn/load nodes.csv')
node_affiliation = pd.read_csv('./static/data/mcn/node affiliation.csv')
road_links = pd.read_csv('./static/data/mcn/road links.csv')
road_nodes = pd.read_csv('./static/data/mcn/road nodes.csv')
mainten_center = pd.read_csv('./static/data/mcn/maintenance center.csv')

# 电网基础隶属关系
aff = node_affiliation
genes = generator_nodes['generator node'].values
loads = load_nodes['load node'].values
roads = road_nodes['road node'].values
gene_loads = {}
for gene in genes:
    gene_load = aff[aff['generator node']==gene] 
    gene_loads[gene] = list(set(gene_load['load node'].values))
load_roads = {}
for load in loads:
    load_road = aff[aff['load node']==load]
    load_roads[load] = list(set(load_road['road node'].values))

# 创建多层耦合网络
mcn = nx.Graph()
#+++++++++++++++++添加网络点元素+++++++++++++++++++++
# 添加电源节点(绿色),并设置sector标签
aff = node_affiliation.copy()
sector_dict = {}
for index, row in generator_nodes.iterrows():
    prob = zone_probs[row['zone']]['generator']
    # 添加节点到耦合网络
    mcn.add_node(row['generator node'], type='generator', zone=row['zone'], dmg_prob = prob,state=row['state'],
        sector=index,color='Green',srv_popu=row['service population'],ns=gene_node_size,loss = 0,\
        baseload=row['basic load capacity'])
    sector_dict[row['generator node']] = index
# 添加荷载节点(紫色),并设置所属分区
for index, row in load_nodes.iterrows():
    gen_label = aff[aff['load node']==row['load node']]
    sector_index = sector_dict[gen_label['generator node'].values[0]]
    prob = zone_probs[row['zone']]['load']
    # 添加节点到耦合网络
    mcn.add_node(row['load node'], type='load', zone=row['zone'], state=row['state'], dmg_prob = prob,\
                sector=sector_index,color='purple',srv_popu=row['service population'],ns=load_node_size,
                loss = 0)
    
# 添加道路节点(橙色),并设置所属分区
for index, row in road_nodes.iterrows():
    gen_label = aff[aff['road node']==row['road node']]
    sector_index = sector_dict[gen_label['generator node'].values[0]]
    mcn.add_node(row['road node'], type='road', zone=row['zone'], state=row['state'],loss = 0,\
                sector=sector_index,color='orange',ns=road_node_size)
# 添加维修中心节点(黄色)
for index, row in mainten_center.iterrows():
    if row['source'] not in mcn.nodes:
        mcn.add_node(row['source'],type='mainten',state=row['state'],loss = 0,color='Yellow',\
                sector=0,ns=fixcenter_size)     

#+++++++++++++++++添加网络边元素+++++++++++++++++++++
# 在电源节点和负载节点之间增加主连接(蓝色)
for index, row in power_links.iterrows():
    if 'G' in row['source'] and 'L' in row['target']:
        # 添加连边到耦合网络
        mcn.add_edge(row['source'], row['target'], type='main',state='normal',\
            color='blue',ls='solid',wd=2,ifc=0,afc=0)

# 荷载节点之间的相互连接(浅蓝色)
for index, row in power_links.iterrows():
    if 'L' in row['source'] and 'L' in row['target']:
        # 添加连边到耦合网络
        mcn.add_edge(row['source'], row['target'], type='backup',state='normal',\
            color='cyan',ls='solid',wd=1,ifc=0,afc=0)
# 添加路网(浅灰色)
for index, row in road_links.iterrows():
    if (row['source'],row['target']) not in mcn.edges:
        mcn.add_edge(row['source'], row['target'], type='roadline',state=row['state'],color='lightgray',\
                ifc=row['initial flow capacity'],afc = row['initial flow capacity'],ls='dashed',wd=0.5)
                
# 添加维修中心到路网(黄色)
for index, row in mainten_center.iterrows():
    #if (row['source'],row['target']) not in mcn.edges:
    mcn.add_edge(row['source'], row['target'], type='roadline',state=row['state'],color='brown',\
            ifc=row['initial flow capacity'],afc = row['initial flow capacity'],ls='solid',wd=2)

# 动态添加耦合网络服务人数和顺坏概率
pwlks = power_links.copy()
for source, target, attr in mcn.edges(data=True):
    try:
        if 'R' not in source:
            s_prob = mcn.nodes[source]['dmg_prob']
        if 'R' not in target:
            t_prob = mcn.nodes[target]['dmg_prob']
        edge_prob = max([s_prob,t_prob])
    except:
        edge_prob = 0
    mcn.edges[(source, target)]['dmg_prob'] = edge_prob
    pwlk = pwlks[(pwlks['source']==source)&(pwlks['target']==target)]
    if pwlk.shape[0] > 0:
        mcn.edges[(target, source)]['srv_popu'] = pwlk['service population'].values[0]

# 创建一个电力网络（用于故障评估和中心性评估）
pw = nx.Graph()
# 添加电源节点(绿色),并设置sector标签
sector_dict = {}
for index, row in generator_nodes.iterrows():
    prob = zone_probs[row['zone']]['generator']
    # 添加节点到耦合网络
    pw.add_node(row['generator node'], type='generator', zone=row['zone'], dmg_prob = prob,state=row['state'],
        sector=index,color='Green',srv_popu=row['service population'],ns=gene_node_size,baseload=row['basic load capacity'])
    sector_dict[row['generator node']] = index
# 添加荷载节点(紫色),并设置所属分区
for index, row in load_nodes.iterrows():
    gen_label = aff[aff['load node']==row['load node']]
    sector_index = sector_dict[gen_label['generator node'].values[0]]
    prob = zone_probs[row['zone']]['load']
    # 添加节点到耦合网络
    pw.add_node(row['load node'], type='load', zone=row['zone'], state=row['state'], dmg_prob = prob,\
                sector=sector_index,color='purple',srv_popu=row['service population'],ns=load_node_size)

#+++++++++++++++++添加网络边元素+++++++++++++++++++++
# 在电源节点和负载节点之间增加主连接(蓝色)
for index, row in power_links.iterrows():
    if 'G' in row['source'] and 'L' in row['target']:
        # 添加连边到耦合网络
        pw.add_edge(row['source'], row['target'], type='main',state='normal',\
            color='blue',ls='solid',wd=2,ifc=0,afc=0)

# 荷载节点之间的相互连接(浅蓝色)
for index, row in power_links.iterrows():
    if 'L' in row['source'] and 'L' in row['target']:
        # 添加连边到耦合网络
        pw.add_edge(row['source'], row['target'], type='backup',state='normal',\
            color='cyan',ls='solid',wd=1,ifc=0,afc=0)

# 动态添加耦合网络服务人数和顺坏概率
pwlks = power_links.copy()
for source, target, attr in pw.edges(data=True):
    try:
        if 'R' not in source:
            s_prob = pw.nodes[source]['dmg_prob']
        if 'R' not in target:
            t_prob = pw.nodes[target]['dmg_prob']
        edge_prob = max([s_prob,t_prob])
    except:
        edge_prob = 0
    pw.edges[(source, target)]['dmg_prob'] = edge_prob
    pwlk = pwlks[(pwlks['source']==source)&(pwlks['target']==target)]
    if pwlk.shape[0] > 0:
        pw.edges[(target, source)]['srv_popu'] = pwlk['service population'].values[0]

# 创建一个道路网络（用于路网可达性评估与流量损失）
rw = nx.Graph()
# 道路网络包含到电力的连边
road_s = road_links['source'].values.tolist()
road_t = road_links['target'].values.tolist()
total_road_nodes = list(set(road_s+road_t))
# 添加道路节点(橙色),并设置所属分区
for node in total_road_nodes:
    if 'R' in node:
        row = road_nodes[road_nodes['road node']==node]
        rw.add_node(node, type='road', zone=row['zone'], state=row['state'],\
                        color='orange',ns=road_node_size)
    else:
        attr = pw.nodes[node]
        rw.add_node(node, type='power', zone=attr['zone'], state=attr['state'],\
                        color=attr['color'],ns=attr['ns'])
# 添加维修中心节点(黄色)
for index, row in mainten_center.iterrows():
    if row['source'] not in rw.nodes:
        rw.add_node(row['source'],type='mainten',state=row['state'],\
                color='Yellow',sector=0,ns=fixcenter_size)

#+++++++++++++++++添加网络边元素+++++++++++++++++++++
# 添加路网(浅灰色)
for index, row in road_links.iterrows():
    if (row['source'],row['target']) not in rw.edges:
        rw.add_edge(row['source'], row['target'], type='roadline',state=row['state'],color='lightgray',\
                ifc=row['initial flow capacity'],afc = row['initial flow capacity'],ls='dashed',wd=0.5)
# 添加维修中心到路网(黄色)
for index, row in mainten_center.iterrows():
    rw.add_edge(row['source'], row['target'], type='roadline',state=row['state'],color='brown',\
            ifc=row['initial flow capacity'],afc = row['initial flow capacity'],ls='solid',wd=2)

def main(i):
    try:
        # 根据各区域的故障概率产生故障
        #np.random.seed(i)  # 确保结果的可重现性
        dmg_genes,dmg_edges,dmg_loads = get_damges_by_zone(mcn,zone_probs)
        sel_edges = [edge for edge in dmg_edges if 'R' not in edge[1]]
        sel_lists = pd.Series(list(range(len(sel_edges))))
        bug_edges = [sel_edges[i] for i in sel_lists.sample(frac=1)]
        dmg_edges = bug_edges[0:np.random.randint(6,10)]
        dmg_loads = dmg_loads[0:np.random.randint(6,10)]
        all_dmgs = [dmg_genes, dmg_edges, dmg_loads]

        # ++++++++++++++++++++++整体模拟执行模块
        # 时间离散格式
        fix = 'S1'
        dt = 0.1
        pop = True
        Nt = 100
        steps = int(Nt/dt) # 迭代步数
        groups = 10 
        at = 3
        order = 8 # 最大的阶数
        threshold = 1000
        t_seq = [round(t, 1) for t in np.linspace(0,Nt,int(Nt/dt))]

        # 初始化一个用于运行的多层耦合网络副本
        mcn_r = mcn.copy()
        pw_r = pw.copy() 
        dmgs = all_dmgs.copy()
        # 初始通行总流量
        total_flows = sum([attr['ifc'] for _,_,attr in mcn_r.edges(data=True)])
        total_flows
        vuls = []
        fixprocess = 0
        overload_dict = {} # 定义一个超负荷状态控制字典
        dmgs_seq,_,_ = seq_dmgs_by_weight_occur(pw_r,t_seq,Nt,dt,at,groups,dmgs,pop=pop)
        acc_fixeds = [] # 计算每次迭代的
        # 定义路网损失传递频率
        freqs = list(range(order))
        fix_ranges = []
        fix_dmgs = []
        maxflows = []
        print('-'*20+'开始执行计算过程'+'-'*20)
        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format('STEP','故障','总消耗','已完成','最大流','脆弱性'))
        for step in range(steps):
            pw_r = pw.copy()
            # 设置时间产生的时间点
            # 1.当前开展修复的节点
            dmg = dmgs_seq[fixprocess] 
            # 2.当前修复节点所需的修复资源(资源+运输)
            if fixprocess not in fix_ranges: 
                fixcost,maxflow = max_cost_cal(mcn_r,dmg,fix,trans_phi,fix_costs)
                if maxflow<threshold:
                    maxflow = threshold + np.random.random()*100
            # 3.根据当前的电网损坏状态评估超负荷运行状态
            bugloads,overload_genes,overload_times = overload_computed(pw_r,dmgs,gene_loads)
            overload_index = sorted(overload_times) 
            for idx in overload_index:   # 初始化超负荷字典
                overload_dict[idx] = 0                   
            # 4.判断当前是否存在超负荷消失的情况
            for gene in overload_index:
                if overload_dict[gene] > overload_times[gene]: # 如果当前的超负荷运行时间
                    # 超负荷小时的荷载节点释放
                    print('>>>>>>>>>>>>>>>>>超负荷状态消失>>>>>>>>>>>>>>>>>>')
                    ov_dmgs = overload_genes[gene] 
                    pw_r = pw.copy()
                    dmgs[2] += ov_dmgs # 将新失效的节点加入到损坏列表中
                    # 重新评估当前路网的超负荷现状
                    bugloads,overload_genes,overload_times = overload_computed(pw_r,dmgs,gene_loads)
                    # 判定是否有新的超负荷电源产生
                    for ov_load in overload_genes.keys():
                        if ov_load not in overload_index:
                            overload_dict[ov_load] = 0

            # 5.根据实际的电网损失计算路网的流量损失传播路径
            _,_,loss_paths = flow_loss_path(mcn_r,bugloads,road_links,order) 
            freq = int(dt*step/2)  # 定义传播刷新的频率
            if freq  in freqs: 
                # 根据传播频率对节点进行流量损失传递计算
                road_afc_loss(mcn_r,freq,loss_paths,flow_rate)
                # 损失传递后的节点最大流动态刷新
                _,maxflow = max_cost_cal(mcn_r,dmg,fix,trans_phi,fix_costs)
                if maxflow < threshold:
                    maxflow = threshold + np.random.random()*100
            
            # 6.整个网络的脆弱性评估
            actul_flows = sum([attr['afc'] for _,_,attr in mcn_r.edges(data=True)])
            vul = (total_flows-actul_flows)/total_flows
            vuls.append(vul)
            
            # 7.修复判定的控制条件
            fixed = sum(acc_fixeds) # 当前损坏点的修复进度
            if fixed <= fixcost: # 损坏修复中
                for ov_load in overload_dict.keys(): # 动态关联超负荷时间
                    overload_dict[ov_load] += dt
                acc_fixeds.append(dt*maxflow) # 动态关联超修复进度
                road_afc_gain(mcn_r,dt,1,fix_dmgs,aff,dmg_genes)
            else:
                fix_ranges.append(fixprocess)
                fixprocess += 1
                # 判定是否计算结束
                if fixprocess > len(dmgs_seq)-1:
                    if vul>0:
                        total_fix_gain(mcn_r,dt,0.5)
                    break
                acc_fixeds = [] # 累计修复资源
                # 修复完成了，刷新dmg
                newdmgs = []
                for dmggs in dmgs:
                    ndmg = []
                    for dmgggs in dmggs:
                        if dmg != dmgggs:
                            ndmg.append(dmgggs)
                        else:
                            fix_dmgs.append(dmgggs)
                    newdmgs.append(ndmg)
                dmgs = newdmgs
                # 网络总体故障的
                total_fix_gain(mcn_r,dt,0.5)
            maxflows.append(maxflows)
            # if dt*step in list(range(Nt)):
            res = (step,dmg,fixcost,fixed,maxflow,vul)
            print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(*res))

        df_vuls = pd.DataFrame(vuls)
        df_vuls.to_csv('./static/res/sim_vuls_{}.csv'.format(i),index=False)
        df_dmgs = pd.DataFrame(all_dmgs)
        df_dmgs.to_csv('./static/res/sim_dmgs_{}.csv'.format(i),index=False)
        return vuls,all_dmgs
    except:
        pass


if __name__ == '__main__':

    # 运行多进程
    N = 1000
    pool = mp.Pool(processes=20)
    result = pool.map(main, range(N))
    res_vuls,res_dmgs = [], []
    for vul,dmg in result:
        res_vuls.append(vul)
        res_dmgs.append(dmg)
    df_vuls = pd.DataFrame(res_vuls)
    df_dmgs = pd.DataFrame(res_dmgs)
    df_vuls.to_csv('./static/res/df_vuls_{}.csv'.format(N),index=False)
    df_dmgs.to_csv('./static/res/df_dmgs_{}.csv'.format(N),index=False)

       

