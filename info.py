gene_node_size = 240
load_node_size = 120
road_node_size = 80
fixcenter_size = 400
# 电网故障概率
zone_probs = {'Z1': {'generator': 0.0340, 'load': 0.0926},
              'Z2': {'generator': 0.0392, 'load': 0.1231},
              'Z3': {'generator': 0.0327, 'load': 0.1054}}

# 电网故障概率（平均）
zone_probs_ave = {'Z1': {'generator': 0.0349, 'load': 0.1077},
                  'Z2': {'generator': 0.0349, 'load': 0.1077},
                  'Z3': {'generator': 0.0349, 'load': 0.1077}}

# 初始化分区路网流量降低系数，维修资源消耗和维修运输时间
flow_rate = {'Z1': 0.2293,'Z2': 0.3029,'Z3': 0.3590}
flow_rate_ave = {'Z1': 0.3153,'Z2': 0.3153,'Z3': 0.3153}

trans_phi = {'Z1': 0.2759,'Z2': 0.3679,'Z3': 0.2453}
trans_phi_ave = {'Z1': 0.2759,'Z2': 0.3679,'Z3': 0.2453}

fix_costs = {'gene': 8000,'edge': 6000, 'load': 4000}
