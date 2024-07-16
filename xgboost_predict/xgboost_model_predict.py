import re
import os
import numpy as np
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score
import pandas as pd
import argparse


def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_features_path', type=str, required=True,
                        help='input_features_path')
    parser.add_argument('-g', action='store', dest='graph', type=str, required=True, 
                        help='graph_id_path')
    # parser.add_argument('-r', action='store', dest='refmap', type=str, required=True, help='refmap_path')
    # parser.add_argument('-gtf', action='store', dest='gtf_file', type=str, required=True,
    #                     help='gtf_file')
    # parser.add_argument('-num', action='store', dest='ground_truth_num', type=int, required=True,
    #                     help='ground_truth_num,5M_2:46647,mouse:141156,human:190850')
    args = parser.parse_args()
    return args


def fit_xgboost_model(labels_train, data_train, n_rounds=60, show_cv=True, max_size_cv=10000):
    if show_cv:
        my_sample = np.random.choice(len(labels_train), size=min(int(np.floor(len(labels_train)/2)), max_size_cv), replace=False)

        data_train_cv = data_train.iloc[my_sample]  # 使用 iloc 按位置选择行
        # print("data_train_cv"+"\n",data_train_cv)
        labels_train_cv = labels_train.iloc[my_sample]  # 使用 iloc 按位置选择行
        unique_sample = np.unique(my_sample)  # 获取唯一的索引值
        data_train_cv_test = data_train.drop(data_train.index[unique_sample])  # 删除指定索引的行
        # print("data_train_cv_test"+"\n",data_train_cv_test)
        labels_train_cv_test = labels_train.drop(labels_train.index[unique_sample])  # 删除指定索引的行
        
        cv_fit = XGBClassifier(n_estimators=n_rounds, objective='binary:logistic',max_depth=9, n_jobs=40, eta = 0.1)
        cv_fit.fit(data_train_cv, labels_train_cv)
        
        predictions = cv_fit.predict(data_train_cv_test)
        accuracy = accuracy_score(labels_train_cv_test, predictions)
        print('Prediction accuracy (CV):', accuracy)
        
    # 使用全部训练数据训练模型
    xgb_model = XGBClassifier(n_estimators=n_rounds, objective='binary:logistic', max_depth=9, n_jobs=40, eta = 0.1)
    xgb_model.fit(data_train, labels_train)  
    return xgb_model


def check_path(path_edge_list):
    # print("path_edge_list,edges",path_edge_list,edges)
    error_count=0
    for spe in path_edge_list:
        if spe not in edges:
            error_count=error_count+1
    print("error_count",error_count)
    if error_count==0:
        return True
    else:
        return False    


def calculate_variance(data):
    # 计算数据集的均值
    n = len(data)
    mean = sum(data) / n 
    # 计算每个数据点与均值的差的平方
    squared_diff = [(x - mean) ** 2 for x in data]
    variance = sum(squared_diff) / n
    return variance


def median(lst):
    # 对列表进行排序
    sorted_lst = sorted(lst)
    n = len(sorted_lst)  
    # 计算中位数值
    if n % 2 == 0:  # 如果列表长度为偶数
        median_index = n // 2
        median_value = (float(sorted_lst[median_index - 1]) + float(sorted_lst[median_index])) / 2
    else:  # 如果列表长度为奇数
        median_index = n // 2
        median_value = float(sorted_lst[median_index])
    return median_value


def split_path(path_sequence):
    # 分割路径字符串并去除首尾空格
    path_info_list = path_sequence.replace(" ", "")
    # path_label = path_info_list.split("->")[-1]
    # print("path_label",path_label)
    segments_with_space = path_info_list.split("->")[:-1]
    node_string=''.join(n for n in segments_with_space)
    # print("node_string",node_string)
    # print("segments_with_space",segments_with_space)
    # 去除空白字符并连接边
    p_edges = [segments_with_space[i].strip() + " -> " + segments_with_space[i + 1].strip() for i in range(len(segments_with_space) - 1)]
    # print("path_sequence,path_edges",p_edges)
    return p_edges,node_string,segments_with_space


# 如果返回值为None,那么这条路是错误的
def get_weight(edges_set):
    path_edge_weigths=[]
    # pattern_to_match = r"^0 -> 1(?:：\d+)?$"
    # print("edges_set",edges_set,type(edges_set[0]))
    for path_edge in edges_set:
        # print("path_edge",path_edge)
        for edge in edges:
            if edge.startswith(path_edge):
                path_edge_weigths.append(float(edge.split()[-2]))
                break
    # print("len(edges_set),len(path_edge_weigths)",len(edges_set),len(path_edge_weigths))
    if len(edges_set)!=len(path_edge_weigths):
        # print("edges",edges,"\n","edges_set","\n",edges_set,"\n")
        # print("error")
        return None  
    else:
        return path_edge_weigths

def calculate_average(data):
    total = sum(data)
    count = len(data)
    average = total / count
    return average


def min_max_median_var(weight_list):
    if weight_list!=[]:
    # print(path_edge_weigths)    
        max_feature=max(weight_list)
        min_feature=min(weight_list)
        median_feature = calculate_average(weight_list)
        # print("median_feature",median_feature)
        # var_feature = np.var(weight_list)
        var_feature = calculate_variance(weight_list)
        # print("var_feature",var_feature)
    else:
        max_feature=0
        min_feature=0
        median_feature=0
        var_feature=0
    return  max_feature,min_feature,median_feature,var_feature    


if __name__ == '__main__':

    args = initialization_parameters()

    graph_id_path = args.graph
    # refmap_path = args.refmap
    input_features_path = args.input_features_path
    # gtf_file = args.gtf_file
    # ground_truth_num = args.ground_truth_num
    
    # print(gtf_file)
    correct=0
    error=0
    data = pd.read_csv(input_features_path)
    # 使用全部数据进行xgboost的模型训练
    # 分离特征和标签
    data_train = data.iloc[:, :-1]
    labels_train = data.iloc[:, -1]
    # 使用 fit_xgboost_model 函数训练模型
    xgb_model = fit_xgboost_model(labels_train, data_train, show_cv=False, max_size_cv=10000)


    # 预测
    # 读取文件
    with open(graph_id_path, "r") as file:
        data = file.read()
    # 利用正则表达式匹配每幅图中的边和点信息（确定一下提取图的个数是否正确？）
    pattern = r"# Graph.*?Assembled_Paths:.*?(?=# Graph|\Z)"
    matches = re.findall(pattern, data, re.DOTALL)
    # print("matches",matches)

    # error_path=0
    # error_id=[]
    # correct_path=0
    # correct_id=[]
    # all_path=0
    # muti_exon_all_id=[]
    # predict_positive_path_ids=[]
    # ids=[]
    # all_ids=[]
    # 处理每个图（match中以字符串的形式存储了每个图的Edges,Nodes,LRP,PathInformation信息）


    # # 存储含有 "=" 的行的列表
    # lines_with_equal = []
    # assembled_correct_ids=[]
    # with open(refmap_path, 'r') as file:
    #     for line in file:
    #         # 去除行尾的换行符
    #         line = line.strip()
    #         # 检查行中是否包含 "="
    #         if '=' in line:
    #             # print(line)
    #             # print(line.split("|")[-1])
    #             line_list=line.split(",")
    #             for i in line_list:
    #                 assembled_correct_ids.append(i.split("|")[-1])
    #                 lines_with_equal.append(i)

    # tong_list=[]
    # all_list=[]
    # predict_pos=[]
    # threshold_list = [0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9]
    # # threshold_list = [0.3,0.325,0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,0.975,1.1]
    # # threshold_list=[0.35]
    # for threshold in threshold_list:
        # error_path=0
        # error_id=[]
        # correct_path=0
        # correct_id=[]
        # all_path=0
        # muti_exon_all_id=[]
        # predict_positive_path_ids=[]
        # ids=[]
        # all_ids=[]
    with open('./xgboost_predict_pro.txt', 'w') as file_pro:
        for match in matches:
            lines = match.strip().split('\n')
            # print("lines",lines[:2])

            # 找到索引位置
            edges_index = lines.index('Edges')
            nodes_index = lines.index('Nodes')
            lrp_index = lines.index('LRP')
            path_info_index = lines.index('Assembled_Paths:')

            # 分隔列表
            edges = lines[edges_index + 1:nodes_index]
            nodes = lines[nodes_index + 1:lrp_index]
            lrp = lines[lrp_index + 1:path_info_index]
            path_info = lines[path_info_index + 1:]
            # print("edges",edges)
            # print("lrp",lrp)
            # print("len(path_info)",len(path_info))
            # print("node",nodes)
            
            # 获取每个match所对应图中的nodes信息
            graph_node_list=[]
            graph_node_weight_list=[]
            graph_node_position_list=[]
            for i in nodes:
                i_info=i.split(":")[0]
                graph_node_weight=i.split(":")[1]
                # print("node_weight",graph_node_weight)
                # print("i_info",i_info)
                graph_node_list.append(i_info.split(" ")[0])
                graph_node_weight_list.append(float(graph_node_weight))
                graph_node_position_list.append(i_info.split(" ")[1:])

            graph_node_weight_var=calculate_average(graph_node_weight_list) 

            # print("graph_node",graph_node,graph_node_position)
            # print("graph_node_list,graph_node_position_list",graph_node_list,graph_node_position_list)
            # 获取每个match对应图中read及其weight信息
            read_node_list=[]
            read_weight_list=[]
            for read_info in lrp:
                # print("read_info",read_info)
                read,read_weight=read_info.split(":")[0],read_info.split(":")[-1]
                # print("read,read_weight",read,read_weight)
                read_node=read.replace(" ","")
                read_weight=float(read_weight.replace(" ",""))
                # print("read_node,read_weight",read_node,read_weight)
                # 去除长度为1的read
                if len(read_node)>=1:
                    read_node_list.append(read_node)
                    read_weight_list.append(read_weight)


            graph_edge_weight_list=[]
            graph_edge_nodes=[]
            if len(edges)!=0:
                for i in edges:
                    graph_edge_nodes.append(i.split(" ")[0])
                    graph_edge_nodes.append(i.split(" ")[2])
                    edge_weight=i.split(" ")[-2]
                    # print("edge_weight",edge_weight)
                    graph_edge_weight_list.append(float(edge_weight))
                graph_edge_weight_var=calculate_average(graph_edge_weight_list)
            else:
                # print("match",match)
                continue

            graph_edge_num_feature=len(edges)
            graph_node_num_feature=len(nodes)

            # print("path_info",len(path_info))
            if len(path_info)==0:
                continue
            else:
                graph_edge_list=[]
                for p in path_info:
                    p = p.split('; ')[-1]
                    e=split_path(p)[0]
                    # print(e)
                    graph_edge_list=graph_edge_list+e
                # print("graph_edge_list",graph_edge_list)
                
                # python对空列表不会执行循环，所以直接跳过了长度为0的路。是不是不存在长度为0的路
                for single_path_id_info in path_info:
                    # all_path=all_path+1

                    single_path_info = single_path_id_info.split('; ')[-1]
                    single_path_edges,path_node_string,path_nodes = split_path(single_path_info)
                    id = single_path_id_info.split('; ')[0]
                    id = id.replace('"', '')
                    # all_ids.append(id)
                    if single_path_edges==[]:
                        # print(single_path_id_info)
                        continue  
                    # ids.append(id)
                    # print("id",id)
                    # if check_path(single_path_edges):
                    # print("a")
                    one_path_features=[]

                    path_edge_num_feature=len(single_path_edges)

                    # 判断这条path中是否所有node都是首尾相接
                    # 获取path中所有node对应的position
                    # print("path_nodes",path_nodes)
                    path_node_weight_list=[]
                    path_nodes_position_list=[]
                    branch_point_num=0
                    for i in path_nodes:
                        if graph_edge_nodes.count(i)>1:
                            branch_point_num=branch_point_num+1
                        path_node_weight_list.append(graph_node_weight_list[int(i)])
                        index_graph_node=graph_node_list.index(i)
                        path_nodes_position_list=path_nodes_position_list+graph_node_position_list[index_graph_node]
                    no_first_and_end=path_nodes_position_list[1:-1]
                    # print("branch_point_num",branch_point_num)
                    # print("no_first_and_end",no_first_and_end)
                    # print("len(no_first_and_end)",len(no_first_and_end))
                    # print("len(path_nodes)*2-1",(len(path_nodes)-1)*2-1)
                    path_node_adjacent_count=0
                    for i in range(len(path_nodes)-1):
                        # print("i*2-1",i*2)
                        if int(no_first_and_end[i*2])==int(no_first_and_end[i*2+1])-1:
                            path_node_adjacent_count=path_node_adjacent_count+1
                    if path_node_adjacent_count==len(path_nodes)-1:
                        continue
                    # ids.append(id)

                    single_path_edges_weights=get_weight(single_path_edges)
                    if single_path_edges_weights:
                        # ids.append(id)
                        path_edge_weigth_max,path_edge_weigth_min,path_edge_weigth_median,path_edge_weigth_var=min_max_median_var(single_path_edges_weights)
                        path_node_weigth_max,path_node_weigth_min,path_node_weigth_median,path_node_weigth_var=min_max_median_var(path_node_weight_list)

                        share_edges=[]
                        for e in single_path_edges:
                            if graph_edge_list.count(e)>1:
                                share_edges.append(e)
                        nonshare_edges = [x for x in single_path_edges if x not in share_edges]
                        # print("share_edges",len(share_edges))        
                        
                        # 获得one_path_info中share边及其权重的信息
                        share_edges_weights=get_weight(share_edges)
                        nonshare_edges_weights=get_weight(nonshare_edges)
                        # print("match",match)
                        # print("share_edges",share_edges)
                        # print("share_edge_weights",share_edge_weights)
                        share_weight_max,share_weight_min,share_weight_medium,share_weight_var=min_max_median_var(share_edges_weights)
                        nonshare_weight_max,nonshare_weight_min,nonshare_weight_medium,nonshare_weight_var=min_max_median_var(nonshare_edges_weights)
                        
                        # 获取read对path的支持情况
                        read_support_path_weight=0
                        for i in range(len(read_node_list)):
                            if read_node_list[i] in path_node_string:
                            # if path_node_string in read_node_list[i]:
                                # print("path_node_string,read_node_list[i]",path_node_string,read_node_list[i])
                                read_support_path_weight=read_support_path_weight+read_weight_list[i]
                        # print("read_support_path_weight",read_support_path_weight)        
                                
                        # one_path_features.append(graph_edge_num_feature)
                        # one_path_features.append(graph_node_num_feature)
                        one_path_features.append(path_edge_num_feature)
                        

                        one_path_features.append(path_edge_weigth_max)
                        one_path_features.append(path_edge_weigth_min)
                        one_path_features.append(path_edge_weigth_median)
                        one_path_features.append(path_edge_weigth_var)
                        
                        one_path_features.append(len(share_edges))
                        one_path_features.append(share_weight_max)
                        one_path_features.append(share_weight_min)
                        one_path_features.append(share_weight_medium)
                        one_path_features.append(share_weight_var)

                        one_path_features.append(len(nonshare_edges))
                        one_path_features.append(nonshare_weight_max)
                        one_path_features.append(nonshare_weight_min)
                        one_path_features.append(nonshare_weight_medium)
                        one_path_features.append(nonshare_weight_var)

                        one_path_features.append(read_support_path_weight)

                        one_path_features.append(graph_edge_weight_var)
                        one_path_features.append(graph_node_weight_var)

                        one_path_features.append(path_node_weigth_max)
                        one_path_features.append(path_node_weigth_min)
                        one_path_features.append(path_node_weigth_median)
                        one_path_features.append(path_node_weigth_var)

                        one_path_features.append(branch_point_num)
                        
                    # 列名定义
                        feature_names = [
                            'path_edge_number',
                            'path_edge_weigth_max', 'path_edge_weigth_min', 'path_edge_weigth_median',
                            'path_edge_weigth_var', 'share_edge_number', 'share_weight_max',
                            'share_weight_min', 'share_weight_medium', 'share_weight_var',
                            'nonshare_path_num', 'nonshare_weight_max', 'nonshare_weight_min',
                            'nonshare_weight_medium', 'nonshare_weight_var','read_support_path_weight',
                            'graph_edge_weight_var','graph_node_weight_var','path_node_weigth_max',
                            'path_node_weigth_min','path_node_weigth_median','path_node_weigth_var','branch_point_num'
                        ]

                        # 将数据列表转换成DataFrame
                        # one_path_features_df = pd.DataFrame([one_path_features])
                        one_path_features_dict = dict(zip(feature_names, one_path_features))
                        # print(one_path_features)
                        one_path_features_df = pd.DataFrame([one_path_features_dict])
                        # print("one_path_features_df")
                        # print(one_path_features_df)
                        

                        # 用xgboost进行预测
                        # one_path_predict_label = xgb_model.predict(one_path_features_df)
                        # 在二分类情况下，你会得到两列，分别对应样本属于类别 0 和类别 1 的概率
                        one_path_predict_label_pro = xgb_model.predict_proba(one_path_features_df)
                        # print(id,one_path_predict_label_pro,one_path_predict_label_pro[0][1])
                        file_pro.write(str(id)+":"+str(one_path_predict_label_pro[0][1])+"\n")
    # threshold_list = [0.35,0.375,0.4,0.425,0.45,0.475,0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9]
