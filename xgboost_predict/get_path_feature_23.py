import re
import numpy as np
import pandas as pd
import argparse


def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', action='store', dest='graph', type=str, required=True, 
                        help='graph_id_path')
    parser.add_argument('-r', action='store', dest='refmap', type=str, required=True, help='refmap_path')
    parser.add_argument('-o', action='store', dest='features_output_path', type=str, required=True,
                        help='features_output_path')
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


def calculate_average(data):
    total = sum(data)
    count = len(data)
    average = total / count
    return average


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
    refmap_path = args.refmap
    features_output_path = args.features_output_path


    with open(graph_id_path, "r") as file:
        data = file.read()
    # 利用正则表达式匹配每幅图中的边和点信息（确定一下提取图的个数是否正确？）
    pattern = r"# Graph.*?Assembled_Paths:.*?(?=# Graph|\Z)"
    matches = re.findall(pattern, data, re.DOTALL)

    # 处理每个图（match中以字符串的形式存储了每个图的Edges,Nodes,LRP,PathInformation信息）


    # 存储含有 "=" 的行的列表
    lines_with_equal = []
    assembled_correct_ids=[]
    with open(refmap_path, 'r') as file:
        for line in file:
            # 去除行尾的换行符
            line = line.strip()
            # 检查行中是否包含 "="
            if '=' in line:
                # print(line)
                # print(line.split("|")[-1])
                line_list=line.split(",")
                for i in line_list:
                    assembled_correct_ids.append(i.split("|")[-1])
                    lines_with_equal.append(i)


    with open(features_output_path+"/23_features_stringtie3.csv", 'w') as f:
        # 向文件中写入数据
        f.write('path_edge_number,path_edge_weigth_max,path_edge_weigth_min,path_edge_weigth_median,path_edge_weigth_var,share_edge_number,share_weight_max,share_weight_min,share_weight_medium,share_weight_var,nonshare_path_num,nonshare_weight_max,nonshare_weight_min,nonshare_weight_medium,nonshare_weight_var,read_support_path_weight,graph_edge_weight_var,graph_node_weight_var,path_node_weigth_max,path_node_weigth_min,path_node_weigth_median,path_node_weigth_var,branch_point_num,path_label'+"\n")
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

                    one_path_features=[]

                    path_edge_num_feature=len(single_path_edges)

                    # 判断这条path中是否所有node都是首尾相接
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

                        
                        for i in one_path_features:
                            f.write(str(i)+",")     
                        
                        if id in assembled_correct_ids:
                            f.write("1")
                        else:
                            f.write("0") 
                        f.write("\n")
                        # print("one_path_features",len(one_path_features), one_path_features)   
