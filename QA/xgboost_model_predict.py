import re
import os
import fcntl
import multiprocessing
import numpy as np
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.metrics import accuracy_score
import pandas as pd
import argparse
import time


def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_features_path', type=str, required=True,
                        help='input_features_path')
    parser.add_argument('-g', action='store', dest='graph', type=str, required=True, 
                        help='graph_id_path')
    parser.add_argument('-o', action='store', dest='output_path', type=str, required=True, help='output_path')
    args = parser.parse_args()
    return args


def fit_xgboost_model(labels_train, data_train, n_rounds=60, show_cv=True, max_size_cv=10000):
    if show_cv:
        my_sample = np.random.choice(len(labels_train), size=min(int(np.floor(len(labels_train)/2)), max_size_cv), replace=False)

        data_train_cv = data_train.iloc[my_sample]
        labels_train_cv = labels_train.iloc[my_sample]
        unique_sample = np.unique(my_sample)
        data_train_cv_test = data_train.drop(data_train.index[unique_sample])
        labels_train_cv_test = labels_train.drop(labels_train.index[unique_sample]) 
        
        cv_fit = XGBClassifier(n_estimators=n_rounds, objective='binary:logistic',max_depth=9, n_jobs=40, eta = 0.1)
        cv_fit.fit(data_train_cv, labels_train_cv)
        
        predictions = cv_fit.predict(data_train_cv_test)
        accuracy = accuracy_score(labels_train_cv_test, predictions)
        print('Prediction accuracy (CV):', accuracy)
        
    xgb_model = XGBClassifier(n_estimators=n_rounds, objective='binary:logistic', max_depth=9, n_jobs=40, eta = 0.1)
    xgb_model.fit(data_train, labels_train)  
    return xgb_model


def check_path(path_edge_list):
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
    n = len(data)
    mean = sum(data) / n 
    squared_diff = [(x - mean) ** 2 for x in data]
    variance = sum(squared_diff) / n
    return variance


def median(lst):
    sorted_lst = sorted(lst)
    n = len(sorted_lst)  
    if n % 2 == 0:
        median_index = n // 2
        median_value = (float(sorted_lst[median_index - 1]) + float(sorted_lst[median_index])) / 2
    else:
        median_index = n // 2
        median_value = float(sorted_lst[median_index])
    return median_value


def split_path(path_sequence):
    path_info_list = path_sequence.replace(" ", "")
    segments_with_space = path_info_list.split("->")[:-1]
    node_string=''.join(n for n in segments_with_space)
    p_edges = [segments_with_space[i].strip() + " -> " + segments_with_space[i + 1].strip() for i in range(len(segments_with_space) - 1)]

    return p_edges,node_string,segments_with_space


def get_weight(edges_set,edges):
    path_edge_weigths=[]
    for path_edge in edges_set:
        for edge in edges:
            if edge.startswith(path_edge):
                path_edge_weigths.append(float(edge.split()[-2]))
                break
    if len(edges_set)!=len(path_edge_weigths):
        return None  
    else:
        return path_edge_weigths

def calculate_average(data):
    total = sum(data)
    count = len(data)
    average = total / count
    return average


def min_max_average_var(weight_list):
    if weight_list!=[]:   
        max_feature=max(weight_list)
        min_feature=min(weight_list)
        average_feature = calculate_average(weight_list)
        var_feature = calculate_variance(weight_list)
    else:
        max_feature=0
        min_feature=0
        average_feature=0
        var_feature=0
    return  max_feature,min_feature,average_feature,var_feature    


def get_transcript_read_info_for_normalization(graphs):
    read_support_path_weight_inf=[]
    for match in graphs:
        lines = match.strip().split('\n')
        path_info_index = lines.index('Assembled_Paths:')

        path_info = lines[path_info_index + 1:]
        for sing_path in path_info:
            sing_path_info = sing_path.split('; ')[-1].split(' ')
            # print(sing_path_info,sing_path_info[-4])
            # p_read_info=sing_path_info[-4:]
            read_support_path_weight_inf.append(sing_path_info[-4])
            # print(len(read_support_path_weight_inf))
    read_support_path_weight_inf = np.array(read_support_path_weight_inf).astype(float)
    read_support_path_weight_sum = np.sum(read_support_path_weight_inf)
    count_non_zero = np.count_nonzero(read_support_path_weight_inf)
    print(np.size(read_support_path_weight_inf) ,read_support_path_weight_sum )
    return read_support_path_weight_sum, count_non_zero



def get_label(match):
    result=[]
    lines = match.strip().split('\n')
    edges_index = lines.index('Edges')
    nodes_index = lines.index('Nodes')
    lrp_index = lines.index('LRP')
    path_info_index = lines.index('Assembled_Paths:')

    edges = lines[edges_index + 1:nodes_index]
    nodes = lines[nodes_index + 1:lrp_index]
    lrp = lines[lrp_index + 1:path_info_index]
    path_info = lines[path_info_index + 1:]

    graph_node_list=[]
    graph_node_weight_list=[]
    graph_node_position_list=[]
    for i in nodes:
        i_info=i.split(":")[0]
        graph_node_weight=i.split(":")[1]
        graph_node_list.append(i_info.split(" ")[0])
        graph_node_weight_list.append(float(graph_node_weight))
        graph_node_position_list.append(i_info.split(" ")[1:])

    graph_node_weight_average=calculate_average(graph_node_weight_list) 

    read_node_list=[]
    read_weight_list=[]
    for read_info in lrp:
        read,read_weight=read_info.split(":")[0],read_info.split(":")[-1]
        read_node=read.replace(" ","")
        read_weight=float(read_weight.replace(" ",""))
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
            graph_edge_weight_list.append(float(edge_weight))
        graph_edge_weight_average=calculate_average(graph_edge_weight_list)
    else:
        return

    graph_edge_num_feature=len(edges)
    graph_node_num_feature=len(nodes)

    if len(path_info)==0:
        return
    else:
        graph_edge_list=[]
        for p in path_info:
            p = p.split('; ')[-1]
            e=split_path(p)[0]
            graph_edge_list=graph_edge_list+e

        for single_path_id_info in path_info:
            single_path_info = single_path_id_info.split('; ')[-1].split(' ')
            # print(single_path_info)
            p_other_info=single_path_info[-4:]
            single_path_info=single_path_info[0]
            single_path_edges,path_node_string,path_nodes = split_path(single_path_info)
            id = single_path_id_info.split('; ')[0]
            id = id.replace('"', '')
            if single_path_edges==[]:
                continue 
            one_path_features=[]

            path_edge_num_feature=len(single_path_edges)

            path_node_weight_list=[]
            path_nodes_position_list=[]
            branch_point_num=0
            for i in path_nodes:
                if graph_edge_nodes.count(i)>2:
                    branch_point_num=branch_point_num+1
                path_node_weight_list.append(graph_node_weight_list[int(i)])
                index_graph_node=graph_node_list.index(i)
                path_nodes_position_list=path_nodes_position_list+graph_node_position_list[index_graph_node]
            no_first_and_end=path_nodes_position_list[1:-1]

            path_node_adjacent_count=0
            for i in range(len(path_nodes)-1):
                if int(no_first_and_end[i*2])==int(no_first_and_end[i*2+1])-1:
                    path_node_adjacent_count=path_node_adjacent_count+1
            if path_node_adjacent_count==len(path_nodes)-1:
                continue

            single_path_edges_weights=get_weight(single_path_edges,edges)
            if single_path_edges_weights:
                path_edge_weigth_max,path_edge_weigth_min,path_edge_weigth_average,path_edge_weigth_var=min_max_average_var(single_path_edges_weights)
                path_node_weigth_max,path_node_weigth_min,path_node_weigth_average,path_node_weigth_var=min_max_average_var(path_node_weight_list)

                share_edges=[]
                for e in single_path_edges:
                    if graph_edge_list.count(e)>1:
                        share_edges.append(e)
                nonshare_edges = [x for x in single_path_edges if x not in share_edges]
                share_edges_weights=get_weight(share_edges,edges)
                nonshare_edges_weights=get_weight(nonshare_edges,edges)

                share_weight_max,share_weight_min,share_weight_average,share_weight_var=min_max_average_var(share_edges_weights)
                nonshare_weight_max,nonshare_weight_min,nonshare_weight_average,nonshare_weight_var=min_max_average_var(nonshare_edges_weights)
                
                read_support_path_weight=0
                for i in range(len(read_node_list)):
                    if read_node_list[i] in path_node_string:
                        read_support_path_weight=read_support_path_weight+read_weight_list[i]

                one_path_features.append(path_edge_num_feature)
                one_path_features.append(path_edge_weigth_max)
                one_path_features.append(path_edge_weigth_min)
                one_path_features.append(path_edge_weigth_average)
                one_path_features.append(path_edge_weigth_var)
                
                one_path_features.append(len(share_edges))
                one_path_features.append(share_weight_max)
                one_path_features.append(share_weight_min)
                one_path_features.append(share_weight_average)
                one_path_features.append(share_weight_var)

                one_path_features.append(len(nonshare_edges))
                one_path_features.append(nonshare_weight_max)
                one_path_features.append(nonshare_weight_min)
                one_path_features.append(nonshare_weight_average)
                one_path_features.append(nonshare_weight_var)

                one_path_features.append(read_support_path_weight)

                one_path_features.append(graph_edge_weight_average)
                one_path_features.append(graph_node_weight_average)

                one_path_features.append(path_node_weigth_max)
                one_path_features.append(path_node_weigth_min)
                one_path_features.append(path_node_weigth_average)
                one_path_features.append(path_node_weigth_var)

                one_path_features.append(branch_point_num)
                # print("p_other_info",p_other_info)
                # print("p_other_info",p_other_info[0])
                # one_path_features.append((float(p_other_info[0])*count_non_zero)/read_support_path_weight_sum)
                        
                one_path_features.append(float(p_other_info[0]))
                one_path_features.append(float(p_other_info[1])) #yuting
                one_path_features.append(float(p_other_info[2])) #yuting

                one_path_features.append(float(p_other_info[-1]))
                #if p_other_info[0] == "0":
                #    one_path_features.append(0)
                #else:
                #    one_path_features.append(1)
                # print(one_path_features[-3:])





                result.append((id,one_path_features))
    return result
 

if __name__ == '__main__':

    args = initialization_parameters()

    graph_id_path = args.graph
    input_features_path = args.input_features_path
    output_path = args.output_path

    correct=0
    error=0
    data = pd.read_csv(input_features_path)

    data_train = data.iloc[:, :-1]
    labels_train = data.iloc[:, -1]
    xgb_model = fit_xgboost_model(labels_train, data_train, show_cv=False, max_size_cv=10000)

    with open(graph_id_path, "r") as file:
        data_raw = file.read()
    pattern = r"# Graph.*?Assembled_Paths:.*?(?=# Graph|\Z)"
    data = re.findall(pattern, data_raw, re.DOTALL)
    # read_support_path_weight_sum, count_non_zero = get_transcript_read_info_for_normalization(data)
    

    pool = multiprocessing.Pool(processes=10)


    results = pool.map(get_label,data)
    pool.close()
    pool.join()

    final_result = []
    for res in results:
        final_result.append(res)
    
    id_list=[]
    feature_list=[]
    with open(output_path +'/TransGram.score', 'w') as file_pro:
        for f_list in final_result:
            if f_list != None:
                for f in f_list:
                    feature_list.append(f[-1])
                    id_list.append(f[0])

        feature_names = [
                'path_edge_number',
                'path_edge_weigth_max', 'path_edge_weigth_min', 'path_edge_weigth_average',
                'path_edge_weigth_var', 'share_edge_number', 'share_weight_max',
                'share_weight_min', 'share_weight_average', 'share_weight_var',
                'nonshare_path_num', 'nonshare_weight_max', 'nonshare_weight_min',
                'nonshare_weight_average', 'nonshare_weight_var','read_support_path_weight',
                'graph_edge_weight_average','graph_node_weight_average','path_node_weigth_max',
                'path_node_weigth_min','path_node_weigth_average','path_node_weigth_var','branch_point_num',
                'raw_read_support_path_num','unreabile_edge_num_in_splicing_graph','edge_num_in_splicing_graph',
                'unreabile_edge_ratio_in_splicing_graph'
            ]


        one_path_features_df = pd.DataFrame(feature_list, columns=feature_names)
        one_path_predict_label_pro = xgb_model.predict_proba(one_path_features_df)
        for i in range(len(id_list)):
            file_pro.write(str(id_list[i])+":"+str(one_path_predict_label_pro[i][1])+"\n")
