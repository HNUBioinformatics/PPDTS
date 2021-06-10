import numpy
import numpy as np
import pandas as pd
import scipy
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.model_selection import KFold
from pylab import *
import random
random.seed(8)#Fixed randomly selected seeds

yuanshi_lianxi='../data/zhang_cir-disease/cicRNA-disease-association.xlsx'
test_yuanshi='../data/test2.xlsx'
lnc_pro_path='../data/lncBOOK/lncBOOK_disease_association.xlsx'
mator_d_p_01='../data/MATADOR/chem_protein_01_1.xlsx'
d_p_1_01='../data/d_p_1_01.xlsx'



''' 5 fold cross validation'''
def get_train_set(interaction_matrix):
    cv = 5
    (rows_index, cols_index) = np.where(interaction_matrix == 1)
    link_num = np.sum(np.sum(interaction_matrix)) #Calculate the number of 1
    random_index=random.sample(range(0,link_num),link_num) #Randomly generated, link_ Random number num without repetition
    size_of_cv=round(link_num//cv) #The number of substitutions for one in each cross validation
    print('cv=%d,size_of_cv=%d'%(cv,size_of_cv))

    result=[]
    for k in range(cv):
        print('The k-th cross validation',k+1)
        if (k!=cv):
            test_row_index = rows_index[random_index[size_of_cv*k:size_of_cv*(k+1)]] #[0:size_of_cv],[size_of_cv:size_of_cv+size_of_cv].....
            test_col_index = cols_index[random_index[size_of_cv*k:size_of_cv*(k+1)]]
        else:
            test_row_index=rows_index[random_index[size_of_cv*k:]]
            test_col_index=cols_index[random_index[size_of_cv*k:]]

        train_set=np.copy(interaction_matrix)
        for i in range(test_row_index.shape[0]):
            train_set[test_row_index[i],test_col_index[i]]=0


        predict_matrix=cf_model(train_set)
        test_index = np.where(train_set == 0)  # Each training matrix element is equal to the element position corresponding to 0
        real_score=interaction_matrix[test_index]
        predict_score=predict_matrix[test_index]

        num_auc=l5_model_evaluate(real_score,predict_score)
        result.append(num_auc)

    print("The last quintic AUC set：",result)
    print("Five times average AUC：\n",np.sum(np.array(result),axis=0)/cv)

''' cos similarity'''
def cos_similarity(inteartion_matrix):
    m1 = np.mat(inteartion_matrix)  # Array to matrix
    cos_similarity = cosine_similarity(m1)  # Calculate the cosine similarity between the rows of the matrix
    return cos_similarity


'''Calculation of drug-target similarity in PPDTS model'''
def ra_similarity(x):

    '''RA part1'''
    result_1 = x / sum(x)  
    y = np.array(result_1)   

    '''RA part2.1'''
    result_2 = []
    y = np.squeeze(y)  # Remove the first dimension
    for i in range(y.shape[0]):
        if x[i, :].sum() == 0:  # If the column of the test set is empty, the division error is avoided
            continue
        result_2.append(y[i, :] / np.count_nonzero(y[i, :]))
    z = np.array(result_2)

    '''RA part2.2'''
    z_result = np.zeros((z.shape[1], z.shape[1]))
    for i in range(z.shape[1]):  # !!!!
        ii = np.where(z[:, i] > 0)
        z_result[i] = z[list(ii)[0], :].sum(axis=0)  # list(ii)[0] Index with non-zero per column 
    return z_result

'''Calculation of drug-target similarity score matrix in PPDTS model'''
def cf_model(inteartion_matrix): # inteartion_matrix 01 matrix；
    simi_maxrix=cos_similarity(inteartion_matrix)

    row_sum_matrix=np.sum(simi_maxrix,axis=1)
    sum_diagonal_matrix=np.linalg.pinv(np.diag(row_sum_matrix))
    row_normalized_protein_similarity_matrix=np.dot(sum_diagonal_matrix,simi_maxrix)
    score_matrix=np.dot(row_normalized_protein_similarity_matrix,inteartion_matrix)
    return score_matrix

''' Baseline methods- RWR model'''
def rwr_model(inteartion_matrix,c):
    cos_maxrix = cos_similarity(inteartion_matrix)
    cos_maxrix=np.array(cos_maxrix)
    num=cos_maxrix.shape[1]
    for i in range(num):
        cos_maxrix[i,i]=0
    column_sum_martrix=np.sum(cos_maxrix,axis=0)
    sum_diagonal_matrix = np.linalg.pinv(np.diag(column_sum_martrix))
    transformmation_matrix=np.dot(cos_maxrix,sum_diagonal_matrix)

    # expand_dims(,axis=1)
    row_sum_interation_matrix=scipy.linalg.pinv(np.expand_dims(np.sum(inteartion_matrix,axis=1),axis=1))
    # aa=reduce(operator.add, row_sum_interation_matrix)
    # print(np.diag(aa).shape)
    # aaa=np.diagflat(row_sum_interation_matrix)
    initial_state_matrix=np.dot(np.diagflat(row_sum_interation_matrix),inteartion_matrix)
    score_matrix=np.dot(scipy.linalg.pinv(np.eye(num)-c*transformmation_matrix),initial_state_matrix)
    return score_matrix

''' Baseline methods- NR model'''
def neighbor_model(train_drug_drug_matrix):
    simi_matrix=cos_similarity(train_drug_drug_matrix)
    return_matrix=np.matrix(simi_matrix)*np.matrix(train_drug_drug_matrix)
    D=np.diag(simi_matrix.sum(axis=1))
    return_matrix=pinv(D)*return_matrix
    return_matrix=return_matrix+np.transpose(return_matrix).T
    return return_matrix

''' Baseline methods- LP model'''
def Label_Propagation(train_drug_drug_matrix):
    simi_matrix = cos_similarity(train_drug_drug_matrix)
    alpha = 0.9
    simi_matrix = np.matrix(simi_matrix)
    train_drug_drug_matrix = np.matrix(train_drug_drug_matrix)
    D = np.diag(((simi_matrix.sum(axis=1)).getA1()))
    N = pinv(D) * simi_matrix

    transform_matrix = (1 - alpha) * pinv(np.identity(len(simi_matrix)) - alpha * N)
    return_matrix = transform_matrix*train_drug_drug_matrix
    return_matrix = return_matrix + np.transpose(return_matrix).T
    return return_matrix

def l5_model_evaluate(interaction_matrix, predict_matrix):

    real_score = np.matrix(np.array(interaction_matrix).flatten())#A two-dimensional matrix becomes one-dimensional by row
    predict_score = np.matrix(np.array(predict_matrix).flatten())
    metrics = get_metrics(real_score, predict_score)
    return metrics

def get_metrics(real_score, predict_score):
    sorted_predict_score = np.array(sorted(list(set(np.array(predict_score).flatten()))))  # Remove duplication
    sorted_predict_score_num = len(sorted_predict_score)
    thresholds = sorted_predict_score[np.int32(sorted_predict_score_num * np.array(range(1, 1000)) / 1000)]
    thresholds = np.mat(thresholds)

    thresholds_num = thresholds.shape[1]

    predict_score_matrix = np.tile(predict_score, (thresholds_num, 1))#Row invariant; Copy by column 999  （n，） --》（999,n）
    negative_index = np.where(predict_score_matrix < thresholds.T)
    positive_index = np.where(predict_score_matrix >= thresholds.T)

    predict_score_matrix[negative_index] = 0
    predict_score_matrix[positive_index] = 1


    TP = predict_score_matrix * real_score.T
    FP = predict_score_matrix.sum(axis=1) - TP
    FN = real_score.sum() - TP
    TN = len(real_score.T) - TP - FP - FN

    fpr = FP / (FP + TN)
    tpr = TP / (TP + FN)

    ROC_dot_matrix = np.mat(sorted(np.column_stack((fpr, tpr)).tolist())).T
    ROC_dot_matrix.T[0] = [0, 0]

    ROC_dot_matrix = np.c_[ROC_dot_matrix, [1, 1]]

    x_ROC = ROC_dot_matrix[0].T
    y_ROC = ROC_dot_matrix[1].T
    auc = 0.5 * (x_ROC[1:] - x_ROC[:-1]).T * (y_ROC[:-1] + y_ROC[1:])


    recall_list = tpr
    precision_list = TP / (TP + FP)
    PR_dot_matrix = np.mat(sorted(np.column_stack((recall_list, precision_list)).tolist())).T
    PR_dot_matrix.T[0] = [0, 1]
    PR_dot_matrix = np.c_[PR_dot_matrix, [1, 0]]
    x_PR = PR_dot_matrix[0].T
    y_PR = PR_dot_matrix[1].T
    aupr = 0.5 * (x_PR[1:] - x_PR[:-1]).T * (y_PR[:-1] + y_PR[1:])

    f1_score_list = 2 * TP / (len(real_score.T) + TP - TN)
    accuracy_list = (TP + TN) / len(real_score.T)
    specificity_list = TN / (TN + FP)

    max_index = np.argmax(f1_score_list)
    f1_score = f1_score_list[max_index, 0]
    accuracy = accuracy_list[max_index, 0]
    specificity = specificity_list[max_index, 0]
    recall = recall_list[max_index, 0]
    precision = precision_list[max_index, 0]
    print("aupr:%f,auc:%f,f1_score:%f,acuuracy:%f,recall:%f,specificity:%f,precision:%f" % (
    aupr[0, 0], auc[0, 0], f1_score, accuracy, recall, specificity, precision))
    return [aupr[0, 0], auc[0, 0], f1_score, accuracy, recall, specificity, precision]


if __name__ == '__main__':
    print('奥力给！')

    interaction_matrix = pd.read_excel(mator_d_p_01, header=None)
    interaction_matrix = numpy.array(interaction_matrix)
    get_train_set(interaction_matrix)



