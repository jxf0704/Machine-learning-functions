# -*- coding: utf-8 -*-

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

%matplotlib inline

from sklearn import set_config
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OrdinalEncoder

from sksurv.datasets import load_gbsg2
from sksurv.preprocessing import OneHotEncoder
from sksurv.ensemble import RandomSurvivalForest
from sksurv.datasets import load_breast_cancer

set_config(display="text")  # displays text representation of estimators
X, y = load_gbsg2()
X, y = load_breast_cancer()
grade_str = X.loc[:, "tgrade"].astype(object).values[:, np.newaxis]
grade_num = OrdinalEncoder(categories=[["I", "II", "III"]]).fit_transform(grade_str)

X_no_grade = X.drop("tgrade", axis=1)
Xt = OneHotEncoder().fit_transform(X_no_grade)
Xt.loc[:, "tgrade"] = grade_num

random_state = 20

X_train, X_test, y_train, y_test = train_test_split(Xt, y, test_size=0.25, random_state=random_state)

frag_arm_traing=pd.read_csv("frag_arm_feature.csv")
frag_arm_testing=pd.read_csv("frag_arm_feature_testing.csv")
rfs_data=pd.read_csv("RFS_data_training_for_rf_python.csv")
#rfs_data_array = np.array(tuple(x) for x in rfs_data.to_numpy())
#rfs_data_array = np.array([(row['RFS_status'], row['RFS_days']) for _, row in rfs_data.iterrows()])

dtype = [('event', bool), ('time', float)]
numpy_array = np.empty(len(rfs_data), dtype=dtype)
numpy_array['event'] = rfs_data['RFS_status'].values
numpy_array['time'] = rfs_data['RFS_days'].values


rsf = RandomSurvivalForest(n_estimators=1000, min_samples_split=10, min_samples_leaf=15, n_jobs=-1, random_state=random_state)
rsf.fit(frag_arm_traing, numpy_array)

rsf.predict(frag_arm_testing.iloc[[0]])

df = pd.DataFrame(rsf.predict(frag_arm_testing), columns=['column_name'])
df.to_csv("frag_arm_testing_risk_score.csv")

def rfs_response_array(rfs_dataframe=rfs_data):
    dtype = [('event', bool), ('time', float)]
    rfs_array = np.empty(len(rfs_dataframe), dtype=dtype)
    rfs_array['event'] = rfs_dataframe['RFS_status'].values
    rfs_array['time'] = rfs_dataframe['RFS_days'].values
    return rfs_array

def loocv_random_forest_one_feature(feature_table=frag_arm_traing,rfs_data_table=rfs_data):
    testing_score=[]
    for i in range(feature_table.shape[0]):
        test_sample=feature_table.iloc[[i]]
        train_sample=feature_table.drop(i)
        train_y_data=rfs_data.drop(i)
        train_y_array=rfs_response_array(rfs_dataframe=train_y_data)
        rsf = RandomSurvivalForest(n_estimators=1000, min_samples_split=10, min_samples_leaf=15, n_jobs=-1, random_state=20)
        rsf.fit(train_sample,train_y_array)
        test_score=rsf.predict(test_sample)[0]
        testing_score.append(test_score)
    risk_score = pd.DataFrame({'risk_score': testing_score})    
    return risk_score
 
rfs_loocv_score = pd.concat([rfs_data, risk_score], axis=1)       
rfs_loocv_score.to_csv("rfs_loocv_score_frag_arm.csv")

frag_training=pd.read_csv("frag_feature.csv")
cnv_training=pd.read_csv("cnv_feature.csv")
griffin_training=pd.read_csv("griffin_feature.csv")
mcms_training=pd.read_csv("mcms_feature.csv")

frag_loocv=loocv_random_forest_one_feature(feature_table=frag_training,rfs_data_table=rfs_data)
cnv_loocv=loocv_random_forest_one_feature(feature_table=cnv_training,rfs_data_table=rfs_data)
griffin_loocv=loocv_random_forest_one_feature(feature_table=griffin_training,rfs_data_table=rfs_data)
mcms_loocv=loocv_random_forest_one_feature(feature_table=mcms_training,rfs_data_table=rfs_data)

rfs_loocv_score_frag = pd.concat([rfs_data, frag_loocv], axis=1)       
rfs_loocv_score_frag.to_csv("rfs_loocv_score_frag.csv")

rfs_loocv_score_cnv = pd.concat([rfs_data, cnv_loocv], axis=1)       
rfs_loocv_score_cnv.to_csv("rfs_loocv_score_cnv.csv")

rfs_loocv_score_griffin = pd.concat([rfs_data, griffin_loocv], axis=1)       
rfs_loocv_score_griffin.to_csv("rfs_loocv_score_griffin.csv")
