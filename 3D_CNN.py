#Import libraries
import tensorflow as tf
import tensorflow.keras
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten, Conv3D, MaxPooling3D,BatchNormalization, SpatialDropout3D

import numpy as np
import pandas as pd
import warnings
import copy
from sklearn.model_selection import StratifiedKFold # For cross-validation training

#Set GPU settings
warnings.filterwarnings(action='ignore')
tf.debugging.set_log_device_placement(True)
tf.test.is_built_with_cuda()
tf.test.is_gpu_available('GPU') 
config=tf.compat.v1.ConfigProto()
config.gpu_options.allow_growth = True


#Download the list of IDs
ii=np.load('/media/leelabsg_storage01/jangho/Aging/ID_128_full.npy',allow_pickle=True)
a=ii.tolist()
ii_pd=pd.DataFrame(a)
cap=ii_pd.rename(columns={0:'id', 1:'age'})

disease=pd.read_csv('/home/leelabsg/jangho/UKB_ALL/PEDMASTER_ALL_20180514.txt', sep = '\t')
mapping = pd.read_csv('/home/lee7801/DATA/UKBB/Mapping/mapping.csv')
mapping_mri = mapping[mapping['Shawn'].isin(cap['id'].copy())]  #37792
print(disease.shape)
print(cap.shape)
print(mapping_mri.shape)

mapping_mri.columns = ['FID', 'Shawn'] 
mri_disease = pd.merge(disease, mapping_mri, on='FID') #37724
print(disease.shape)



# Explore what diseases the individuals have
mri_disease_count= np.sum(mri_disease.iloc[:, 70:1756], axis = 0)
# top ten disease
mri_disease_count[mri_disease_count.isin(np.sort(mri_disease_count)[::-1][:10])]
# exploration diagnosed
np.sum(mri_disease.iloc[:, 70:1756], axis = 1).value_counts()[:20]

#Filtering out individuals with diseases we chose
def preprocess_filt_disease(data, disease):
    hi_disease = data[np.sum(data[disease] == 1, axis = 1) >= 1].index
    data.drop(hi_disease, axis = 0, inplace = True)
    
    print(len(hi_disease), 'dropped')


#processed = preprocess_make_train(mri_disease, 9)
disease_list = ['X191','X198','X225','X249', 'X250', 'X290','X290.11', 'X291', 'X292', 'X303', 'X306', 'X315', 'X430','X433', 'X752',
                'X753']+list(mri_disease.columns[128:260])+list(mri_disease.columns[547:550])
# 191 Manlignant and unknown neoplasms of brain and nervous system
# 198 Secondary malignant neoplasm
# 225 Benign neoplasm of brain and other parts of nervous system
# 249, 250 diabetes
# 290~292 dimensia
# 303 Psychogenic and somatoform disorders
# 306 Other mental disorder
# 315 Develomental delays and disorders
# 331 Other cerebral degenerations
# 346 Abnormal findings on study of brain and/or nervous system
# 348 Other conditions of brain
# 430~433 cerebrovascular disease
# 752 congenital anomalies of nervous system, spine
# 753 congenital anomalies of nervous system, spine
# 781 Symptoms involving nervous and musculoskeletal systems
# -ing

# Cancers
# 145, 149, 150, 151 153, 155, 157, 158, 159, 164, 165, 170, 172, 173, 174, 175, 180, 182, 184, 185, 187, 189, 190, 191
# 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 204, 208, 210, 211, 212, 213, 214, 215, 216, 217, 218, 220, 221, 222, 
# 223, 224, 225, 226, 227, 228, 229, 230
processed=copy.deepcopy(mri_disease)
preprocess_filt_disease(processed, disease_list)
processed

seeya = mri_disease.drop(processed.index, axis = 0)
print(seeya.shape)

#MRI data 
final_array=np.load('/media/leelabsg_storage01/jangho/FINAL_128.npy')    

ID=np.array(cap)
ID_dat=pd.DataFrame(ID)
ID_dat=ID_dat.rename(columns={0:'id',1:'age'})
ID_dat['label'] = ID_dat.index
ID_dat=ID_dat.iloc[0:35012,0:3]


processed.rename(columns = {'Shawn' : 'id'}, inplace = True)
ok = pd.merge(ID_dat, processed[['id']], on='id')  
#Filtered out 8,473 individuals with diseases, 25,656 left 
print(ok.shape)

y = ok['age']  # Age of the samples
X = ok.drop(['age'], axis = 1) # ID of the samples

#Cross-validation
skf = StratifiedKFold(n_splits=4, random_state=7, shuffle=True) 
cv_index = {}

i = 0
for train, test in skf.split(X, y):
    cv_index[i] = [train, test]
    i += 1


#CV number
# n = 3
n = 0 
# n = 1
# n = 0

y_train, y_test = y.loc[cv_index[n][0]], y.loc[cv_index[n][1]] 
train_index = X.loc[cv_index[n][0]]['label']
test_index = X.loc[cv_index[n][1]]['label']
train_in = cv_index[n][0]
test_in = cv_index[n][1]

#ID list for getting result later
ID_result_train=X.loc[cv_index[n][0]]['id']
ID_result_test=X.loc[cv_index[n][1]]['id']

#Get MRI data of training set and test set
x_train = final_array[0:128,0:128,0:128,train_index]
# x_valid = final_array[0:128,0:128,0:20,valid_index]
x_test = final_array[0:128,0:128,0:128,test_index]

#Adjust coordinates
x_train = np.expand_dims(np.transpose(x_train, (3, 0, 1, 2)), axis = 4)
# x_valid = np.expand_dims(np.transpose(x_valid, (3, 0, 1, 2)), axis = 4)
x_test = np.expand_dims(np.transpose(x_test, (3, 0, 1, 2)), axis = 4)

#Get age of training set and test set
y_train = np.asarray(y_train).reshape(len(train_index))
# y_valid = np.asarray(y_valid).reshape(len(valid_index))
y_test = np.asarray(y_test).reshape(len(test_index))   

train_dataset = tf.data.Dataset.from_tensor_slices((x_train, y_train))
# valid_dataset = tf.data.Dataset.from_tensor_slices((x_valid, y_valid))
test_dataset = tf.data.Dataset.from_tensor_slices((x_test, y_test))


input_shape=(128,128,128, 1)
#Making empty model
def just_make():
    model=Sequential()
    model.add(Conv3D(32,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=input_shape))
    model.add(BatchNormalization())
    model.add(Conv3D(32,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=input_shape))
 #   model.add(BatchNormalization())
    model.add(SpatialDropout3D(0.2))
    model.add(MaxPooling3D(pool_size=(2,2,2)))
    model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(64,64,64)))
    model.add(BatchNormalization())
    model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(64,64,64)))
 #   model.add(BatchNormalization())
    model.add(SpatialDropout3D(0.2))
    model.add(MaxPooling3D(pool_size=(2,2,2)))
    model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(32,32,32)))
    model.add(BatchNormalization())
    model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(32,32,32)))
    model.add(BatchNormalization())
#     model.add(SpatialDropout3D(0.2))
    model.add(MaxPooling3D(pool_size=(2,2,2)))
    model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(16,16,16)))
    model.add(BatchNormalization())
    model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(16,16,16)))
    model.add(BatchNormalization())
    model.add(Conv3D(64,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(16,16,16)))
    model.add(BatchNormalization())
#     model.add(SpatialDropout3D(0.2))
    model.add(MaxPooling3D(pool_size=(2,2,2)))
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(8,8,8)))
    model.add(BatchNormalization())
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(8,8,8)))
    model.add(BatchNormalization())
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(8,8,8)))
    model.add(BatchNormalization())
#     model.add(SpatialDropout3D(0.2))
    model.add(MaxPooling3D(pool_size=(2,2,2)))
    
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(4,4,4)))
    model.add(BatchNormalization())
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(4,4,4)))
    model.add(BatchNormalization())
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(4,4,4)))
    model.add(BatchNormalization())
    
    model.add(MaxPooling3D(pool_size=(2,2,2)))
    
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(2,2,2)))
    model.add(BatchNormalization())
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(2,2,2)))
    model.add(BatchNormalization())
    model.add(Conv3D(96,kernel_size=(3,3,3),activation='relu',kernel_initializer='he_uniform',padding='same',input_shape=(2,2,2)))
    model.add(BatchNormalization())
    model.add(Flatten())
    model.add(Dense(96, activation='relu', kernel_initializer='he_uniform'))
    #     model.add(Dropout(0.2))
    model.add(Dense(32, activation='relu',kernel_initializer='he_uniform'))
    model.add(Dense(1))
    
    return model




def knot(train_dataset, model, directory_withname):
    
    optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
    
    model.compile(optimizer=optimizer,loss='mse',metrics=['mse', 'mae'])
    BATCH_SIZE = 32
    SHUFFLE_BUFFER_SIZE = 100
    train_dataset = train_dataset.shuffle(SHUFFLE_BUFFER_SIZE).batch(BATCH_SIZE)
    
    rep=1
    for i in range(rep):
        model.fit(train_dataset, epochs=40)
    model.save_weights(directory_withname, overwrite=True)
    
    return model

ensemble1=just_make()
ensemble2=just_make()
ensemble3=just_make()

ensemble1 = knot(train_dataset, ensemble1, '/home/leelabsg/jangho/FINAL_oneset_CV1_m1')
ensemble2 = knot(train_dataset, ensemble2, '/home/leelabsg/jangho/FINAL_oneset_CV1_m2')
ensemble3 = knot(train_dataset, ensemble3, '/home/leelabsg/jangho/FINAL_oneset_CV1_m3')

BATCH_SIZE = 1
train_dataset = train_dataset.batch(BATCH_SIZE)
test_dataset = test_dataset.batch(BATCH_SIZE)

RESULT_train=np.zeros(shape=(25656*3/4,3)) #19,242
RESULT=np.zeros(shape=(25656/4,3)) #6,414
RES_INDEX_train=0
RES_INDEX=0

#Get result from the test set
ind=0
for a,b in test_dataset:
    res1=ID_result_test.iloc[ind]   #ID
    res2_1=float(ensemble1(a))
    res2_2=float(ensemble2(a))
    res2_3=float(ensemble3(a))
    res2=(res2_1+res2_2+res2_3)/3  #Predicted Age
    res3=float(b) #True Age
    ind+=1
    RESULT[RES_INDEX,0]=res1
    RESULT[RES_INDEX,1]=res2
    RESULT[RES_INDEX,2]=res3
    RES_INDEX+=1


#Get result from the training set
ind=0
for a,b in train_dataset:
    res1=ID_result_train.iloc[ind]
    res2_1=float(ensemble1(a))
    res2_2=float(ensemble2(a))
    res2_3=float(ensemble3(a))
    res2=(res2_1+res2_2+res2_3)/3
    res3=float(b)
                                                                        
    ind+=1
    RESULT_train[RES_INDEX_train,0]=res1
    RESULT_train[RES_INDEX_train,1]=res2
    RESULT_train[RES_INDEX_train,2]=res3
    RES_INDEX_train+=1


np.save('CV1_RESULT_TEST',RESULT)
np.save('CV1_RESULT_TRAINING',RESULT_train)


