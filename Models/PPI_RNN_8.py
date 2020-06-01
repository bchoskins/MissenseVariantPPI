#PPI Prediction with a Siamese RCNN 
import pandas as pd 
import sys
import keras
from keras.models import Sequential, Model
from keras.layers import Dense, Activation, Dropout, Embedding, LSTM, Bidirectional, BatchNormalization, merge, add
from keras.layers.core import Flatten, Reshape
from keras.layers.merge import Concatenate, concatenate, subtract, multiply
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import MaxPooling1D, AveragePooling1D, GlobalAveragePooling1D
from keras.optimizers import Adam,  RMSprop
import os
import tensorflow as tf
import keras.backend.tensorflow_backend as KTF
import numpy as np
from tqdm import tqdm
from keras.layers import Input, CuDNNGRU, GRU, CuDNNLSTM
from numpy import linalg as LA
import scipy

#data = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/UniprotModelDataWithNegatives.csv')
#proteins = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/AllProteins.csv', index_col=0)
#protein_good = np.load('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/JustProteinEmbeddings.npy')

#import os
#os.system('top')

#run on UI-GPU or UI-GPU-HM
binary = pd.read_csv('/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/UniprotModelDataBinary.csv',index_col=0)
seq_A = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqA_pad.npy")
seq_B = np.load("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqB_pad.npy")


#char_count = proteins['Sequence_A'].apply(lambda x: len(x))
#char_count.describe()
#count    17497.000000
#mean       591.633366
#std        619.547231
#min         24.000000
#25%        284.000000
#50%        441.000000
#75%        706.000000
#max      34350.000000

""" import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sns.boxplot(char_count)
plt.title(f'Sequence char count')
plt.grid(True)
plt.xticks(np.arange(0, 15000, 2000.0))
plt.savefig("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SequenceDistribution.png")
 """

#creating dictionary for embedding
""" embDict = {}
embDict = {'A': 1, 'G':2, 'V':3, 'I':4, 'L':5, 'F':6, 'P':7, 'Y':8, 'M':9, 'T':10, 'S':11, 'H':12, 'N':13, 'Q':14, 'W':15, 'R':16, 'K':17, 'D':18, 'E':19, 'C':20, 'U':21}

#loop through seq_index to get the sequences, split them, embed them with above dict, one hot encode the values, store into vector matrix 
def aminoEncoding(data):
  encode_list = []
  for row in data['Sequence_A'].values:
    row_encode = []
    for code in row:
      print(code)
      row_encode.append(embDict.get(code, 0))
    encode_list.append(np.array(row_encode))
  
  return encode_list

Protein_Encode = aminoEncoding(proteins)

from keras.preprocessing.sequence import pad_sequences
#a little more than one standard deviation above the third quartile lengths
max_length = 1000
protein_pad = pad_sequences(Protein_Encode, maxlen=max_length, padding='post', truncating='post')
protein_pad.shape

from keras.utils import to_categorical
protein_good = to_categorical(protein_pad)
protein_good.shape
type(protein_good)

#np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/JustProteinEmbeddings.npy", protein_good)

seq_A_pad = np.array([protein_good[np.hstack(np.where(x == proteins['Protein_A']))] for x in binary['Protein_A']])
seq_B_pad = np.array([protein_good[np.hstack(np.where(x == proteins['Protein_A']))] for x in binary['Protein_B']])

seq_A = np.array(seq_A_pad)
seq_B = np.array(seq_B_pad)

seq_A = seq_A[:,0,:,:]
seq_B = seq_B[:,0,:,:]

np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqA_pad.npy", seq_A)
np.save("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqB_pad.npy", seq_B)
  """
###############Look into remaking model################ 
def build_model():
    seq_input1 = Input(shape=(1000, 22), name='seq1')
    seq_input2 = Input(shape=(1000, 22), name='seq2')
    l1=Conv1D(hidden_dim, 3)
    r1=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l2=Conv1D(hidden_dim, 3)
    r2=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l3=Conv1D(hidden_dim, 3)
    r3=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l4=Conv1D(hidden_dim, 3)
    r4=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l5=Conv1D(hidden_dim, 3)
    r5=Bidirectional(CuDNNGRU(hidden_dim, return_sequences=True))
    l6=Conv1D(hidden_dim, 3)
    s1=MaxPooling1D(3)(l1(seq_input1))
    s1=concatenate([r1(s1), s1])
    s1=MaxPooling1D(3)(l2(s1))
    s1=concatenate([r2(s1), s1])
    s1=MaxPooling1D(3)(l3(s1))
    s1=concatenate([r3(s1), s1])
    s1=MaxPooling1D(3)(l4(s1))
    s1=concatenate([r4(s1), s1])
    s1=MaxPooling1D(3)(l5(s1))
    s1=concatenate([r5(s1), s1])
    s1=l6(s1)
    s1=GlobalAveragePooling1D()(s1)
    s2=MaxPooling1D(3)(l1(seq_input2))
    s2=concatenate([r1(s2), s2])
    s2=MaxPooling1D(3)(l2(s2))
    s2=concatenate([r2(s2), s2])
    s2=MaxPooling1D(3)(l3(s2))
    s2=concatenate([r3(s2), s2])
    s2=MaxPooling1D(3)(l4(s2))
    s2=concatenate([r4(s2), s2])
    s2=MaxPooling1D(3)(l5(s2))
    s2=concatenate([r5(s2), s2])
    s2=l6(s2)
    s2=GlobalAveragePooling1D()(s2)
    merge_text = multiply([s1, s2])
    x = Dense(100, activation='linear')(merge_text)
    x = keras.layers.LeakyReLU(alpha=0.3)(x)
    x = Dense(int((hidden_dim+7)/2), activation='linear')(x)
    x = keras.layers.LeakyReLU(alpha=0.3)(x)
    main_output = Dense(1, activation='sigmoid')(x)
    merge_model = Model(inputs=[seq_input1, seq_input2], outputs=[main_output])
    return merge_model


hidden_dim = 25
n_epochs=50
batch_size = 512
adam = Adam(lr=0.001, amsgrad=True, epsilon=1e-6) 

#class label train, validate, test split
class_lbl =  binary.iloc[:, 9].as_matrix()

pos = np.where(class_lbl == 1)
neg = np.where(class_lbl == 0)

import random
random.seed(4)
np.random.shuffle(pos[0])
np.random.shuffle(neg[0])

#train, test indeces
pos = np.array((pos[0][0:125000], pos[0][125000:]))
neg = np.array((neg[0][0:520000], neg[0][520000:]))

train = list(pos[0]) + list(neg[0])
test =  list(pos[1]) + list(neg[1])

np.random.shuffle(train)
np.random.shuffle(test)

train_test = []
train_test.append((train, test))

es = keras.callbacks.EarlyStopping(monitor='val_loss', min_delta=0,patience=5,verbose=0, mode='auto')

cp = keras.callbacks.ModelCheckpoint(filepath="/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqModel_8.h5",
        verbose=1, save_best_only=True)

from tensorflow.python.client import device_lib
print(device_lib.list_local_devices())

for train, test in train_test:
  merge_model = None
  merge_model = build_model()
  adam = Adam(lr=0.001, amsgrad=True, epsilon=1e-6)
  merge_model.compile(optimizer=adam, loss='mse', metrics=['accuracy'])
  history = merge_model.fit([seq_A[train], seq_B[train]], class_lbl[train], callbacks=[es,cp], batch_size=batch_size, epochs=n_epochs, validation_split=0.2) 


hist_df = pd.DataFrame(history.history) 
hist_df.to_csv("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/history8.csv")

""" import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
#plt.grid(True)
#plt.xticks(np.arange(0, 15000, 2000.0))
plt.savefig("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/train_val_loss.png")

#post 33 epochs, stopped running manually to check
#Epoch 33/50
#516000/516000 [==============================] - 248s 480us/step - loss: 0.0383 - acc: 0.9530 - val_loss: 0.0438 - val_acc: 0.9457
pred = merge_model.predict([seq_A[test], seq_B[test]])

true = class_lbl[test]

from sklearn.metrics import roc_auc_score
roc = roc_auc_score(true, pred)
#Model1 = 0.9571048531995481

#model = keras.models.load_model("/Dedicated/jmichaelson-wdata/rotating_students/bhoskins/ProteinData/SeqModel.h5") """
