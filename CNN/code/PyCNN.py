import keras
from keras.utils import to_categorical
from keras.models import Model, Sequential
from keras.layers import Input, Dense, Dropout, Flatten, Conv2D, MaxPooling2D, concatenate, BatchNormalization, \
    Activation

import argparse
import os,sys

import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import pkg_resources
from keras import backend as K
from keras.layers import Conv2D
from keras.layers import GlobalAveragePooling2D
from keras.utils import to_categorical
from keras.layers import Maximum
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import confusion_matrix

from janggu import Janggu
from janggu import Scorer
from janggu import inputlayer
from janggu import outputdense
from janggu.data import Array
from janggu.data import Bioseq
from janggu.layers import Complement
from janggu.layers import DnaConv2D
from janggu.layers import Reverse
from janggu.utils import ExportClustermap
from janggu.utils import ExportJson
from janggu.utils import ExportScorePlot
from janggu.utils import ExportTsne
from janggu.utils import ExportTsv

# Visualize training history
from keras.models import Sequential
from keras.layers import Dense
import matplotlib.pyplot as plt

# fix random seed for reproducibility
seed = 7
np.random.seed(seed)

for i in range(len(sys.argv)):
    "sys.argv[%d] = %s" % (i, sys.argv[i])

if len(sys.argv) != 5:
    print('''python PyCNN.py positive negative num_eachSet windowsize''')
    exit(0)

positive = sys.argv[1]
negative = sys.argv[2]
num = int(sys.argv[3])
windowsize = int(sys.argv[4])

def nseqs(filename):
    return sum((1 for line in open(filename) if line[0] == '>'))


# load the dataset
SAMPLE_1 = positive
SAMPLE_2 = negative

# DNA sequences in one-hot encoding will be used as input
DNA = Bioseq.create_from_seq('dna', fastafile=[SAMPLE_1, SAMPLE_2],
                            order=1, datatags=['train'], cache=True)
DNA = np.resize(DNA, (DNA.shape[0], DNA.shape[1], DNA.shape[3], DNA.shape[2]))

# An array of 1/0 will be used as labels for training
Y = np.asarray([1 for line in range(nseqs(SAMPLE_1))] +
               [0 for line in range(nseqs(SAMPLE_2))])


LABELS = Array('y', Y, conditions=['Methylation'])
annot = pd.DataFrame(Y[:], columns=LABELS.conditions).applymap(
    lambda x: 'CTCF binding' if x == 1 else 'not binding').to_dict(orient='list')

num_train = num
num_test = num
num_train_pos = int(num_train*0.6)
num_total = 2*num

x_train = np.append(DNA[:num_train_pos], DNA[num_train:(num_train+num_train_pos)], axis = 0)
y_train = np.append(Y[:num_train_pos], Y[num_train:(num_train+num_train_pos)], axis = 0)
x_test = np.append(DNA[num_train_pos:num_train], DNA[(num_train+num_train_pos):num_total], axis=0)
y_test = np.append(Y[num_train_pos:num_train], Y[(num_train+num_train_pos):num_total], axis=0)

# suffle function
def shufflelists(lists):
    
    ri = np.random.permutation(len(lists[0]))
    out = []
    for l in lists:
        out.append(l[ri])
    return out

# suffle the datasets
[x_Train4D, y_Train] = shufflelists([x_train, y_train])
[x_Test4D, y_Test] = shufflelists([x_test, y_test])

print x_Train4D.shape, y_Train.shape
print x_Test4D.shape, y_Test.shape
y_Train = to_categorical(y_Train, num_classes=2)
y_Test = to_categorical(y_Test, num_classes=2)

def DNA(seq_length=1000, num_filters=256,filter_sizes=3,dropout_rate=0.3, num_classes=2, num_hidden=512):
    # initialization
    in_shape = (seq_length, 4, 1)
    input_shape = Input(shape=in_shape)
    conv = Conv2D(num_filters, (filter_sizes, 4), padding='valid', activation='relu')(input_shape)
    pool = MaxPooling2D((2, 1), padding='valid')(conv)
    x = Flatten()(pool)
    x = Dropout(dropout_rate)(x)
    x = Dense(num_hidden, activation='relu')(x)
    out = Dense(num_classes, activation='softmax')(x)
    model = Model(input_shape, out)
    return model

model = DNA(windowsize-1)
print model.summary()

model.compile(optimizer='sgd', loss='mean_squared_error', metrics=['acc'])

#model.fit(x_Train4D,y_Train,batch_size=32,epochs=100,verbose=1,shuffle='batch',validation_data=(x_test, y_test),callbacks=[SensitivitySpecificityCallback()])

train_his = model.fit(x_Train4D, y_Train, epochs=100, batch_size=32, validation_split=0.2)
print(train_his.history.keys())
# summarize history for accuracy
plt.figure()
plt.plot(train_his.history['acc'])
plt.plot(train_his.history['val_acc'])
plt.title('model accuracy')
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()
plt.savefig("modelaccuracy.png")
# summarize history for loss
plt.figure()
plt.plot(train_his.history['loss'])
plt.plot(train_his.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.show()
plt.savefig("modelloss.png")

predictions = model.predict(x_test)
#y_test = np.argmax(y_test, axis=-1)

predictions = np.argmax(predictions, axis=-1)
c = confusion_matrix(y_test, predictions)
print('Confusion matrix:\n', c)
print('sensitivity', float(c[0, 0])/(float(c[0,0])+float(c[0,1])))
print('specificity', float(c[1, 1])/(float(c[1,0])+float(c[1,1])))



score = model.evaluate(x_Test4D, y_Test)
print(score)
    
