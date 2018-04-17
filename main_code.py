from rdkit import Chem
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn import cross_validation
from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.optimizers import SGD
import matplotlib.pyplot as plt
data = pd.read_table('output.smi', sep='\t', header=None, names=['smiles','mol_id'])
data['Mol'] = data['smiles'].apply(Chem.MolFromSmiles)
def fps_plus_mw(mol):
    return np.append(Fingerprinter.FingerprintMol(mol)[0],Descriptors.MolWt(mol))
data['Descriptors'] = data['Mol'].apply(fps_plus_mw)

datart = pd.read_table('rtout.csv', sep=' ', header=None, names=['mol_id','rt1','rt2'])
result = pd.merge(data, datart[['mol_id', 'rt1', 'rt2']], on='mol_id')
#print result['Descriptors']

X = np.array(list(result['Descriptors']))
y = result['rt2'].values
st = StandardScaler()
X= st.fit_transform(X)
X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, y, test_size=0.25, random_state=42)
model = Sequential()
model.add(Dense(output_dim=50, input_dim=X.shape[1]))
model.add(Activation("sigmoid"))
model.add(Dense(output_dim=1))
model.add(Activation("linear"))
model.compile(loss='mean_squared_error', optimizer=SGD(lr=0.009, momentum=0.9, clipnorm=1.,nesterov=True))

history = model.fit(X_train, y_train, nb_epoch=10000, batch_size=20)
y_pred = model.predict(X_test)
# print "Accuracy", binary_accuracy(y_pred, y_test)
rms = (np.mean((y_test.reshape(-1,1) - y_pred)**2))**0.5
#s = np.std(y_test -y_pred)
print "Neural Network RMS", rms
plt.scatter(y_train,model.predict(X_train), label = 'Train', c='blue')
plt.title('Neural Network Predictor')
plt.xlabel('Measured Retention Time')
plt.ylabel('Predicted Retention Time')
plt.scatter(y_test,model.predict(X_test),c='lightgreen', label='Test', alpha = 0.8)
plt.legend(loc=4)
plt.show()
