from rdkit import Chem
from rdkit.Chem.EState import Fingerprinter
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn import cross_validation
from sklearn.neural_network import MLPRegressor
import matplotlib.pyplot as plt
data = pd.read_table('output.smi', sep='\t', header=None, names=['smiles','mol_id'])
data['Mol'] = data['smiles'].apply(Chem.MolFromSmiles)
def fps_plus_mw(mol):
    return np.append(Fingerprinter.FingerprintMol(mol)[0],Descriptors.MolWt(mol))
data['Descriptors'] = data['Mol'].apply(fps_plus_mw)

datart = pd.read_table('rtout.csv', sep=' ', header=None, names=['mol_id','rt1','rt2'])
result = pd.merge(data, datart[['mol_id', 'rt1', 'rt2']], on='mol_id')

X = result['Descriptors']
y = result['rt2']
X_list = []
Y_list = []
for x,_y in zip(X,y):
    X_list.append(x)
    Y_list.append(_y)
    print _y

X_list=np.asarray(X_list)
Y_list=np.asarray(Y_list)

clf = MLPRegressor(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(100, 100, 20), activation="identity",verbose=True)
clf.fit(X_list, Y_list)
y_pred = clf.predict(X_list)
rms = (np.mean((Y_list.reshape(-1,1) - y_pred)**2))**0.5
#s = np.std(y_test -y_pred)
print "Neural Network RMS", rms
x_range = range(1, len(X_list)+1)
# print len(y_pred)

plt.scatter(x_range, Y_list, label = 'Train', c='blue')
plt.title('Neural Network Predictor')
plt.xlabel('Measured Retention Time')
plt.ylabel('Predicted Retention Time')
plt.scatter(x_range, y_pred, c='lightgreen', label='Test', alpha = 0.8)
plt.legend(loc=4)
plt.show()
