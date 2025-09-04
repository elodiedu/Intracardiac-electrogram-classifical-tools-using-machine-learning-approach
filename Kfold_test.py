import os
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.model_selection import StratifiedKFold, GridSearchCV, train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.utils import shuffle
from matplotlib.colors import ListedColormap
from skimage.measure import marching_cubes
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# --- Load and Prepare Data ---
directory = 'C:/01_Study/09_2023_Fall/data'

def load_and_flatten_mat_variable(file, key):
    cell_array = sio.loadmat(file)[key]
    reshaped = np.array([[cell_array[0, i][0, j] for j in range(15)] for i in range(102)])
    return reshaped.ravel()

RAKT = load_and_flatten_mat_variable(os.path.join(directory, 'RA1KT.mat'), 'RA1KT')
RAMSE = load_and_flatten_mat_variable(os.path.join(directory, 'RA1MSE.mat'), 'RA1MSE')
RAMSF = load_and_flatten_mat_variable(os.path.join(directory, 'RA1MSF.mat'), 'RA1MSF')

mat_contents = sio.loadmat(os.path.join(directory, 'iEGMsInformationSinus.mat'))
DF3 = np.vstack(mat_contents['BipolarFE_Hz'][0]).ravel()
KT3 = np.vstack(mat_contents['BipolarKT'][0]).ravel()
MI3 = np.vstack(mat_contents['BiSampEn'][0]).ravel()

X_orange = np.column_stack((RAKT, RAMSE, RAMSF))
X_blue = np.column_stack((KT3, MI3, DF3))
X_final = np.vstack((X_orange, X_blue))
y_final = np.concatenate((np.zeros(len(RAKT)), np.ones(len(KT3))))
X_final, y_final = shuffle(X_final, y_final, random_state=42)

# --- Performance Metric Function ---
def evaluate_performance(y_true, y_pred):
    acc = accuracy_score(y_true, y_pred)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    Se = tp / (tp + fn)
    Sp = tn / (tn + fp)
    f1 = f1_score(y_true, y_pred)
    return acc, Se, Sp, f1

# --- Plotting Helpers ---




# --- Enhanced run_classifier ---
def run_classifier(model_cls, param_grid, X, y, title):
    # train/validation/test split
    X_tv, X_test, y_tv, y_test = train_test_split(X, y, test_size=0.2, stratify=y, random_state=42)

    grid = GridSearchCV(model_cls(), param_grid, cv=5, return_train_score=True)
    grid.fit(X_tv, y_tv)

    # Top 5 param sets
    res = grid.cv_results_
    idx = np.argsort(-res['mean_test_score'])[:5]
    print(f"\n=== {title} Top 5 Params ===")
    for i in idx:
        print(f"params={res['params'][i]}, train_acc={res['mean_train_score'][i]:.3f}, val_acc={res['mean_test_score'][i]:.3f}")

    # 5-fold train vs val metrics
    print(f"\n=== {title} 5-Fold Metrics (Train vs Val) ===")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    for fold, (ti, vi) in enumerate(skf.split(X_tv, y_tv),1):
        X_tr, X_val = X_tv[ti], X_tv[vi]; y_tr, y_val = y_tv[ti], y_tv[vi]
        mdl = model_cls(**grid.best_params_); mdl.fit(X_tr, y_tr)
        tr_acc, _, _, tr_f1 = evaluate_performance(y_tr, mdl.predict(X_tr))
        vl_acc, _, _, vl_f1 = evaluate_performance(y_val, mdl.predict(X_val))
        print(f"Fold{fold}: Train Acc={tr_acc:.3f}, F1={tr_f1:.3f} | Val Acc={vl_acc:.3f}, F1={vl_f1:.3f}")

    # final test
    best = model_cls(**grid.best_params_); best.fit(X_tv, y_tv)
    t_acc, _, _, t_f1 = evaluate_performance(y_test, best.predict(X_test))
    print(f"\n=== {title} Test Performance ===\nTest Acc={t_acc:.3f}, F1={t_f1:.3f}")



# --- Execute classifiers ---

run_classifier(SVC, {'C':[0.1,1,10], 'gamma':['scale','auto'], 'kernel':['linear']}, X_final, y_final, 'SVM')
run_classifier(RandomForestClassifier, {'n_estimators':[10,50,100], 'max_depth':[None,3,5,10]}, X_final, y_final, 'Random Forest')

# --- KMeans Evaluation (unchanged) ---
def run_kmeans_with_evaluation(X, y, title):
    scaler = StandardScaler(); X_s = scaler.fit_transform(X)
    accs, ses, sps, f1s = [],[],[],[]
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    for ti, vi in skf.split(X_s,y):
        X_tr, X_te = X_s[ti], X_s[vi]; y_tr, y_te = y[ti], y[vi]
        km = KMeans(n_clusters=2, random_state=42).fit(X_tr)
        pl = km.predict(X_te)
        if accuracy_score(y_te, pl)<0.5: pl=1-pl
        acc, Se, Sp, f1 = evaluate_performance(y_te, pl)
        accs.append(acc); ses.append(Se); sps.append(Sp); f1s.append(f1)
    print(f"\n{title} Performance: Acc={np.mean(accs):.2f}, Se={np.mean(ses):.2f}, Sp={np.mean(sps):.2f}, F1={np.mean(f1s):.2f}")
    class KMeansWrapper:
        def __init__(self, km): self.km=km
        def predict_proba(self, X): d=np.exp(-self.km.transform(X)); return d/d.sum(axis=1,keepdims=True)
    km_final = KMeans(n_clusters=2, random_state=42).fit(X_s)
    wrap = KMeansWrapper(km_final)


run_kmeans_with_evaluation(X_final, y_final, 'KMeans Clustering')
