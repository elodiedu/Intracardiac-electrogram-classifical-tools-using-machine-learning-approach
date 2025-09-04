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

# Load Data Directory, replace that with the data directory folder
directory = 'C:/01_Study/09_2023_Fall/data'

# Load RA1 features
def load_and_flatten_mat_variable(file, key):
    cell_array = sio.loadmat(file)[key]
    reshaped = np.array([[cell_array[0, i][0, j] for j in range(15)] for i in range(102)])
    return reshaped.ravel()

RAKT = load_and_flatten_mat_variable(os.path.join(directory, 'RA1KT.mat'), 'RA1KT')
RAMSE = load_and_flatten_mat_variable(os.path.join(directory, 'RA1MSE.mat'), 'RA1MSE')
RAMSF = load_and_flatten_mat_variable(os.path.join(directory, 'RA1MSF.mat'), 'RA1MSF')

# Load Sinus features
mat_contents = sio.loadmat(os.path.join(directory, 'iEGMsInformationSinus.mat'))
DF3 = np.vstack(mat_contents['BipolarFE_Hz'][0]).ravel()
KT3 = np.vstack(mat_contents['BipolarKT'][0]).ravel()
MI3 = np.vstack(mat_contents['BiSampEn'][0]).ravel()

# Combine into feature matrix and labels
X_orange = np.column_stack((RAKT, RAMSE, RAMSF))
X_blue = np.column_stack((KT3, MI3, DF3))
X_final = np.vstack((X_orange, X_blue))
y_final = np.concatenate((np.zeros(len(RAKT)), np.ones(len(KT3))))

X_final, y_final = shuffle(X_final, y_final, random_state=42)

# Performance Metrics
def evaluate_performance(y_true, y_pred):
    acc = accuracy_score(y_true, y_pred)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    Se = tp / (tp + fn)
    Sp = tn / (tn + fp)
    f1 = f1_score(y_true, y_pred)
    return acc, Se, Sp, f1

# Plotting Functions
from skimage.measure import marching_cubes
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection

# Updated 3D Decision Boundary Plotter with border lines
def plot_3d_decision_boundary(X, y, model, title):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter data
    ax.scatter(X[y == 0][:, 0], X[y == 0][:, 1], X[y == 0][:, 2], c='orange', alpha=0.5, label='RA Sinus')
    ax.scatter(X[y == 1][:, 0], X[y == 1][:, 1], X[y == 1][:, 2], c='blue', alpha=0.5, label='LA Sinus')

    # Mesh grid for marching cubes
    x_range = np.linspace(X[:, 0].min(), X[:, 0].max(), 30)
    y_range = np.linspace(X[:, 1].min(), X[:, 1].max(), 30)
    z_range = np.linspace(X[:, 2].min(), X[:, 2].max(), 30)
    xx, yy, zz = np.meshgrid(x_range, y_range, z_range)
    grid_points = np.c_[xx.ravel(), yy.ravel(), zz.ravel()]

    try:
        values = model.decision_function(grid_points)
    except:
        try:
            values = model.predict_proba(grid_points)[:, 1] - 0.5
        except:
            return  # No decision surface can be computed (like KMeans)

    values = values.reshape(xx.shape)

    # Marching cubes to extract boundary surface
    verts, faces, _, _ = marching_cubes(values, level=0, spacing=(
        x_range[1] - x_range[0], y_range[1] - y_range[0], z_range[1] - z_range[0]))
    verts[:, 0] += x_range[0]
    verts[:, 1] += y_range[0]
    verts[:, 2] += z_range[0]

    # Mesh and wireframe
    mesh = Poly3DCollection(verts[faces], alpha=0.2, facecolor='green', edgecolor='black', linewidth=0.2)
    ax.add_collection3d(mesh)

    ax.set_xlabel('KT')
    ax.set_ylabel('MSE')
    ax.set_zlabel('MSF')
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()
    plt.show()


def plot_certainty_histogram(model, X, y, title):
    try:
        probs = model.predict_proba(X)[:, 1]
    except AttributeError:
        probs = model.decision_function(X)
        probs = (probs - probs.min()) / (probs.max() - probs.min())
    class0 = probs[y == 0]
    class1 = probs[y == 1]
    plt.hist(class0, bins=30, alpha=0.5, label='RA Sinus', color='lightgrey')
    plt.hist(class1, bins=30, alpha=0.5, label='LA Sinus', color='dimgray', hatch='\\')
    plt.title(title)
    plt.xlabel('Predicted Probability')
    plt.ylabel('Count')
    plt.legend()
    plt.grid(True)
    plt.show()

# Function to run classifier with GridSearchCV and evaluation
def run_classifier(model_instance, param_grid, X, y, title):
    grid_search = GridSearchCV(model_instance, param_grid, cv=5)
    grid_search.fit(X, y)
    print(f"\nBest {title} Parameters:", grid_search.best_params_)

    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    acc_list, se_list, sp_list, f1_list = [], [], [], []

    for train_idx, test_idx in kf.split(X, y):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        model = type(model_instance)(**grid_search.best_params_)
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        acc, Se, Sp, f1 = evaluate_performance(y_test, y_pred)
        acc_list.append(acc)
        se_list.append(Se)
        sp_list.append(Sp)
        f1_list.append(f1)

    print(f"{title} Performance:")
    print(f"Accuracy: {np.mean(acc_list):.2f}, Sensitivity: {np.mean(se_list):.2f}, Specificity: {np.mean(sp_list):.2f}, F1 Score: {np.mean(f1_list):.2f}")
    plot_3d_decision_boundary(X_test, y_test, model, f'{title} 3D Boundary')
    plot_certainty_histogram(model, X_test, y_test, f'{title} Prediction Certainty')

# Run classifiers
run_classifier(DecisionTreeClassifier(), {'max_depth': [3, 5, 10, None], 'min_samples_split': [2, 5, 10]}, X_final, y_final, 'Decision Tree')
run_classifier(SVC(), {'C': [0.1, 1, 10], 'gamma': ['scale', 'auto'], 'kernel': ['rbf','polynomial','linear'], 'probability': [True]}, X_final, y_final, 'SVM')
run_classifier(RandomForestClassifier(), {'n_estimators': [10, 50, 100], 'max_depth': [None, 3,5,10]}, X_final, y_final, 'Random Forest')


# KMeans Clustering with decision surface and evaluation
def run_kmeans_with_evaluation(X, y, title):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Run KMeans multiple times for evaluation
    acc_list, se_list, sp_list, f1_list = [], [], [], []
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    for train_idx, test_idx in kf.split(X_scaled, y):
        X_train, X_test = X_scaled[train_idx], X_scaled[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        kmeans = KMeans(n_clusters=2, random_state=42)
        kmeans.fit(X_train)

        # Align labels
        pred_labels = kmeans.predict(X_test)
        if accuracy_score(y_test, pred_labels) < 0.5:
            pred_labels = 1 - pred_labels

        acc, Se, Sp, f1 = evaluate_performance(y_test, pred_labels)
        acc_list.append(acc)
        se_list.append(Se)
        sp_list.append(Sp)
        f1_list.append(f1)

    print(f"{title} Performance:")
    print(
        f"Accuracy: {np.mean(acc_list):.2f}, Sensitivity: {np.mean(se_list):.2f}, Specificity: {np.mean(sp_list):.2f}, F1 Score: {np.mean(f1_list):.2f}")

    # Train final model for visualization
    kmeans = KMeans(n_clusters=2, random_state=42).fit(X_scaled)
    final_labels = kmeans.labels_
    if accuracy_score(y, final_labels) < 0.5:
        final_labels = 1 - final_labels
        kmeans.labels_ = final_labels

    # Use a dummy class wrapper to mimic classifier
    class KMeansWrapper:
        def __init__(self, km):
            self.km = km

        def predict_proba(self, X):
            dists = self.km.transform(X)
            probs = np.exp(-dists)
            return (probs / probs.sum(axis=1, keepdims=True))

    wrapped_kmeans = KMeansWrapper(kmeans)

    plot_3d_decision_boundary(X_scaled, y, wrapped_kmeans, f'{title} 3D View')
    plot_certainty_histogram(wrapped_kmeans, X_scaled, y, f'{title} Certainty')


# Run the KMeans clustering as classifier-like evaluation
run_kmeans_with_evaluation(X_final, y_final, 'KMeans Clustering')
