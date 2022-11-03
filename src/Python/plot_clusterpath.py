import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def extract(data):
    n_lambdas = len(np.unique(data[:,1]))
    X = np.empty([n_lambdas, int(data.shape[0]/n_lambdas)])
    Y = np.empty_like(X)
    
    for i in range(X.shape[1]):
        X[:,i] = data[data[:,0]==(i+1),2]
        Y[:,i] = data[data[:,0]==(i+1),3]
    
    return X, Y


def find_start(arr):
    start = 0
    
    for i in range(1, len(arr)):
        if arr[i]!=arr[i-1]+1:
            start=i
    
    arr = arr[start:]
    
    return arr


def no_duplicate_lines(X, Y, r):
    Xs = []
    Ys = []
    Xs.append(X[:,0])
    Ys.append(Y[:,0])
    
    
    for i in range(1, X.shape[1]):
        ps = []
        for j in range(len(Xs)):
            l1 = len(Xs[j])
            r1 = np.argwhere(Xs[j].round(r)==X[:l1,i].round(r))[:,0]
            r2 = np.argwhere(Ys[j].round(r)==Y[:l1,i].round(r))[:,0]
            
            r1 = find_start(r1)
            r2 = find_start(r2)
            
            if len(r1)>0 and len(r2)>0:
                p = max(r1[0], r2[0])
                if p in r1 and p in r2:
                    ps.append(p)
        
        if len(ps)==0:
            Xs.append(X[:,i])
            Ys.append(Y[:,i])
        else:
            Xs.append(X[:(min(ps)+1),i])
            Ys.append(Y[:(min(ps)+1),i])
    
    return Xs, Ys


def plot_clusterpath(data, y, colors=None, r_tol=3, s=3, lw=0.4, alpha=1, 
                     lc="darkgrey", title=None, fs_label=7, fs_title=7,
                     labels=None, lab_y_off=None, lab_x_off=None,
                     as_background=False, mark_end_points=True, Roman=False):
    if colors is None:
        # n_colors = max(3, np.unique(np.array(y, dtype=int)).size)
        # colors = sns.color_palette("Set1", n_colors)
        colors = sns.color_palette("muted")
    elif colors == "Greys":
        n_colors = np.unique(np.array(y, dtype=int)).size
        colors = sns.color_palette("Greys", n_colors)
    
    X, Y = extract(data)
    Xs, Ys = no_duplicate_lines(X, Y, r_tol)
    
    unique_idx = np.unique(np.round(X[-1,:], 3), return_index=True)[1]
    
    y_cols = []
    for i in range(len(y)):
        y_cols.append(colors[int(y[i])])

    plt.scatter(X[0,:], Y[0,:], s=s, c=y_cols, zorder=1)
    
    if not (labels == None):
        for x, y, l, l1, l2 in zip(X[0,:], Y[0,:], labels, lab_x_off, 
                                   lab_y_off):
            plt.text(x+l1, y+l2, l, fontsize=fs_label, 
                     horizontalalignment='center')
    
    for j in range(len(Xs)):
        plt.plot(Xs[j], Ys[j], lw=lw, c=lc, zorder=0, alpha=alpha)
    
    if not as_background and mark_end_points:
        plt.scatter(X[-1,unique_idx], Y[-1,unique_idx], s=s, c="black")
    
    plt.xticks([])
    plt.yticks([])
    
    if not as_background:
        if not Roman:
            plt.xlabel("$x_1$", fontsize=fs_label)
            plt.ylabel("$x_2$", fontsize=fs_label, rotation=0, labelpad=9)
        else:
            plt.xlabel("$\mathrm{x_1}$", fontsize=fs_label)
            plt.ylabel("$\mathrm{x_2}$", fontsize=fs_label, rotation=0,
                       labelpad=9)
    else:
        plt.axis("off")
    
    if title is not None:
        plt.title(title, fontsize=fs_title)
    
    
    return None
