import numpy as np
import matplotlib.pyplot as plt

def load_dendr(path):
    merge = np.genfromtxt(path + "merge.csv")
    height = np.genfromtxt(path + "height.csv")
    order = np.genfromtxt(path + "order.csv")
    
    return merge, height, order


def plot_dendr(dendr, leaf_len, lc, lw, text_offset, fs, draw_labels=True,
               leaf_points=None, linecols=None):
    merge, height, order = dendr
    
    order_idx = np.argsort(order)
    mid_lines = []
    
    leafs_drawn = 0
    
    for i in range(merge.shape[0]):
        if linecols is not None:
            lc = linecols[i]
        if sum(merge[i, :] < 0) == 2:
            x1_idx = int(abs(merge[i, 0]) - 1)
            x2_idx = int(abs(merge[i, 1]) - 1)
        
            x1 = order_idx[[x1_idx, x2_idx]]
            if leaf_points is not None:
                x1 = np.array([leaf_points[x1[0]], leaf_points[x1[1]]])
            y1 = [height[i], height[i]]
            plt.plot(x1, y1, c=lc, lw=lw)
            
            x2 = order_idx[[x1_idx, x1_idx]]
            if leaf_points is not None:
                x2 = np.array([leaf_points[x2[0]], leaf_points[x2[1]]])
            y2 = [height[i], 0] #height[i] - leaf_len]
            plt.plot(x2, y2, c=lc, lw=lw)
            
            x3 = order_idx[[x2_idx, x2_idx]]
            if leaf_points is not None:
                x3 = np.array([leaf_points[x3[0]], leaf_points[x3[1]]])
            y3 = [height[i], 0] #height[i] - leaf_len]
            plt.plot(x3, y3, c=lc, lw=lw)
            
            l1 = "$" + f"{int(abs(merge[i, 0]))}" + "$"
            l2 = "$" + f"{int(abs(merge[i, 1]))}" + "$"
            
            if draw_labels:
                plt.text(x2[0], y2[1] + text_offset, l1, fontsize=fs,
                         horizontalalignment='center')
                plt.text(x3[0], y3[1] + text_offset, l2, fontsize=fs, 
                         horizontalalignment='center')
            
            mid_lines.append(x1.mean())
            leafs_drawn += 2
        
        if sum(merge[i, :] < 0) == 0:
            x1_idx = int(merge[i, 0] - 1)
            x2_idx = int(merge[i, 1] - 1)
            
            x1 = [mid_lines[x1_idx], mid_lines[x2_idx]]
            y1 = [height[i], height[i]]
            plt.plot(x1, y1, c=lc, lw=lw)
            
            x2 = [mid_lines[x1_idx]]*2
            y2 = [height[x1_idx], height[i]]
            plt.plot(x2, y2, c=lc, lw=lw)                
            
            x3 = [mid_lines[x2_idx]]*2
            y3 = [height[x2_idx], height[i]]
            plt.plot(x3, y3, c=lc, lw=lw)
            
            mid_lines.append(sum(x1)/len(x1))
        
        if sum(merge[i, :] < 0) == 1:
            pos_idx = np.argwhere(merge[i, :] > 0)[0][0]
            neg_idx = np.argwhere(merge[i, :] < 0)[0][0]
            
            x1_idx = int(abs(merge[i, neg_idx]) - 1)
            x2_idx = int(merge[i, pos_idx] - 1)
            
            x1 = [order_idx[x1_idx], mid_lines[x2_idx]]
            if leaf_points is not None:
                x1 = np.array([leaf_points[x1[0]], x1[1]])
            y1 = [height[i], height[i]]
            plt.plot(x1, y1, c=lc, lw=lw)
            
            x2 = [mid_lines[x2_idx]]*2
            y2 = [height[x2_idx], height[i]]
            plt.plot(x2, y2, c=lc, lw=lw)
            
            x3 = order_idx[[x1_idx, x1_idx]]
            if leaf_points is not None:
                x3 = np.array([leaf_points[x3[0]], leaf_points[x3[1]]])
            y3 = [height[i], 0] #height[i] - leaf_len]
            plt.plot(x3, y3, c=lc, lw=lw)
            
            l1 = "$" + f"{int(abs(merge[i, neg_idx]))}" + "$"

            if draw_labels:
                plt.text(x3[0], y3[1] + text_offset, l1, fontsize=fs,
                         horizontalalignment='center')
            
            mid_lines.append(sum(x1)/len(x1))
            leafs_drawn += 1
    
    if leafs_drawn < len(order):
        leafs_left = order[leafs_drawn:]
        
        for i in range(leafs_left.size):
            x1 = [leafs_drawn + i, leafs_drawn + i]
            y1 = [leaf_len, 0]
            
            l1 = "$" + f"{int(leafs_left[i])}" + "$"
            
            if draw_labels:
                plt.text(x1[0], y1[1] + text_offset, l1, fontsize=fs,
                         horizontalalignment='center')
            plt.plot(x1, y1, c=lc, lw=lw)
    
    return None
