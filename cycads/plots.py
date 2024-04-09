import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns

plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['figure.dpi'] = 300
palette = dict(A="tab:red", T="tab:green", C="tab:blue", G="tab:purple", S='tab:black')
figure_kw = dict(figsize = (5, 4), constrained_layout=True)
hist_kw = dict(facecolor='tab:blue', edgecolor='k', linewidth=0.5)
grid_kw = dict(color='k', alpha=0.1)
title_kw = dict(fontsize=16)
labelsize_kw = dict(fontsize=14)
ticksize_kw = dict(labelsize=12)

def post_process_ax(ax):
    ax.spines[['right', 'top']].set_visible(False)

