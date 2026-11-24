import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'figure.dpi': 300
})

base_dir = "path/to/project_root"
fig_dir = os.path.join(base_dir, "results/figures")
os.makedirs(fig_dir, exist_ok=True)

# Run same code to generate PDF and PNG directly on cluster
# ...
