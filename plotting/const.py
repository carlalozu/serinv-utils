import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update({
    'axes.labelsize': 16,  # X and Y label font size
    'axes.titlesize': 16,  # Title font size
    'xtick.labelsize': 14,  # X-tick labels font size
    'ytick.labelsize': 14,  # Y-tick labels font size
    'legend.fontsize': 14,   # Legend font size
    'lines.linewidth': 3,
})

    
PEAK_PERFORMANCE = {
    'fritz': 2765, # GFLOPS
    'alex':  37400,  # 37.4 TFLOPS
}