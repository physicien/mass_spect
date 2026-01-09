#!/usr/bin/python3

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

#GLOBAL CONSTANTS
DATAPATH0 = "./data/MALDI/C130-D5h_MALDI.csv"
DATAPATH1 = "./data/MALDI/C130_isotopic_distribution.csv"
df0 = pd.read_csv(DATAPATH0)
df1 = pd.read_csv(DATAPATH1)

#plot parameters
#xmin=0.0053
#xmax=0.0173
ymin=0.00
#ymax=3.23
inset_h=0.07
inset_w=0.20
acs_w=390
acs_h=100
figprefix="graph_SI_"
minor_ticks=True
show_grid=False
save_plot=False
show_plot=True
acs_format=True
x_label=r'm/z'
y_label=r'Relative Abundance'

#Functions
def fig_name(path,prefix):
    basename = os.path.basename(path)
    filename, file_extension = os.path.splitext(basename)
    return prefix + filename

def inset_range(df):
    xmin=df["mz"].min()
    xmax=df["mz"].max()
    return xmin-1,xmax+1

def pos_ratio(pos,xmin,xmax):
    return (pos-xmin)/(xmax-xmin)

def pos_inset(df):
    x_data = df["mz"]
    peak_val = df["abundance"].argmax()
    peak_pos = x_data[peak_val]
    xmin=x_data.min()
    xmax=x_data.max()
    mid=xmin+(xmax-xmin)/2
    ratio=pos_ratio(peak_pos,xmin,xmax)
    if peak_pos <= mid: 
        pos=ratio+0.07
    else:
        pos=ratio-0.22
    return pos

#Data processing
main_x_data = df0["mz"]
main_y_data = df0["abundance"]
calc_x_data = df1["mz"]
calc_y_data = df1["abundance"]

pos_y_max = df0['abundance'].argmax()
mass_peak = main_x_data[pos_y_max]

x_inset=pos_inset(df0)

#set font size and plot parameters
if acs_format:
    width=acs_w/72
    height=acs_h/72
    s1,s2,s3,FS=5,6,7,7
    ms=8
    lw=0.50
    capsize=2
    major_tick,minor_tick=4,2
else:
    width=10.0
    height=3.0
    s1,s2,s3,FS=14,16,18,18
    ms=10
    lw=1.25
    capsize=4
    major_tick,minor_tick=8,4

plt.rcParams['figure.figsize'] = [width, height]
plt.rc('font', size=s2)              # Set the default text font size
plt.rc('axes', titlesize=s2)         # Set the axes title font size
plt.rc('axes', labelsize=s2)         # Set the axes labels font size
plt.rc('xtick', labelsize=s1)        # Set the font size for x tick labels
plt.rc('ytick', labelsize=s1)        # Set the font size for y tick labels
plt.rc('legend', fontsize=s1)        # Set the legend font size
plt.rc('figure', titlesize=s3)       # Set the font size of the figure title

#prepare plot
palette = sns.color_palette("colorblind")
sns.set_palette(palette)
fig, main_ax = plt.subplots()

main_ax.plot(main_x_data,main_y_data,color='k',lw=lw)

exp_inset_ax = fig.add_axes([x_inset,.65,inset_h,inset_w])
exp_inset_ax.plot(main_x_data,main_y_data,color='k',lw=lw)

calc_inset_ax = fig.add_axes([x_inset+.08,.65,inset_h,inset_w])
calc_inset_ax.bar(calc_x_data,calc_y_data,width=0.5,color='k')

#annotations
main_ax.annotate(f'{mass_peak:.0f}',
        xy=(main_x_data[pos_y_max],main_y_data[pos_y_max]),
        xycoords='data',fontsize=FS,
        xytext=(0,4),textcoords='offset points',
        ha='center')
exp_inset_ax.set_title('Expt.',y=-0.7)
calc_inset_ax.set_title('Calc.',y=-0.7)

#legend
#main_ax.legend(loc='lower right')

#set xlim
main_ax.set_xlim(min(main_x_data),max(main_x_data))
xmin_exp,xmax_exp = inset_range(df1)
exp_inset_ax.set_xlim(xmin_exp,xmax_exp)
calc_inset_ax.set_xlim(xmin_exp,xmax_exp)

#set ylim
ymax = max(main_y_data)
main_ax.set_ylim(ymin,ymax)

#set labels
main_ax.set_xlabel(x_label)
main_ax.set_ylabel(y_label)

#show minor ticks
if minor_ticks:
    main_ax.minorticks_on()

#hide axis
main_ax.get_yaxis().set_visible(False)
exp_inset_ax.axis('off')
calc_inset_ax.axis('off')

#hide spines
main_ax.spines[['top','left','right']].set_visible(False)

#tick parameters
main_ax.tick_params(which='major',length=major_tick)
main_ax.tick_params(which='minor',length=minor_tick)

#show grid
if show_grid:
    main_ax.grid(True,which='major',axis='x',color='black',linestyle='dotted',
            linewidth=0.5)

#tight layout
plt.tight_layout()

if save_plot:
    figname = fig_name(DATAPATH0,figprefix)
    plt.savefig(f"{figname}.svg")

if show_plot:
    plt.show()
