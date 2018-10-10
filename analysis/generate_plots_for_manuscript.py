# -*- coding: utf-8 -*-

import os
import copy
import pickle

import numpy
import scipy.stats
from scipy.stats import linregress, mannwhitneyu, spearmanr, ttest_ind
from scipy.stats.mstats import theilslopes, kendalltau
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable


# # # # #
# CONSTANTS

# Variable by which the data was sorted, e.g. 'capacity' or 'apathy'.
SORTBY = 'capacity'

# Plotting settings.
# 'Points of view' colour palette (Wong, 2011, Nature Methods)
COLS = { \
    'orange':       (230, 159,   0), \
    'sky':          ( 86, 180, 233), \
    'green':        (  0, 158, 115), \
    'yellow':       (240, 228,  66), \
    'blue':         (  0, 114, 178), \
    'vermillion':   (213,  94,   0), \
    'purple':       (204, 121, 167), \
    }
# Transform colours to Matplotlib-compatible colours.
for key in COLS.keys():
    COLS[key] = (COLS[key][0]/255.0, COLS[key][1]/255.0, COLS[key][2]/255.0)
# Colours for different aspects of plotting.
PLOTCOLS = { \
    'samples':          COLS['orange'], \
    'regression':       COLS['vermillion'], \
    'high':             COLS['purple'], \
    'low':              COLS['sky'], \
    'pairline':         '#000000', \
    }
# Font sizes for different elements.
FONTSIZE = { \
    'title':            30, \
    'legend':           14, \
    'bar':              14, \
    'label':            20, \
    'ticklabels':       14, \
    'annotation':       14, \
    }

# Files and folders.
# Folder that contains this script file.
CODEDIR = os.path.dirname(os.path.abspath(__file__))
# Highest-level directory that we need.
DIR = os.path.dirname(CODEDIR)
# Directory for pre-processed data (intermediate steps with numbers that are
# generated through different analysis steps.)
PPDIR = os.path.join(DIR, 'processed_data_sorted-by-%s' % (SORTBY))
if not os.path.isdir(PPDIR):
    raise Exception("Could not find pre-processed data folder at '%s'" % (PPDIR))
# Directory for generated output, such as graphs.
OUTDIR = os.path.join(DIR, 'output_sorted-by-%s' % (SORTBY))
if not os.path.isdir(OUTDIR):
    os.mkdir(OUTDIR)
POUTDIR = os.path.join(OUTDIR, 'pretty_plots')
if not os.path.isdir(POUTDIR):
    os.mkdir(POUTDIR)


# # # # #
# LOAD DATA

# Load the pre-processed data form pickle files.
with open(os.path.join(PPDIR, 'idata.pickle'), 'rb') as f:
    idata = pickle.load(f)
with open(os.path.join(PPDIR, 'pdata.pickle'), 'rb') as f:
    pdata = pickle.load(f)
with open(os.path.join(PPDIR, 'sesdata.pickle'), 'rb') as f:
    sesdata = pickle.load(f)
with open(os.path.join(PPDIR, 'wmdata.pickle'), 'rb') as f:
    wmdata = pickle.load(f)
with open(os.path.join(PPDIR, 'qdata.pickle'), 'rb') as f:
    qdata = pickle.load(f)

# Count the number of trials in which one participant claimed no items.
nonclaimed = pdata['cont']['d_nclaimed'] == 8
print("A participant claimed no items in a total of %d trials from %d different sessions." % \
    (numpy.sum(nonclaimed.astype(int)), \
    numpy.sum( (numpy.sum(nonclaimed.astype(int), axis=1) > 0).astype(int))))


# # # # #
# PLOTTING

def draw_complot(x, y, xname, yname, varname, save_to_file=True, ax=None, \
    stats_title=False, stats_annotate=True, stats_parametric=True):
    
    # Compute descriptive statistics.
    m = [numpy.nanmean(x), numpy.nanmean(y)]
    sd = [numpy.nanstd(x), numpy.nanstd(y)]
    sem = [sd[0]/numpy.sqrt(len(x)), sd[1]/numpy.sqrt(len(y))]
    
    # Statistical test to determine whether the groups are different.
    t, tp = ttest_ind(x, y)
    u, up = mannwhitneyu(x, y)
    
    # Create a new plot.
    if ax  is None:
        fig, ax = pyplot.subplots(nrows=1, ncols=1)
    xloc = [0.5, 1.0]
    boxloc = [0.3, 1.2]
    lbl = [xname, yname]
    # Plot lines to connect the pairs.
    for i in range(len(x)):
        ax.plot(xloc, [x[i],y[i]], '-', color=PLOTCOLS['pairline'], alpha=0.1)
    # Plot the groups.
    for i, _y in enumerate([x,y]):
        _x = numpy.ones(len(x)) * xloc[i]
        # Draw the individual points.
        ax.plot(_x, _y, 'o', label=lbl[i], color=PLOTCOLS[lbl[i]], alpha=0.3)
        # Draw an error bar.
#        ax.errorbar(xloc[i], m[i], yerr=sem[i], fmt='-', capsize=10, \
#            ecolor=PLOTCOLS[lbl[i]], elinewidth=5, markeredgewidth=5)
        # Draw a violin plot.
        vp = ax.violinplot(_y, vert=True, \
            showmeans=False, showmedians=False, showextrema=False, \
            positions=[boxloc[i]], widths=0.15)
        vp['bodies'][0].set_facecolor(PLOTCOLS[lbl[i]])
        vp['bodies'][0].set_edgecolor(PLOTCOLS[lbl[i]])
#        vp['cbars'].set_color(PLOTCOLS[lbl[i]])
#        vp['cbars'].set_linewidth(3)
#        vp['cmins'].set_color(PLOTCOLS[lbl[i]])
#        vp['cmaxes'].set_color(PLOTCOLS[lbl[i]])
#        vp['cmedians'].set_color(PLOTCOLS[lbl[i]])
#        vp['cmedians'].set_linewidth(3)
        # Draw a boxplot within the violin plot.
        bp = ax.boxplot(_y, notch=False, sym='ko', vert=True, \
            positions=[boxloc[i]], widths=0.03, \
            patch_artist=True, \
            boxprops={'edgecolor':'black', 'facecolor':'white', 'linewidth':3}, \
            capprops={'color':'black', 'linewidth':3}, \
            whiskerprops={'color':'black', 'linewidth':3}, \
            flierprops={'color':'black', 'linewidth':3}, \
            medianprops={'color':PLOTCOLS[lbl[i]], 'linewidth':3}, \
            )
    # Annotate the statistics.
    if stats_annotate:
        if stats_parametric:
            stat_str = r"t = %.2f" % (t)
            if tp < 0.001:
                p_str = r"p<0.001"
            else:
                p_str = r"p=%.3f" % (tp)
        else:
            stat_str = r"u = %.2f" % (u)
            if up < 0.001:
                p_str = r"p<0.001"
            else:
                p_str = r"p=%.3f" % (up)
        xpos = xloc[0] + 0.02*(xloc[1]-xloc[0])
        ypos = max(numpy.nanmax(x), numpy.nanmax(y)) * 1.03
        ax.annotate(r"$%s, %s$" % (stat_str, p_str), (xpos,ypos), \
            fontsize=FONTSIZE['annotation'])
    # Finish the plot.
    ax.set_ylim(top=max(numpy.nanmax(x), numpy.nanmax(y)) * 1.1)
    ax.set_xlim([0, 1.5])
    ax.set_xticks(xloc)
    ax.set_xticklabels([xname.capitalize(), yname.capitalize()], \
        fontsize=FONTSIZE['label'])
    ax.set_ylabel(varname, fontsize=FONTSIZE['ticklabels'])
    if stats_title:
        ax.set_title("t=%.2f, p=%.3f; U=%.2f, p=%.3f" % (t, tp, u, up))
    # Save the plot.
    if save_to_file:
        fig.savefig(os.path.join(OUTDIR, "complot_%s.png" % (varname)))
    if ax is None:
        pyplot.close(fig)
    

def draw_corplot(x, y, xname, yname, add_robust=False, save_to_file=True, \
    ax=None, stats_title=True, stats_legend=False, customcol=None, \
    legendprefix=''):
    # Choose the right colour for the plot.
    if customcol is None:
        regress_col = PLOTCOLS['regression']
        sample_col = PLOTCOLS['samples']
        sample_alpha = 0.75
    else:
        regress_col = customcol
        sample_col = customcol
        sample_alpha = 0.75
    # Create a new plot.
    if ax is None:
        fig, ax = pyplot.subplots(nrows=1, ncols=1)
    # Plot a scatter plot of the x and y values.
    ax.plot(x, y, 'o', color=sample_col, alpha=sample_alpha)
    # Plot the regression line.
    if add_robust:
        # Perform a linear regression.
        slope, intercept, lo_slope, up_slope = theilslopes(y, x, alpha=0.95)
        # Plot the regression line.
        x_pred = numpy.array([numpy.min(x), numpy.max(x)])
        y_pred = slope * x_pred + intercept
        y_lo = lo_slope * x_pred + intercept
        y_up = up_slope * x_pred + intercept
        ax.plot(x_pred, y_pred, '-', color=regress_col)
        ax.fill_between(x_pred, y_lo, y_up, linewidth=3, alpha=0.2, \
            color=PLOTCOLS['regression'])
    # Perform a linear regression.
    model = linregress(x, y)
    try:
        r = model.rvalue
        p = model.pvalue
        slope = model.slope
        intercept = model.intercept
    except:
        slope, intercept, r, p, stderr = model
    # Perform a Spearman correlation.
    spearman = spearmanr(x, y)
    try:
        spearman_rho = spearman.correlation
        spearman_p = spearman.pvalue
    except:
        spearman_rho, spearman_p = spearman
    # Compute Kendall's Tau.
    kendall = kendalltau(x, y)
    try:
        kendall_tau = kendall.correlation
        kendall_p = kendall.pvalue
    except:
        kendall_tau, kendall_p = kendall
    # Set the regression line's label.
    if stats_legend:
        # Uncomment if you'd like to see both parametric and non-parametric
        # test results.
        #lbl = r"$R=%.2f, p=%.2f$" % (r, p)
        #lbl = lbl + "\n" + r"$\tau=%.2f, p=%.2f$" % (kendall_tau, kendall_p)
        # Show Kendall's tau, as we're using a lowish N.
        if kendall_p < 0.001:
            kendall_pstr = r"p<0.001"
        else:
            kendall_pstr = r"p=%.3f" % (kendall_p)
        lbl = r"%s$\tau=%.2f, %s$" % (legendprefix, kendall_tau, kendall_pstr)
    else:
        lbl = None
    # Plot the regression line.
    x_pred = numpy.array([numpy.min(x), numpy.max(x)])
    y_pred = slope * x_pred + intercept
    ax.plot(x_pred, y_pred, '-', color=regress_col, linewidth=3, label=lbl)
    # Finish the plot.
    ax.set_xlabel(xname.capitalize(), fontsize=FONTSIZE['label'])
    ax.set_ylabel(yname.capitalize(), fontsize=FONTSIZE['label'])
    if stats_title:
        ax.set_title("R=%.2f, p=%.3f; Rho=%.2f, p=%.3f; Tau=%.3f, p=%.3f" % \
            (r, p, spearman_rho, spearman_p, kendall_tau, kendall_p))
    if stats_legend:
        ax.legend(loc="best", fontsize=FONTSIZE['legend'])
    # Save the plot.
    if save_to_file:
        fig.savefig(os.path.join(OUTDIR, "corplot_%sx%s.png" % (xname, yname)))
    if ax is None:
        pyplot.close(fig)


# # # # #
# FIGURE 4

# Individual error distributions, and intra-dyadic differences in capacity,
# in-game error, and reward contribution.
fig, axes = pyplot.subplots(nrows=2, ncols=2, figsize=(16.0, 12.0), dpi=300.0)
fig.subplots_adjust(left=0.07, bottom=0.10, right=0.95, top=0.95,
    wspace=0.25, hspace=0.25)
axes = [axes[0,0], axes[0,1], axes[1,0], axes[1,1]]

# Plot all participants' individual error.
# Choose the first axis.
ax = axes[0]
# Create a colour map to colour participant's performance with.
cmap = matplotlib.cm.get_cmap("inferno")
vmin = 1 #numpy.nanmin(idata['capacity'])
vmax = 6 #numpy.nanmax(idata['capacity'])
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
# Sort by capacity.
ppnames = numpy.array(idata['ppnames'])
ppnames = ppnames[numpy.argsort(idata['capacity'])]
# Loop through all participants.
for i, ppname in enumerate(ppnames):
    # Find the participant's capacity.
    pi = idata['ppnames'].index(ppname)
    c = idata['capacity'][pi]
    # Label the lowest, middle, and highest.
    if i == 0:
        lbl = r"Lowest: $\kappa$=%.2f" % (c)
    elif i == len(ppnames)//2:
        lbl = r"Middle: $\kappa$=%.2f" % (c)
    elif i == len(ppnames)-1:
        lbl = r"Highest: $\kappa$=%.2f" % (c)
    else:
        lbl = None
    # Select the appropriate colour.
    col = cmap(norm(c))
    # Compute the error histogram.
    notnan = numpy.isnan(wmdata[ppname]['E']) == False
    hist, bin_edges = numpy.histogram(wmdata[ppname]['E'][notnan], \
        bins=11, density=True)
    # Plot the histogram.
    bin_centres = bin_edges[:-1] + numpy.diff(bin_edges)/2.0
    ax.plot(bin_centres, hist, '-', linewidth=3, color=col, alpha=0.4, \
        label=lbl)
# Clarify the plot.
ax.legend(loc="upper right", fontsize=FONTSIZE['legend'])
ax.set_xlim([-90, 90])
ax.set_xlabel("Recall error in individual assessment (deg)", fontsize=FONTSIZE['label'])
ax.set_ylim(bottom=0, top=0.05)
ax.set_ylabel("Probability density", fontsize=FONTSIZE['label'])
# Add a colour bar.
divider = make_axes_locatable(ax)
bax = divider.append_axes("right", size="5%", pad=0.05)
cbar = matplotlib.colorbar.ColorbarBase(bax, cmap=cmap, norm=norm, \
    ticks=range(vmin, vmax+1, 1), orientation='vertical')
cbar.set_label(r"Short-term memory capacity $\kappa$", fontsize=FONTSIZE['bar'])

# CAPACITY, ERROR, AND REWARD
for i, varname in enumerate(["capacity", "abserr", "earned"]):
    # The first ax is used by a previous plot, so we need to update the
    # iterator.
    i += 1
    # Choose the right data.
    x = pdata['player_low'][varname]
    y = pdata['player_high'][varname]
    xname = 'low'
    yname = 'high'
    # Draw a plot.
    draw_complot(x, y, xname, yname, varname, save_to_file=False, ax=axes[i], \
        stats_title=False)
    # Clarify the axis labels.
    axes[i].set_xlabel(SORTBY.capitalize(), fontsize=FONTSIZE['label'])
    if varname == 'capacity':
        ylbl = r"Capacity $\kappa$ (bits)"
    elif varname == 'abserr':
        ylbl = "Absolute error in game (deg)"
    elif varname == 'earned':
        ylbl = "Reward contributed (points)"
    elif varname == 'nclaimed':
        ylbl = "Number of claimed items"
    axes[i].set_ylabel(ylbl, fontsize=FONTSIZE['label'])
fig.savefig(os.path.join(POUTDIR, "fig_04.png"))
pyplot.close(fig)


# # # # #
# FIGURE 5

# NCLAIMED AND CAPACITY
fig, axes = pyplot.subplots(nrows=1, ncols=2, figsize=(16.0, 6.0), dpi=300.0)
fig.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.93,
    wspace=0.2, hspace=0.2)
for i, varname in enumerate(['nclaimed', 'ratios']):
    # Create a regular comparison plot for nclaimed.
    if varname == 'nclaimed':
        # Choose the right data.
        x = pdata['player_low'][varname]
        y = pdata['player_high'][varname]
        xname = 'low'
        yname = 'high'
        # Draw a plot.
        draw_complot(x, y, xname, yname, varname, save_to_file=False, ax=axes[i], \
            stats_title=False)
        # Clarify the axis labels.
        axes[i].set_xlabel(SORTBY.capitalize(), fontsize=FONTSIZE['label'])
        ylbl = "Number of claimed items"
        axes[i].set_ylabel(ylbl, fontsize=FONTSIZE['label'])

    # Create a custom plot for the ratios.
    elif varname == 'ratios':
        # Compute the capacity-ratios for each pair.
        k_ratio = pdata['player_high']['capacity'] / pdata['player_low']['capacity']
        claim_ratio = pdata['player_high']['nclaimed'] / pdata['player_low']['nclaimed']
        # Convenience renaming of the data.
        x = k_ratio
        y = claim_ratio
        xname = "Player capacity ratio (high/low)"
        yname = "Player claim ratio (high/low)"
        # Create two lines: one for equal and one for proportional division.
        max_val = max(numpy.nanmax(x), numpy.nanmax(y)) * 1.02
        fair_line = numpy.arange(0, max_val+0.1, 0.1)
        equal_line = numpy.ones(fair_line.shape)
        # Plot the division lines.
        axes[i].plot(fair_line, equal_line, '--', lw=3, color='#000000', label="Equal division")
        axes[i].plot(fair_line, fair_line, ':', lw=3, color='#000000', label="Proportional division")
        axes[i].fill_between(fair_line, fair_line, numpy.ones(fair_line.shape)*max_val, \
            color='#000000', alpha=0.1)
        axes[i].annotate("Higher-capacity player over-claims", (0.6, 0.9*max_val), \
            fontsize=FONTSIZE['annotation'], color=PLOTCOLS['high'], alpha=1.0)
        axes[i].annotate("Lower-capacity player over-claims", (0.53*max_val, 0.5), \
            fontsize=FONTSIZE['annotation'], color=PLOTCOLS['low'], alpha=1.0)
        # Plot the actual data.
        draw_corplot(x, y, xname, yname, add_robust=False, save_to_file=False, \
            ax=axes[i], stats_title=False, stats_legend=True)
        # Reset the legend.
        axes[i].legend(loc="center right", fontsize=FONTSIZE['legend'])
        # Set the limits on the axes.
        axes[i].set_xlim([0.5, max_val])
        axes[i].set_ylim([0.4, max_val])
# Save and close the figure.
fig.savefig(os.path.join(POUTDIR, "fig_05.png"))
pyplot.close(fig)


# # # # #
# FIGURE 6

# RATIOS AND COLLABORATION
fig, axes = pyplot.subplots(nrows=1, ncols=3, figsize=(24.0, 6.0), dpi=300.0)
fig.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.85,
    wspace=0.2, hspace=0.2)
for i, varname in enumerate(['good', 'fair', 'likeable']):
    # Get the collaboration scores.
    y1 = pdata['player_high'][varname]
    y2 = pdata['player_low'][varname]
    min_y = 0 # min(numpy.nanmin(y1), numpy.nanmin(y2)) * 1.1
    max_y = 1.4 #min(numpy.nanmax(y1), numpy.nanmax(y2)) * 1.1
    yname = "Individual player rating"
    # Compute the capacity-ratios for each pair.
    k_ratio = pdata['player_high']['capacity'] / pdata['player_low']['capacity']
    claim_ratio = pdata['player_high']['nclaimed'] / pdata['player_low']['nclaimed']
    x = claim_ratio - k_ratio
    min_x = numpy.nanmin(x) * 1.1
    max_x = numpy.nanmax(x) * 1.1
    xname = "Claim-capacity disproportionality"
    # Plot a proportional division line.
    prop_y = numpy.arange(min_y, max_y+0.1, 0.1)
    prop_division = numpy.zeros(prop_y.shape)
    axes[i].plot(prop_division, prop_y, ':', lw=3, label="Proportional division", \
        color='#000000')
    fill_x = numpy.arange(min_x, 0, 0.01)
    fill_ymin = numpy.ones(fill_x.shape) * min_y
    fill_ymax = numpy.ones(fill_x.shape) * max_y
    axes[i].fill_between(fill_x, fill_ymin, fill_ymax, color='#000000', alpha=0.1)
    axes[i].annotate("Higher-capacity player\nover-claims", (0.1, 0.05), \
        fontsize=FONTSIZE['annotation'], color=PLOTCOLS['high'], alpha=1.0)
    axes[i].annotate("Lower-capacity player\nover-claims", (min_x+0.1, 0.05), \
        fontsize=FONTSIZE['annotation'], color=PLOTCOLS['low'], alpha=1.0)
    # Plot a correlation between the division disproportionality and players'
    # collaboration assessments.
    draw_corplot(x, y1, xname, yname, add_robust=False, save_to_file=False, \
        ax=axes[i], stats_title=False, stats_legend=True, \
        customcol=PLOTCOLS['high'], legendprefix='Higher-capacity player: ')
    draw_corplot(x, y2, xname, yname, add_robust=False, save_to_file=False, \
        ax=axes[i], stats_title=False, stats_legend=True, \
        customcol=PLOTCOLS['low'], legendprefix='Lower-capacity player: ')
    # Set an axis title.
    if varname == 'good':
        tlbl = "Collaboration was 'good'"
    elif varname == 'fair':
        tlbl = "Collaborator was 'fair'"
    elif varname == 'likeable':
        tlbl = "Collaborator was 'likeable'"
    tit = axes[i].set_title(tlbl, fontsize=FONTSIZE['title'])
    tit.set_position([.5, 1.05])
    # Reset the legend.
    axes[i].legend(loc="upper left", fontsize=FONTSIZE['legend'])
    # Set the limits on the axes.
    axes[i].set_xlim([min_x, max_x])
    axes[i].set_ylim([min_y, max_y])
# Save and close the figure.
fig.savefig(os.path.join(POUTDIR, "fig_06.png"))
pyplot.close(fig)


# # # # #
# FIGURE 7

# INTER-DISTANCE AND COLLABORATION
fig, axes = pyplot.subplots(nrows=2, ncols=3, figsize=(24.0, 12.0), dpi=300.0)
fig.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.93,
    wspace=0.2, hspace=0.2)
for j, predictorname in enumerate(['d_avg', 'd_var']):
    # Get the predictor values.
    x = pdata[predictorname]
    min_x = numpy.nanmin(x) * 0.9
    max_x = numpy.nanmax(x) * 1.1
    if predictorname == 'd_avg':
        xname = "Average inter-player distance"
    elif predictorname == 'd_var':
        xname = "Variability in inter-player distance"
    for i, varname in enumerate(['good', 'fair', 'likeable']):
        # Get the collaboration scores.
        y1 = pdata['player_high'][varname]
        y2 = pdata['player_low'][varname]
        yname = "Individual player rating"
        # Plot a correlation between the division disproportionality and players'
        # collaboration assessments.
        draw_corplot(x, y1, xname, yname, add_robust=False, save_to_file=False, \
            ax=axes[j,i], stats_title=False, stats_legend=True, \
            customcol=PLOTCOLS['high'], legendprefix='Higher-capacity player: ')
        draw_corplot(x, y2, xname, yname, add_robust=False, save_to_file=False, \
            ax=axes[j,i], stats_title=False, stats_legend=True, \
            customcol=PLOTCOLS['low'], legendprefix='Lower-capacity player: ')
        # Set an axis title.
        if j == 0:
            if varname == 'good':
                tlbl = "Collaboration was 'good'"
            elif varname == 'fair':
                tlbl = "Collaborator was 'fair'"
            elif varname == 'likeable':
                tlbl = "Collaborator was 'likeable'"
            tit = axes[j,i].set_title(tlbl, fontsize=FONTSIZE['title'])
            tit.set_position([.5, 1.05])
        # Reset the legend.
        axes[j,i].legend(loc="upper left", fontsize=FONTSIZE['legend'])
        # Set the limits on the axes.
        axes[j,i].set_xlim([min_x, max_x])
        axes[j,i].set_ylim([0, 1.3])
# Save and close the figure.
fig.savefig(os.path.join(POUTDIR, "fig_07.png"))
pyplot.close(fig)


# # # # #
# TEXT DATA FOR JASP

# We need the data on each participants' questionnaire outcomes, as well as
# the number of stimuli that they claimed, and their performance.
jasp_vars = ['nclaimed', 'capacity', \
    'apathy', 'apathy_ba', 'apathy_sm', 'apathy_es', \
    'B5_extraversion', 'B5_agreeableness', 'B5_conscientiousness', 'B5_emotional_stability', 'B5_openness', \
    ]

with open(os.path.join(OUTDIR, 'values_for_jasp.csv'), 'w') as f:
    header = []
    for playertype in ['high', 'low', 'd']:
        for varname in jasp_vars:
            header.append("%s_%s" % (playertype, varname))
    f.write(','.join(header))
    for i in range(len(pdata["player_high"]["nclaimed"])):
        line = []
        for playertype in ['high', 'low', 'd']:
            for varname in jasp_vars:
                if playertype == 'd':
                    line.append(pdata["player_high"][varname][i] - pdata["player_low"][varname][i])
                else:
                    line.append(pdata["player_%s" % (playertype)][varname][i])
        f.write('\n' + ','.join(map(str, line)))


# # # # #
# OTHER FIGURES

# DISTANCE AND REWARD
fig, axes = pyplot.subplots(nrows=1, ncols=2, figsize=(16.0, 6.0), dpi=300.0)
fig.subplots_adjust(left=0.05, bottom=0.1, right=0.98, top=0.98,
    wspace=0.2, hspace=0.2)
for j, predictorname in enumerate(['d_avg', 'd_var']):
    # Get the predictor values.
    x = pdata[predictorname]
    min_x = numpy.nanmin(x) * 0.9
    max_x = numpy.nanmax(x) * 1.1
    if predictorname == 'd_avg':
        xname = "Average inter-player distance"
    elif predictorname == 'd_var':
        xname = "Variability in inter-player distance"
    # Get the reward scores.
    y1 = pdata['player_high']['earned']
    y2 = pdata['player_low']['earned']
    yname = "Individual player rating"
    # Plot a correlation between the division disproportionality and players'
    # collaboration assessments.
    draw_corplot(x, y1, xname, yname, add_robust=False, save_to_file=False, \
        ax=axes[j], stats_title=False, stats_legend=True, \
        customcol=PLOTCOLS['high'], legendprefix='Higher-capacity player: ')
    draw_corplot(x, y2, xname, yname, add_robust=False, save_to_file=False, \
        ax=axes[j], stats_title=False, stats_legend=True, \
        customcol=PLOTCOLS['low'], legendprefix='Lower-capacity player: ')
    draw_corplot(x, y1+y2, xname, "Total reward", add_robust=False, save_to_file=False, \
        ax=axes[j], stats_title=False, stats_legend=True, \
        legendprefix='Both players: ')
    # Reset the legend.
    axes[j].legend(loc="upper left", fontsize=FONTSIZE['legend'])
    # Set the limits on the axes.
    axes[j].set_xlim([min_x, max_x])
# Save and close the figure.
fig.savefig(os.path.join(POUTDIR, "fig_x0.png"))
pyplot.close(fig)


# NCLAIMED, AND COLLABORATION
for figtype in ['collaboration-differences', 'inter-player-agreement']:
    if figtype == 'collaboration-differences':
        fig, axes = pyplot.subplots(nrows=1, ncols=3, figsize=(24.0, 6.0), dpi=300.0)
        fig.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.93,
            wspace=0.2, hspace=0.2)
        axes = axes.reshape((1,len(axes)))
    elif figtype == 'inter-player-agreement':
        fig, axes = pyplot.subplots(nrows=2, ncols=3, figsize=(24.0, 12.0), dpi=300.0)
        fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95,
            wspace=0.2, hspace=0.2)
    for i, varname in enumerate(["good", "fair", "likeable"]):
        # Choose the right data.
        x = pdata['player_low'][varname]
        y = pdata['player_high'][varname]
        xname = 'low'
        yname = 'high'
        # Draw a comparison plot.
        draw_complot(x, y, xname, yname, varname, save_to_file=False, ax=axes[0,i], \
            stats_title=False)
        # Clarify the axis labels.
        axes[0,i].set_xlabel(SORTBY.capitalize(), fontsize=FONTSIZE['label'])
        axes[0,i].set_ylim([0,1.1])
        axes[0,i].set_ylabel("Player rating (0-1)", fontsize=FONTSIZE['label'])
        if varname == 'good':
            tlbl = "Collaboration was 'good'"
        elif varname == 'fair':
            tlbl = "Collaborator was 'fair'"
        elif varname == 'likeable':
            tlbl = "Collaborator was 'likeable'"
        axes[0,i].set_title(tlbl, fontsize=FONTSIZE['title'])
        # Draw a correlation plot.
        if figtype == 'inter-player-agreement':
            draw_corplot(x, y, xname, yname, add_robust=False, save_to_file=False, \
                ax=axes[1,i], stats_title=False, stats_legend=True)
            # Clarify the axis labels.
            axes[1,i].set_xlim([0,1.02])
            axes[1,i].set_xlabel("Lower-capacity player", fontsize=FONTSIZE['label'])
            axes[1,i].set_ylim([0,1.02])
            axes[1,i].set_ylabel("Higher-capacity player", fontsize=FONTSIZE['label'])
    fig.savefig(os.path.join(POUTDIR, "fig_x1_%s" % (figtype)))
    pyplot.close(fig)

# COLLABORATION, AND DISTANCE (AVERAGE AND DEVIATION)
fig, axes = pyplot.subplots(nrows=2, ncols=3, figsize=(24.0, 12.0), dpi=300.0)
fig.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95,
    wspace=0.25, hspace=0.25)
for i, collabname in enumerate(['good', 'fair', 'likeable']):
    for j, varname in enumerate(["d_avg", "d_var"]):
        # Choose the right data.
        x = pdata['player_high'][collabname] + pdata['player_low'][collabname]
        y = pdata[varname]
        xname = "Total '%s' collaboration rating" % (collabname)
        if varname == 'd_avg':
            yname = "Average inter-player distance"
        elif varname == 'd_var':
            yname = "Variability in inter-player distance"
        # Draw a correlation plot.
        draw_corplot(x, y, xname, yname, add_robust=False, save_to_file=False, \
            ax=axes[j,i], stats_title=False, stats_legend=True)
        # Limit the axes.
        axes[j,i].set_xlim([0,2.02])
        if varname == 'd_avg':
            axes[j,i].set_ylim([450,700])
        elif varname == 'd_var':
            axes[j,i].set_ylim([100,220])
        # Remove the axis labels from some of the plots.
        if i > 0:
            axes[j,i].set_ylabel("")
        if j == 0:
            axes[j,i].set_xlabel("")
        # Set the axis title.
        if j == 0:
            if collabname == 'good':
                tlbl = "Collaboration was 'good'"
            elif collabname == 'fair':
                tlbl = "Collaborator was 'fair'"
            elif collabname == 'likeable':
                tlbl = "Collaborator was 'likeable'"
            axes[j,i].set_title(tlbl, fontsize=FONTSIZE['title'])
fig.savefig(os.path.join(POUTDIR, "fig_x2.png"))
pyplot.close(fig)
