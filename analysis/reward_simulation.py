# -*- coding: utf-8 -*-

import numpy
from scipy import stats
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable

import vonmises_toolbox as vt
from entropy import KL_divergence_continuous


# # # # #
# CONSTANTS

# Simulation properties.
NSTIM = 8
KMIN = 0.0
KMAX = 6.0
KSTEP = 0.5
KRANGE = numpy.arange(KMIN, KMAX+KSTEP, KSTEP)
# Reset the minimum to avoid a K of 0.
KMIN = 0.01
KRANGE[0] = KMIN

# Reward scale properties.
MINREW = 1.0
MAXREW = 100.0
REWSD = [5, 15, 25, 'lin']

# Plotting settings
FONTSIZE = { \
    'title':18, \
    'axtitle':16, \
    'bigaxis':14, \
    'bigticks':12, \
    'bigbartitle':14, \
    'bigbarticks':12, \
    'axis':12, \
    'ticks':10, \
    'bartitle':10, \
    'barticks':8, \
    }


# # # # #
# DKL LOOKUP TABLE

# Inversing the DKL computation in an exact way is non-triavial, but generating
# a lookup table is super easy. Guess which way we'll implement here?
print("Generating DKL lookup table.")

# Create a reasonably fine-grained space of possible DKL values. We do this
# by going through a reasonable number of potential standard deviations of
# error distributions, and by computing the corresponding DKL.
# NEVER SET sd_min TO 0! THIS IS IMPOSSIBLE, AND WILL RESULT IN A FASLE DKL OF 0!
sd_min = 0.001
sd_max = 3.0
sd_step = 0.001
sd_range = numpy.arange(sd_min, sd_max+sd_step, sd_step)
dkl_lookup = numpy.zeros(sd_range.shape, dtype=float) * numpy.NaN

# Create a probability density function for random guessing.
dx = 0.01
x = numpy.arange(-numpy.pi, numpy.pi+dx, dx)
q = numpy.ones(x.shape, dtype=float) / (2*numpy.pi)
# Loop through all SD values and compute the correspoding DKL.
for i, sd in enumerate(sd_range):
    # Compute the corresponding spread parameter, and then construct the
    # probability density function.
    p = vt.vonmisespdf(x, 0, vt.sd2kappa(sd))
    # For very low standard deviations, the Von Mises distribution generates
    # NaNs. For those, step in with a regular normal distribution.
    nans = numpy.isnan(p)
    if numpy.sum(nans.astype(int)) > 0:
        normpdf = stats.norm.pdf(x, loc=0, scale=sd)
        p[nans] = normpdf[nans]
    # Compute the KL divergence.
    dkl_lookup[i] = KL_divergence_continuous(p, q, dx, mode='bits')

# Compute the theoretical DKL minimum and maximum within this simulation, and
# make sure they are in the lookup table.
dkl_min = KMIN / float(NSTIM)
dkl_max = KMAX / 1.0

all_included = (numpy.nanmin(dkl_lookup) <= dkl_min) and (numpy.nanmax(dkl_lookup) >= dkl_max)
if not all_included:
    raise Exception( \
        "ERROR: The DKL lookup table does not include all values! Current min=%.3f and max=%.3f; required are min=%.3f and max=%.3f" \
        % (numpy.nanmin(dkl_lookup), numpy.nanmax(dkl_lookup), dkl_min, dkl_max))


# # # # #
# REWARD SIMULATION

print("Running expected reward simulation.")

# Construct the reward scales.
x = numpy.arange(0.0, 90.1, 1.0)
rewardscale = {}
for sd in REWSD:
    if type(sd) in [int, float]:
        rewardscale[sd] = stats.norm.pdf(x, loc=0, scale=sd)
        rewardscale[sd] *= 1.0 / rewardscale[sd].max()
    elif sd == 'lin':
        rewardscale[sd] = 1.0 - (numpy.array(x, dtype=float) / numpy.max(x))

# Construct a claim range.
claim_step = 1
claim_range = numpy.arange(0, NSTIM+claim_step, claim_step)

# Construct empty vectors to hold the potential reward outcomes.
optimal_claim_difference = {}
proportional_claim_difference = {}
rew_equal = {}
rew_prop = {}
rew_optimal = {}
for sd in rewardscale.keys():
    optimal_claim_difference[sd] = numpy.zeros((len(KRANGE),len(KRANGE)), \
        dtype=float) * numpy.NaN
    proportional_claim_difference[sd] = numpy.zeros((len(KRANGE),len(KRANGE)), \
        dtype=float) * numpy.NaN
    rew_equal[sd] = numpy.zeros((len(KRANGE),len(KRANGE)), dtype=float) \
        * numpy.NaN
    rew_prop[sd] = numpy.zeros((len(KRANGE),len(KRANGE)), dtype=float) \
        * numpy.NaN
    rew_optimal[sd] = numpy.zeros((len(KRANGE),len(KRANGE)), dtype=float) \
        * numpy.NaN

# Loop through all possible differences in capacity.
for i, k_one in enumerate(KRANGE):
    for j, k_two in enumerate(KRANGE):
        # Create a vector to contain all possible differences in the number of
        # claimed items, and one for all associated rewards.
        d_claimed = numpy.zeros(len(claim_range), dtype=float) * numpy.NaN
        t_reward = {}
        for sd in rewardscale.keys():
            t_reward[sd] = numpy.zeros(len(claim_range), dtype=float) * numpy.NaN
        # Loop through all possible numbers of items the first participant
        # could have claimed.
        for a, claimed_one in enumerate(claim_range):

            # Compute how much the second participant claimed.
            claimed_two = NSTIM - claimed_one
            
            # Store the claim difference.
            d_claimed[a] = claimed_one - claimed_two

            # Only compute the following if the participant actually claimed
            # an item. If they didn't, the reward is 0.
            if claimed_one > 0:
                # Compute the DKL per stimulus.
                dkl_one = k_one / float(claimed_one)
                # Find the closest error SD value.
                sd_one = sd_range[numpy.argmin(numpy.abs(dkl_lookup-dkl_one))]
                # Compute the expected error (divide by two as the DKL
                # estimates are based on error*2 distributions).
                err_one = numpy.rad2deg(stats.norm.expect(loc=0, scale=sd_one, lb=0)) / 2.0
                # Compute the reward for each individual.
                rew_one = {}
                for sd in rewardscale.keys():
                    rew_one[sd] = (rewardscale[sd][int(abs(err_one))] * (MAXREW - MINREW)) + MINREW
            else:
                rew_one = {}
                for sd in rewardscale.keys():
                    rew_one[sd] = 0.0
            
            # Only compute the following if the participant actually claimed
            # an item. If they didn't, the reward is 0.
            if claimed_two > 0:
                # Compute the DKL per stimulus.
                dkl_two = k_two / float(claimed_two)
                # Find the closest error SD value.
                sd_two = sd_range[numpy.argmin(numpy.abs(dkl_lookup-dkl_two))]
                # Compute the expected error (divide by two as the DKL
                # estimates are based on error*2 distributions).
                err_two = numpy.rad2deg(stats.norm.expect(loc=0, scale=sd_two, lb=0)) / 2.0
                # Compute the reward for each individual.
                rew_two = {}
                for sd in rewardscale.keys():
                    rew_two[sd] = (rewardscale[sd][int(abs(err_two))] * (MAXREW - MINREW)) + MINREW
            else:
                rew_two = {}
                for sd in rewardscale.keys():
                    rew_two[sd] = 0.0

            # Compute the combined reward.
            for sd in rewardscale.keys():
                t_reward[sd][a] = rew_one[sd] + rew_two[sd]
        
        for sd in rewardscale.keys():
            # Store the difference in claimed items with the highest expected
            # total reward.
            highest = numpy.argmax(t_reward[sd])
            optimal_claim_difference[sd][j,i] = d_claimed[highest]
            # Find the index at which the items are equally divided.
            equal_index = numpy.argmin(numpy.abs(d_claimed))
            # Store the rewards under equal division.
            rew_equal[sd][j,i] = t_reward[sd][equal_index]
            # Store the reward under proportional claimed item difference.
            # For this, we only count the rewards in situations where the better
            # player claimed more items. The exact ratio of stimulus division
            # depends on the ratio between the players' capacities.
            k_ratio = float(k_one) / float(k_one+k_two)
            ideal_dclaimed = k_ratio*NSTIM - (1.0-k_ratio)*NSTIM
            prop_index = numpy.argmin(numpy.abs(d_claimed-ideal_dclaimed))
            if (k_one > k_two) and (prop_index == equal_index):
                prop_index = equal_index+1
            elif (k_one < k_two) and (prop_index == equal_index):
                prop_index = equal_index-1
            rew_prop[sd][j,i] = t_reward[sd][prop_index]
            proportional_claim_difference[sd][j,i] = d_claimed[prop_index]
            # Store the reward under optimal claimed item difference. Note that
            # optimal here does not mean equal or proportional. Instead, an optimal
            # strategy could be to give the best player only 1 item, and the worst
            # player all other items. This does not always work, though, and is
            # hard to compute without knowing the exact payoff matrix.
            rew_optimal[sd][j,i] = t_reward[sd][highest]

# Compute the proportional increases from using proportional or optimal division.
prop_increase = {}
optimal_increase = {}
for sd in rewardscale.keys():
    prop_increase[sd] = (rew_prop[sd] - rew_equal[sd]) / rew_equal[sd]
    optimal_increase[sd] = (rew_optimal[sd] - rew_equal[sd]) / rew_equal[sd]


# # # # #
# PLOT OPTIMALS

print("Generating plots of optimal item division and rewards.")

# Loop through all reward schemes and create a plot for each.
for sd in rewardscale.keys():

    # Establish what the capacity ranged and tickslabels are.
    sel = KRANGE % 1 == 0
    sel[0] = True
    k_ticks = numpy.arange(0, len(KRANGE))[sel]
    k_ticklabels = map(str, KRANGE[sel].astype(int))
    
    # Create a new figure.
    fig = pyplot.figure(figsize=(8.0,8.0), dpi=600.0)
    fig.subplots_adjust(left=0.1, bottom=0.06, right=0.9, top=0.96,
        wspace=0.5, hspace=0.5)
    
    # Top plot: potential gains from proportional division.
    vmin = -0.1
    vmax = 0.1
    ax = pyplot.subplot2grid((3, 3), (0, 0), colspan=3, rowspan=2)
    #fig.suptitle("Potential gain from propotional labour division", fontsize=FONTSIZE['title'])
    ax.set_title("Benefit of proportional labour division", fontsize=FONTSIZE['title'])
    cax = ax.imshow(prop_increase[sd], vmin=vmin, vmax=vmax, cmap='coolwarm', \
        aspect='equal', interpolation='none', origin='lower')
    cbar = fig.colorbar(cax, ticks=[vmin, vmax], orientation='vertical')
    cbar.ax.set_yticklabels(map(str, (cbar.get_ticks()*100).astype(int)), \
        fontsize=FONTSIZE['bigbarticks'])
    cbar.ax.set_ylabel("Potential reward change (%)", fontsize=FONTSIZE['bigbartitle'])
    ax.set_xlabel(r"Capacity $\kappa$ player 1 (bits)", fontsize=FONTSIZE['bigaxis'])
    ax.set_xticks(k_ticks)
    ax.set_xticklabels(k_ticklabels, fontsize=FONTSIZE['bigticks'])
    ax.set_ylabel("Capacity $\kappa$ player 2 (bits)", fontsize=FONTSIZE['bigaxis'])
    ax.set_yticks(k_ticks)
    ax.set_yticklabels(k_ticklabels, fontsize=FONTSIZE['bigticks'])
    
    # Bottom-left
    ax = fig.add_subplot(337)
    ax.set_title("Reward function", fontsize=FONTSIZE['axtitle'])
    x = range(0, len(rewardscale[sd]))
    ax.plot(x, (rewardscale[sd]*(MAXREW-MINREW))+MINREW, '-', lw=3, color='#FF69B4')
    ax.set_xlabel("Error (deg)", fontsize=FONTSIZE['axis'])
    ax.set_xticks([0, 45, 90, 180])
    ax.set_ylabel("Reward (points)", fontsize=FONTSIZE['axis'])
    ax.set_yticks(range(0, 101, 20))
    
    # Bottom-centre
    ax = fig.add_subplot(338)
    ax.set_title("Equal division", fontsize=FONTSIZE['axtitle'])
    cax = ax.imshow(rew_equal[sd], vmin=0, vmax=200, cmap='inferno', \
        aspect='equal', interpolation='none', origin='lower')
    divider = make_axes_locatable(ax)
    bax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(cax, cax=bax, ticks=range(0, 201, 50), orientation='vertical')
    cbar.ax.set_yticklabels(map(str, range(0, 201, 50)), fontsize=FONTSIZE['barticks'])
    #cbar = fig.colorbar(cax, ticks=range(0, 201, 50), orientation='vertical')
    #cbar.ax.set_ylabel(r"Expected joint reward", fontsize=FONTSIZE['bartitle'])
    ax.set_xlabel(r"$\kappa$ player 1", fontsize=FONTSIZE['axis'])
    ax.set_xticks(k_ticks)
    ax.set_xticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])
    ax.set_ylabel(r"$\kappa$ player 2", fontsize=FONTSIZE['axis'])
    ax.set_yticks(k_ticks)
    ax.set_yticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])
    
    # Bottom-right
    ax = fig.add_subplot(339)
    ax.set_title("Proportional division", fontsize=FONTSIZE['axtitle'])
    divider = make_axes_locatable(ax)
    cax = ax.imshow(rew_prop[sd], vmin=0, vmax=200, cmap='inferno', \
        aspect='equal', interpolation='none', origin='lower')
    bax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(cax, cax=bax, ticks=range(0, 201, 50), orientation='vertical')
    cbar.ax.set_yticklabels(map(str, range(0, 201, 50)), fontsize=FONTSIZE['barticks'])
    #cbar = fig.colorbar(cax, ticks=range(0, 201, 50), orientation='vertical')
    cbar.ax.set_ylabel(r"Expected joint reward", fontsize=FONTSIZE['bartitle'])
    ax.set_xlabel(r"$\kappa$ player 1", fontsize=FONTSIZE['axis'])
    ax.set_xticks(k_ticks)
    ax.set_xticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])
    #ax.set_ylabel(r"$\kappa$ player 2", fontsize=FONTSIZE['axis'])
    ax.set_yticks(k_ticks)
    ax.set_yticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])
    
    # SAVE FULL FIGURE
    fig.savefig("equal_vs_proportional_sd-%s.png" % (sd))
    pyplot.close(fig)


# Create a plot to compare the different reward structures.
fig, axes = pyplot.subplots(nrows=len(REWSD), ncols=6, \
    figsize=(16.0,12.0), dpi=300.0)
fig.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.90,
    wspace=0.3, hspace=0.3)

# Establish what the capacity ranged and tickslabels are.
sel = KRANGE % 1 == 0
sel[0] = True
k_ticks = numpy.arange(0, len(KRANGE))[sel]
k_ticklabels = map(str, KRANGE[sel].astype(int))

for i, sd in enumerate(REWSD):

    # Reward scale.
    if type(sd) is int:
        label = r"$\sigma$ = %d" % (sd)
    elif type(sd) is float:
        label = r"$\sigma$ = %.2f" % (sd)
    elif sd == 'lin':
        label = "Linear"
    if i == 0:
        axes[i,0].set_title("Reward function", fontsize=FONTSIZE['axtitle'])
    x = range(0, len(rewardscale[sd]))
    axes[i,0].plot(x, (rewardscale[sd]*(MAXREW-MINREW))+MINREW, '-', lw=3, \
        color='#FF69B4', label=label)
    if i == len(REWSD)-1:
        axes[i,0].set_xlabel("Error (deg)", fontsize=FONTSIZE['axis'])
    axes[i,0].set_xticks([0, 15, 45, 90])
    axes[i,0].set_ylabel("Reward (points)", fontsize=FONTSIZE['axis'])
    axes[i,0].set_yticks(range(0, 101, 20))
    axes[i,0].legend(loc="upper right", fontsize=FONTSIZE['axis'])
    
    # Winnings.
    for j, divide_type in enumerate(["Equal", "Proportional", "Paradoxical"]):
        if divide_type == "Equal":
            val = rew_equal[sd]
        elif divide_type == "Proportional":
            val = rew_prop[sd]
        elif divide_type == "Paradoxical":
            val = rew_optimal[sd]
        if i == 0:
            axes[i,j+1].set_title("%s" % (divide_type), fontsize=FONTSIZE['axtitle'])
        cax = axes[i,j+1].imshow(val, vmin=0, vmax=200, cmap='inferno', \
            aspect='equal', interpolation='none', origin='lower')
        divider = make_axes_locatable(axes[i,j+1])
        bax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(cax, cax=bax, ticks=range(0, 201, 50), orientation='vertical')
        cbar.ax.set_yticklabels(map(str, range(0, 201, 50)), fontsize=FONTSIZE['barticks'])
        #cbar = fig.colorbar(cax, ticks=range(0, 201, 50), orientation='vertical')
#        if j == 2:
#            cbar.ax.set_ylabel(r"Expected joint reward", fontsize=FONTSIZE['bartitle'])
        if i == len(REWSD)-1:
            axes[i,j+1].set_xlabel(r"$\kappa$ player 1 (bits)", fontsize=FONTSIZE['axis'])
        axes[i,j+1].set_xticks(k_ticks)
        axes[i,j+1].set_xticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])
        if j == 0:
            axes[i,j+1].set_ylabel(r"$\kappa$ player 2 (bits)", fontsize=FONTSIZE['axis'])
        axes[i,j+1].set_yticks(k_ticks)
        axes[i,j+1].set_yticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])
    
    # Relative benefit.
    for j, divide_type in enumerate(["Proportional", "Paradoxical"]):
        if divide_type == "Proportional":
            vmin = -0.1
            vmax = 0.1
            val = prop_increase[sd]
        elif divide_type == "Paradoxical":
            vmin = -1.0
            vmax = 1.0
            val = optimal_increase[sd]
        if i == 0:
            axes[i,j+4].set_title("%s" % (divide_type), fontsize=FONTSIZE['axtitle'])
        cax = axes[i,j+4].imshow(val, vmin=vmin, vmax=vmax, cmap='coolwarm', \
            aspect='equal', interpolation='none', origin='lower')
        divider = make_axes_locatable(axes[i,j+4])
        bax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = fig.colorbar(cax, cax=bax, ticks=[vmin, 0, vmax], orientation='vertical')
        yticklabels = map(str, (cbar.get_ticks()*100).astype(int))
        for k in range(len(yticklabels)):
            yticklabels[k] = yticklabels[k] + '%'
        cbar.ax.set_yticklabels(yticklabels, fontsize=FONTSIZE['barticks'])
#        if j == 1:
#            cbar.ax.set_ylabel("Potential reward change (%)", fontsize=FONTSIZE['bartitle'])
        if i == len(REWSD)-1:
            axes[i,j+4].set_xlabel(r"$\kappa$ player 1 (bits)", fontsize=FONTSIZE['axis'])
        axes[i,j+4].set_xticks(k_ticks)
        axes[i,j+4].set_xticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])
#        if j == 0:
#            axes[i,j+4].set_ylabel("$\kappa$ player 2 (bits)", fontsize=FONTSIZE['axis'])
        axes[i,j+4].set_yticks(k_ticks)
        axes[i,j+4].set_yticklabels(k_ticklabels, fontsize=FONTSIZE['ticks'])

# Save full figure.
fig.savefig("reward_scheme_comparison.png")
pyplot.close(fig)

