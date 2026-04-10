import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import LogLocator
from matplotlib.transforms import Affine2D

import math
import csv
import re

def get_bg_cross_section(cms_energy):
    if cms_energy == 250:
        qqbar = 127965.53 + 70416.743
        ww_had = 14866.42 + 136.82153
        w_sl = 10264.016 + 86.696149
        ww_sl = 18779.145 + 173.46829
        Zee_sl = 1423.3098 + 1219.3967
        Znunu_sl = 453.86976 + 131.21958
        ZZ_had = 1405.06 + 606.70978
        ZZWWMix = 12389.292 + 225.56868 # ??
    elif cms_energy == 500:
        qqbar = 31686.082 + 17581.895
        ww_had = 7680.69 + 33.52
        w_sl = 7805.53 + 22.83
        ww_sl = 9521.45 + 45.58
        Zee_sl = 1961.13 + 1726.55
        Znunu_sl = 951.72 + 58.98
        ZZ_had = 680.22 + 271.88
        ZZWWMix = 6400.11 + 78.70 # ??
    else:
        print("WARNING: get_bg_cross_section(): wrong CMS energy!")
    
    tot_bg = qqbar + ww_had + w_sl + ww_sl + Zee_sl + Znunu_sl + ZZ_had + ZZWWMix
    print("qqbar: ", qqbar, "tot. bg.: ", tot_bg)
    return tot_bg

def get_bg_events(cms_energy, eff, bg_err):
    if cms_energy == 250:
        tot_lumi=2000
        eLpR_lumi = 526.5
        eRpL_lumi = 526.5
        eLpL_lumi = 58.5
        eRpR_lumi = 58.5
        aa_lumi_BB = tot_lumi*0.4337 # 862.972721997
        aa_lumi_BW = tot_lumi*0.5257
        aa_lumi_WW = tot_lumi


        qqbar = 127965.53*eLpR_lumi + 70416.743*eRpL_lumi
        ww_had = 14866.42*eLpR_lumi + 136.82153*eRpL_lumi
        w_sl = 10264.016*eLpR_lumi + 86.696149*eRpL_lumi + 190.53144*eLpL_lumi + 190.63749*eRpR_lumi
        ww_sl = 18779.145*eLpR_lumi + 173.46829*eRpL_lumi
        Zee_sl = 1423.3098*eLpR_lumi + 1219.3967*eRpL_lumi + 1155.8334*eLpL_lumi + 1157.2006*eRpR_lumi
        Znunu_sl = 453.86976*eLpR_lumi + 131.21958*eRpL_lumi
        ZZ_had = 1405.06*eLpR_lumi + 606.70978*eRpL_lumi
        ZZ_sl = 838.07949*eLpR_lumi + 466.81644*eRpL_lumi
        ZZWWMix = 12389.292*eLpR_lumi + 225.56868*eRpL_lumi
        llbar = 21214.001*eLpR_lumi + 16363.043*eRpL_lumi
        aa_had = (90338.374*aa_lumi_BW + 90119.741*aa_lumi_BW + 71505.71*aa_lumi_WW + 42149.6384*aa_lumi_BB) #* aa_lumi
    elif cms_energy == 500:
        eLpR_lumi = 936
        eRpL_lumi = 936
        eLpL_lumi = 234
        eRpR_lumi = 234

        qqbar = 31686.082*eLpR_lumi + 17581.895*eRpL_lumi
        ww_had = 7680.69*eLpR_lumi + 33.52*eRpL_lumi
        w_sl = 7805.53*eLpR_lumi + 22.83*eRpL_lumi + 753.07*eLpL_lumi + 750.07*eRpR_lumi
        ww_sl = 9521.45*eLpR_lumi + 45.58*eRpL_lumi
        Zee_sl = 1961.13*eLpR_lumi + 1726.55*eRpL_lumi + 1775.49*eLpL_lumi + 1778.01*eRpR_lumi
        Znunu_sl = 951.72*eLpR_lumi + 58.98*eRpL_lumi
        ZZ_had = 680.22*eLpR_lumi + 271.88*eRpL_lumi
        ZZ_sl = 608.57*eLpR_lumi + 288.36*eRpL_lumi
        ZZWWMix = 6400.11*eLpR_lumi + 78.70*eRpL_lumi
    else:
        print("WARNING: get_bg_cross_section(): wrong CMS energy!")
    
    tot_bg = qqbar + ww_had + w_sl + ww_sl + Zee_sl + Znunu_sl + ZZ_had + ZZ_sl + ZZWWMix + aa_had + llbar
    tot_bg_after = qqbar*eff["qqbar"] + ww_had*eff["WWhad"] + (w_sl + ww_sl)*(eff["qqbar"]+eff["llbar"]/2.) \
        + Zee_sl*eff["qqbar"] + Znunu_sl*eff["qqbar"] + ZZ_had*eff["WWhad"] + ZZ_sl*(eff["qqbar"]+eff["llbar"]) + ZZWWMix*eff["WWhad"] \
        + aa_had*eff["aahad"] \
        + llbar*eff["llbar"] 
    
    N_dijet = qqbar + w_sl + ww_sl + Zee_sl + Znunu_sl + ZZ_sl
    N_fjet = ww_had + ZZ_had + ZZWWMix
    err_bg_after = math.sqrt( (N_dijet*bg_err["qqbar"])**2 + (N_fjet*bg_err["WWhad"])**2 \
                             + (aa_had*bg_err["aahad"])**2 + (llbar*bg_err["llbar"])**2 )
    
    print("efficiency; qqbar: ", eff["qqbar"], " WWhad: ", eff["WWhad"])
    print('')
    print("qqbar: ", qqbar, "ww_had: ",ww_had, "w_sl: ",w_sl, "ww_sl: ",ww_sl, "Zee_sl: ",Zee_sl, "Znunu_sl: ",Znunu_sl, \
          "ZZ_had: ",ZZ_had, "ZZ_sl: ",ZZ_sl, "ZZWWMix: ",ZZWWMix, "aa_had: ", aa_had, "llbar: ", llbar)
    print("tot. bg.: ", tot_bg)
    print('')
    print("after cuts:")
    print("qqbar: ", qqbar*eff["qqbar"], "+/-",  "ww_had: ",ww_had*eff["WWhad"], "w_sl: ",w_sl*eff["qqbar"], \
          "ww_sl: ",ww_sl*eff["qqbar"], "Zee_sl: ",Zee_sl*eff["qqbar"], "Znunu_sl: ",Znunu_sl*eff["qqbar"], \
          "ZZ_had: ",ZZ_had*eff["WWhad"], "ZZ_sl: ",ZZ_sl*eff["qqbar"], "ZZWWMix: ",ZZWWMix*eff["WWhad"], \
            "aa_had: ", aa_had*eff["aahad"], "llbar: ", llbar*eff["llbar"])
    print("tot. bg. after cuts: ", tot_bg_after, "+/-", err_bg_after, '\n')
    
    return tot_bg_after, err_bg_after



def get_limits_plot(filename, ax, marker_style, line_style, color_cycle, eff, bg_err, branching_ratio, save_to_csv=False):

    n_ltimes = 9 # Number of lifetimes in file
    i_ltime = 0 # initial lifetime to plot
    f_ltime = n_ltimes  # final lifetime to plot

    
    cmsen = 250
    
    if cmsen == 250:
        lumi = 2000
        # bxs = 4.14e11 * 5 # 5 years of running effectively
        bxs = 9.7e11 # 7.4E+7 sec. * 5 Hz * 2625 to get 2 ab-1
        ev_per_bx = 1.55
    if cmsen == 500:
        bxs = 4.14e11 * 8.5
        ev_per_bx = 1.05
        lumi = 4000
    if filename.find("fimp") != -1:
        leg_label = "$m_{F} - m_H=$"
        signal = 'fimp'
        xfactor = 0.1 # IDM mass splitting is in 10 GeV and ALP masses in MeV

    hnunu_cr = (60.35+21.46+67.11+42.93)

                # Z to nunu
    hz_cr = (60.35+21.46+67.11+42.93 \
                # Z to qq
                + 343.03023 + 219.48615 \
                # Z to ll
                + 17.671491 + 11.138876 + 16.970655 + 10.869108 + 16.940726 + 10.843428 \
                ) + (0.62348544 + 0.62348544)
                            # Z to ee
                

    # ggtohad_eff = 1.e-9
    # eepairs_eff = 1.e-10
    overlay_eff = 1.26e-10
    if filename.find("fimp") != -1:
        overlay_eff = 0.

    ggtohad = ev_per_bx * bxs
    eepairs = 1. * bxs
    # tot_bg = get_bg_cross_section(cmsen) * lumi
    tot_sm_after, sm_bg_err = get_bg_events(cmsen, eff, bg_err)

    print('overlay before cuts: ', (ggtohad + eepairs), '\n')
    print('overlay after cuts: ', (ggtohad + eepairs)*overlay_eff, '\n')

    # bg_after_sel = ggtohad*ggtohad_eff + eepairs*eepairs_eff + qqbar*qqbar_eff
    bg_after_sel = (ggtohad + eepairs)*overlay_eff + tot_sm_after #tot_bg*qqbar_eff
    tot_bg_err = sm_bg_err # assuming no error on overlay
    # n_obs = 1.64 * math.sqrt(bg_after_sel)
    n_obs = 1.96 * math.sqrt(bg_after_sel)
    if bg_after_sel < 1:
        n_obs = 4.64
    if bg_after_sel == 0:
        n_obs = 2.63

    # Read data from the file using NumPy
    data = np.loadtxt(filename, dtype=str, skiprows=1)

    # Parse the data and convert columns to appropriate data types
    scenarios = []
    labels = []
    for entry in data[:, 0]:
        # Use regex to extract all numerical values in the string
        numbers = list(map(float, re.findall(r'\d+\.?\d*', entry)))
        unit = 'GeV'
        if entry.find('kev') != -1:
            unit = 'keV'
        if len(numbers) == 2:
            scenarios.append(tuple(numbers))
            labels.append('$m_{F} =$ '+str(numbers[1])+' GeV, $m_{s} =$ '+str(numbers[0])+' ' + unit)

    efficiencies = data[:, 8].astype(float) / data[:, 7].astype(float) 
    frac_in_tpc = data[:, 7].astype(float) / data[:, 6].astype(float) 
    eff_pass = data[:, 8].astype(float) / data[:, 6].astype(float) 
    n_rej = data[:, 6].astype(float)  - data[:, 8].astype(float)
    n_pass = data[:, 8].astype(float)
    err_rej = data[:, 9].astype(float)
    err_pass = data[:, 10].astype(float)
    # err_acc = err_pass / data[:, 6].astype(float) # approx.!!
    err_acc = np.sqrt( (err_pass**2 * n_rej**2 + err_rej**2 * n_pass**2) ) \
                / (n_pass+n_rej)**2
    decay_len = data[:, 1].astype(float)

    # Plot the data
    # plt.figure(figsize=(8, 6))
    # fig, ax = plt.subplots()
    # ax.locator_params(nbins=10, axis='y')
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=range(100)))

    # Set x-axis range
    # ax.set_xlim(decay_len[0]/2, decay_len[-1]*1.5)
    ax.set_xlim(1.e1, decay_len[-1]*1.5)
    ax.set_ylim(9.e-2, 2.e3)
    if branching_ratio:
        ax.set_ylim(1.e-6, 1.)

    handles = []  # To store legend handles


    t = 0
    # while t < scenarios.size:
    while t < len(scenarios):

        n_lim = n_obs / (eff_pass[(t+i_ltime):(t+f_ltime)])
        cl_limits = n_lim / lumi

        # errors = n_obs * np.sqrt( (err_rej[(t+i_ltime):(t+f_ltime)]/n_pass[(t+i_ltime):(t+f_ltime)])**2 \
        #                             + (n_rej[(t+i_ltime):(t+f_ltime)]*err_pass[(t+i_ltime):(t+f_ltime)]/(n_pass[(t+i_ltime):(t+f_ltime)]**2))**2 ) \
        #                     / lumi
        if bg_after_sel < 1:
            errors = n_obs * err_acc[(t+i_ltime):(t+f_ltime)] / (lumi * eff_pass[(t+i_ltime):(t+f_ltime)]**2) # signal err only
        else:
            errors = np.sqrt(
                     ( n_obs * err_acc[(t+i_ltime):(t+f_ltime)] / (lumi * eff_pass[(t+i_ltime):(t+f_ltime)]**2) )**2 \
                    +( 1.96*tot_bg_err / (lumi*eff_pass[(t+i_ltime):(t+f_ltime)]*np.sqrt(bg_after_sel))       )**2 \
            )

        if branching_ratio == True:
            # cl_limits /= hz_cr
            # errors    /= hz_cr
            cl_limits /= hnunu_cr
            errors    /= hnunu_cr

        # print("eff**2",eff_pass[(t+i_ltime):(t+f_ltime)]**2)
        # print("err",err_acc[(t+i_ltime):(t+f_ltime)])
        # print("n_obs",n_obs)
        # print("lumi",lumi,'\n')
        print(cl_limits, '\n')
        print('Signal error part: ', n_obs * err_acc[(t+i_ltime):(t+f_ltime)] / (lumi * eff_pass[(t+i_ltime):(t+f_ltime)]**2))
        print('Bckg. error part: ', 1.96*tot_bg_err / (lumi*eff_pass[(t+i_ltime):(t+f_ltime)]*np.sqrt(bg_after_sel)))
        
        # trans = Affine2D().translate(float(t)/5, 0.0) + ax.transData
        trans = Affine2D().translate(0.0, 0.0) + ax.transData

        color = color_cycle[int(t / n_ltimes)]
        # print(t/n_ltimes, color)

        er = ax.errorbar(decay_len[(t+i_ltime):(t+f_ltime)], cl_limits, yerr=errors, \
                         marker=marker_style, linestyle=line_style, color=color, \
                            label=str(labels[t]), \
                                ms=4, capsize=5.0, transform=trans) #'$c\tau$='+str(scenarios[t]) label='$\sqrt{s}=$'+str(cmsen) \
                                                                        #+' GeV, $\mathcal{L}=$'+str(lumi)+' fb$^{-1}$')
        handles.append(er)  # Add the line handle to the list

        # if save_to_csv:
        #     # Save cl_limits and errors to CSV file
        #     with open(f'./Limits/limits_{signal}_all_scenarios.csv', mode='a', newline='') as file:
        #         writer = csv.writer(file)
        #         if t == 0:
        #             writer.writerow(['Scenario (GeV)', 'Decay Length (mm)', 'CL Limits', 'Errors'])
        #         for dl, cl, err in zip(decay_len[(t+i_ltime):(t+f_ltime)], cl_limits, errors):
        #             writer.writerow([scenarios[t], dl, cl, err])
        
        t += n_ltimes

    plt.xlabel(r"$c\tau$ [mm]", fontsize=16)
    if branching_ratio == True:
        plt.ylabel(r"$\sigma_{\mathrm{95\% C.L.}}/\sigma_{h\nu\nu}$", fontsize=16)
    else:
        plt.ylabel('$\sigma_{\mathrm{95\% C.L.}}$ [fb]', fontsize=16)
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=12)

    if branching_ratio:
        plt.yticks([1.e-6,1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.],fontsize=12) 
    elif filename.find("fimp") != -1:
        plt.yticks([0.1, 1., 1.e1, 1.e2, 1.e3, 1.e4],fontsize=12)


    return handles


# Check if the correct number of arguments is provided
if len(sys.argv) != 3 and len(sys.argv) != 2:
    print("Usage: python plot_limits.py <standard-selection-data> <tight-selection-data>")
    sys.exit(1)

# Get the filename from the command-line argument
filename1 = sys.argv[1]
filename2 = ""
if len(sys.argv) > 2:
    filename2 = sys.argv[2]

branching_ratio = False

# plot1 = get_limits_plot(filename1)
# plot2 = get_limits_plot(filename2)

# # qqbar_eff_loose = 0.00113498389907
# qqbar_eff_loose = 0.000901803607214 # loose mass with RCut
# # qqbar_eff_tight = 6.0120240481e-05
# qqbar_eff_tight = 3.00601202405e-05 # tight mass with RCut
signal = 'fimp'
name = ''

# 250 gev only
n_evts = {"qqbar": 171034, "WWhad": 179600, "aahad": 2973773, "llbar": 479000}
# eff_loose = {"qqbar": 0.00026310558134639898, "WWhad": 0.0004064587973273942, "aahad": 0.1681365726301234e-05} # loose mass, RCut, refDist50
eff_loose = {"qqbar": 0.00004092753487610651, "WWhad": 0.00010579064587973273, 
             "aahad": 6.725462905204937e-07, "llbar": 0.00020041753653444676} # loose mass, RCut, refDist50
eff_tight = {"qqbar": 2.2979525243e-05, "WWhad": 3.56730917869e-05, "aahad": 1.06398410266e-06} # tight mass, RCut, refDist50, iso < 1
err_loose = {"qqbar": math.sqrt(eff_loose["qqbar"]*(1.-eff_loose["qqbar"])/n_evts["qqbar"]), 
             "WWhad": math.sqrt(eff_loose["WWhad"]*(1.-eff_loose["WWhad"])/n_evts["WWhad"]), 
             "aahad": math.sqrt(eff_loose["aahad"]*(1.-eff_loose["aahad"])/n_evts["aahad"]), 
             "llbar": math.sqrt(eff_loose["llbar"]*(1.-eff_loose["llbar"])/n_evts["llbar"]), 
            } # loose mass, RCut, refDist50
err_tight = {"qqbar": 1.15e-05, "WWhad": 1.14e-05, "aahad": 0.61e-06} # tight mass, RCut, refDist50, iso < 1

print("eff_loose: ", eff_loose)
print("err_loose: ", err_loose)


fig, ax = plt.subplots()

colors_loose = ['tomato', 'yellowgreen', 'orange', 'dodgerblue', 'violet', 'cyan']
colors_tight = ['tomato', 'yellowgreen', 'orange', 'dodgerblue']

save_to_csv = True

if len(sys.argv) > 2:
    handles2 = get_limits_plot(filename2, ax, marker_style='s', line_style='--', color_cycle=colors_tight, 
                            eff=eff_tight, bg_err=err_tight, branching_ratio=branching_ratio, 
                            save_to_csv=save_to_csv)
handles1 = get_limits_plot(filename1, ax, marker_style='o', line_style='-', color_cycle=colors_loose, 
                           eff=eff_loose, bg_err=err_loose, branching_ratio=branching_ratio, 
                           save_to_csv=save_to_csv)

ax.legend(fontsize=11,handles=handles1)

ax.text(0.01, 1.055, 'ILD Preliminary', transform=ax.transAxes, verticalalignment='top', horizontalalignment='left', fontsize=12, weight='bold')
ax.text(0.02, 0.01, name, transform=ax.transAxes, verticalalignment='bottom', horizontalalignment='left', fontsize=12)
ax.text(0.999, 1.07, '$\sqrt{s}=250$ GeV, $\int\mathcal{L}\mathrm{d}t=2\,\mathrm{ab}^{-1}$', transform=ax.transAxes, verticalalignment='top', \
        horizontalalignment='right', fontsize=12)

if branching_ratio:
    plt.savefig('./Limits/' + signal + '_allCuts_BR.pdf')
else:
    plt.savefig('./Limits/' + signal + '_allCuts_CS.pdf')
plt.show()
