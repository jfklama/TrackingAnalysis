import sys
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import LogLocator
from matplotlib.transforms import Affine2D

from labellines import labelLine, labelLines


from lhc_limits import plot_lhc_limits

import math
import csv

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
    
    tot_bg = qqbar + ww_had + w_sl + ww_sl + Zee_sl + Znunu_sl + ZZ_had + ZZ_sl + ZZWWMix + llbar + aa_had
    tot_bg_after = qqbar*eff["qqbar"] + ww_had*eff["WWhad"] + w_sl*eff["qqbar"] + ww_sl*(eff["qqbar"]+0.5*eff["llbar"]) \
        + Zee_sl*eff["qqbar"] + Znunu_sl*eff["qqbar"] + ZZ_had*eff["WWhad"] + ZZ_sl*(eff["qqbar"]+eff["llbar"]) \
        + ZZWWMix*eff["WWhad"] + aa_had*eff["aahad"] + llbar*eff["llbar"]
    
    N_dijet = qqbar + w_sl + ww_sl + Zee_sl + Znunu_sl + ZZ_sl
    N_fjet = ww_had + ZZ_had + ZZWWMix
    N_lep = llbar + 0.5*ww_sl + ZZ_sl
    err_bg_after = math.sqrt( (N_dijet*bg_err["qqbar"])**2 + (N_fjet*bg_err["WWhad"])**2  \
                             + (N_lep*bg_err["llbar"])**2 + (aa_had*bg_err["aahad"])**2 )
    # err_bg_after = math.sqrt( ((qqbar+w_sl)*bg_err["qqbar"])**2 )
    
    print("efficiency; qqbar: ", eff["qqbar"], "+/-", bg_err["qqbar"], ", WWhad: ", eff["WWhad"], "+/-", bg_err["WWhad"], \
          ", aahad: ", eff["aahad"], "+/-", bg_err["aahad"], ", llbar: ", eff["llbar"], "+/-", bg_err["llbar"])
    print('')
    print("qqbar: ", qqbar, "ww_had: ",ww_had, "w_sl: ",w_sl, "ww_sl: ",ww_sl, "Zee_sl: ",Zee_sl, "Znunu_sl: ",Znunu_sl, \
          "ZZ_had: ",ZZ_had, "ZZ_sl: ",ZZ_sl, "ZZWWMix: ",ZZWWMix, "aa_had: ", aa_had, "llbar: ", llbar)
    print("tot. bg.: ", tot_bg)
    print('')
    print("after cuts:")
    print("qqbar: ", qqbar*eff["qqbar"], "ww_had: ",ww_had*eff["WWhad"], "w_sl: ",w_sl*eff["qqbar"], \
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
    if filename.find("idm") != -1:
        x_label = "$m_{A} - m_H$ [GeV]"
        leg_label = "$m_{A} - m_H=$"
        signal = 'idm'
        xfactor = 0.1 # IDM mass splitting is in 10 GeV and ALP masses in MeV
    elif filename.find("alp") != -1:
        signal = 'alp'
        x_label = "$m_a$ [GeV]"
        leg_label = "$m_a=$"
        xfactor = 0.001
    elif filename.find("trsm") != -1:
        leg_label = "$m_s=$"
        xfactor = 1.
        lumi = 526.5*2
        if filename.find("MassScenarios") != -1:
            leg_label = "ILD, 250 GeV, $2\,\mathrm{ab}^{-1}$, " + leg_label

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
    if filename.find("trsm") != -1:
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
    if bg_after_sel < 10 and bg_after_sel > 7.: # FIXME: hardcored for bg_after_sel = 8.36 +/- 4.11
        n_obs = 7.604


    # Read data from the file using NumPy
    data = np.loadtxt(filename, dtype=str, skiprows=1)

    # Parse the data and convert columns to appropriate data types
    scenarios = data[:, 0].astype(float) * xfactor 
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
    ax.set_xlim(decay_len[0]/2, decay_len[-1]*1.5)
    if branching_ratio:
        ax.set_ylim(1.e-6, 1.)

    handles = []  # To store legend handles


    t = 0
    while t < scenarios.size:

        n_lim = n_obs / (eff_pass[(t+i_ltime):(t+f_ltime)])
        cl_limits = n_lim / lumi

        # errors = n_obs * np.sqrt( (err_rej[(t+i_ltime):(t+f_ltime)]/n_pass[(t+i_ltime):(t+f_ltime)])**2 \
        #                             + (n_rej[(t+i_ltime):(t+f_ltime)]*err_pass[(t+i_ltime):(t+f_ltime)]/(n_pass[(t+i_ltime):(t+f_ltime)]**2))**2 ) \
        #                     / lumi
        if bg_after_sel < 1:
            errors = n_obs * err_acc[(t+i_ltime):(t+f_ltime)] / (lumi * eff_pass[(t+i_ltime):(t+f_ltime)]**2) # signal err only
        elif bg_after_sel < 10 and bg_after_sel > 7.:
            # lower error on limit assuming 4 observed events
            errors = np.array([5.970 / (eff_pass[(t+i_ltime):(t+f_ltime)])]) / lumi 
            # upper error on limit assuming 11 observed events
            upper_error = (8.832 / (eff_pass[(t+i_ltime):(t+f_ltime)])) / lumi
            errors = np.vstack([errors, upper_error])
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
        if bg_after_sel > 0:
            print('Bckg. error part: ', 1.96*tot_bg_err / (lumi*eff_pass[(t+i_ltime):(t+f_ltime)]*np.sqrt(bg_after_sel)))
        print('Total error: ', errors)
        print()
        
        # trans = Affine2D().translate(float(t)/5, 0.0) + ax.transData
        trans = Affine2D().translate(0.0, 0.0) + ax.transData

        color = color_cycle[int(t / n_ltimes)]
        # print(t/n_ltimes, color)

        if filename.find("trsm") != -1: #and int(t / n_ltimes) == 0:
            # Plot the central line
            er, = ax.plot(decay_len[(t+i_ltime):(t+f_ltime)], cl_limits, 
                    marker=marker_style, linestyle=line_style, color=color, 
                    label=leg_label+str(scenarios[t])+' GeV', ms=4, transform=trans)

            # Plot the error band
            ax.fill_between(decay_len[(t+i_ltime):(t+f_ltime)], 
                    errors[0,:], errors[1,:], 
                    color=color, alpha=0.2, linewidth=0, transform=trans)
            
        else:
            er = ax.errorbar(decay_len[(t+i_ltime):(t+f_ltime)], cl_limits, yerr=errors, \
                         marker=marker_style, linestyle=line_style, color=color, fillstyle='full', \
                            label=leg_label+str(scenarios[t])+' GeV', ms=4, capsize=5.0, transform=trans) #'$c\tau$='+str(scenarios[t]) label='$\sqrt{s}=$'+str(cmsen) \
                                                                        # +' GeV, $\mathcal{L}=$'+str(lumi)+' fb$^{-1}$')

        # plt.fill_between(decay_len[(t+i_ltime):(t+f_ltime)], cl_limits-errors, cl_limits+errors)
        handles.append(er)  # Add the line handle to the list

        if save_to_csv:
            # Save cl_limits and errors to CSV file
            with open(f'./Limits/limits_{signal}_all_scenarios.csv', mode='a', newline='') as file:
                writer = csv.writer(file)
                if t == 0:
                    writer.writerow(['Scenario (GeV)', 'Decay Length (mm)', 'CL Limits', 'Errors'])
                for dl, cl, err in zip(decay_len[(t+i_ltime):(t+f_ltime)], cl_limits, errors):
                    writer.writerow([scenarios[t], dl, cl, err])
        
        t += n_ltimes

    plt.xlabel(r"$c\tau$ [mm]", fontsize=16)
    if branching_ratio == True:
        plt.ylabel(r"$\sigma_{\mathrm{95\% C.L.}}/\sigma_{h\nu\nu}$", fontsize=16)
    else:
        plt.ylabel('$\sigma_{\mathrm{95\% C.L.}}$ [fb]', fontsize=16)
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=14)

    if branching_ratio:
        plt.yticks([1.e-6,1.e-5,1.e-4,1.e-3,1.e-2,1.e-1,1.],fontsize=14) 
    elif filename.find("trsm") != -1:
        plt.yticks([1.e-2,1.e-1,1.,1.e1,1.e2,1.e3,1.e4,1.e5],fontsize=14)
    else:
        plt.yticks([1.e-1,1.,1.e1,1.e2,1.e3,1.e4,1.e5,1.e6],fontsize=14)

    return handles


# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python plot_limits.py <standard-selection-data> <tight-selection-data>")
    sys.exit(1)

# Get the filename from the command-line argument
filename1 = sys.argv[1]
filename2 = sys.argv[2]

branching_ratio = False

# plot1 = get_limits_plot(filename1)
# plot2 = get_limits_plot(filename2)

# # qqbar_eff_loose = 0.00113498389907
# qqbar_eff_loose = 0.000901803607214 # loose mass with RCut
# # qqbar_eff_tight = 6.0120240481e-05
# qqbar_eff_tight = 3.00601202405e-05 # tight mass with RCut
signal = 'idm'
name = 'heavy scalars'
if filename1.find("alp") != -1:
    signal = 'alp'
    name = 'light pseudoscalar'
if filename1.find("trsm") != -1:
    signal = 'trsm'
    name = ''
    branching_ratio = True
    # branching_ratio = False

# 250 gev only
eff_loose = {"qqbar": 0.000798538502195, "WWhad": 0.00148637882445, "aahad": 2.12796820532e-05, "llbar": 0.000127348643006} # loose mass, RCut, refDist50
err_loose = {"qqbar": 0.68e-04, "WWhad": 0.094e-03, "aahad": 0.28e-05, "llbar": 1.63042846508e-05} # loose mass, RCut, refDist50
eff_tight = {"qqbar": 2.2979525243e-05, "WWhad": 3.56730917869e-05, "aahad": 1.06398410266e-06, "llbar": 4.17536534447e-06} # tight mass, RCut, refDist50, iso < 1
err_tight = {"qqbar": 1.15e-05, "WWhad": 1.14e-05, "aahad": 0.61e-06, "llbar": 2.95242298526e-06} # tight mass, RCut, refDist50, iso < 1

if filename1.find("trsm") != -1:
    eff_loose = {"qqbar": 0.0002814991842268539, "WWhad": 0.000428077101442382, "aahad": 0.0, "llbar": 0.000108559498956} # loose mass, RCut, refDist50, vtxPt > 10
    eff_tight = {"qqbar": 0.0, "WWhad": 5.945515297810861e-06, "aahad": 0.0} # tight mass, RCut, refDist50, iso > 1, vtxPt > 10
    err_loose = {"qqbar": 0.68e-04, "WWhad": 0.094e-03, "aahad": 0.0, "llbar": 1.50536766669e-05} # loose mass, RCut, refDist50
    err_tight = {"qqbar": 1.15e-05, "WWhad": 1.14e-05, "aahad": 0.61e-06} # tight mass, RCut, refDist50, iso < 1

    eff_prompt = {"qqbar": 0.000134273, "WWhad": 0.0, "aahad": 0.208834356, "llbar": 0.001348921} # efficiency of prompt track cut
    err_prompt = {"qqbar": math.sqrt(eff_prompt["qqbar"]*(1.-eff_prompt["qqbar"])/14895), 
             "WWhad": math.sqrt(eff_prompt["WWhad"]*(1.-eff_prompt["WWhad"])/26747), 
             "aahad": math.sqrt(eff_prompt["aahad"]*(1.-eff_prompt["aahad"])/8150), 
                "llbar": math.sqrt(eff_prompt["llbar"]*(1.-eff_prompt["llbar"])/2224)
            }
    print("eff_prompt: ", eff_prompt, "err_prompt: ", err_prompt)
    
    # combined errors
    # err_loose = {"qqbar": math.sqrt(eff_loose["qqbar"]*(1.-eff_loose["qqbar"])/174068), 
    #          "WWhad": math.sqrt(eff_loose["WWhad"]*(1.-eff_loose["WWhad"])/168194), 
    #          "aahad": math.sqrt(eff_loose["aahad"]*(1.-eff_loose["aahad"])/2789591), 
    #         } # loose mass, RCut, refDist50
    err_loose = {"qqbar": math.sqrt( (eff_loose["qqbar"]*err_prompt["qqbar"])**2 + (err_loose["qqbar"]*eff_prompt["qqbar"])**2 ), 
                "WWhad": math.sqrt( (eff_loose["WWhad"]*err_prompt["WWhad"])**2 + (err_loose["WWhad"]*eff_prompt["WWhad"])**2 ),
                "aahad": math.sqrt( (eff_loose["aahad"]*err_prompt["aahad"])**2 + (err_loose["aahad"]*eff_prompt["aahad"])**2 ),
                "llbar": math.sqrt( (eff_loose["llbar"]*err_prompt["llbar"])**2 + (err_loose["llbar"]*eff_prompt["llbar"])**2 )
            } # loose mass, RCut, refDist50

    # factors from cut on prompt tracks
    # eff_loose["qqbar"] *= 1.14897e-05 # 0.004237288 
    eff_loose["qqbar"] *= eff_prompt["qqbar"] 
    eff_loose["WWhad"] *= eff_prompt["WWhad"]
    # eff_loose["aahad"] *= 0.00061 # 0.008823529
    eff_loose["aahad"] *= eff_prompt["aahad"] 
    eff_loose["llbar"] *= eff_prompt["llbar"]    

    eff_tight["qqbar"] *= eff_prompt["qqbar"] 
    eff_tight["WWhad"] *= eff_prompt["WWhad"]
    eff_tight["aahad"] *= eff_prompt["aahad"] 



fig, ax = plt.subplots()

colors_loose = ['tomato', 'yellowgreen', 'orange', 'dodgerblue']
colors_tight = ['tomato', 'yellowgreen', 'orange', 'dodgerblue']
# if signal == 'alp' or signal == 'trsm':
if signal == 'alp':
    colors_tight = ['yellowgreen', 'orange', 'dodgerblue']
if signal == 'trsm' and branching_ratio and filename1.find("MassScenarios") != -1:
    colors_loose = ['tomato', 'yellowgreen']
    if filename1.find("highMass") != -1:
        colors_loose = ['orange', 'dodgerblue']

save_to_csv = False

if signal != 'trsm':
    handles2 = get_limits_plot(filename2, ax, marker_style='s', line_style='--', color_cycle=colors_tight, 
                            eff=eff_tight, bg_err=err_tight, branching_ratio=branching_ratio, 
                            save_to_csv=save_to_csv)
handles1 = get_limits_plot(filename1, ax, marker_style='o', line_style='-', color_cycle=colors_loose, 
                           eff=eff_loose, bg_err=err_loose, branching_ratio=branching_ratio, 
                           save_to_csv=save_to_csv)

if signal == 'trsm' and branching_ratio and filename1.find("MassScenarios") != -1:
    # plt.ylabel(r"$\mathcal{B}$(h$\rightarrow$XX)", fontsize=16)
    plt.ylabel(r"BR(h$\rightarrow$XX)", fontsize=16)
    ax.set_ylim(1.e-7, 1.)
    if filename1.find("lowMass") != -1:
        ax.set_xlim(0.1, 1.e7)
        # files = {0.4: 'dielectron', 2: 'dimuon'}
        files = [[0.4, 'dielectron'], [2, 'dimuon']]
    elif filename1.find("highMass") != -1:
        ax.set_xlim(1., 1.e7)
        ax.set_ylim(1.e-8, 1.)
        # files = {50: 'dimuon', 60: 'recast'}
        files = [[60, 'dimuon'], [60, 'recast']]
        colors_loose = ['dodgerblue', 'dodgerblue']

    for f, col in zip(files, colors_loose):
        h = plot_lhc_limits(ax, f[0], f[1], col)
        handles1.append(h[0])

    # labelLines(ax.get_lines(), yoffsets=-1.e-3, xvals=[1,100])
    # print([l.get_label() for l in handles1])

ax.legend(fontsize=11,handles=handles1)

if filename1.find("MassScenarios") == -1:
    ax.text(0.01, 1.055, 'ILD Simulation', transform=ax.transAxes, verticalalignment='top', horizontalalignment='left', fontsize=12, weight='bold')
    ax.text(0.02, 0.01, name, transform=ax.transAxes, verticalalignment='bottom', horizontalalignment='left', fontsize=12)
    ax.text(0.999, 1.07, '$\sqrt{s}=250$ GeV, $\int\mathcal{L}\mathrm{d}t=2\,\mathrm{ab}^{-1}$', transform=ax.transAxes, verticalalignment='top', \
            horizontalalignment='right', fontsize=12)

plt.tight_layout()

if branching_ratio:
    plt.savefig('./Limits/' + signal + '_allCuts_BR.pdf')
else:
    plt.savefig('./Limits/' + signal + '_allCuts_CS.pdf')
plt.show()
