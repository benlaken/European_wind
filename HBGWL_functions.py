import numpy as np
import pandas as pd
import random
import sys
from scipy import stats
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import seaborn as sns


def calc_sem(data):
    """standard error of mean = sample standard deviation / square root(sample size - 1)"""
    return np.std(data)/np.sqrt(len(data)-1.)


def extractCompositeStatsAllDir(data, epoch=0):
    """
    From a dictionary object produced by gen_composite() create
    extract a list of mean and sem values for all directions 
    for a given epoch period (default 0).
    Input: the output of gen_composite()
    Output: a list of means and sem values for a specified epoch.
    """
    keys=["N","NE","E","SE","S","SW","W","NW"]
    index = data[keys[0]]["epoch"] == epoch
    means = [data[key]["means"][index] for key in keys]
    sem = [data[key]["sem"][index] for key in keys]
    return means, sem


def extractConfsByDirection(data):
    """
    Converts the dictionary output of a Monte Carlo call like:
    > confs[key]=hbgwl.monteCarlo(df=df, its=1000, key=key)
    into an easy to work with dataframe.
    Input: list of dictionaries
    Output: Pandas dataframe of confidence intervals by direction
    """
    percentiles=[0.5,2.5,5,50,95,97.5,99.5]
    keys=["N","NE","E","SE","S","SW","W","NW"]
    tmp_hold=[]
    for key in keys:
        tmp = [data[key][per] for per in percentiles]
        tmp_hold.append(tmp)
    return pd.DataFrame(data=tmp_hold,index=keys,columns=percentiles)


def find_non_overlapping_sample(datelist, prd = 365):
    """
    After feeding the function a filtered subset of dates,
    e.g. dataframe.sort("column",ascending=True).head(100).index
    this function will find the a non-overlapping dates within a
    specified period (default of 365). 
    User also specifies the minimum number of samples they wish
    to get from the input (i.e. the number of unique dates).
    Input: a list of dates ordered by magnitude.
    Output: A list of dates that do not have overlap within a
    user-specified period.
    
    """
    start_size = len(datelist)
    idx = 0
    while(idx < len(datelist)):
        #print(idx,len(datelist))
        timeDiff = [abs((datelist[idx] - nn).days) for nn in datelist]
        timeDiff = np.array(timeDiff)
        mask = ((timeDiff == 0) | (timeDiff > 365))
        datelist = datelist[mask]
        idx += 1
    print("Reduced {0} to {1} non-overlapping dates in ±{2} dy prd".format(
            start_size,len(datelist),prd,idx))
    return datelist


def gen_composite(data, keyDates, months = range(-12,13)):
    """
    Create composites.
    Inputs:
        data = the data to composite from
        key_dates = a list of key dates in pandas datetime format.
        months = a range of integers (e.g. range(-12,13)). These integers
        will be subtracted from a base time in steps of months. If not
        specified it will default from -12 to 12.
    Outputs: A dictionary object with keys of months (as epoch),
        means, and standard error of means (sem).
    """
    composite ={}
    mn_tmp = []
    sem_tmp = []
    for n in months:
        index = [(date-relativedelta(months=n)) for date in keyDates] 
        mn_tmp.append(np.mean(data[index]))
        sem_tmp.append(calc_sem(data[index]))
    composite["means"]=np.array(mn_tmp)
    composite["sem"]=np.array(sem_tmp)
    composite["epoch"]=np.array(months)
    return composite


def gen_keydates(df,key):
    """
    Identifies the largest/smallest values associated with a given key. 
    Then calls a function to filter out overlapping dates within 365 days.
    Finally, it selects the top 11, sorts them to ascending order, 
    and returns them as a dictionary with keys of 'max' and 'min'.
    """
    tmpmax = find_non_overlapping_sample(df.sort(key,ascending=False
                                                      ).head(150).index)
    tmpmin = find_non_overlapping_sample(df.sort(key,ascending=True
                                                      ).head(150).index)
    compPhase ={}
    compPhase['max'] = tmpmax[0:11].order()
    compPhase['min'] = tmpmin[0:11].order()
    print("Average for max/min sample of {0}: {1:2.3f} and {2:2.3f}".format(
            key,np.mean(df[key][compPhase['max']]),
            np.mean(df[key][compPhase['min']])))
    return compPhase


def get_p_from_kdes(df, solarPhase, naoPhase, ensoPhase, its=1000, epoch=0):
    """
    Use kernel density estimates to Monte Carlo outputs to identify p values for the
    mean and mean uncertainty range values. Print the output to the screen.
    """
    keys=["N","NE","E","SE","S","SW","W","NW"]
    date_groups = [solarPhase["max"],solarPhase["min"],
                   naoPhase["max"],naoPhase["min"],
                   ensoPhase["max"],ensoPhase["min"]]
    comp_groups = ["Solar maximum", "Solar minimum", "positive NAO",
                  "negative NAO", "El Nino", "La Nina"]
    # Populate a dataframe with monte carlo outputs for each direction
    mc = pd.DataFrame()
    print("Running Monte Carlo...")
    for n,key in enumerate(keys):
        mc[key] = monteCarlo(df=df, key=key, its=its, give_array=True)
        status(current=n,end_val=len(keys)-1,key=key)
        
    print("\n\nDisplaying mean (uncertainty) and p-value (with range covered by uncertainty):\n")
    for i, date_group in enumerate(date_groups):
        comp = {key : gen_composite(data=df[key],keyDates=date_group,
                                        months=[epoch]) for key in keys}
        means, sem = extractCompositeStatsAllDir(comp, epoch=epoch)
        print("{0}, epoch {1}, δ wind".format(comp_groups[i], epoch))
        for n, key in enumerate(keys):
            kde = stats.gaussian_kde(mc[key])
            mval = means[n][0]
            err = sem[n][0]
            print("{0} {1:3.2f}(±{2:3.2f}) p:{3:1.2} ({4:1.2}--{5:1.2})".format(
                    key,mval, err, kde(mval)[0],kde(mval-err)[0],kde(mval+err)[0]))
        print(end='\n')
    return


def getWindDic(hb_to_direction = False):
    """
    Return a dictionary of weather system origin data, either for HB direction,
    or direction by HB code.
    """
    if hb_to_direction:
        tmp = {'NA':'N','NZ':'N','HNA':'N','HNZ':'N','HB':'N',
               'TRM':'N','NEA':'NE','NEZ':'NE','HFA':'E',
               'HFZ':'E','HNFA':'E','HNFZ':'E','SEA':'SE','SEZ':'SE',
               'SA':'S','SZ':'S','TB':'S','TRW':'S','SWA':'SW','SWZ':'SW',
               'WZ':'W','WS':'W','WA':'W','WW':'W','NWA':'NW','NWZ':'NW',
               'HM':0,'TM':0,'U':0,'BM':0}
        return tmp
    else:
        tmp = {"N":["NA","NZ","HNA","HNZ","HB","TRM"],
                 "NE":["NEA","NEZ"],
                 "E":['HFA','HFZ','HNFA','HNFZ'],
                 "SE":['SEA','SEZ'],
                 "S":['SA','SZ','TB','TRW'],
                 "SW":['SWA','SWZ'],
                 "W":['WZ','WS','WA','WW'],
                 "NW":['NWA','NWZ']}
        return tmp

    
def monteCarlo(df, key, its=1000, size = 11,show_kde = False, give_array=False,
               percentiles=[0.5,2.5,5,50,95,97.5,99.5]):
    """
    For a given number of iterations (its) select sample of n (size)
    from a population given in in dataframe (df[key]). Identify the
    percentiles at specified intervals, and return them in a dictionary.
    Optionally return a sns.kdeplot matplotlib ax instance for plotting.
    This returns percentile values by default, but if give_array is true
    this behaviour is over-ridden, and it provides the whole MC-array.
    Input: df (Dataframe)
        key, iteration, size, percentiles
    Output: Either a dictionary of percentile scores, or matplotlib ax object.
    """
    rnd_means = np.array([np.mean(df[key][random.sample(list(df.index), size)])
                          for n in range(its)])
    if give_array:
        return rnd_means
    if show_kde:
        tmp = pd.DataFrame(data=rnd_means,columns=[key])
        return sns.kdeplot(tmp[key],legend=True,alpha=0.75)
    else: 
        return {pcn: stats.scoreatpercentile(rnd_means, pcn) for pcn in percentiles}
    
    
def season_climatology(data,chatty=False):
    """
    DJF, FMA, JJA, SON months are subset, and stats calculated.
    Input: data, must be a PD Dataframe object with a pd.datetime index
    Output: A dictionary object, of Dic[seazon][direction][mean, sem]
    """
    directions = ["N","NE","E","SE","S","SW","W","NW"]
    seasons = {"DJF":[12,1,2],"MAM":[3,4,5],"JJA":[6,7,8],"SON":[9,10,11]}
    output ={}
    if chatty:
        print("Mean (μ) and SEM frequency (days/month) by Season")
    for season in seasons:
        if chatty:
            print("For {0}".format(season))
        mnmask = [indx in seasons[season] for indx in data.index.month]
        tmp = {}
        for direction in directions:
            mu = np.mean(data[direction][mnmask])
            sem = calc_sem(data[direction][mnmask])
            if chatty:
                print("|------->{0:3s} {1:3.2f}μ, {2:3.2f}sem".format(
                        direction,mu,sem))
            tmp[direction]=(mu,sem)
        output[season]=tmp
    return output


def status(current, end_val, key, bar_length=20):
    ''' Creates a small status bar so users can track the MC easily.
    '''
    percent = float(current) / end_val
    hashes = '#' * int(round(percent * bar_length))
    spaces = ' ' * (bar_length - len(hashes))
    sys.stdout.write("\rMC progress: [{1}] {2}% ({0})".format(
            key ,hashes + spaces, int(round(percent * 100))))
    sys.stdout.flush()



def figure_composite_complex(df, conf_df, naoPhase, ensoPhase, solarPhase):
    """
    Plot of the composite means at specific epochs with Conf Intervals.
    """
    props = dict(boxstyle='round', facecolor='w', alpha=1.0)
    majorLocator   = plt.MultipleLocator(1)
    ticks = np.arange(0, 8, 1)
    fsize=11
    cfilled = ['','N','NE','E','SE','S','SW','W','NW','']
    fig_results = plt.figure()
    fig_results, axs = plt.subplots(2, 4, sharex=True, sharey=True)
    fig_results.set_size_inches(14,8)
    big_ax = fig_results.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', 
                       bottom='off',left='off', right='off')
    big_ax.set_frame_on(False)
    big_ax.grid(False)
    big_ax.set_xlabel(r"$\delta$ weather system origin (days/month)", fontsize=fsize)
    big_ax.set_ylabel(r"Wind direction", fontsize=fsize)
    axs[0, 0].set_ylim([-0.1,7.1])
    axs[0, 0].set_xlim([-10,10])
    epoch = [-20, -10, -5, 0, 5, 10, 20, 19]
    rstart = min(epoch)
    rfin = max(epoch)+1
    keys=["N","NE","E","SE","S","SW","W","NW"]
    for n in range(8):
        comp_smax = {key : gen_composite(data=df[key],keyDates=solarPhase["max"],
                                        months=range(rstart,rfin)) for key in keys}
        comp_smin = {key : gen_composite(data=df[key],keyDates=solarPhase["min"],
                                        months=range(rstart,rfin)) for key in keys}
        comp_elnino = {key : gen_composite(data=df[key],keyDates=ensoPhase["max"],
                                          months=range(rstart,rfin)) for key in keys}
        comp_lanina = {key : gen_composite(data=df[key], keyDates=ensoPhase["min"],
                                          months=range(rstart,rfin)) for key in keys}
        comp_posnao = {key : gen_composite(data=df[key],keyDates=naoPhase["max"],
                                          months=range(rstart,rfin)) for key in keys}
        comp_negnao = {key : gen_composite(data=df[key], keyDates=naoPhase["min"],
                                          months=range(rstart,rfin)) for key in keys}
        smax_means, smax_sem = extractCompositeStatsAllDir(comp_smax, epoch=epoch[n])
        smin_means, smin_sem = extractCompositeStatsAllDir(comp_smin, epoch=epoch[n])
        elnino_means, elnino_sem = extractCompositeStatsAllDir(comp_elnino, epoch=epoch[n])
        lanina_means, lanina_sem = extractCompositeStatsAllDir(comp_lanina, epoch=epoch[n])
        posnao_means, posnao_sem = extractCompositeStatsAllDir(comp_posnao, epoch=epoch[n])
        negnao_means, negnao_sem = extractCompositeStatsAllDir(comp_negnao, epoch=epoch[n])
        if n < 4:       
            axs[0,n].errorbar(posnao_means,ticks, xerr=posnao_sem, fmt='o-',
                        color=sns.color_palette()[1], ms=5.0, linewidth=1.0)
            axs[0,n].errorbar(negnao_means,ticks, xerr=negnao_sem, fmt='o-',
                        color=sns.color_palette()[5], ms=5.0, linewidth=1.0)
            axs[0,n].errorbar(elnino_means,ticks, xerr=elnino_sem, fmt='o-',
                        color=sns.color_palette()[4],ms=5.0, linewidth=1.0)
            axs[0,n].errorbar(lanina_means,ticks, xerr=lanina_sem, fmt='o-',
                    color=sns.color_palette()[3], ms=5.0, linewidth=1.0)
            axs[0,n].errorbar(smax_means,ticks, xerr=smax_sem, fmt='o-',
                            color=sns.color_palette()[2], ms=5.0, linewidth=1.0)             
            axs[0,n].errorbar(smin_means,ticks, xerr=smin_sem, fmt='o-',
                    color=sns.color_palette()[0], ms=5.0, linewidth=1.0)
            axs[0,n].fill_betweenx(ticks, conf_df[2.5], conf_df[97.5],
                             color='gray',linewidth=0.5,alpha=0.3)
            axs[0,n].fill_betweenx(ticks, conf_df[0.5], conf_df[99.5],
                             color='gray',linewidth=0.5,alpha=0.3)
            axs[0,n].set_title(r"Epoch "+str(epoch[n]), fontsize=fsize)
        elif n >= 4 and n < 7:
            axs[1,n-4].errorbar(posnao_means,ticks, xerr=posnao_sem, fmt='o-',
                        color=sns.color_palette()[1], ms=5.0, linewidth=1.0)
            axs[1,n-4].errorbar(negnao_means,ticks, xerr=negnao_sem, fmt='o-',
                        color=sns.color_palette()[5], ms=5.0, linewidth=1.0)
            axs[1,n-4].errorbar(elnino_means,ticks, xerr=elnino_sem, fmt='o-',
                        color=sns.color_palette()[4],ms=5.0, linewidth=1.0)
            axs[1,n-4].errorbar(lanina_means,ticks, xerr=lanina_sem, fmt='o-',
                        color=sns.color_palette()[3], ms=5.0, linewidth=1.0)            
            axs[1,n-4].errorbar(smax_means,ticks, xerr=smax_sem, fmt='o-',
                        color=sns.color_palette()[2], ms=5.0, linewidth=1.0)
            axs[1,n-4].errorbar(smin_means,ticks, xerr=smin_sem, fmt='o-',
                    color=sns.color_palette()[0], ms=5.0, linewidth=1.0)        
            axs[1,n-4].fill_betweenx(ticks, conf_df[2.5], conf_df[97.5],
                             color='gray',linewidth=0.5,alpha=0.3)
            axs[1,n-4].fill_betweenx(ticks, conf_df[0.5], conf_df[99.5],
                             color='gray',linewidth=0.5,alpha=0.3)
            axs[1,n-4].set_title(r"Epoch "+str(epoch[n]), fontsize=fsize)
        if n == 7:
            fig_results.delaxes(axs[1,n-4]) # Erase the last figure block
    axs[0,0].legend(["Pos. NAO","Neg. NAO","El Niño","La Niña","S. Max.","S. Min."],
               loc=6, prop={'size':12}, numpoints=1,markerscale=1.0,
               fancybox=True,frameon=True,bbox_to_anchor=(3.8, -0.5))

    axs[0, 0].set_yticklabels(cfilled, fontsize=fsize) # place labels on x-axis
    axs[1, 0].set_yticklabels(cfilled, fontsize=fsize) # place labels on y-axis
    fig_results.savefig('Figs/epoch_complex.pdf', dpi=300)
    fig_results.show()
    return


def figure_composite_perEpoch(df, conf_df, solarPhase, ensoPhase, epoch=0):
    """
    Old version of the composite plot, shows less info, but still nice to keep.
    """
    # Extract the data from the dataframe of monthly frequency and confidence intervals 
    keys=["N","NE","E","SE","S","SW","W","NW"]
    comp_smax = {key : gen_composite(data=df[key], keyDates=solarPhase["max"]) for key in keys}
    comp_smin = {key : gen_composite(data=df[key], keyDates=solarPhase["min"]) for key in keys}
    comp_elnino = {key : gen_composite(data=df[key], keyDates=ensoPhase["max"]) for key in keys}
    comp_lanina = {key : gen_composite(data=df[key], keyDates=ensoPhase["min"]) for key in keys}
    
    smax_means, smax_sem = extractCompositeStatsAllDir(comp_smax, epoch=epoch)
    smin_means, smin_sem = extractCompositeStatsAllDir(comp_smin, epoch=epoch)
    elnino_means, elnino_sem = extractCompositeStatsAllDir(comp_elnino, epoch=epoch)
    lanina_means, lanina_sem = extractCompositeStatsAllDir(comp_lanina, epoch=epoch)
    
    # Create the plot
    props = dict(boxstyle='round', facecolor='w', alpha=0.75)
    majorLocator   = plt.MultipleLocator(1)
    ticks = np.arange(0, 8, 1)
    cfilled = ['','N','NE','E','SE','S','SW','W','NW','']
    with sns.axes_style("whitegrid"):
        fig_results = plt.figure()
        fig_results,(ax1,ax2)=plt.subplots(2, 1, sharex=True, sharey=True)
        ax1,ax2
        fig_results.set_size_inches(8,6)
        ax1.set_xlim([-0.1,7.1])
        ax1.errorbar(ticks, elnino_means, yerr=elnino_sem, fmt='o', capsize=5.0, color='r',
                     ms=5.0, alpha=0.8, linewidth=2)
        ax1.errorbar(ticks, lanina_means, yerr=lanina_sem, fmt='o', capsize=5.0, color='b',
                     ms=5.0, alpha=0.8, linewidth=2)
        legb=ax1.legend(["El Niño","La Niña"], loc=0, prop={'size':11}, numpoints=1,
                        markerscale=1., frameon=True, fancybox=True)

        ax2.errorbar(ticks, smax_means, yerr=smax_sem, fmt='o', capsize=5.0, color='r',
                     ms=5.0, alpha=0.8, linewidth=2)
        ax2.errorbar(ticks, smin_means, yerr=smin_sem, fmt='o', capsize=5.0, color='b',
                     ms=5.0, alpha=0.8, linewidth=2)
        legb=ax2.legend(["Solar Max.","Solar Min."], loc=0, prop={'size':11}, numpoints=1,
                        markerscale=1., frameon=True, fancybox=True)
        for ax in [ax1,ax2]:
            ax.fill_between(ticks, conf_df[5.0], conf_df[95.0],color='gray',linewidth=0.5,alpha=0.6)
            ax.fill_between(ticks, conf_df[2.5], conf_df[97.5],color='gray',linewidth=0.5,alpha=0.3)
            ax.fill_between(ticks, conf_df[0.5], conf_df[99.5],color='gray',linewidth=0.5,alpha=0.3)
        ax1.set_xticklabels(cfilled, fontsize=11) # place labels on x-axis
        ax1.set_title(r"Epoch "+str(epoch), fontsize=11)
        ax1.set_ylabel(r"$\delta$ weather system origin (days/month)", fontsize=11)
        ax2.set_ylabel(r"$\delta$ weather system origin (days/month)", fontsize=11)
        ax2.set_xlabel(r"Direction", fontsize=11)
        fig_results.savefig('Figs/epoch'+str(epoch)+'_composite.pdf', dpi=300)
        fig_results.show()
        return

    
def fig_distribution(data):
    """
    KDE and CDF of wind frequency for main compass directions
    Input: Pandas dataframe of monthly weather system origin frequency
    Output: figure
    """
    fsize=11 # <-- Change font-size here
    fig_kde, (ax1, ax2) = plt.subplots(2, sharex=True)
    fig_kde.set_size_inches(4,8) # <-- Change size of plot here 
    for key in ["N","E","S","W"]:
        ax1 = sns.kdeplot(data[key], ax=ax1)
    ax1.set_xlim(0,30)
    for key in ["N","E","S","W"]:
        ax2 = sns.kdeplot(data[key],cumulative=True,
                            ax=ax2)
    ax2.set_xlim(0,30)
    ax1.set_ylabel(r"KDE", fontsize=fsize)
    ax2.set_ylabel(r"CDF", fontsize=fsize)
    ax2.set_xlabel(r"Frequency (days/month)", fontsize=fsize)
    ax1.set_title("Frequency distribution", fontsize=fsize) 
    fig_kde.savefig("Figs/freq_dist.pdf",dpi=300)
    fig_kde.show()
    return


def figure_forcing_TS(data):
    """
    Show NAO index, MEI and Wolf sunspot number as monthly
    time series data.
    Input: Pandas dataframe containing the NAO, MEI & Wolf data
    Output: Plot
    """
    fsize = 11 # <-- Font size kwarg
    fig_TS,(ax1, ax2, ax3)=plt.subplots(3,1,sharex=True)
    fig_TS.set_size_inches(9,6)
    
    ax1.plot(data.index,data.NAO,
             lw=1.,color=sns.color_palette()[1])    
    ax2.plot(data.index,data.MEI,
             lw=1.,color=sns.color_palette()[0])
    ax3.plot(data.index,data.Wolf,
             lw=1.,color=sns.color_palette()[2])
    ax1.set_ylabel("NAO index",size=fsize)    
    ax2.set_ylabel("ENSO index",size=fsize)
    ax3.set_ylabel("Sunspot number",size=fsize)
    ax3.set_xlabel("Year",size=fsize)
    fig_TS.savefig('Figs/Forcing_TS.pdf', dpi=300)
    fig_TS.show()
    return


def figure_seasons(data):
    """
    Create a violin plot of the wearth system origin during the 
    DJF, FMA, JJA, SON months. This replaces an polar version of this
    figure which remains in the figure_SeasonalClimo() function.
    Input: data, must be a PD Dataframe object with a pd.datetime index
    Output: None. Only saves a plot to pdf, and also displays to screen.
    """
    fsize = 11 # <-- specify font size
    fig_TS,(ax1,ax2,ax3,ax4)=plt.subplots(4,1, sharey=True,sharex=True)
    axobs =[ax1,ax2,ax3,ax4]
    fig_TS.set_size_inches(8,6)
    big_ax = fig_TS.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', 
                       bottom='off',left='off', right='off')
    big_ax.set_frame_on(False)
    big_ax.set_ylabel(r"Frequency (days/month)", labelpad=20, fontsize=fsize)
    big_ax.set_xlabel(r"Season", fontsize=fsize)
    big_ax.set_title("Weather system origin by season", fontsize=fsize)
    big_ax.grid(False)

    directions = ["N","E","S","W"]
    seasons = {"DJF":[12,1,2],"MAM":[3,4,5],"JJA":[6,7,8],"SON":[9,10,11]}
    for j, direction in enumerate(directions):
        szn_frame = pd.DataFrame()
        for season in ["DJF","MAM","JJA","SON"]:
            mnmask = [indx in seasons[season] for indx in data.index.month]
            szn_frame[season] = data[direction][mnmask].values
        clr=sns.color_palette()[j]
        sns.violinplot(szn_frame, cut=0, linewidth=1.0, inner="quartiles",
                       color=clr,ax=axobs[j])
        axobs[j].set_ylabel(direction)
    fig_TS.savefig('Figs/Seasonal_violin.pdf', dpi=300)
    fig_TS.show()
    return


def fig_forcing_composite(df,naoPhase,ensoPhase,solarPhase):
    """
    Create a composite figure of the NAO, ENSO and solar forcing.
    Use the composite function to create the values within 
    this routine.
    Output: figure 
    """
    fsize = 11 # <-- Font size kwarg
    fig_comp,(ax1, ax2, ax3)=plt.subplots(3,1,sharex=True)
    ax1,ax2
    fig_comp.set_size_inches(9,6)
    mrange = range(-25,26)  # <-- Set composite window here
    cols = [sns.color_palette()[0],sns.color_palette()[2]]
    
    for n, key in enumerate(["min","max"]):
        composite = gen_composite(data=df.NAO, keyDates=naoPhase[key],
                                  months=mrange)
        ax1.errorbar(composite["epoch"], composite["means"],
                 xerr=None,yerr=composite["sem"],color=cols[n])
    
    for n, key in enumerate(["min","max"]):
        composite = gen_composite(data=df.MEI, keyDates=ensoPhase[key],
                                  months=mrange)
        ax2.errorbar(composite["epoch"], composite["means"],
                 xerr=None,yerr=composite["sem"],color=cols[n])
        
    for n, key in enumerate(["min","max"]):
        composite = gen_composite(data=df.Wolf, keyDates=solarPhase[key],
                                  months=mrange)
        ax3.errorbar(composite["epoch"], composite["means"],
                 xerr=None,yerr=composite["sem"],color=cols[n])
    ax1.set_xlim(min(mrange),max(mrange))
    ax1.set_ylabel("NAO index",size=fsize)
    ax2.set_ylabel("ENSO index",size=fsize)
    ax3.set_ylabel("Sunspot number",size=fsize)
    ax3.set_xlabel("Months since peak values",size=fsize)
    fig_comp.savefig('Figs/composite_forcing.pdf', dpi=300)
    fig_comp.show()
    return


def figure_montecarlo_kde(df,its=1000):
    """
    Create a KDE figure of the Monte Carlo-generated samples
    for wind direction.
    """
    fig_kde = plt.figure()
    fig_kde.set_size_inches(4,4)
    ax1 = fig_kde.add_subplot(111)
    # nb. order of keys are intentional so colors
    # of cardinal directions are consistent.
    keys=["N","E","S","W","NE","SE","SW","NW"]
    for n,key in enumerate(keys):
        status(current=n, end_val=len(keys)-1, key=key)
        #print("\rCalculating {0}...".format(key),end="")
        ax1 = monteCarlo(df=df, its=its, show_kde=True, key=key)
    print("\rFinished MC calculation")
    ax1.set_title("Monte Carlo distributions")
    ax1.set_ylabel("KDE")
    ax1.set_xlabel(r"$\delta$ weather system origin (days/month)")
    fig_kde.savefig('Figs/MC_kde.pdf', dpi=300)
    fig_kde.show()
    return


def figure_MonthlyTS(df):
    """
    Create a multi-panel simple time-series figure for monthly data.
    Input: Pandas Dataframe of Monthly data with directions as columns.
    Output: Time series figure
    """
    fsize = 11 # <-- specify font size
    fig_TS,(ax1,ax2,ax3,ax4)=plt.subplots(4,1,sharey=True,sharex=True)
    axobs =[ax1,ax2,ax3,ax4]
    fig_TS.set_size_inches(8,6)
    big_ax = fig_TS.add_subplot(111)
    big_ax.set_axis_bgcolor('none')
    big_ax.tick_params(labelcolor='none', top='off', 
                       bottom='off',left='off', right='off')
    big_ax.set_frame_on(False)
    big_ax.set_ylabel(r"Frequency (days/month)", fontsize=fsize)
    big_ax.set_xlabel(r"Year", fontsize=fsize)
    big_ax.set_title("Monthly weather system origin", fontsize=fsize)
    big_ax.grid(False)
    #props = ['b.','g.','r.','m.'] # <- some plot color and marker
    for n, wind in enumerate(["N","E","S","W"]):
        axobs[n].plot(df.index, df[wind],".",
                      color=sns.color_palette()[n],ms=3.0)
        leg1=axobs[n].legend([wind], loc=2,prop={'size':fsize},
                                   numpoints=1,markerscale=2.0)
        leg1.get_frame().set_alpha(0.9)
    fig_TS.savefig('Figs/Monthly_freq.pdf', dpi=300)
    fig_TS.show()
    return


def figure_SeasonalClimo(data):
    """
    Uses the season_climatology() function of this file to calculate a
    dictionary object of mean wind frequency by direction by season.
    Then plots the results in a polar plot.
    
    Input: Monthly wind data in pandas dataframe format (with columns
        as different compass directions.
    Output: A polar climatology plot.
    """
    # Get seasonal data
    szn_dat = season_climatology(data, chatty=False) 

    # Set up plot
    fig_pl = plt.figure()
    fig_pl.set_size_inches(4,4)
    ax1 = fig_pl.add_subplot(111, polar=True)

    # Set properties and lists
    fsize = 11    # <-- Font size
    ylabs = ['','4','6','8','10']
    ticks = [0.0,45.,90.,135.,180.,225.,270.,315.,360.]
    seasons = ["DJF","MAM","JJA","SON"]
    directions = ['E','NE','N','NW','W','SW','S','SE','']    # Blank at end, as will be plot label
    rad_ticks = [n*(np.pi/180.) for n in ticks]
    for n, season in enumerate(seasons):
        sz_means = [szn_dat[season][direction][0] for direction in directions[0:-1]]
        sz_means.append(sz_means[0])                        # Close loop for the polar plot

        ax1.plot(rad_ticks,sz_means,"-",color=sns.color_palette()[n],
                 linewidth=1.5,alpha=0.9)                   # Add info to plot in a loop

    leg1=ax1.legend(["DJF","MAM","JJA","SON"], loc=2,prop={'size':fsize},
                    numpoints=1,markerscale=5.,frameon=True,fancybox=True)
    leg1.get_frame().set_alpha(1.0)
    ax1.set_xticklabels(directions, fontsize=fsize)
    ax1.set_yticklabels(ylabs, fontsize=fsize)
    ax1.set_xlabel("Weather system origin (days/month)", fontsize=fsize)
    ax1.grid(True)
    fig_pl.subplots_adjust(left=0.2, bottom=0.13, right=0.95, 
                           top=0.92, wspace=0, hspace=0)
    plt.savefig('Figs/Seasonal_climo.pdf', dpi=300)
    fig_pl.show()
    return