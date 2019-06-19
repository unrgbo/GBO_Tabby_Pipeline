import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Tabby_pipe import get_paths


dpath, opath, calpath, execpath, npath, home = get_paths()
simJD = 2457800

def ref(filters=['g','r','i','V'], average = False, write=False):

    if not average:
        for band in filters:
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.set_title("Reference stars: {}' filter\n 2018 March - current".format(band), fontsize=16)
            ax.set_xlabel('Time (JD-2457800)', fontsize=14)
            ax.set_ylabel('Relative Flux', fontsize=14)
            ax.axhline(y=1, linestyle='--', linewidth=1.5, color='black')
            ax.minorticks_on()
            ax.tick_params(which='major', labelsize=12, length=6, width=1.5, direction='in', right=True, top=True)
            ax.tick_params(which='minor', length=3, width=1.5, direction='in', right=True, top=True)
            fig.tight_layout(pad=1)
            ax.set_ylim(0.95, 1.04)
            ax.set_xlim(390, 685)


            sum = pd.read_csv(opath + band + '_summary.csv')
            refmean = sum[['ref1mean', 'ref2mean', 'ref3mean']].copy()
            refmean['mean'] = refmean.mean(axis = 1)
            totalnorm = refmean['mean'].mean()
            refstd = sum[['ref1err', 'ref2err', 'ref3err']].copy()
            refstd['error'] = refstd.mean(axis = 1)

            refnorm = []
            n = 0

            for i in range(3):
                refnorm.append(sum['ref{}mean'.format(i + 1)])
                refnorm[i] = np.mean(refnorm[i])
                n += .01
                ax.errorbar(sum.loc[:, 'median_JD'] - simJD, (sum.loc[:, 'ref{}mean'.format(i + 1)] / refnorm[i]) - n,
                            yerr=sum.loc[:, 'ref{}err'.format(i + 1)] / refnorm[i],
                            markersize=4, linewidth=1, fmt='o', label="ref{}".format(i + 1))
            ax.errorbar(sum.loc[:, 'median_JD'] - simJD, refmean.loc[:, 'mean'] / totalnorm,
                        yerr=refstd.loc[:, 'error'] / totalnorm,
                        markersize=4, linewidth=1, fmt='o', label="Average")

            ax.legend(loc=1)

            if write:
                plt.savefig(opath + 'plots/References_{}.png'.format(band), dpi=300)

    else:
        for band in filters:
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.set_title("Reference stars: {}' filter\n 2018 March - current".format(band), fontsize=16)
            ax.set_xlabel('Time (JD-2457800)', fontsize=14)
            ax.set_ylabel('Relative Flux', fontsize=14)
            ax.axhline(y=1, linestyle='--', linewidth=1.5, color='black')
            ax.minorticks_on()
            ax.tick_params(which='major', labelsize=12, length=6, width=1.5, direction='in', right=True, top=True)
            ax.tick_params(which='minor', length=3, width=1.5, direction='in', right=True, top=True)
            fig.tight_layout(pad=1)
            ax.set_ylim(0.95, 1.04)
            ax.set_xlim(390, 685)


            sum = pd.read_csv(opath + band + '_summary.csv')
            refmean = sum[['ref1mean', 'ref2mean', 'ref3mean']].copy()
            refmean['mean'] = refmean.mean(axis=1)
            totalnorm = refmean['mean'].mean()
            refstd = sum[['ref1err', 'ref2err', 'ref3err']].copy()
            refstd['error'] = refstd.mean(axis=1)

            refnorm = []
            n = 0

            for i in range(3):
                refnorm.append(sum['ref{}mean'.format(i + 1)])
                refnorm[i] = np.mean(refnorm[i])
            ax.errorbar(sum.loc[:, 'median_JD'] - simJD, refmean.loc[:, 'mean'] / totalnorm,
                        yerr=refstd.loc[:, 'error'] / totalnorm,
                        markersize=4, linewidth=1, fmt='o', label="{} Filter".format(band))

            ax.legend(loc=1)

            if write:
                plt.savefig(opath + 'plots/RefAvg_{}.png'.format(band), dpi=300)

def best_fit(filters=['g','r','i','V'], individual=False, write=False):
    icolor = '#860000'
    Vcolor = '#99ff00'

    gcolor = 'blue'
    rcolor = 'orange'
    # rcolor = 'green'
    # icolor = 'purple'
    Vcolor = 'green'

    for band in filters:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.set_xlabel('Time (JD-2457800)', fontsize=14)
        ax.set_ylabel('Relative Flux', fontsize=14)
        ax.set_title("2018 June - December".format(band), fontsize=16)
        ax.axhline(y=1, linestyle='--', linewidth=1.5, color='black')
        ax.minorticks_on()
        ax.tick_params(which='major', labelsize=12, length=6, width=1.5, direction='in', right=True, top=True)
        ax.tick_params(which='minor', length=3, width=1.5, direction='in', right=True, top=True)
        fig.tight_layout(pad=1.5)
        ax.set_ylim(0.96, 1.04)
        ax.set_xlim(470, 680)

        sum = pd.read_csv(opath + band + '_summary.csv')
        tabby = sum[['mean', 'error']].copy()
        tabbynorm = sum['mean'].mean()

        if band == 'g' or band == 'r':
            x = sum.loc[11:, 'median_JD'] - simJD
            y = tabby.loc[11:, 'mean']/tabbynorm
            yerr = tabby.loc[11:, 'error']/tabbynorm
        elif band == 'V':
            x = sum.loc[4:, 'median_JD'] - simJD
            y = tabby.loc[4:, 'mean']/tabbynorm
            yerr = tabby.loc[4:, 'error'] / tabbynorm
        else:
            x = sum.loc[:, 'median_JD'] - simJD
            y = tabby.loc[:, 'mean'] / tabbynorm
            yerr = tabby.loc[:, 'error'] / tabbynorm

        linex = np.linspace(400,700, 300)

        m,b = np.polyfit(x, y, 1)
        print m, b

        percentday = m*100
        slope = percentday *365

        ax.errorbar(x, y, yerr=yerr,
                    markersize=4, linewidth=1, fmt='o', label="{}' Filter".format(band),
                    color = '#ff7500')
        ax.plot(linex, (m*linex) + b, color='black', label='Slope = {} percent/year'.format(str(slope)[:4]))

        plt.show()

        ax.legend(loc=4)

        if write:
            plt.savefig(opath + 'plots/Best_Fit_{}.png'.format(band), dpi=300)


