# helper.py
# Defines helper functions

import numpy as np # 1.13.3
from scipy.integrate import odeint # 1.0.0
import scipy.optimize as op
import scipy.stats as stats
import matplotlib.pyplot as plt # 2.1.1
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import pandas as pd # 0.22.0
import corner # 2.0.1
import progressbar # 3.34.3
import seaborn as sns # 0.8.1
from cycler import cycler # 0.10.0

DIR_DATA = './data/'
DIR_PLOTS = './plots/'
DIR_OUT = './output/'

# This helper file defines plotting functions for the jupyter notebooks


######### 1. PLOT FORMATTING #########

def formatplot(ax,xlabel,ylabel,xlim,ylim,logx=False,logy=False,logxy=False,symlogx=False):

    # Set titles and labels
    #ax.set_title('Plot title')
    if xlabel!=False:
        ax.set_xlabel(xlabel, labelpad=12)
    if ylabel!=False:    
        ax.set_ylabel(ylabel, labelpad=12)

    # Set axis limits
    if xlim!=False:
        ax.set_xlim(xlim)
    if ylim!=False:
        ax.set_ylim(ylim)

    # Set tick values
    # ax.set_xticks([0,0.5,1])
    # ax.set_yticks([0,2,4,6,8])

    # Set line thicknesses
    #ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%1.e"))
    #ax.axhline(linewidth=2, color='k')      
    #ax.axvline(linewidth=2, color='k')
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2) 

    # Set ticks
    if logx==True:
        ax.set_xscale("log")

    elif logy==True:
        ax.set_yscale("log")

    elif logxy==True:
        ax.set_xscale("log")
        ax.set_yscale("log")
    
    elif symlogx==True:
        ax.set_xscale("symlog",linthreshx=1e-4)
        ax.set_yscale("log")

    else:
        minorLocatorx=AutoMinorLocator(2) # Number of minor intervals per major interval
        minorLocatory=AutoMinorLocator(2)
        ax.xaxis.set_minor_locator(minorLocatorx)
        ax.yaxis.set_minor_locator(minorLocatory)

    ax.tick_params(which='major', width=2, length=8, pad=9,direction='in',top='on',right='on')
    ax.tick_params(which='minor', width=2, length=4, pad=9,direction='in',top='on',right='on')


######### 2. PLOT TRACES #########


def plottraces(y1,yerr1,y2,yerr2,y3,yerr3,y4,yerr4,
               titration,name,model_single,model_dual,DIR_PLOTS,
               numberofmodeltraces=50,samples=0,quant=None):
	# Plot dose response traces (Figure 5)

    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8,6)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, axarr=plt.subplots()
    plt.subplots_adjust(left=0.3,bottom=0.2,right=0.95,top=0.9)
    
    # Plot data
    axarr.errorbar(titration,y1,yerr=yerr1,fmt='o',ms=12,label='BCB(+)ADD',color='#C4122C')
    axarr.errorbar(titration,y2,yerr=yerr2,fmt='^',ms=12,label='BCB(-)ADD',color='#228863')
    if (y3.any()!=None and y4.any()!=None and yerr3.any()!=None and yerr4.any()!=None):
        axarr.errorbar(titration/2,y3,yerr=yerr3,fmt='s',ms=12,label='BCB',color='#28A0A3')
        axarr.errorbar(titration/2,y4,yerr=yerr4,fmt='D',ms=12,label='ADD',color='#0A719F')
    else:
        pass

    modelscale=np.linspace(0,1,100)

    # Plot model
    
    if quant!=None:
        
        quant1=quant[0]
        quant2=quant[1]
        quant3=quant[2]
        quant4=quant[3]
        
        axarr.fill_between(modelscale,quant1[0],quant1[2],color='#C4122C',alpha=0.1)
        axarr.fill_between(modelscale,quant2[0],quant2[2],color='#228863',alpha=0.1)
        axarr.fill_between(modelscale,quant3[0],quant3[2],color='#28A0A3',alpha=0.1)
        axarr.fill_between(modelscale,quant4[0],quant4[2],color='#0A719F',alpha=0.1)
        
        axarr.plot(modelscale,quant1[1],'-',color='#C4122C',alpha=1,lw=1.5)
        axarr.plot(modelscale,quant2[1],'-',color='#228863',alpha=1,lw=1.5)
        axarr.plot(modelscale,quant3[1],'-',color='#28A0A3',alpha=1,lw=1.5)
        axarr.plot(modelscale,quant4[1],'-',color='#0A719F',alpha=1,lw=1.5)

        axarr.set_yticks([0,2e3,4e3])
        axarr.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)

        formatplot(axarr,'ZF DNA ratio','deGFP (RFU)', xlim=([-0.05,1.05]),ylim=[0,4600])

        axarr.legend(loc='best', fontsize=18,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)
        
    elif samples.any()!=None:
        for A,C0,K1,K2,Erp1,Erp2,Er1r2 in samples[np.random.randint(len(samples), 
                                                  size=numberofmodeltraces)]:
            ypred4=model_single(A,C0,K1,modelscale,Erp1)
            ypred3=model_single(A,C0,K2,modelscale,Erp2)
            ypred2=model_dual(A,C0,K1,K2,modelscale/2,modelscale/2,Erp1,Erp2,0)
            ypred1=model_dual(A,C0,K1,K2,modelscale/2,modelscale/2,Erp1,Erp2,Er1r2)            
            axarr.plot(modelscale,ypred1,'-',color='#C4122C',alpha=0.05)
            axarr.plot(modelscale,ypred2,'-',color='#228863',alpha=0.05)
            axarr.plot(modelscale,ypred3,'-',color='#28A0A3',alpha=0.05)
            axarr.plot(modelscale,ypred4,'-',color='#0A719F',alpha=0.05)
 
        formatplot(axarr,'ZF DNA ratio','deGFP (RFU)', xlim=([-0.05,1.05]),ylim=[0,4600])

        axarr.legend(loc='best', fontsize=15,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)

def plottracesAA(y1,yerr1,y2,yerr2,titration,name,model_dual,
                 DIR_PLOTS,numberofmodeltraces=50,samples=0,quant=None):

    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8,6)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, axarr=plt.subplots()
    plt.subplots_adjust(left=0.3,bottom=0.2,right=0.95,top=0.9)
    
    # Plot data
    axarr.errorbar(titration,y1,yerr=yerr1,fmt='o',ms=12,label='AA-LZ(+)',color='#C4122C')
    axarr.errorbar(titration,y2,yerr=yerr2,fmt='^',ms=12,label='AA-LZ(-)',color='#228863')

    modelscale=np.linspace(0,1,100)

    # Plot model
    
    if quant!=None:
        
        quant1=quant[0]
        quant2=quant[1]
     
        axarr.fill_between(modelscale,quant1[0],quant1[2],color='#C4122C',alpha=0.1)
        axarr.fill_between(modelscale,quant2[0],quant2[2],color='#228863',alpha=0.1)      
        axarr.plot(modelscale,quant1[1],'-',color='#C4122C',alpha=1,lw=1.5)
        axarr.plot(modelscale,quant2[1],'-',color='#228863',alpha=1,lw=1.5)

        axarr.set_yticks([0,2e3,4e3])
        axarr.ticklabel_format(style='sci', axis='y', scilimits=(0,0),useMathText=True)

        formatplot(axarr,'ZF DNA ratio','deGFP (RFU)', xlim=([-0.05,1.05]),ylim=[0,4600])

        axarr.legend(loc='best', fontsize=18,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)
        
    elif samples.any()!=None:
        for A,C0,K,Erp,Er1r2 in samples[np.random.randint(len(samples), size=numberofmodeltraces)]:
            ypred2=model_dual(A,C0,K,K,modelscale/2,modelscale/2,Erp,Erp,0)
            ypred1=model_dual(A,C0,K,K,modelscale/2,modelscale/2,Erp,Erp,Er1r2)            
            axarr.plot(modelscale,ypred1,'-',color='#C4122C',alpha=0.05)
            axarr.plot(modelscale,ypred2,'-',color='#228863',alpha=0.05)
 
        formatplot(axarr,'ZF DNA ratio','deGFP (RFU)', xlim=([-0.05,1.05]),ylim=[0,4600])

        axarr.legend(loc='best', fontsize=15,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)
        
def plottraces2(data,parameternames,nwalkers,niterations,ZFname,DIR_PLOTS):

    numberofplots=data.shape[1]

    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8.27,5.83)} #figure size in inches. A4=11.7x8.3. 
    font_options={'size':'18','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)


    # Call plots
    if numberofplots>1:
        f, axarr=plt.subplots(numberofplots)
        for i in range(0,numberofplots):
            for j in range(1,nwalkers+1):
                axarr[i].plot(np.arange(niterations),data[niterations*j-niterations:niterations*j,i],'k-',lw=0.5)

            formatplot(axarr[i],False,parameternames[i],xlim=False,ylim=False)
    else:
        f, axarr=plt.subplots()
        for i in range(1,nwalkers+1):
            axarr.plot(np.arange(niterations),data[niterations*i-niterations:niterations*i],'k-',lw=0.5)

        formatplot(axarr,False,parameternames[0],xlim=False,ylim=False)
    axarr[numberofplots-1].set_xlabel('Iterations', labelpad=12)
    plt.savefig(DIR_PLOTS+ZFname+'trace.pdf',dpi=my_dpi,bbox_inches='tight')

def plottracesHill(y1,yerr1,y2,yerr2,y3,yerr3,y4,yerr4,
                   titration,name,Hill,DIR_PLOTS,
                   numberofmodeltraces=50,samples=0,quant=None):

    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8.27,5.83)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, axarr=plt.subplots()
    plt.subplots_adjust(left=0.25,bottom=0.2,right=0.95,top=0.95)
    
    # Plot data
    axarr.errorbar(titration,y1,yerr=yerr1,fmt='o',ms=7,label='BCB+ADD coop',color='#C4122C')
    axarr.errorbar(titration,y2,yerr=yerr2,fmt='^',ms=7,label='BCB+ADD non-coop',color='#228863')
    if (y3.any()!=None and y4.any()!=None and yerr3.any()!=None and yerr4.any()!=None):
        axarr.errorbar(titration/2,y3,yerr=yerr3,fmt='s',ms=7,label='BCB',color='#28A0A3')
        axarr.errorbar(titration/2,y4,yerr=yerr4,fmt='D',ms=7,label='ADD',color='#0A719F')
    else:
        pass

    modelscale=np.linspace(0,1,100)

    # Plot model
    
    if quant!=None:
        
        quant1=quant[0]

        axarr.fill_between(modelscale,quant1[0],quant1[2],color='k',alpha=0.1)        
        axarr.plot(modelscale,quant1[1],'-',color='k',alpha=1,lw=1.5)


        formatplot(axarr,'ZF DNA ratio','deGFP (RFU)', xlim=([-0.05,1.05]),ylim=[0,4600])

        axarr.legend(loc='best', fontsize=15,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)
        
    elif samples.any()!=None:
        for y_0,y_1,K,n in samples[np.random.randint(len(samples), size=numberofmodeltraces)]:
            ypred1=Hill(y_0,y_1,modelscale,K,n)           
            axarr.plot(modelscale,ypred1,'-',color='k',alpha=0.05)
 
        formatplot(axarr,'ZF DNA ratio','deGFP (RFU)', xlim=([-0.05,1.05]),ylim=[0,4600])

        axarr.legend(loc='best', fontsize=15,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)
    
def plottraces3L(df3,x,y0,y1,yerr1,y2,yerr2,name,DIR_PLOTS,
                 helixdist,modelhelix,
                 numberofmodeltraces=50,samples=0,quant=None):

    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8,6)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, axarr=plt.subplots()
    plt.subplots_adjust(left=0.3,bottom=0.2,right=0.95,top=0.9)
    
    # Plot data    
    axarr.errorbar(df3['Spacing'],df3['Control'],yerr=df3['Cont_err'],fmt='s',ms=12,color='k',label='Unrepressed')    
    axarr.errorbar(df3['Spacing'],df3['Non-cognate'],yerr=df3['NC_err'],fmt='^',ms=12,color='#0A719F',label='Non-cooperative')
    axarr.errorbar(df3['Spacing'],df3['Cognate'],yerr=df3['C_err'],fmt='o',ms=12,color='#C4122C',label='Cooperative')

    # Plot model
    
    if quant!=None:
        
        quant4=quant[3]
        quant5=quant[4]
        
        axarr.fill_between(x,quant4[0],quant4[2],color='#C4122C',alpha=0.1,label='__nolegend__')
        axarr.fill_between(x,quant5[0],quant5[2],color='#228863',alpha=0.1,label='__nolegend__')

        axarr.plot(x,quant4[1],'-',color='#C4122C',alpha=1,lw=1.5,label='__nolegend__')
        axarr.plot(x,quant5[1],'-',color='#228863',alpha=1,lw=1.5,label='__nolegend__')


        formatplot(axarr,'spacing (bp)','deGFP (RFU)', xlim=([5,32]),ylim=False)

        axarr.legend(loc='best', fontsize=15,numpoints=1)
        plt.savefig('plots/'+name+'.pdf',dpi=my_dpi,transparent=True)
        
    elif samples.any()!=None:
        for lnlamb,phi,lnR0 in samples[np.random.randint(len(samples), size=numberofmodeltraces)]:
 
            lamb=np.exp(lnlamb)
            R0=np.exp(lnR0)
    
            # Fixed parameters from dose response experiments
            A=y0
            C0=0.66070446353476231
            r1=R0/0.16746036268850761
            r2=R0/0.0082708295083317035
            E10=1.4349132332094823
            E2=1.3110190282670602
            E120=-3.4804403863425248
    
            xmod=x
            E1,E12=helixdist(xmod,E10,E120,lamb,phi)
            p_c,p_nc,p_0=modelhelix(y0,C0,r1,r2,E1,E2,E12)
        
            axarr.plot(x,p_0,'-',color='k',alpha=0.05,label='__nolegend__')
            axarr.plot(x,p_nc,'-',color='#228863',alpha=0.05,label='__nolegend__')
            axarr.plot(x,p_c,'-',color='#C4122C',alpha=0.05,label='__nolegend__')
            
 
        formatplot(axarr,'spacing (bp)','deGFP (RFU)', xlim=([5,32]),ylim=False)
        axarr.legend(loc='best', fontsize=15,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)

def plottraces3FR(x,y0,y1,yerr1,y2,yerr2,name,DIR_PLOTS,
                  helixdist,modelhelix,
                  numberofmodeltraces=50,samples=0,quant=None):

    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8,6)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, axarr=plt.subplots()
    plt.subplots_adjust(left=0.3,bottom=0.2,right=0.95,top=0.9)
    
    # Plot data    
    axarr.errorbar(x,y1,yerr=yerr1,fmt='o',ms=12,color='#C4122C',label='(+)')
    axarr.errorbar(x,y2,yerr=yerr2,fmt='^',ms=12,color='#228863',label='(-)')

    # Plot model
    
    if quant!=None:
        
        quant1=quant[0]
        quant2=quant[1]
        
        axarr.fill_between(x,quant1[0],quant1[2],color='#C4122C',alpha=0.1,label='__nolegend__')
        axarr.fill_between(x,quant2[0],quant2[2],color='#228863',alpha=0.1,label='__nolegend__')
        
        axarr.plot(x,quant1[1],'-',color='#C4122C',alpha=1,lw=1.5,label='__nolegend__')
        axarr.plot(x,quant2[1],'-',color='#228863',alpha=1,lw=1.5,label='__nolegend__')


        formatplot(axarr,'spacing (bp)','fold-repression', xlim=([5,32]),ylim=([0,15]))
        axarr.legend(loc='best', fontsize=18,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)
        
    elif samples.any()!=None:
        for lnlamb,phi,lnR0 in samples[np.random.randint(len(samples), size=numberofmodeltraces)]:
 
            lamb=np.exp(lnlamb)
            R0=np.exp(lnR0)
    
            # Fixed parameters from dose response experiments
            A=y0
            C0=0.66070446353476231
            r1=R0/0.16746036268850761
            r2=R0/0.0082708295083317035
            E10=1.4349132332094823
            E2=1.3110190282670602
            E120=-3.4804403863425248
    
            xmod=x
            E1,E12=helixdist(xmod,E10,E120,lamb,phi)
            p_c,p_nc,p_0=modelhelix(y0,C0,r1,r2,E1,E2,E12)
        
            axarr.plot(x,p_0/p_c,'-',color='#C4122C',alpha=0.05,label='__nolegend__')
            axarr.plot(x,p_0/p_nc,'-',color='#228863',alpha=0.05,label='__nolegend__')
            
        formatplot(axarr,'spacing (bp)','fold-repression', xlim=([5,32]),ylim=False)
        axarr.legend(loc='best', fontsize=15,numpoints=1)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)    
        
def plottraces3CR(x,y0,y1,yerr1,y2,yerr2,name,DIR_PLOTS,
                  helixdist,modelhelix,
                  numberofmodeltraces=50,samples=0,quant=None):

    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8,6)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, axarr=plt.subplots()
    plt.subplots_adjust(left=0.3,bottom=0.2,right=0.95,top=0.9)
    
    # Plot data    
    axarr.errorbar(x,y1/y2,yerr=y1/y2*np.sqrt((yerr1/y1)**2+(yerr2/y2)**2),fmt='o',ms=12,color='k')

    # Plot model
    
    if quant!=None:
        
        quant3=quant[2]

        
        axarr.fill_between(x,quant3[0],quant3[2],color='k',alpha=0.1,label='__nolegend__')
        
        axarr.plot(x,quant3[1],'-',color='k',alpha=1,lw=1.5,label='__nolegend__')


        formatplot(axarr,'spacing (bp)','cooperativity ratio', xlim=([5,32]),ylim=False)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)
        
    elif samples.any()!=None:
        for lnlamb,phi,lnR0 in samples[np.random.randint(len(samples), size=numberofmodeltraces)]:
 
            lamb=np.exp(lnlamb)
            R0=np.exp(lnR0)
    
            # Fixed parameters from dose response experiments
            A=y0
            C0=0.66070446353476231
            r1=R0/0.16746036268850761
            r2=R0/0.0082708295083317035
            E10=1.4349132332094823
            E2=1.3110190282670602
            E120=-3.4804403863425248
    
            xmod=x
            E1,E12=helixdist(xmod,E10,E120,lamb,phi)
            p_c,p_nc,p_0=modelhelix(y0,C0,r1,r2,E1,E2,E12)
        
            axarr.plot(x,p_nc/p_c,'-',color='k',alpha=0.05,label='__nolegend__')
            
        formatplot(axarr,'spacing (bp)','cooperativity ratio', xlim=([5,32]),ylim=False)
        plt.savefig(DIR_PLOTS+name+'.pdf',dpi=my_dpi,transparent=True)      
        

######### 3. BOXPLOTS #########


def boxplots(ZFname,parameternames,nwalkers,iterations,tburn,DIR_PLOTS,DIR_OUT):
    
    # Read data
    df1=pd.read_csv(DIR_OUT+'samplesout_'+ZFname+'.csv',delimiter=',')
    data=np.zeros(df1.shape[0]*(df1.shape[1]-1)).reshape(df1.shape[0],(df1.shape[1]-1))
    for i in range(0,int(df1.shape[1]-1)):
        data[:,i]=np.array(df1.iloc[:,i+1]) # Put dataframe into array. Dataframe has no. columns = no. parameters.

    # Burn-in time correction

    data2=np.zeros((df1.shape[0]-tburn*nwalkers)*(df1.shape[1]-1)).reshape((df1.shape[0]-(tburn*nwalkers)),(df1.shape[1]-1))
    for i in range(0,int(df1.shape[1]-1)):
        for j in range(1,nwalkers+1):
            data2[(iterations-tburn)*(j-1):(iterations-tburn)*(j),i]=np.array(df1.iloc[iterations*j-iterations+tburn:iterations*j,i+1])
    df=pd.DataFrame(data2,columns=parameternames)
    
    # Plot data
    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8.27,5.83)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, ax=plt.subplots()
    plt.subplots_adjust(left=0.25,bottom=0.2,right=0.95,top=0.95)
    
    
    ax=sns.boxplot(data=np.abs(df),linewidth=2)
    formatplot(ax,False,'Parameter values', xlim=False,ylim=False,logy=True)
    
    plt.savefig(DIR_PLOTS+ZFname+'box.pdf',dpi=my_dpi,transparent=True)


def boxplots2(nwalkers,iterations,tburn,dataset,DIR_PLOTS,DIR_OUT):
    if dataset=='GCN':
        ZFnames=['GCN_Hill_ADD','GCN_Hill_BCB','GCN_Hill_nc','GCN_Hill_c']
        filename='GCN'
    elif dataset=='PDZ':
        ZFnames=['PDZ_Hill_ADD','PDZ_Hill_BCB','PDZ_Hill_nc','PDZ_Hill_c']
        filename='PDZ'
        
    parameternames=['$n_{ADD}$','$n_{BCB}$','$n_{noncoop}$','$n_{coop}$']
    
    df=pd.DataFrame(columns=parameternames)
    k=0
    for ZFname in ZFnames:
        # Read data
        df1=pd.read_csv(DIR_OUT+'samplesout_'+ZFname+'.csv',delimiter=',')
        data=np.zeros(df1.shape[0]*(df1.shape[1]-1)).reshape(df1.shape[0],(df1.shape[1]-1))
        for i in range(0,int(df1.shape[1]-1)):
            data[:,i]=np.array(df1.iloc[:,i+1]) # Put dataframe into array. Dataframe has no. columns = no. parameters.

        # Burn-in time correction

        data2=np.zeros((df1.shape[0]-tburn*nwalkers)*(df1.shape[1]-1)).reshape((df1.shape[0]-(tburn*nwalkers)),(df1.shape[1]-1))
        for i in range(0,int(df1.shape[1]-1)):
            for j in range(1,nwalkers+1):
                data2[(iterations-tburn)*(j-1):(iterations-tburn)*(j),i]=np.array(df1.iloc[iterations*j-iterations+tburn:iterations*j,i+1])
        df[parameternames[k]]=data2[:,-1]
        k+=1
        
    # Plot data
    plt.close("all")

    my_dpi=150

    figure_options={'figsize':(8.27,5.83)} #figure size in inches. A4=11.7x8.3. A5=8.27,5.83
    font_options={'size':'28','family':'sans-serif','sans-serif':'Arial'}
    plt.rc('figure', **figure_options)
    plt.rc('font', **font_options)

    current_palette=sns.color_palette("deep", 4)
    plt.rc('axes',prop_cycle=(cycler('color',current_palette)))
    f, ax=plt.subplots()
    plt.subplots_adjust(left=0.25,bottom=0.2,right=0.95,top=0.95)   
    
    ax=sns.boxplot(data=np.abs(df),linewidth=2)
    formatplot(ax,False,'Parameter values', xlim=False,ylim=([0,4]),logy=False)
    
    plt.savefig(DIR_PLOTS+'Hill_n_box_'+filename+'.pdf',dpi=my_dpi,transparent=True)

