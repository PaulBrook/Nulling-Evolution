import numpy as np
from numpy import log as ln
import matplotlib
import matplotlib.pylab as plt
from matplotlib.patches import Rectangle
import sys

isq2pi = 1/(2*np.pi)**(0.5)


def histogram(on_data,off_data,bin_master):
    fig = plt.figure(figsize=(20,6))
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    
    hy,hx,_=plt.hist(on_data,bins=bin_master,alpha=.3,label='on-pulse data',color='green')
    hy,hx,_=plt.hist(off_data,bins=bin_master,alpha=.3,label='off-pulse data',color='darkorange')

    fig.text(0.07, 0.5, 'Count', ha='center', va='center', rotation='vertical', size=22)
    fig.text(0.504, 0.02, 'Flux Density in On-Pulse Window', ha='center', va='center', rotation='horizontal', size=22)
    plt.legend(prop={'size': 20},loc='upper left')
    plt.grid()

    #fig.suptitle('Histogram of Flux Density in On- and Off-Windows', fontsize=30)

    return fig

def histogram_neg(on_data,off_data,bin_master):
    fig = plt.figure(figsize=(20,6))
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    
    hy_on,hx_on,_=plt.hist(on_data,bins=bin_master,alpha=.3,label='on-pulse data',color='green')
    hy_off,hx_off,_=plt.hist(off_data,bins=bin_master,alpha=.3,label='off-pulse data',color='darkorange')

    fig.text(0.07, 0.5, 'Count', ha='center', va='center', rotation='vertical', size=22)
    fig.text(0.504, 0.02, 'Flux Density in Window', ha='center', va='center', rotation='horizontal', size=22)
    plt.legend(prop={'size': 20},loc='upper left')
    plt.xlim(right=0)
    plt.grid()

    fig.suptitle('Histogram with Negative Values of Flux Density Only', fontsize=30)    

    return fig, hy_on, hy_off, 
    
def scaled_hist(x, y_off, y_on, best_scaling):
    fig, ax = plt.subplots(figsize =(20, 6))
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    plt.fill_between(x, y_off, step="post", alpha=0.3, color='darkorange',label='off-pulse data')
    plt.fill_between(x, y_on, step="post", alpha=0.3, color='green',label='on-pulse data')
    plt.grid()
    fig.text(0.07, 0.5, 'Count', ha='center', va='center', rotation='vertical', size=22)
    fig.text(0.504, 0.02, 'Bins', ha='center', va='center', rotation='horizontal', size=22)
    plt.legend(prop={'size': 20},loc='upper left')

    plt.ylim(bottom=0)

    fig.suptitle('Histogram After On-Window Data Has Been Scaled By {:.3f}'.format(best_scaling), fontsize=30)    

def bayes_hist(hist_data, x_func, tot_func, on_func, null_func):
    fig = plt.figure(figsize=(20,6))
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    
    hy,hx,_=plt.hist(hist_data,50,alpha=.3,label='data',color='gray')
    peak_hist = np.max(hy)
    peak_function = np.max(tot_func)
    peak_ratio = peak_function/peak_hist
    plt.plot(x_func, tot_func/peak_ratio,'k--',linewidth=3)
    plt.plot(x_func, null_func/peak_ratio,'b-',linewidth=1)
    plt.plot(x_func, on_func/peak_ratio,'r-',linewidth=1)
    hx=(hx[1:]+hx[:-1])/2 # for len(x)==len(y)
    
    plt.grid()
    fig.text(0.07, 0.5, 'Count', ha='center', va='center', rotation='vertical', size=22)
    fig.text(0.504, 0.02, 'Flux Density in Window', ha='center', va='center', rotation='horizontal', size=22)
    fig.suptitle('Histogram of On-Window Flux Density with Recovered Bayesian Parameters', fontsize=30)    

    return fig

def waterfall_and_average(data,mean,windows):
    fig = plt.figure(figsize=(20,8))
    plt.rcParams['xtick.labelsize'] = 16
    plt.rcParams['ytick.labelsize'] = 16
    ax1 = plt.subplot2grid((3,1),(0,0),rowspan = 2)                                                                                                                                                          
    ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False, top=False, labeltop=False)                                                                                                                                              
    ax1.imshow(data[:,:], aspect='auto',cmap='binary',origin='lower')
    fig.text(0.07, 0.645, 'Pulse Number', ha='center', va='center', rotation='vertical', size=18)
    fig.text(0.07, 0.25, 'Normalized\nFlux Density', ha='center', va='center', rotation='vertical', size=18)
    fig.text(0.504, 0.03, 'Phase Bin', ha='center', va='center', rotation='horizontal', size=18)
    ax2 = plt.subplot2grid((3,1),(2,0))
    ax2.plot(mean,color='k')
    ax2.set_xlim(0,data.shape[1])
    ax2.set_ylim(-0.1,1.1)
    ax2.grid()
    plt.vlines(windows[0],-0.1,1.1,linestyles='dashed',lw=4,color='green')
    plt.vlines(windows[1],-0.1,1.1,linestyles='dashed',lw=4,color='green')
    ax2.add_patch(Rectangle((windows[0],-0.1),(windows[1]-windows[0]),1.2,facecolor='green',alpha = 0.1, zorder = 39))
    plt.vlines(windows[2],-0.1,1.1,linestyles='dashed',lw=4,color='darkorange')
    plt.vlines(windows[3],-0.1,1.1,linestyles='dashed',lw=4,color='darkorange')
    ax2.add_patch(Rectangle((windows[2],-0.1),(windows[3]-windows[2]),1.2,facecolor='darkorange',alpha = 0.1, zorder = 39))                                                                    

    plt.vlines(1,-0.1,1.1,linestyles='dashed',lw=4,color='darkorange')
    plt.vlines(windows[3]-windows[2],-0.1,1.1,linestyles='dashed',lw=4,color='darkorange')
    ax2.add_patch(Rectangle((0,-0.1),(windows[3]-windows[2]),1.2,facecolor='darkorange',alpha = 0.1, zorder = 39))


    fig.subplots_adjust(hspace=0.0)
    fig.subplots_adjust(wspace=0.0)

    #fig.suptitle('All Profiles and Average Profile', fontsize=30)
    return fig,ax2;
    
def logL_lognormal_conv(params, data, mean_off_window_nulls, std_off_window_nulls):
    
    """likelihood"""
    
    # do the convolution

    #print(f'DATA!{data}')
    #print(f'PARAMS at top: {params}')
    
    total_function_x, total_function_y, _, _ = convolve(params, data, mean_off_window_nulls, std_off_window_nulls)

    #if np.any(np.isnan(total_function_x)) or np.any(np.isnan(total_function_y)):
    #    print("!!!!!!!!!!!!!!NaN detected in total_function_x or total_function_y")
    #    sys.exit()
    
    #print(f'TOTAL FUNCTION X{total_function_x}')
    
    #print(f'size of total func x and y {total_function_x.shape} {total_function_y.shape}')
    #print(f'y array {total_function_y}')
    
    N = len(data)    
    
    # now put the data values into the function
    
    likelihood = []

    is_nan = np.isnan(total_function_y)
    #print(f'IS NAN??? {np.any(is_nan)}')
    
    #print(f'PARAMS further down: {params}')

    #if np.any(np.isnan(total_function_x)) or np.any(np.isnan(total_function_y)):
    #    print("!!NaN detected in total_function_x or total_function_y")
    #    sys.exit()

    
    if params[0] > (N*1.0) or params[0] < 0 or params[1] < 0 or params[1] > 10.0 or params[2] < 0 or params[2] > 2.0 or np.any(is_nan):
    # I used (N*0.9) because the sampler sometimes tends to simplify and make NF=1 when it's clearly not. I'm keeping it just as N for now and see how it goes.  
        sum_out = -1e50
    
    else:

        ###function_y_is_neg = 0        
        
        for each in data:
            #print(f'eACCCCCCHHHH {each}')
            difference_array = np.absolute(total_function_x-each)
            #print(f'diffference array {difference_array}')
            index = difference_array.argmin()
            #print(f'INDEXXXXX {index}')
            ###if total_function_y[index] > 0:
            likelihood.append(total_function_y[index])
            ###else:
                # I added this because if the total_function_y is not >0 then I think the LogL and the sum_out will yeild NaN
                # I should probably stop the total_function_y from being negative in the first place.
            ###    function_y_is_neg = 1
            ###    likelihood.append(total_function_y[index])

        ###if function_y_is_neg == 0:

        #print(f'likelihood!!!! {likelihood}')
        
        logLs = np.log(likelihood)

        #print(f'logLs!!!! {logLs}')
        
        sum_out = np.sum(logLs)

        ###elif function_y_is_neg == 1:

        ###    sum_out = -1e50

        #print(f'SUM OUT!!!!!: {sum_out}')
        
    return sum_out

def convolve(params, data,  mean_off_window_nulls, std_off_window_nulls):
    """
    Convolves the lognormal function with the noise in the data.
    
    params = 
    
    data = 
    
    mean_off_windows_nulls = 
    
    std_off_window_nulls = 
    """
        
    Noff, u, sig = params
    
    N = len(data)    

    # the resolution of the of the functions
    
    function_resolution = 0.1
    
    # the range of the functions
    
    if abs(np.min(data)) > abs(np.max(data)):
        set_range = abs(np.min(data))
    else:
        set_range = abs(np.max(data))

    #print(f'set range is: {set_range}')    
        
    #function_resolution = (set_range*2.0)/100.0
    #function_resolution = 0.1
    
    # the off on and conv x-range
    
    off_func_x = np.arange(-set_range,set_range,function_resolution)
    null_func_x = np.arange(-set_range,set_range,function_resolution)
    on_func_x = np.arange(function_resolution,set_range,function_resolution)

    conv_func_x = np.arange(-set_range,(2.0*set_range)-function_resolution,function_resolution)
    
    # the y values after the function has acted
    
    
    func_arg_off = -0.5 * ((off_func_x-mean_off_window_nulls)/std_off_window_nulls)**2
    off_func_y_for_conv = 1*np.exp(func_arg_off)
    #off_func_y_for_conv = np.ones((func_arg_off.shape[0]))
    
    func_norm_null = Noff/std_off_window_nulls
    func_arg_null = -0.5 * ((null_func_x-mean_off_window_nulls)/std_off_window_nulls)**2
    #func_arg_null = -0.5 * ((null_func_x-uoff)/std_off_window_nulls)**2
    func_norm_on = (1/on_func_x) * (N-Noff)/sig
    func_arg_on = -0.5 * ((ln(on_func_x)-u)/sig)**2
    
    
    null_func_y = func_norm_null*np.exp(func_arg_null)
    on_func_y_for_conv = func_norm_on*np.exp(func_arg_on)
    
    # do the convolution

    conv_func_y = np.convolve(on_func_y_for_conv,off_func_y_for_conv,mode='full')/(np.sum(on_func_y_for_conv))
    
    area_on = np.sum(on_func_y_for_conv)
    area_conv = np.sum(conv_func_y)
    
    conv_func_y *= (area_on/area_conv)
    
    # pad the on function now so that it can be plotted nicely with off function
    
    on_func_y_expanded = np.zeros((off_func_y_for_conv.shape))
    
    if on_func_y_expanded.shape[0] - (2*on_func_y_for_conv.shape[0]) == 2: # to stop it sometimes not having same length as other arrays
        on_func_y_expanded[on_func_y_for_conv.shape[0]+2:] = on_func_y_for_conv
    elif on_func_y_expanded.shape[0] - (2*on_func_y_for_conv.shape[0]) == 1:
        on_func_y_expanded[on_func_y_for_conv.shape[0]+1:] = on_func_y_for_conv
    else:
        sys.exit()

    # now null function must be expaned to be added to the conv function
    
    null_func_y_expanded = np.zeros((conv_func_y.shape))
    null_func_y_expanded[:null_func_y.shape[0]] = null_func_y
    
    # add the conv and the off:
    
    total_function_y = null_func_y_expanded + conv_func_y
    total_function_x = conv_func_x

    return total_function_x, total_function_y, null_func_y_expanded, conv_func_y

def wang(hist_list_on,hist_list_off):
    scaling = 1.0
    if len(hist_list_on) == len(hist_list_off):
        pass
    else:
        print('Different number of bins!')
    smallest_diff = 1e9 # just to make the while loop work
    if np.mean(hist_list_on) == 0.0:
        best_scaling = 'There are no on-pulse windows with negative flex density value. NF will be assumed to be 0.0'
    else:
        while np.mean(hist_list_on)*scaling/np.mean(hist_list_off) < 1.5 :
            diff = 0.0
            for i in range(len(hist_list_on)):
                diff += hist_list_off[i]-hist_list_on[i]*scaling
            if abs(diff) < abs(smallest_diff):
                smallest_diff = diff
                best_scaling = scaling
            scaling += 0.001
            
    return best_scaling
