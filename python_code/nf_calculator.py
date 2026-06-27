#!/usr/bin/env python
# coding: utf-8

### PULSAR NULLING FRACTION CALCULATOR
### 

import numpy as np
import matplotlib.pylab as plt
import os
import argparse
from scipy import signal
import emcee
from corner import corner
import nf_calculator_functions as nff

### INITIALIZE THE PARSER TO ACCEPT USER INPUTS

parser = argparse.ArgumentParser()

parser.add_argument("-pulsar", required = True, help = "Name of pulsar")
parser.add_argument("-fb", required = True, type = int, help = "First bin of pulse window")
parser.add_argument("-lb", required = True, type = int, help = "Last bin of pulse window")
parser.add_argument("-outdir", required = True, help = "Directory for output plots and files")
parser.add_argument("-datadir", required = True, help = "Directory where the pulsar data is found")

## LOAD THE ARGUMENTS

args = parser.parse_args()

### PULSAR NAME

p_name = args.pulsar

output_dir = f'{args.outdir}{p_name}'

if not os.path.exists(args.outdir):
    os.makedirs(args.outdir)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

### DEFINING THE OUTPUT LISTS

mjd_list = []
nf_list = []
hs_nf_list = []
nf_upper_unc = []
nf_lower_unc = []
sn_proxy_list = []
max_nulls = []
all_nulls = []
all_non_nulls = []
n_prof_list = []

### NUMBER OF OBSERVATIONS

datadir = f'{args.datadir}{p_name}/'
file_tot = len(os.listdir(datadir))
print('There are {} observations in this data set.'.format(file_tot))

### DEFINING THE ON-PULSE WINDOW AND OFF PULSE WINDOW

on_window_begin = args.fb
on_window_end = args.lb

### CALCULATE NULLING FRACTION

file_num = 0

for file in os.listdir(datadir):
        
        null_count_obs = []
        non_null_count_obs = []
        in_null = 0
        null_count = 0
        non_null_count = 0
        file_num += 1
        obs_file = datadir+file

        print('----------------------------------------------------------------------------------------')
        print('----------------------------------------------------------------------------------------\n')
        print(f'Now working on file: {file}')
        print(f'File number: {file_num} of {file_tot} for pulsar {p_name}')

        ########################################################
        ################# DATA PREPARATION #####################
        ########################################################
                
        ### ADD MJDS TO FILE
        
        f = open(obs_file, 'r')
        mjd_list.append(int(float(f.readline().lstrip('#'))))
        print('MJD:',mjd_list[-1],'\n')
        f.close()

        ### LOAD DATA AND INFER NUMBER OF BINS FROM BIN COUNTER RESET
        rawdata = np.loadtxt(obs_file)
        bin_nums = rawdata[:, 0]
        resets = np.where((bin_nums[1:] == 1) & (bin_nums[:-1] > 1))[0]
        no_bins = int(resets[0] + 1) if len(resets) > 0 else len(bin_nums)
        print(f'Number of phase bins inferred from data: {no_bins}')

        off_window_begin = on_window_begin + int(no_bins/4.)
        off_window_end = on_window_end + int(no_bins/4.)
        windows = [on_window_begin, on_window_end, off_window_begin, off_window_end]

        ### RESHAPE INTO (n_pulses, n_bins)
        print('Data shape before reshape:',rawdata.shape)
        data = np.reshape(rawdata[:,1],(int(rawdata.shape[0]/no_bins),no_bins))
        print('Data shape after reshape:',data.shape)
            
        n_prof = float(data.shape[0])
        n_prof_list.append(n_prof)

        ### FIND MEAN PROFILE                                                                                                                                                                                         
        mean = np.zeros((data.shape[1]))
        for b in range(no_bins):
            mean[b] = np.mean(data[:,b])
        mean = mean/np.max(mean)

        ### MOVE PEAK TO BIN 128                                                                                                                                                             
        peak_bin = np.argmax(mean)
        data = np.roll(data,int(no_bins/4)-peak_bin,axis=1)
        mean = np.roll(mean,int(no_bins/4)-peak_bin,axis=0)
            
        ### BASELINE BASED ON OFF WINDOW

        all_off_bins = []
            
        for i in range(int(data.shape[0])):
            all_off_bins.append(data[i,off_window_begin:off_window_end])
            
        off_baseline = np.mean(all_off_bins)
        print('Subtracting a baseline of', off_baseline)

        data = data - off_baseline

        ### FIND MEAN PROFILE                                                                                                                                                                                         
        mean = np.zeros((data.shape[1]))
        for b in range(no_bins):
            mean[b] = np.mean(data[:,b])
        mean = mean/np.max(mean)

        ### PLOT WATERFALL AND THE MEAN PROFILE

        peak_window_begin = on_window_begin
        peak_window_end = on_window_end

        nff.waterfall_and_average(data,mean,windows);
        plt.savefig(f'{output_dir}/{file}_waterfall_mean.png')
        plt.close()

        ### TAKE PEAK/(STD IN THE OFF-PULSE WINDOW) FOR EACH SINGLE PULSE
        ### TAKE THE LARGEST 5% AND ASSUME THAT THESE ARE PULSES (RATHER THAN NULLS)
        ### THIS IS NOT VALID IF THE NULLING FRACTION OF THE OBSERVATION IS > 95%
        
        find_sn_threshold = []
        for i in range(data.shape[0]):
            find_sn_threshold.append((np.max(data[i,peak_window_begin:peak_window_end])/np.std(data[i,off_window_begin:off_window_end])))
        sn_threshold = np.percentile(find_sn_threshold,95)

        ### NOW TAKE THESE SINGLE PULSES THAT WE HAVE DEEMED TO BE REAL PULSES (RATHER THAN NULLS)
        ### AND CALCULATE THEIR MEDIAN VALUE. THIS WILL SERVE AS A PROXY FOR THE S/N OF THE OBSERVATION
        ### THIS VALUE IS NOT USED IN THE CALCULATION OF NF BUT IS JUST A DIAGNOSTIC AND PLOTTED ON THE FINAL NULLING FRACTION PLOT
          
        temp_sn_list = []
            
        for i in range(data.shape[0]):
            if np.max(data[i,peak_window_begin:peak_window_end])/np.std(data[i,off_window_begin:off_window_end]) > sn_threshold:
                temp_sn_list.append(np.max(data[i,peak_window_begin:peak_window_end])/np.std(data[i,off_window_begin:off_window_end]))
        print('( Number of profiles that made it into the s/n calc:',len(temp_sn_list),'out of',int(n_prof),')')
        sn_proxy_list.append(np.median(temp_sn_list))
            
        ### FIND THE FLUX IN ALL THE ON WINDOWS ND NORMALISE BY THE MEAN FLUX DENSITY

        flux_density = []
        for i in range(data.shape[0]):
            flux_density.append(np.sum(data[i,on_window_begin:on_window_end]))

        norm_flux_density = [j/1.0 for j in flux_density]

        ### FIND THE FLUX IN ALL THE OFF WINDOWS AND NORMALISE BY THE MEAN FLUX DENSITY.                                                                                                                                                                       
        flux_density_off = []
        std_flux_density_off = []

        for i in range(int(data.shape[0])):
            flux_density_off.append(np.sum(data[i,off_window_begin:off_window_end]))
            std_flux_density_off.append(np.std(data[i,off_window_begin:off_window_end]))

        # NORMALISE

        norm_flux_density_off = [j/1.0 for j in flux_density_off]
        std_norm_flux_density_off = [j/1.0 for j in std_flux_density_off]

        ### FIND THE STATS TO EITHER FIX OR HAVE AS STARTING POINT (FOR THE MCMC SAMPLING) ###
        new_data_on_and_off_zipped = list(sorted(zip(norm_flux_density,norm_flux_density_off,std_norm_flux_density_off)))
        new_data_on_and_off_unzipped = list(zip(*new_data_on_and_off_zipped))
        new_data = np.array(new_data_on_and_off_unzipped[0])
        new_data_off = np.array(new_data_on_and_off_unzipped[1])
        std_new_data_off = np.array(new_data_on_and_off_unzipped[2])

        mean_of_nulls = np.mean(new_data_off)
        print(f'mean of nulls: {mean_of_nulls}')
        std_of_nulls = np.std(new_data_off)

        null_bar = mean_of_nulls + std_of_nulls
        pulse_bar = mean_of_nulls + (2.0*std_of_nulls)
        print(f'The null bar is set at {null_bar}')
        print(f'The pulse bar is set at {pulse_bar}')
                       
            
        ### PLOT FLUX DENSITY OF EACH PULSE WITH NULL/PULSE THRESHOLD LINES

        fig = plt.figure(figsize=(20,6))
        resamp = signal.resample(norm_flux_density,num=(int(n_prof/1)))
        plt.plot(resamp,'kx')
        plt.plot(resamp,'k-',alpha=0.1)
        plt.hlines(null_bar,0,n_prof,linestyle='dashed',linewidth=4,color='C0')
        plt.hlines(pulse_bar,0,n_prof,linestyle='dashed',linewidth=4,color='C1')
        plt.grid()
        fig.suptitle('Flux Density of Each Pulse', fontsize=30)
        plt.savefig(f'{output_dir}/{file}_flux.png')
        plt.close()
            
        ### FIND THE MAXIMUM NUMBER OF NULLS IN A ROW TO FACILITATE THE CALCULATION OF BINOMIAL UNCERTAINTY
            
        for p in range(len(resamp)):
            if resamp[p] < null_bar and in_null == 0 and p != len(resamp)-1:
                null_count += 1
                in_null = 1
                if p != 0:
                    non_null_count_obs.append(non_null_count)
                    non_null_count = 0
            elif resamp[p] < null_bar and in_null == 1 and p != len(resamp)-1:
                null_count += 1
            elif resamp[p] > null_bar and in_null == 0 and p != len(resamp)-1:
                non_null_count += 1
            elif resamp[p] > null_bar and in_null == 1 and p != len(resamp)-1:
                non_null_count += 1
                null_count_obs.append(null_count)
                in_null = 0
                null_count = 0
            elif p == len(resamp)-1: # at the end of the pulses it should append the number of nulls too
                if resamp[p] < null_bar and in_null == 1:
                    null_count+=1
                    null_count_obs.append(null_count)
                elif resamp[p] > null_bar and in_null == 0:
                    non_null_count+=1
                    non_null_count_obs.append(non_null_count)
                elif resamp[p] < null_bar and in_null == 0:
                    null_count_obs.append(1)
                    non_null_count_obs.append(non_null_count)
                elif resamp[p] > null_bar and in_null == 1:
                    non_null_count_obs.append(1)
                    null_count_obs.append(null_count)
                
            else:
                pass

        #max_nulls.append(np.max(null_count_obs))
        max_nulls.append(np.max(null_count_obs) if null_count_obs else 0)
        all_nulls.append(null_count_obs) # trying to say all the null trains. Maybe an average of all will give a better idea of trials?
        all_non_nulls.append(non_null_count_obs) # trying to say all the null trains. Maybe an average of all will give a better idea of trials?
        print(f'Max nulls: {max_nulls[-1]}')
        print('\n\n\n')
            
        ########################################################
        ############## HISTOGRAM SCALING METHOD ################
        ########################################################

        ### SET UP THE HISTOGRAM
        
        num_bins = 100
        
        upper = np.max((new_data[-1],np.max(new_data_off))) # highest flux density value
        lower = np.min((new_data[0],np.min(new_data_off))) # lowest flux density value
        
        bin_master_0 = (np.linspace(lower,upper,num_bins))
        
        bin_master_1 = np.append(bin_master_0,bin_master_0[1]-bin_master_0[0]+bin_master_0[-1])
        
        first_positive_arg = np.argmax(bin_master_1>0)
        
        bin_master = bin_master_1 - (bin_master_1[first_positive_arg]-0.0)

        ### PLOT THE HISTOGRAM BEFORE SCALING
        
        nff.histogram(norm_flux_density,norm_flux_density_off,bin_master)
        plt.savefig(f'{output_dir}/{file}_hist_before.png')
        plt.close()
        
        ### GET THE NEGATIVE VALUES
        
        neg_on = [x for x in norm_flux_density if x <= 0]
        neg_off = [x for x in norm_flux_density_off if x <= 0]
        
        ### PLOT HISTOGRAM OF NEGATIVE FLUX ONLY
        
        _, hy_on, hy_off = nff.histogram_neg(neg_on,neg_off,bin_master)
        plt.savefig(f'{output_dir}/{file}_hist_neg.png')
        plt.close()
        
        ### CALCULATE SCALING
        
        best_scaling = nff.histogram_scaling(hy_on,hy_off)

        ### PRODUCE A WARNING IF THERE WERE NO NEGATIVE ON-PULSE_WINDOW VALUES AND SET THE NULLING FRACTION TO 0.0

        if isinstance(best_scaling, str):
            print(best_scaling)
            hs_nf = 0.0
        else:
            hs_nf = 1./best_scaling

        hs_nf_list.append(hs_nf)
            
        ### PRODUCE A POST-SCALING HISTOGRAM (UNLESS THERE ARE NO NEGATIVE ON-PULSE-WINDOW VALUES)
        
        if isinstance(best_scaling, str):
            print('There will be no scaled plot because the nulling fraction is 0.0')
        else:
            
            while hy_off[-1] == 0 and hy_on[-1] == 0:
                hy_off = hy_off[:-1]
                hy_on = hy_on[:-1]
                
            hs_x = np.arange(0,hy_off.shape[0]+1)
            hy_off = np.append(hy_off,0)
            hy_on = np.append(hy_on,0)

            # Show scaled plot
            nff.scaled_hist(hs_x, hy_off, hy_on*best_scaling, best_scaling)
            plt.savefig(f'{output_dir}/{file}_hist_scaled.png')
            plt.close()
            
        print('\n')
            
        #########################################################
        ################# BAYESIAN METHOD #######################
        #########################################################

        ### STARTING POINT FOR THE SAMPLER IS TAKEN FROM THE HS NF RESULTS.
        ### BUT WE DON'T WANT TO START RIGHT AT 0.0 OR 1.0

        if hs_nf == 1.0:
            start_from_hs_nf = 0.9
        elif hs_nf == 0.0:
            start_from_hs_nf = 0.1
        else:
            start_from_hs_nf = hs_nf

        Nexp = start_from_hs_nf * n_prof

        print(f'Number of nulls from HS results: {Nexp}')
        print('\n')
        
        ### ESTIMATE OF THE MEAN AND STD OF THE ON-PULSE DATA
        
        mean_of_on = np.mean(new_data[int(Nexp):])
        std_of_on = np.std(new_data[int(Nexp):])

        ### WE CAN ESTIMATE THE MU AND SIGMA PARAMETERS FOR THE LOGNORMAL DISTRIBUTION FROM THE MEAN AND STD ESTIMATES WE JUST CALCULATED
        
        estimated_sigma = np.sqrt(np.log((np.square(std_of_on) / np.square(mean_of_on)) + 1))
        estimated_mu = np.log(mean_of_on) - (np.square(estimated_sigma) / 2)

        ### STARTING POINTS FOR THE SAMPLER

        pTrue = [Nexp, estimated_mu, estimated_sigma]
        
        ### SAMPLER SET UP
        
        ndim = 3
        nwalkers = 10
        n_burn_steps = 100
        n_sample_steps = 1000
            
        p0 = [pTrue + np.random.randn(ndim)*[pTrue[0]/1000., pTrue[1]/20., pTrue[2]/20.] for i in range(nwalkers)]

        sampler = emcee.EnsembleSampler(nwalkers, ndim, nff.logL_lognormal_conv, args=[new_data, mean_of_nulls, std_of_nulls])
        
        ### BURN IN
        
        pos, prob, state = sampler.run_mcmc(p0, n_burn_steps)  # burn in
        sampler.reset()

        ### SAMPLE
        
        pox, prob, state = sampler.run_mcmc(pos, n_sample_steps, rstate0=state)  # sample

        ### PLOT THE POSTERIORS
        
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10
        plt.rcParams['axes.labelsize'] = 14
        corner(sampler.flatchain, 100, title='' ,labels=["Number of Nulls", "$\\mu$\nOn-Pulses", "$\\sigma$\nOn-Pulses"])
        plt.savefig(f'{output_dir}/{file}_corner_plot.png')            
        plt.close()

        ### CALCULATE THE MEDIAN VALUES OF THE POSTERIORS
        
        med_of_no_nulls = np.percentile(sampler.flatchain[:,0], 50)
        med_of_on_mean = np.percentile(sampler.flatchain[:,1], 50)
        med_of_on_std = np.percentile(sampler.flatchain[:,2], 50)
        
        print(f'Median of number of nulls from Bayesian analysis: {med_of_no_nulls}')
        print(f'Median of mu value from Bayesian analysis: {med_of_on_mean}')
        print(f'Median of sigma value from Bayesian analysis: {med_of_on_std}')
        
        ### CALCULATE THE 1 SIGMA ERROR BARS
        
        low, med, up = np.percentile(sampler.flatchain[:,0], [15.865,50,84.135])
        bayes_nf = med/n_prof
        
        ### ADD CALCULATED NF AND UNCERTAINTY TO A LIST

        nf_list.append(bayes_nf)
        nf_upper_unc.append((up-med)/n_prof)
        nf_lower_unc.append((med-low)/n_prof)
        
        ### PLOT HISTOGRAM OF THE FLUX DENSITY DATA WITH THE BEST FITTING FUNCTION

        total_function_x, total_function_y, null_func, on_func = nff.convolve([med_of_no_nulls, med_of_on_mean, med_of_on_std], new_data, mean_of_nulls, std_of_nulls)
        
        if total_function_x.shape[0] == total_function_y.shape[0]:
            pass
        else:
            total_function_x = total_function_x[:-1]
            
        nff.bayes_hist(new_data,total_function_x,total_function_y,on_func,null_func)

        plt.savefig(f'{output_dir}/{file}_bayes_fit.png')
        plt.close()
            
        print('\n\n')
            
        print("                                   +{:.5f}".format(np.sqrt(hs_nf_list[-1]*n_prof_list[-1]) / n_prof_list[-1]))
        print('The HS nulling fraction: {:.5f}'.format(hs_nf))
        print("                                   +{:.5f}".format(np.sqrt(hs_nf_list[-1]*n_prof_list[-1]) / n_prof_list[-1]))
        print('\n')
        print("                                   +{:.5f}".format((up-med)/n_prof))
        print("Bayesian nulling fraction: {:.5f}".format(med/n_prof))
        print("                                   -{:.5f}\n".format((med-low)/n_prof))

        print('\n')
        print('Above uncertainties do not include binomial uncertainy from finite observation length')
        print('\n\n')

### ZIP LISTS OF RESULTS TOGETHER AND ORDER THEM CHRONOLOGICALLY:

zipped = list(sorted(zip(mjd_list, nf_list, nf_lower_unc, nf_upper_unc, sn_proxy_list, n_prof_list, hs_nf_list, max_nulls)))

if len(mjd_list) != len(nf_list) or len(mjd_list) != len(nf_lower_unc) or len(mjd_list) != len(nf_upper_unc) or len(mjd_list) != len(sn_proxy_list) or len(mjd_list) != len(n_prof_list) or len(mjd_list) != len(hs_nf_list) or len(mjd_list) != len(max_nulls) or len(mjd_list) != len(all_nulls) or len(mjd_list) != len(all_non_nulls):
    raise RuntimeError('Output lists are not all of equal length — one or more observations may have failed to process.')

print(f'All {file_tot} observations processed successfully.')

### SAVE RESULTS

unzipped_to_save = list(zip(*zipped))
np.savetxt(f'{output_dir}/{p_name}.txt',unzipped_to_save)

# Save the consecutive null lists to a text file
with open(f'{output_dir}/{p_name}_consecutive_nulls_and_non.txt', 'w') as f:
    f.write(f"{repr(all_nulls)}\n")
    f.write(f"{repr(all_non_nulls)}\n")
