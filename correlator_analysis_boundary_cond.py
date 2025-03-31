# Fernando Alvarado, March 2025

# Program utility: Application of the antiperiodic boundary condition in time to the correlator data.
# For the ensembles with antiperiodic boundary conditions, the quark fields must obey q(T)= - q(t1).
# Requirements: numpy and matplot libraries.
# In this program "@user" written where the user has to input information.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#############################################################################
# Importing correlator data
# @user:  names and path should be changed for the specific cases.
# In this case we have a nucleon at rest on N451: t00/1_f/bwd_row0/1 possibilities.
# If one has a matrix of coupled correlators < O_i \bar O_j >, one will have several elements {i,j} corresponding to individual correlators.
# If that is the case, just apply the change of sign procedure to all indices.

typestring=["t00_fwd_row0","t00_fwd_row1","t00_bwd_row0","t00_bwd_row1","t01_fwd_row0","t01_fwd_row1","t01_bwd_row0","t01_bwd_row1"];
tottypes=len(typestring);
corrfile = [None] * tottypes;
for i in range(tottypes):
    corrfile[i]="corr_"+typestring[i]+".dat";
print("The correlators are:")
print(corrfile)
print("The program is running...")
corr_list=[[None] * 2 for _ in range(tottypes)];
for i in range(tottypes):
    # just the correlator type name
    corr_list[i][0]=corrfile[i];
    # the rows correspond to time separations organised by config, column 0 = tsep, column 1 = correlator value:
    corr_list[i][1]=np.genfromtxt("chimera_extraction/Ops3_converted/"+ corrfile[i],skip_header=1,usecols=(0,1),delimiter=' ');

#############################################################################

#############################################################################
# @user, input the total number of studied configs, the maximum time separation, "maxtsep",  and the minimum one,"mintsep".
totconfigs=101; maxtsep=20; mintsep=2;
# derived quantities, just for the looping:
tottseps = 1+maxtsep-mintsep;# total number of tseps
tsepsperiod=tottseps+mintsep;# total number of tseps plus the skipped ones

# Removing the first two tseps, which are zero in the last_laph file.
corr_list_no_zero=[[None] * 2 for _ in range(tottypes)];
tsep_by_cnfg=[None]*tottypes; # to store tsep by config

for i in range(tottypes):
    corr_list_no_zero[i][0]=corrfile[i];
    corr_list_no_zero[i][1]=np.zeros([tottseps*totconfigs,2]);# Nrows=tsep*config, Ncolumns=2
    tsep_by_cnfg[i]=np.zeros([tottseps, totconfigs]);# matrix tsep*config
    for tsepi in range(tottseps):
        for configi in range(totconfigs):
            tsep_by_cnfg[i][tsepi,configi]=corr_list[i][1][tsepi + 2 + tsepsperiod*configi,0];
            corr_list_no_zero[i][1][tsepi+tottseps*configi, :]=corr_list[i][1][tsepi + 2 + tsepsperiod*configi,:];
#############################################################################


#############################################################################
# separating by config:

corr_list_by_cnfg=[[None] * 2 for _ in range(tottypes)];# list of different corr types, to be separated by config
tsep_by_cnfg_no_zero=[None]*tottypes;# to store tsep by config


for i in range(tottypes):
    corr_list_by_cnfg[i][0]=corrfile[i];
    corr_list_by_cnfg[i][1]=np.zeros([tottseps, totconfigs]);# to separate by config: it's a ( tsep x config ) matrix
    tsep_by_cnfg_no_zero[i]=np.zeros([tottseps, totconfigs]);

    for tsepi in range(tottseps):
        for configi in range(totconfigs):
            tsep_by_cnfg_no_zero[i][tsepi,configi]=corr_list_no_zero[i][1][tsepi+tottseps*configi,0];
            corr_list_by_cnfg[i][1][tsepi,configi]=corr_list_no_zero[i][1][tsepi+tottseps*configi,1];
#############################################################################

#############################################################################
# Importing the time sources
# @user: edit file and path.
# in this case we import the N451 source times, configs=1,2,...1011 and t=0,1,...,127 (notice the c-counting for t).
# It contains the config and source times: [cnfg_abs,tsrc00,tsrc01]
cnfg_abs_tsrc_list=np.genfromtxt("/home/falvar/Start/Projects/N451_correlators/N451_last_laph/reading_time_sources/source_times_N451.dat",dtype=int,usecols=(0,1,2),delimiter=' ');

# @user: Insert the absolute first config and the config interval:
# in this example,
cnfg_abs_start=10; cnfg_abs_interval=10; #the total number of cnfgs already inserted, totconfigs (=101 in this example).
# cnfg_abs= config_start+(configi)*config_interval


# @user, insert the number of t slices:
#in this example:
tott=128;# last one is t=tott-1 (c-counting)
#############################################################################

#############################################################################
# Transforming to the relative configs , "configi" (configi is just an index to label the studied configs, running from 0 to "totconfigs"), so that one has the list [configi,tsrc00,tsrc01]
configi_tsrc_list=np.zeros([totconfigs,3]);
for configi in range(totconfigs):
    cnfg_abs= cnfg_abs_start+(configi)*cnfg_abs_interval; # cnfg_abs=1,2,...,totconfigs
    configi_tsrc_list[configi,0]=(cnfg_abs_tsrc_list[cnfg_abs-1,0]-cnfg_abs_start)/cnfg_abs_interval;# cnfg_abs-1 because of python
    configi_tsrc_list[configi,1:]=cnfg_abs_tsrc_list[cnfg_abs-1,1:]; # t_src
#############################################################################

#############################################################################
# WARNING: running twice might cause error.

# Here the boundary condition that changes the sign is applied.
corr_list_by_cnfg_sign=corr_list_by_cnfg; # list of different corr types, to be separated by config

for i in range(tottypes):
    if "t00" in corrfile[i]:
        tsrci=0;
    elif "t01" in corrfile[i]:
        tsrci=1;
    else:
        print("could not tell if it is a t00 or a t01 source");


    if "fwd" in corrfile[i]:
        for configi in range(totconfigs):
            for tsepi in range(tottseps):
                tsrc=configi_tsrc_list[configi,1+tsrci]; # source time
                tfwd=tsrc+mintsep+tsepi; # In the example case, before the boundary, 2,3,...,127 or if tserc=127, would start at 127+2=129
                if tfwd >= tott:
                    corr_list_by_cnfg_sign[i][1][tsepi:,configi]= - corr_list_by_cnfg[i][1][tsepi:,configi]; # changes sign
                    break;
        # drawing plots
        pdf_filename = "corr_"+typestring[i]+"_sign.pdf";
        with PdfPages(pdf_filename) as pdf:
            for configi in range(totconfigs):
                plt.figure()
                plt.scatter(tsep_by_cnfg_no_zero[i][:,configi],corr_list_by_cnfg_sign[i][1][:,configi], marker='.')
                plt.title(typestring[i]+"config"+str(configi))
                plt.grid()
                pdf.savefig()
                plt.close()
        print(f"Plots saved to {pdf_filename}")

    elif "bwd" in corrfile[i]:
        corr_list_by_cnfg_sign[i][1][:,:]= - corr_list_by_cnfg_sign[i][1][:,:]; # the bwd corr has a minus sign to be accounted for.
        for configi in range(totconfigs):
            for tsepi in range(tottseps):
                tsrc=configi_tsrc_list[configi,1+tsrci]; # source time
                tbwd=tsrc-mintsep-tsepi;
                if tbwd <= -1:
                    corr_list_by_cnfg_sign[i][1][tsepi:,configi]= - corr_list_by_cnfg[i][1][tsepi:,configi]; # changes sign
                    break;
        # drawing plots:
        pdf_filename = "corr_"+typestring[i]+"_sign.pdf";
        with PdfPages(pdf_filename) as pdf:
            for configi in range(totconfigs):
                plt.figure()
                plt.scatter(tsep_by_cnfg_no_zero[i][:,configi],corr_list_by_cnfg_sign[i][1][:,configi], marker='.')
                plt.title(typestring[i]+"config"+str(configi))
                plt.grid()
                pdf.savefig()
                plt.close()
        print(f"Plots saved to {pdf_filename}")

    else:
        print("could not tell if it is a fwd or a bwd correlator")
        break;

    # Might be good to delete tsrc, tfwd, tbwd, tsrci, configi.
#############################################################################

#############################################################################
# check that the abs values are exactly the same as before processing:
for i in  range(tottypes):
    if np.all(abs(corr_list_by_cnfg_sign[i][1][:,:])==abs(corr_list_by_cnfg[i][1][:,:])):
        print("ok");
    elif np.any(abs(corr_list_by_cnfg_sign[i][1][:,:])!=abs(corr_list_by_cnfg[i][1][:,:])):
        print("ERROR");
    else:
        print("ERROR number 2");

# @user: One can see in the figures whether the sign problem is gone.
#############################################################################
