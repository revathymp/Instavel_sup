######## GNSS instavel Vs strong-motion derived velocity analysis and plots ##################
####################################### Goals ################################################
#     1) Read a GNSS instavel and convert it into an obspy stream
#     2) Filter instavel
#     3) Load acceleration from IRIS client
#     4) Filter and integrate strong-motion acceleration to corresponding velocity
#     5) Loop through different filters
#     6) write out CC and lag values
#     7) comparison plots for GNSS -seismic traces for different filters
#     8) CC lag plot
##############################################################################################

from obspy import UTCDateTime 				
from obspy import read, Trace
from obspy import Stream
from obspy import read_inventory
from obspy.clients.fdsn import Client 		
from obspy.signal.cross_correlation import correlate
import numpy as np
#import pandas as pd
#import noise
import matplotlib.pyplot as plt
client=Client("IRIS")

def compute_shift(cc):
    zero_index = int((len(cc)-1)/2)
    shift = zero_index - np.argmax(cc)
    return shift

# READ INSTAVEL ASCII FILE
infile = 'Chignik_phase_vels/ac34_1hz.L2.igr.phase_vels'
site = np.genfromtxt(infile, dtype=None, encoding=None, names=['date','time','dt','lat','lon','height','vx', 'vy', 'vz', 'sigx', 'sigy', 'sigz', 'vdt', 'vn', 've', 'vu', 'sign', 'sige', 'sigu'])

# CREATE AN EMPTY STREAM AND POPULATE WITH INSTAVEL DATA
gnss = Trace()
gnss.stats.network = '--'
gnss.stats.station = 'AC34'
gnss.stats.channel = 'IVE'       # for [I]nst[V]el [E] component]
gnss.stats.starttime = site["date"][0] + 'T' + site["time"][0]
gnss.stats.sampling_rate = 1
gnss.stats.calib = 1
gnss.data = site['ve']
gnss.detrend('linear')
gnss.detrend('demean')
dt = UTCDateTime("2021-07-29T06:15:49.000000Z")
gnss.trim(dt-40, dt+360)


# LOAD ONE SEISMIC ACCELERATION STREAM FROM IRIS DMC
start = UTCDateTime("2021-07-29T06:15:08.998400Z")
seis_a = client.get_waveforms('AK', 'S19K', "*", 'HNE', (start), (start+400), attach_response=True)

# TIDY UP ACCELERATION TRACE
seis_a[0].remove_response(inventory=None, output="ACC", plot=False)
seis_a[0].detrend('linear')

# TIDY UP GNSS VELOCITY TRACE -- FILTER THE HECK OUTTA THAT THING!!!!!!!!!!!
gnss.filter('bandpass', freqmin=0.001, freqmax= 0.45, corners=4, zerophase=True)
seis_a[0].filter('bandpass', freqmin=0.001, freqmax= 0.45, corners=4, zerophase=True)

    #CREATE A SEISMIC VELOCITY TRACE
seis_v = seis_a[0].copy()
seis_v.integrate(method='cumtrapz')
hz = seis_v.data[::100]   # Brute-force 'downsampling' of strong-motion derived velocity to 1 Hz post-filtering and -integration

    # CREATE AN EMPTY STREAM AND POPULATE WITH 1Hz STRONG-MOTION DERIVED VELOCITY DATA
sm_v = Trace()
sm_v.stats.network = 'AK'
sm_v.stats.station = 'S19K'
sm_v.stats.channel = 'HNE'    
sm_v.stats.starttime = UTCDateTime("2021-07-29T06:15:09.0Z")
sm_v.stats.sampling_rate = 1
sm_v.stats.calib = 1
sm_v.data = hz

#G.append(gnss)
#S.append(sm_v)

# NORMALIZE EACH TRACE WITH ITS ABSOLUTE MAXIMUM
    #st_seis_vz.normalize()
    #st_sm_vz[0].normalize()
    #st_gnss_vz[0].normalize()

# CROSS CORRELATION
lags = np.arange(-400, 401)  
cc = correlate(gnss.data, sm_v.data, 400)
corr = max(cc)
shift = compute_shift(cc)
print(corr, shift)
#print(G)
#print(S)    

# CREATE STREAMS WITH ALL TRACES AND PLOT THEM
plt.rcParams['figure.figsize'] = [12,8]
plt.rcParams.update({'font.size': 12})
#plt.style.use('seaborn')


#fig, axx = plt.subplots(1, 1, sharex=True, sharey=False)
fig,axx = plt.subplots()
axx.plot(gnss.times(), gnss.data, 'k-', linewidth=1)
axx2=axx.twinx()
axx2.plot(sm_v.times(), sm_v.data, 'r-', linewidth=1)
#axx[2, 0].plot(st_seis_vz.times(), st_seis_vz.data, 'k-', linewidth=1, label=st_seis_vz.stats['station'])
axx.set_title(f'GNSS instantaneous velocity (1Hz)')
#axx[1, 0].set_title(f'Strong-motion velocity (1Hz)')
#axx[2, 0].set_title(f'Strong-motion velocity (100Hz)')
#axx[0, 0].set_ylabel('Amplitude')
#axx[1, 0].set_ylabel('Amplitude')
#axx[2, 0].set_ylabel('Amplitude')

axx.set_xlabel('time (s)')
plt.tight_layout()

# TRACE LEGENDS
#for row in axx[:,0]:
#    ll = row.legend(loc=1)
#    plt.setp(ll.get_texts(), color='red') #color legend

#ll = axx[2,1].legend(loc=1)
#plt.setp(ll.get_texts(), color='red') #color legend

plt.show()
