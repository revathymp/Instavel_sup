######## GNSS instavel Vs strong-motion derived velocity analysis and plots ##################
####################################### Goals ################################################
#     1) Read a GNSS instavel and convert it into an obspy stream
#     2) Filter instavel
#     3) Load acceleration from IRIS client
#     4) Filter and integrate strong-motion acceleration to corresponding velocity
#     5) Normalize all streams and plot
#     6) Include spectrograms and cross-correlation plots
##############################################################################################

from obspy import UTCDateTime 				
from obspy import Trace
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
    zero_index = int(len(cc) / 2) - 1
    shift = zero_index - np.argmax(cc)
    return shift

# READ INSTAVEL ASCII FILE
infile = 'Chignik_phase_vels/ac52_1hz.L2.igr.phase_vels'
site = np.genfromtxt(infile, dtype=None, encoding=None, names=['date','time','dt','lat','lon','height','vx', 'vy', 'vz', 'sigx', 'sigy', 'sigz', 'vdt', 'vn', 've', 'vu', 'sign', 'sige', 'sigu'])

# CREATE AN EMPTY STREAM AND POPULATE WITH INSTAVEL DATA
st_gnss_vz = Stream(Trace())
st_gnss_vz[0].stats.network = '--'
st_gnss_vz[0].stats.station = 'AC52'
st_gnss_vz[0].stats.channel = 'IVE'       # for [I]nst[V]el [E] component]
st_gnss_vz[0].stats.starttime = site["date"][0] + 'T' + site["time"][0]
st_gnss_vz[0].stats.sampling_rate = 1
st_gnss_vz[0].stats.calib = 1
st_gnss_vz[0].data = site['ve']
st_gnss_vz[0].detrend('linear')
st_gnss_vz[0].detrend('demean')
dt = UTCDateTime("2021-07-29T06:15:49.000000Z")
st_gnss_vz[0].trim(dt+50, dt+210)
#st_gnss_vz[0].plot()
#print(st_gnss_vz[0].max())

# TIDY UP GNSS VELOCITY TRACE -- FILTER THE HECK OUTTA THAT THING!!!!!!!!!!!
st_gnss_vz[0].filter('bandpass', freqmin=0.01, freqmax=0.2, corners=4, zerophase=True)
#st_gnss_vz[0].plot()
#print(st_gnss_vz[0].max())
#st_gnss_vz[0].filter('lowpass', freq=0.45, corners=2, zerophase=True)

# LOAD ONE SEISMIC ACCELERATION STREAM FROM IRIS DMC
start = UTCDateTime("2021-07-29T06:16:59.000000Z")
st_seis_vz = client.get_waveforms('AK', 'P16K', "*", 'BHE', (start), (start+160), attach_response=True)
#st_seis_az[0].plot()
#print(st_seis_az[0].max())

############ What's-in-there flags ################################
### The nearest strong-motion time-stamps are rounded to the ######
### corresponding GNSS ones. 'print' lines help with this #########
#print(st_seis_az[0].data)
#print(len(st_seis_az[0]))
#print(st_seis_az[0].stats.starttime)
#print(st_seis_az[0].times("utcdatetime"))
####################################################################

# TIDY UP VEL TRACE
st_seis_vz[0].remove_response(inventory=None, output="VEL", plot=False)
st_seis_vz[0].detrend('linear')
st_seis_vz[0].filter('bandpass', freqmin=0.01, freqmax=0.2, corners=4, zerophase=True)
#st_seis_az[0].filter('lowpass', freq=5, corners=2, zerophase=True)
#st_seis_az[0].plot()
#print(st_seis_az[0].max())

# CREATE A SEISMIC VELOCITY TRACE
st_seis_vz = st_seis_vz[0].copy()
#st_seis_vz.integrate(method='cumtrapz')
hz = st_seis_vz.data[::50]   # Brute-force 'downsampling' of strong-motion derived velocity to 1 Hz post-filtering and -integration
#st_seis_vz.plot()
#print(st_seis_vz.max())

# CREATE AN EMPTY STREAM AND POPULATE WITH 1Hz STRONG-MOTION DERIVED VELOCITY DATA
st_sm_vz = Stream(Trace())
st_sm_vz[0].stats.network = 'AK'
st_sm_vz[0].stats.station = 'P16K'
st_sm_vz[0].stats.channel = 'BHE'    
st_sm_vz[0].stats.starttime = UTCDateTime("2021-07-29T06:16:59.0Z")
st_sm_vz[0].stats.sampling_rate = 1
st_sm_vz[0].stats.calib = 1
st_sm_vz[0].data = hz

# NORMALIZE EACH TRACE WITH ITS ABSOLUTE MAXIMUM
#st_seis_vz.normalize()
#st_sm_vz[0].normalize()
#st_gnss_vz[0].normalize()

# CROSS CORRELATION
lags = np.arange(-210, 211)  
cc = correlate(st_gnss_vz[0].data, st_sm_vz[0].data, 210)
print(max(cc))
print(compute_shift(cc))

# CREATE STREAMS WITH ALL TRACES AND PLOT THEM
plt.rcParams['figure.figsize'] = [12,8]
plt.rcParams.update({'font.size': 12})
#plt.style.use('seaborn')

fig, axx = plt.subplots(3,2, sharex=False, sharey=False)
axx[0, 0].plot(st_gnss_vz[0].times(), st_gnss_vz[0].data, 'k-', linewidth=1, label=st_gnss_vz[0].stats['station'])
axx[1, 0].plot(st_sm_vz[0].times(), st_sm_vz[0].data, 'k-', linewidth=1, label=st_sm_vz[0].stats['station'])
axx[2, 0].plot(st_seis_vz.times(), st_seis_vz.data, 'k-', linewidth=1, label=st_seis_vz.stats['station'])
axx[0, 0].set_title(f'GNSS instantaneous velocity (1Hz)')
axx[1, 0].set_title(f'lag-corrected Strong-motion velocity (1Hz)')
axx[2, 0].set_title(f'lag-corrected Strong-motion velocity (50Hz)')
axx[0, 0].set_ylabel('Amplitude (m/s)')
axx[1, 0].set_ylabel('Amplitude (m/s)')
axx[2, 0].set_ylabel('Amplitude (m/s)')

st_gnss_vz[0].spectrogram(log=False, dbscale=False, wlen=30, show=False, axes=axx[0, 1], samp_rate=st_gnss_vz[0].stats.sampling_rate, per_lap=0.99, cmap='nipy_spectral')
mappable = axx[0,1].images[0]
ax33 = fig.add_axes([0.95, 0.731, 0.01, 0.217])
plt.colorbar(mappable=mappable, cax=ax33)
st_sm_vz[0].spectrogram(log=False, dbscale=False, wlen=30, show=False, axes=axx[1, 1], samp_rate=st_sm_vz[0].stats.sampling_rate, per_lap=0.99, cmap='nipy_spectral')
mappable = axx[1,1].images[0]
ax23 = fig.add_axes([0.95, 0.405, 0.01, 0.217])
plt.colorbar(mappable=mappable, cax=ax23)
axx[2, 1].plot(lags, cc, 'k', linewidth=1, label='lag-corrected \n AC52*P16K (1Hz)')
axx[0, 1].set_title(f'GNSS spectrogram')
axx[1, 1].set_title(f'Seismic (1Hz) spectrogram')
axx[2, 1].set_title(f'Cross-correlation with bandpass 0.01-0.2 Hz') #Cross-correlation for bandpass 0.001-0.5 Hz
axx[1, 1].set_xlabel('time 70s after origin-time (s)')
axx[0, 0].set_xlabel('time 50s after origin-time (s)')
axx[2, 0].set_xlabel('time 70s after origin-time (s)')
axx[2, 1].set_xlabel('lags')
axx[0, 1].set_ylabel('Frequency (Hz)')
axx[1, 1].set_ylabel('Frequency (Hz)')
plt.tight_layout()

# TRACE LEGENDS
for row in axx[:,0]:
    ll = row.legend(loc=1)
    plt.setp(ll.get_texts(), color='red') #color legend

ll = axx[2,1].legend(loc=1)
plt.setp(ll.get_texts(), color='red') #color legend

plt.show()
