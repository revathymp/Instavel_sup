# This is an example script with two goal:
#     1) read a GNSS instavel and convert it into an obspy stream
#     2) demonstrate a few basic stream manipulations
#
# Obspy is pretty easy to get started with but can be a bit overwhelming. 
#      https://docs.obspy.org/tutorial/ is a good place to start
# At a minimum it is worth understanding the somewhat subtle difference
# between an obspy Trace and an Obspy Stream. Both are used below.
#
# Obspy can be found in the conda-forge channel. 
# I *think* obspy requires Python <=3.8
# For example:
#
# conda create -n obspy python=3.7
# conda config --add channels conda-forge
# conda activate obspy  
# conda install obspy

from obspy import UTCDateTime 				
from obspy import Trace
from obspy import Stream
from obspy import read_inventory
from obspy.clients.fdsn import Client 		
from obspy.signal.cross_correlation import correlate
import numpy as np
import pandas as pd
#import noise
import matplotlib.pyplot as plt
client=Client("IRIS")
import numpy as np

# READ INSTAVEL ASCII FILE
infile = 'Chignik_phase_vels/ac34_1hz.L2.igr.phase_vels'
site = np.genfromtxt(infile, dtype=None, encoding=None, names=['date','time','dt','lat','lon','height','vx', 'vy', 'vz', 'sigx', 'sigy', 'sigz', 'vdt', 'vn', 've', 'vu', 'sign', 'sige', 'sigu'])

# CREATE AN EMPTY STREAM AND POPULATE WITH INSTAVEL DATA
st_gnss_vz = Stream(Trace())
st_gnss_vz[0].stats.network = '--'
st_gnss_vz[0].stats.station = 'AC34'
st_gnss_vz[0].stats.channel = 'IVE'       # for [I]nst[V]el [E] component]
st_gnss_vz[0].stats.starttime = site["date"][0] + 'T' + site["time"][0]
st_gnss_vz[0].stats.sampling_rate = 1
st_gnss_vz[0].stats.calib = 1
st_gnss_vz[0].data = site['ve']
st_gnss_vz[0].detrend('spline', order = 3, dspline = 100)
st_gnss_vz[0].detrend('demean')
dt = UTCDateTime("2021-07-29T06:15:49.000000Z")
st_gnss_vz[0].trim(dt-40, dt+360)
#print(st_gnss_vz[0].times("utcdatetime"))

##### Frequency spectrum testing flags #############################
#st_gn_vz = st_gnss_vz[0].copy()
#st_gn_vz.trim(dt, dt+120)
#st_gn_vz.spectrogram(log=False, title='AC34', wlen=5)
#st_gn_vz.plot()
#st_gnss_vz[0].plot()
#st_gnss_vz.spectrogram(log=False, title='AC34', wlen=5)
#####################################################################


# TIDY UP GNSS VELOCITY TRACE -- FILTER THE HECK OUTTA THAT THING!
st_gnss_vz[0].filter('bandpass', freqmin=0.001, freqmax=0.05, corners=4, zerophase=True)
#st_gnss_vz[0].filter('lowpass', freq=0.45, corners=2, zerophase=True)


####### Frequency spectrum testing flags ###########################
#st_gn_vz = st_gnss_vz[0].copy()
#st_gn_vz.trim(dt, dt+120)
#st_gn_vz.spectrogram(log=False, title='AC34', wlen=5)
#st_gn_vz.plot()
####################################################################

# LOAD ONE SEISMIC ACCELERATION STREAM FROM IRIS DMC
start = UTCDateTime("2021-07-29T06:15:08.998400Z")
st_seis_az = client.get_waveforms('AK', 'S19K', "*", 'HNE', (start), (start+400), attach_response=True)

############ What's-in-there flags ################################
### The nearest strong-motion time-stamps are rounded to the ######
### corresponding GNSS ones. 'print' lines help with this #########

#print(st_seis_az[0].data)
#print(len(st_seis_az[0]))
#print(st_seis_az[0].stats.starttime)
#print(st_seis_az[0].times("utcdatetime"))
####################################################################


# TIDY UP ACCELERATION TRACE
st_seis_az[0].detrend('linear')
#st_seis_az[0].filter('lowpass', freq=5, corners=2, zerophase=True)
st_seis_az[0].filter('bandpass', freqmin=0.001, freqmax=0.05, corners=4, zerophase=True)


# CREATE A SEISMIC VELOCITY TRACE
st_seis_vz = st_seis_az[0].copy()
st_seis_vz.integrate(method='cumtrapz')
#print(st_seis_vz.data)
hz = st_seis_vz.data[::100]   # Brute-force 'downsampling' of strong-motion velocity to 1 Hz post-filtering and -integration
#print(st_seis_vz.times("utcdatetime")[::100])
#print(hz)
ahz = st_seis_az[0].data[::100]    # Brute-force 'downsampling' of strong motion acceleration to 1 Hz post-filtering

# CREATE AN EMPTY STREAM AND POPULATE WITH 1Hz STRONG-MOTION VELOCITY DATA
st_sm_vz = Stream(Trace())
st_sm_vz[0].stats.network = 'AK'
st_sm_vz[0].stats.station = 'S19K'
st_sm_vz[0].stats.channel = 'HNE'    
st_sm_vz[0].stats.starttime = UTCDateTime("2021-07-29T06:15:09.0Z")
st_sm_vz[0].stats.sampling_rate = 1
st_sm_vz[0].stats.calib = 1
st_sm_vz[0].data = hz
print(st_sm_vz[0].times("utcdatetime"))

# CREATE AN EMPTY STREAM AND POPULATE WITH 1HZ STRONG-MOTION ACCELERATION DATA
st_seis_daz = Stream(Trace())
st_seis_daz[0].stats.network = 'AK'
st_seis_daz[0].stats.station = 'S19K'
st_seis_daz[0].stats.channel = 'HNE' 
st_seis_daz[0].stats.starttime = UTCDateTime("2021-07-29T06:15:09.0Z")
st_seis_daz[0].stats.sampling_rate = 1
st_seis_daz[0].stats.calib = 1
st_seis_daz[0].data = ahz
#print(st_sm_vz[0].times("utcdatetime"))

# CREATE A GNSS ACCELERATION TRACE
st_gnss_az = st_gnss_vz[0].copy()
st_gnss_vz.differentiate(method='gradient')

# NORMALIZE EACH TRACE WITH ITS ABSOLUTE MAXIMUM
st_seis_az[0].normalize()
st_seis_vz.normalize()
st_sm_vz[0].normalize()
st_gnss_az.normalize()
st_seis_daz[0].normalize()
st_gnss_vz[0].normalize()

### Cross Correlation using Pandas
series1, series2 = pd.Series(st_gnss_vz[0].data),  pd.Series(st_sm_vz[0].data)
#window = 1
lags = np.arange(-400, 401)  # contrained
#rs = rs = np.nan_to_num([crosscorr(series1, series2, lag) for lag in lags])
cc = correlate(st_gnss_vz[0].data, st_sm_vz[0].data, 400) 
#cc = correlate(st_gnss_vz[0].data, st_sm_vz[0].data, 400)
print(lags)
print(cc)
#traceCCF2 = Stream(Trace())


#traceCCF2[0].stats.network = '-'
#traceCCF2[0].stats.station = 'Cross-correlated'
#traceCCF2[0].stats.channel = '-'      
#traceCCF2[0].stats.starttime = UTCDateTime("2021-07-29T06:15:09.0Z")
#traceCCF2[0].stats.sampling_rate = 1
#traceCCF2[0].stats.calib = 1
#traceCCF2[0].data = cc
#traceCCF2[0].plot()

#traceCCF2 = Trace()
#traceCCF2.data = cc
#traceCCF2.times(reftime=st_sm_vz[0].stats.starttime)

# CREATE STREAMS WITH ALL TRACES AND PLOT THEM

plt.rcParams['figure.figsize'] = [10,8]
plt.rcParams.update({'font.size': 12})
#plt.style.use('seaborn')

fig, axx = plt.subplots(3,2, sharex=False, sharey=False)
axx[0, 0].plot(st_gnss_vz[0].times(), st_gnss_vz[0].data, 'k-', linewidth=0.2, label=st_gnss_vz[0].stats['station'])
axx[1, 0].plot(st_sm_vz[0].times(), st_sm_vz[0].data, 'k-', linewidth=0.2, label=st_sm_vz[0].stats['station'])
axx[2, 0].plot(lags, cc, 'k', linewidth=0.5)
#axx[2, 0].plot(st_seis_vz.times(), st_seis_vz.data, 'k-', linewidth=0.2, label=st_seis_vz.stats['station'])
axx[0, 0].set_title(f'GNSS instantaneous velocity (1Hz)')
axx[1, 0].set_title(f'Strong-motion velocity (1Hz)')
axx[2, 0].set_title(f'cross-correlation for 0.001-0.05 Hz')
#axx[2, 0].set_title(f'Strong-motion velocity (100Hz)')
axx[0, 0].set_ylabel('Normalized amplitude')
axx[1, 0].set_ylabel('Normalized amplitude')
#axx[2, 0].set_ylabel('Normalized amplitude')

st_gnss_vz[0].spectrogram(log=False, wlen=10,show=False, axes=axx[0, 1], samp_rate=st_gnss_vz[0].stats.sampling_rate)
st_sm_vz[0].spectrogram(log=False, wlen=10,show=False, axes=axx[1, 1], samp_rate=st_sm_vz[0].stats.sampling_rate)
#traceCCF2.spectrogram(log=False, wlen=10,show=False, axes=axx[2, 1])
#st_seis_vz.spectrogram(log=False, wlen=10,show=False, axes=axx[2, 1], samp_rate=st_seis_vz.stats.sampling_rate)
axx[0, 1].set_title('Spectrogram GNSS')
axx[1, 1].set_title('Spectrogram seismic (1Hz)')
#axx[2, 0].set_title('Spectrogram of cross-correlated traces')
#axx[2, 1].set_title('Spectrogram seismic (100Hz)')
axx[2, 0].set_xlabel('time (s)')
axx[2, 1].set_xlabel('time (s)')
axx[0, 1].set_ylabel('Frequency')
axx[1, 1].set_ylabel('Frequency')
#axx[2, 1].set_ylabel('Frequency')
plt.tight_layout()

## put legend
for row in axx[:,0]:
    ll = row.legend(loc=1)
    plt.setp(ll.get_texts(), color='red') #color legend

plt.show()

#st_gnss_vz[0].spectrogram(log=False, title='AB13' + str(st_gnss_vz[0].stats.starttime), wlen=8)
#st = Stream([st_gnss_vz[0], st_sm_vz[0], st_seis_vz])
#ast = Stream([st_gnss_az, st_seis_daz[0], st_seis_az[0]])
#st.plot(equal_scale=True, automerge=False)
#ast.plot(equal_scale=True, automerge=False)
#st.plot(outfile="AC67-S19K_vel.png", format="png", equal_scale=True, automerge=False, ticks="5")
#ast.plot(outfile="AC67-S19K_accl.png", format="png", equal_scale=True, automerge=False, ticks="5")
