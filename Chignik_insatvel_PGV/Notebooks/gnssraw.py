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
client=Client("IRIS")
import numpy as np

# READ INSTAVEL ASCII FILE
infile = 'Chignik_phase_vels/ab13_1hz.L2.igr.phase_vels'
site = np.genfromtxt(infile, dtype=None, encoding=None, names=['date','time','dt','lat','lon','height','vx', 'vy', 'vz', 'sigx', 'sigy', 'sigz', 'vdt', 'vn', 've', 'vu', 'sign', 'sige', 'sigu'])

# CREATE AN EMPTY STREAM AND POPULATE WITH INSTAVEL DATA
st_gnss_vz = Stream(Trace())
st_gnss_vz[0].stats.network = '--'
st_gnss_vz[0].stats.station = 'AB13'
st_gnss_vz[0].stats.channel = 'IVE'       # for [I]nst[V]el [E] component]
st_gnss_vz[0].stats.starttime = site["date"][0] + 'T' + site["time"][0]
st_gnss_vz[0].stats.sampling_rate = 1
st_gnss_vz[0].stats.calib = 1
st_gnss_vz[0].data = site['ve']

# TIDY UP GNSS VELOCITY TRACE -- FILTER THE HECK OUTTA THAT THING!
#st_gnss_vz[0].filter('bandpass', freqmin=0.01, freqmax=0.49, corners=4, zerophase=True)
#st_gnss_vz[0].filter('lowpass', freq=0.499, corners=2, zerophase=True)

# LOAD ONE SEISMIC ACCELERATION STREAM FROM IRIS DMC
origintime = UTCDateTime("2021-07-29T06:15:49")
st_seis_az = client.get_waveforms('AK', 'S15K', "*", 'HNE', (origintime-60), (origintime+540), attach_response=True) 

# NORMALIZE EACH TRACE WITH ITS ABSOLUTE MAXIMUM
#st_gnss_vz[0].normalize()

# TIDY UP ACCELERATION TRACE
st_seis_az.detrend('linear')
#st_seis_az.filter('bandpass', freqmin=0.01, freqmax=0.49, corners=4, zerophase=True)
st_seis_az.filter('lowpass', freq=20, corners=2, zerophase=True)

# CREATE A SEISMIC VELOCITY TRACE
st_seis_vz = st_seis_az.copy()
st_seis_vz.integrate(method='cumtrapz')

# NORMALIZE EACH TRACE WITH ITS ABSOLUTE MAXIMUM
#st_seis_az.normalize()
#st_seis_vz.normalize()

# CREATE A STREAM WITH ALL TRACES AND PLOT THEM
st = Stream([st_gnss_vz[0], st_seis_vz[0], st_seis_az[0]])
#st.plot(equal_scale=False, automerge=False)
st.plot(outfile="AB13-S15K-E-RAW.pdf", format="pdf", equal_scale=False, automerge=False, ticks="5")

# IGNORE
# st_tmp = Stream([st_gnss[0], st_seis[0], st_seis[0], st_seis[0]])
# pre_filt = (0.005, 0.006, 30.0, 35.0)
# st_tmp[2].remove_response(output='VEL', pre_filt=pre_filt)
# st_tmp[3].remove_response(output='DISP', pre_filt=pre_filt)
# st_tmp.plot()
