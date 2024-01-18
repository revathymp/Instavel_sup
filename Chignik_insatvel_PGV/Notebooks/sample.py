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
from obspy import UTCDateTime 				
from obspy import Trace
from obspy import Stream
from obspy import read_inventory
from obspy.clients.fdsn import Client 		
client=Client("IRIS")
import numpy as np

# LOAD ONE SEISMIC ACCELERATION STREAM FROM IRIS DMC
origintime = UTCDateTime("2021-07-29T06:15:49.000000Z")
st_seis_az = client.get_waveforms('AK', 'S15K', "*", 'HNE', (origintime-20), (origintime+180), attach_response=True) 
#times = st_seis_az[0].times('utcdatetime')
#index = times.searchsorted(origintime)
#print(times[index])
#print(st_seis_az[0].data[index])

st_seis_az[0].stats.sampling_rate = 1
print(st_seis_az[0].times("utcdatetime"))
#tr = st_seis_az[0]
#print(tr.stats)

