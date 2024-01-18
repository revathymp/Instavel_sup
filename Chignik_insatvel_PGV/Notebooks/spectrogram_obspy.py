from obspy.clients.fdsn import Client
client_wm = Client("IRIS")
from obspy import read
from obspy import Stream
from obspy import UTCDateTime
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates


# LOAD WAVEFORM DATA. DO A LITTLE PRE-PROCESSING
#tday = UTCDateTime("2021-08-09 07:45") #landslide
#tday = UTCDateTime("2021-10-11 09:10:24.000") #Chignik
tday = UTCDateTime("2015-10-18 05:17:00.000") #taan fjord
event_name = "Taan_fjord_log_temp"
print('grabbing waveforms for ' + tday.strftime("%Y%m%d"))
st = client_wm.get_waveforms("AK", "RKAV", "*", "BHZ", tday, tday+360, attach_response=True)
st.detrend("linear")
st.detrend("demean")


#fig, (ax1, ax2) = plt.subplots(2,1)   
fig, (ax1, ax2, ax3) = plt.subplots(3,1)   
date_format = mdates.DateFormatter('%H:%M:%S')
timevector = st[0].times("matplotlib")

# PLOT WAVEFORM USING BUILT-IN OBSPY PLOT ROUTINE
#st.plot()
ax1.plot(st[0].times("matplotlib"),st[0].data,'-',color='k',lw=.2)
ax1.xaxis.set_major_formatter(date_format)
#ax1.xaxis.set_minor_locator(minutes)
ax1.xaxis_date()
ax1.tick_params(direction='in')
ax1.axes.xaxis.set_ticklabels([])
ax1.set_xlim(timevector[0], timevector[-1]) 
ax1.set_ylabel('counts')
#ax1.set_xlabel("Time [UTC]" )
text = (st[0].stats.station)+' '+(st[0].stats.channel)
ax1.text(0.01, 0.88, text, transform=ax1.transAxes, fontsize=10, fontweight='bold', color='black')


# PLOT SPECTROGRAM USING BUILT-IN OBSPY ROUTINE
# this works, but is hard to control
#Linear scale
st.spectrogram(log=False, axes=ax2, per_lap=0.5, dbscale=True, wlen=5,  cmap='nipy_spectral')
#####THESE SHOULD WORK FOR COLORBAR BUT NEED SOME ADJUSTMENTS: 
###mappable = ax2.images[0]
####ax33 = fig.add_axes([0.83, 0.1, 0.03, 0.6])
####plt.colorbar(mappable=mappable, cax=ax33)

ax2.set_ylabel('freq. (Hz)')
ax2.text(0.01, 0.88, text, transform=ax2.transAxes, fontsize=10, fontweight='bold', color='white')
ax2.axes.xaxis.set_ticklabels([])
#ax2.set_xlabel("Time [sec]" )

#Log scale
st.spectrogram(log=True, axes=ax3, per_lap=0.5, dbscale=True, wlen=5,  cmap='nipy_spectral')

ax3.set_ylabel('freq. (Hz)')
ax3.text(0.01, 0.88, text, transform=ax3.transAxes, fontsize=10, fontweight='bold', color='white')
ax3.set_xlabel("Time [sec]" )

fig.suptitle(event_name + ' ' + tday.strftime("%Y/%m/%d"))
fig.savefig(event_name + '_spectogram_obspy', bbox_inches='tight')

