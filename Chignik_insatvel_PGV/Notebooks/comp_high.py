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
import matplotlib.pyplot as plt
client=Client("IRIS")

def compute_shift(cc):
    zero_index = int((len(cc)-1)/2)
    print(zero_index)
    if np.argmax(cc) > zero_index:
        shift = zero_index - np.argmax(cc)
    else:
        shift = np.argmax(cc) - zero_index
    return shift

# READ INSTAVEL ASCII FILE
infile = 'Chignik_phase_vels/ac12_1hz.L2.igr.phase_vels'
site = np.genfromtxt(infile, dtype=None, encoding=None, names=['date','time','dt','lat','lon','height','vx', 'vy', 'vz', 'sigx', 'sigy', 'sigz', 'vdt', 'vn', 've', 'vu', 'sign', 'sige', 'sigu'])

# CREATE AN EMPTY STREAM AND POPULATE WITH INSTAVEL DATA
gnss = Trace()
gnss.stats.network = '--'
gnss.stats.station = 'AC12'
gnss.stats.channel = 'IVE'       # for [I]nst[V]el [E] component]
gnss.stats.starttime = site["date"][0] + 'T' + site["time"][0]
gnss.stats.sampling_rate = 1
gnss.stats.calib = 1
gnss.data = site['ve']
gnss.detrend('linear')
gnss.detrend('demean')
dt = UTCDateTime("2021-07-29T06:15:49.000000Z")
gnss.trim(dt+10, dt+170) #70 sec after 2021-07-29T06:15:48.998400Z


# LOAD ONE SEISMIC ACCELERATION STREAM FROM IRIS DMC
start = UTCDateTime("2021-07-29T06:15:59.008400Z")
seis_a = client.get_waveforms('AK', 'CHN', "*", 'BNE', (start), (start+160), attach_response=True)

# TIDY UP ACCELERATION TRACE
seis_a[0].remove_response(inventory=None, output="ACC", plot=False)
seis_a[0].detrend('linear')

# TIDY UP GNSS VELOCITY TRACE -- FILTER THE HECK OUTTA THAT THING!!!!!!!!!!!
stg = read()
stg.clear()
sts = read()
sts.clear()
corrs = []

shifts = []
frequencies = [0.499, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05]
for fs in frequencies:
    gn = gnss.copy()
    ss = seis_a[0].copy()
    gn.filter('bandpass', freqmin=0.001, freqmax=fs, corners=4, zerophase=True) 
    ss.filter('bandpass', freqmin=0.001, freqmax=fs, corners=4, zerophase=True)

    #CREATE A SEISMIC VELOCITY TRACE
    seis_v = ss.copy()
    seis_v.integrate(method='cumtrapz')
    hz = seis_v.data[::50]   # Brute-force 'downsampling' of strong-motion derived velocity to 1 Hz post-filtering and -integration

    # CREATE AN EMPTY STREAM AND POPULATE WITH 1Hz STRONG-MOTION DERIVED VELOCITY DATA
    sm_v = Trace()
    sm_v.stats.network = 'AK'
    sm_v.stats.station = 'CHN'
    sm_v.stats.channel = 'BNE'    
    sm_v.stats.starttime = UTCDateTime("2021-07-29T06:15:59.0Z")
    sm_v.stats.sampling_rate = 1
    sm_v.stats.calib = 1
    sm_v.data = hz

    stg.append(gn)
    sts.append(sm_v)

# NORMALIZE EACH TRACE WITH ITS ABSOLUTE MAXIMUM
    #st_seis_vz.normalize()
    #st_sm_vz[0].normalize()
    #st_gnss_vz[0].normalize()

# CROSS CORRELATION
    lags = np.arange(-160, 161)  
    cc = correlate(gn.data, sm_v.data, 160)
    corr = max(cc)
    print(np.argmax(cc))
    shift = compute_shift(cc)
    print(corr, shift)
    corrs.append(corr)
    shifts.append(shift)

###Flags###
#print(corrs)
#print(shifts)
#print(stg)
#print(sts)    

# CREATE STREAMS WITH ALL TRACES AND PLOT THEM
plt.rcParams['figure.figsize'] = [10,8]
plt.rcParams.update({'font.size': 8})
#plt.style.use('seaborn')

fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(nrows=6, sharex=False)
ax1.plot(frequencies, corrs, 'cd')
ax1.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax1.set_title('AC12 (GNSS) vs CHN (strong-motion) (1Hz sampling)', fontsize = 14)
ax12=ax1.twinx()
ax12.plot(frequencies, shifts, 'g1')
ax1.set_xlabel('Frequency, f (Hz) [Bandpass filter 0.001-f Hz]', fontsize = 10)
ax1.set_ylabel('Cross-correlation', fontsize = 10, color = 'c')
ax12.set_yticks([-2, -1, 0, 1, 2])
ax12.set_ylabel('Shifts (s)', fontsize = 10, color = 'g')
 
ax2.plot(stg[0].times(), stg[0].data, 'k-', linewidth=1)
ax22=ax2.twinx()
ax22.plot(sts[0].times(), sts[0].data, 'r-', linewidth=1)
ax2.text(120, -0.022, r"Bandpass 0.001-0.5Hz", color="b", fontsize=10)

ax3.plot(stg[5].times(), stg[5].data, 'k-', linewidth=1)
ax3.set_ylabel('Amplitudes in m/s (GNSS)', fontsize = 10, color = 'k')
ax32=ax3.twinx()
ax32.plot(sts[5].times(), sts[5].data, 'r-', linewidth=1)
ax32.set_ylabel('Amplitudes in m/s (Strong motion)', fontsize = 10, color = 'r')
ax3.text(120, -0.016, r"Bandpass 0.001-0.25Hz", color="b", fontsize=10)

ax4.plot(stg[6].times(), stg[6].data, 'k-', linewidth=1)
ax42=ax4.twinx()
ax42.plot(sts[6].times(), sts[6].data, 'r-', linewidth=1)
ax4.text(120, -0.015, r"Bandpass 0.001-0.2Hz", color="b", fontsize=10)

ax5.plot(stg[8].times(), stg[8].data, 'k-', linewidth=1)
ax52=ax5.twinx()
ax52.plot(sts[8].times(), sts[8].data, 'r-', linewidth=1)
ax5.text(120, -0.006, r"Bandpass 0.001-0.1Hz", color="b", fontsize=10)

ax6.plot(stg[9].times(), stg[9].data, 'k-', linewidth=1)
ax62=ax6.twinx()
ax62.plot(sts[9].times(), sts[9].data, 'r-', linewidth=1)
ax6.set_xlabel('Time 10s after origin-time (s)')
ax6.text(120, -0.005, r"Bandpass 0.001-0.05Hz", color="b", fontsize=10)

plt.tight_layout()

plt.show()
