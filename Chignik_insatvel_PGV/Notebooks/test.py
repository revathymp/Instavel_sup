from obspy import UTCDateTime
from obspy import Trace
from obspy import Stream
from obspy import read_inventory
from obspy.clients.fdsn import Client
client = Client("IRIS")
starttime = UTCDateTime("2021-07-29T06:10:49")
endtime = UTCDateTime("2021-07-29T06:22:49")
inventory = client.get_waveforms(network="AK", station="P17K", location="*", channel="??E",
                                starttime=starttime,
                                endtime=endtime)
print(inventory)
