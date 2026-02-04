import ControlFlowIMERG as cf

## trigger download of IMERG data
### PARAMETER TO CHANGE ###
start = "2014-05-01 00:00:00"
end = "2014-05-31 00:00:00"
#lat_min, lat_max = 14.0, 33.0  
#lon_min, lon_max = -118.0, -86.0
#region = [lat_min, lat_max, lon_min, lon_max]


cf.daily_routine(start, end, nside=128)