# 
# Python script to download GEFS data from the AWS server.
# 
# For PICCOLO period, GEFS data is available on the AWS server at 0.5º resolution.
# Forecasts go out 384 hours (16 days) with initializations every 6 hours.
# 
# This script takes the start date as the desired initialization date, with all
# data coming from that forecast.
# 
# James Ruppert
# 17 Nov 2024

import numpy as np
import requests

##### MAIN SETTINGS #####

download_dir = 'bcsics/gefs/' # Directory to save files

# Start, end dates
YYMMDDHH1 = '2024090100'
YYMMDDHH2 = '2024090400'
forecast_interval = 3 # hours

# Choose number of ensemble members
nens = 5 # e.g., will use perturbation ensemble members 1-5 for nens=5 (neglect 0 = CTL)

##### MAIN PROCS #####

# Example web addresses:
# Variables are contained in three separate paths/file types
# https://noaa-gefs-pds.s3.amazonaws.com/gefs.20240901/00/atmos/pgrb2ap5/gep01.t00z.pgrb2a.0p50.f102
# https://noaa-gefs-pds.s3.amazonaws.com/gefs.20240901/00/atmos/pgrb2bp5/gep01.t00z.pgrb2b.0p50.f111
# File set ['pgrb2sp25', 'pgrb2s'] is 0.5 deg single-layer fields, so not used

server = "https://noaa-gefs-pds.s3.amazonaws.com/"

# Create time array
YY1, MM1, DD1, HH1 = YYMMDDHH1[:4], YYMMDDHH1[4:6], YYMMDDHH1[6:8], YYMMDDHH1[8:10]
YY2, MM2, DD2, HH2 = YYMMDDHH2[:4], YYMMDDHH2[4:6], YYMMDDHH2[6:8], YYMMDDHH2[8:10]
time1 = np.datetime64(YY1+'-'+MM1+'-'+DD1+'T'+HH1+':00:00')
time2 = np.datetime64(YY2+'-'+MM2+'-'+DD2+'T'+HH2+':00:00')
forecast_hours = np.arange(time1, time1+np.timedelta64(385, 'h'), forecast_interval, dtype='datetime64[h]')
forecast_hours = forecast_hours[np.where(forecast_hours <= time2)]
fhh = np.arange(0, forecast_hours.size*forecast_interval, 3)

# Use a single forecast initialization
http_init = server+"gefs."+YY1+MM1+DD1+"/"+HH1

# Variables are contained in three separate paths/file types
dset = [
    ['pgrb2ap5', 'pgrb2a'],
    ['pgrb2bp5', 'pgrb2b'],
    # ['pgrb2sp25', 'pgrb2s']
]

for ifh in fhh:
    ifh_str = str(ifh).zfill(3)
    for iens in range(1, nens+1):
        iens_str = str(iens).zfill(2)
        for d in dset:
            file = "gep"+iens_str+".t"+HH1+"z."+d[1]+".0p50.f"+ifh_str
            http = http_init+"/atmos/"+d[0]+"/"+file
            response = requests.get(http)
            if response.status_code == 200:
                file_path = download_dir + file
                with open(file_path, 'wb') as file:
                    file.write(response.content)
            else:
                print("Failed on: ",http)
