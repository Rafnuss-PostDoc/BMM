#!/usr/bin/env python
import cdsapi
import os

os.chdir('/mnt/c/Users/rnussba1/switchdrive/BMM/windfarm/ECMWF')

# https://cds.climate.copernicus.eu/api-how-to
c = cdsapi.Client()
data = c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': ['10m_u_component_of_wind', '10m_v_component_of_wind','100m_u_component_of_wind', '100m_v_component_of_wind'],
        'year': '2018',
        'month': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'],
        'day': [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'],
        'time': [ '00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
        'format': 'netcdf',
        'area': [56, -5, 42, 17],  # North, West, South, East. Default: global
    },"2018_srf.nc")