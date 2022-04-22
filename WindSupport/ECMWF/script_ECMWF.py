#!/usr/bin/env python
import cdsapi
import os

os.chdir('/mnt/c/Users/rnussba1/switchdrive/BMM/WindSupport/ECMWF')

# https://cds.climate.copernicus.eu/api-how-to
c = cdsapi.Client()
data = c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': ['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_temperature','total_precipitation'],
        'year': '2018',
        'month': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'],
        'day': [ '01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'],
        'time': [ '00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
        'format': 'netcdf',
        'area': [55, -5, 43, 14],  # North, West, South, East. Default: global
    },"2018_srf.nc")


 # 1000/975/950/925/900/875/850/825/800/775/750/700/650/600/550/500/450/400/350/300/250/225/200/175/150/125/100/70/50/30/20/10/7/5/3/2/1
# https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels
data = c.retrieve(
    'reanalysis-era5-pressure-levels',
    {
        'product_type':'reanalysis',
        'variable':[
            'u_component_of_wind', 'v_component_of_wind', 'temperature'
        ],
        'pressure_level':[
            #"1000","975","950","925", #-> 2018_pressure_1
            #"900","875","850","825", #-> 2018_pressure_2
            #"800","775","750","700" #-> 2018_pressure_3
            "650","600","550","500" #-> 2018_pressure_4
        ],
        'year': '2018',
        'month': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'],
        'day': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31'],
        'time': [ '00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00'],
        'format': 'netcdf',
        'area': [55, -5, 43, 14],  # North, West, South, East. Default: global
    }, "2018_pressure_4.nc")



data = c.retrieve(
    'reanalysis-era5-pressure-levels-monthly-means',
    {
        'product_type':'monthly_averaged_reanalysis',
        'variable':[
            'u_component_of_wind', 'v_component_of_wind'
        ],
        'pressure_level':[
            "1000","975","950","925","900","875","850","825","800","775","750","700","650","600","550",'500'
        ],
        'year': ['2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019'],
        'month': ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12'],
         'time': '00:00',
        'format': 'netcdf',
        'area': [55, -5, 43, 14],  # North, West, South, East. Default: global
    }, "2000_2019_monthly_pressure2.nc")