% folder ='Density_estimationMap_ImageOverlay/';
% load('data/Density_estimationMap')
% data = log(g.dens_est);

% folder ='Density_simulationMap_ImageOverlay/';
% load('data/Density_simulationMap')
% data = log(real_dens);

% folder ='FlightSpeed_estimationMap_ImageOverlay/';
% load('data/FlightSpeed_estimationMap')
% data = g.dens_est;

% folder ='FlightSpeed_simulationMap_ImageOverlay/';
% load('data/FlightSpeed_simulationMap_reassemble')
% data = real_dens;

% folder ='FlightDir_estimationMap_ImageOverlay/';
% load('data/FlightDir_estimationMap')
% data = gdir.dd_est;


R = georefcells([g.lat(1) g.lat(end)],[g.lon(1) g.lon(end)], [g.nlat g.nlon])

i=330;

geotiffwrite('test.tif',data(:,:,i),R)