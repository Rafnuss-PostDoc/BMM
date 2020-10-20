# Wind Support to Bird migration

## Method

Cleaning procedure of bird density and speed vector (bu,bv) as in 2nd paper with the exepction of no interpolation to ground level. 

Wind data 
- from ERA5 reanalysis at pressure level (from 1000hPa to 550hPa). 
- Hourly resolution and 0.25°x0.25°. 
- Linearly interpolated (time-space 4D) to each WR datapoint.

ground speed: sqrt(bu^2 + bv^2)
air speed: sqrt( (bu-wu)^2 + (bv-wv)^2) 


How do bird compensate for wind?
- strengh and direction of wind
- prefered direction
- altitude
- time/season
- timing of deparute: are bird compensating more early in the night?
(- specie of bird, size)

The overwhelming evidence in support of the prediction that birds increase their air speed in a headwind and decrease their air speed in a tailwind.

NEED TO READ:
10.1016/j.jtbi.2012.05.026
10.1098/rsif.2014.0588
10.1093/beheco/ars078