function score = autoCleanInterp_calibration(di,x)

    
    tmp_interp = autoCleanInterp(di.dens3,x(1),x(2));
    
    di.dens3(tmp_interp)
    di.dens2(tmp_interp)
    
    
    tmp = sum(isnan(di.dens),2)>10;
    tmp(di.day)=true;

    errt1 = sum(tmp_rain(~tmp) & ~true_rainday(~tmp)); % false positive: remove while no rain
    errt2 = sum(~tmp_rain(~tmp) & true_rainday(~tmp)); % false negative: leave rain
    score = (errt1+errt2)/sum(~tmp);

end