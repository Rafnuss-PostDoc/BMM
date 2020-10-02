function score = autoCleanRain_calibration(di,x)

    init_nan = sum(isnan(di.dens),2)>10;
    true_rainday = sum(isnan(di.dens2),2)>10;
    tmp_rain = autoCleanRain(di,x(1),x(2),x(3),x(4),x(5));
    
    % set as nan if no value available. (rain is not detected when there is
    % not data)
    tmp_rain(init_nan)=true;
    
    % Day is anyway removed
    init_nan(di.day)=true;
    true_rainday(di.day)=true;
    tmp_rain(di.day)=true;

    % figure;imagesc([init_nan true_rainday tmp_rain]')

    errt1 = sum(tmp_rain(~init_nan) & ~true_rainday(~init_nan)); % false positive: remove while no rain
    errt2 = sum(~tmp_rain(~init_nan) & true_rainday(~init_nan)); % false negative: leave rain
    score = (errt1+errt2)/sum(~init_nan);

end