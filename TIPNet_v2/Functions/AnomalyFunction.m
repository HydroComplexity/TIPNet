function [data_out] = AnomalyFunction(data,dt, dt_units, anom,anom_units)
%AnomalyFunction: calculate anomaly of data (i.e. 5-day anomaly)

% data: vector of input timeseries data
% dt and dt_units: timestep and associated unit (e.g. 5 minutes)
% anom and anom_units: anomaly and unit (e.g. 10 day anomaly)

    tsteps = length(data);
    anomaly = zeros(1,tsteps);
    
    %dt_units: 'mins', 'hrs','days', 'months'
    %anom_units: 'days','years'
    
%determine number of timesteps per unit in terms of anom_units
if strcmp(anom_units,'days')==1 && strcmp(dt_units,'mins')==1
  steps_per_unit = 24*60/dt;
elseif strcmp(anom_units,'days')==1 && strcmp(dt_units,'hrs')==1
  steps_per_unit = 24/dt;
elseif strcmp(anom_units,'days')==1 && strcmp(dt_units,'days')==1
  steps_per_unit = 1;
elseif strcmp(anom_units,'years')==1 && strcmp(dt_units,'mins')==1
  steps_per_unit = 24*60*365/dt;
elseif strcmp(anom_units,'years')==1 && strcmp(dt_units,'hrs')==1
  steps_per_unit = 24*365/dt;
elseif strcmp(anom_units,'years')==1 && strcmp(dt_units,'days')==1
  steps_per_unit = 365/dt;
elseif strcmp(anom_units,'years')==1 && strcmp(dt_units,'months')==1
  steps_per_unit = 12/dt;
elseif strcmp(anom_units,'years')==1 && strcmp(dt_units,'years')==1
  steps_per_unit = 1;
end



%should already be integer, but in case not..
steps_per_unit = round(steps_per_unit);

%if anom is an odd number (e.g. 5 days) --> use 2 days on each side
%if anom is an even number (e.g. 10 days) --> use 5 days on each side
div2 = floor(anom./2);



% compute anomaly
for j = 1:tsteps
    
    if (j-div2*steps_per_unit)>0 && (j+div2*steps_per_unit)<(tsteps+1)
        anomaly(j) = mean(data((j-div2*steps_per_unit):steps_per_unit:(j+div2*steps_per_unit)));
    elseif (j-div2*steps_per_unit)<1
        %fprintf('index = %d, average from %d to %d\n',j,j,j+div2*steps_per_unit)
        anomaly(j) = mean(data(j:steps_per_unit:(j+div2*steps_per_unit)));
    else
        anomaly(j) = mean(data((j-div2*steps_per_unit):steps_per_unit:j));
    end
    
end

    
data_out = data - anomaly';




