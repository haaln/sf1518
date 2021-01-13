function [CONSUMPTION] = consumption(VELOCITY)
%CONSUMPTION    Power consumption of an electric vehicle
%  CONSUMPTION(VELOCITY) computes power consumption of a Tesla roadster by
%  interpolating using cubic splines inside a pre-defined boundary
%  with its corresponding function values imported from a data structure.
%
%  The function will, by default, return a value of the power consumption
%  and optionally return the life-expectancy of the battery measured in
%  travel distance for that constant velocity.
%
%  The function will also accept matrices which enables for computation of
%  several consumption values of their corresponding velocity
%
%  Example: 
%       consumption(50)
%
%       94.6112
%
%  Example:
%       consumption = consumption([10 20 30])
%
%       consumption = [108.6916   86.2297   83.4635]
%
%
%  See also VELOCITY, TIME_TO_DESTINATION, TOTAL_CONSUMPTION


%tries to load roadster returns error if it cant find it
try load roadster.mat speed_kmph consumption_Whpkm
catch
    error('file not found')
end

%checks if the velocity is bounded to the data
%checks if input argument is a number
if ~ismatrix(VELOCITY)
    if VELOCITY < min(speed_kmph) || VELOCITY > max(speed_kmph)
        error('Velocity is outside defined boundary for interpolated function');
    elseif isnumeric(VELOCITY) ~= 1
        error('invalid argument, see example');
    end
end

%checks if the velocity (vector) is bounded to the data
%checks if input argument is a number
if ismatrix(VELOCITY)
    if sum(VELOCITY < min(speed_kmph)) >= 1 || sum(VELOCITY > max(speed_kmph)) >= 1
        error('Velocity vector is outside defined boundary for interpolated function');
    elseif isnumeric(VELOCITY) ~= 1
        error('invalid argument, see example');
    end
end

% create Piece-wise Polynomial and evaluate at VELOCITY
pp_consumption = spline(speed_kmph,consumption_Whpkm);
CONSUMPTION = ppval(pp_consumption, VELOCITY);

% compute potential mileage with full battery at current consumption rate
MILEAGE = round(55000./CONSUMPTION,1);

% ~ plot sub-routine ~
% if nargin == 1
%     speed = min(speed_kmph):0.01:max(speed_kmph);
%     consumption = spline(speed_kmph,consumption_Whpkm,speed);
%     min_consumption = min(consumption_Whpkm);
%     optimal_speed = speed(min_consumption == consumption);
%     
%     clf
%     subplot(2,1,1)
%     plot(speed,consumption)
%     hold on
%     p2 = plot(optimal_speed,min_consumption,'xr');
%     hold off
%     legend(p2, {['Optimal speed = ', num2str(optimal_speed) 'km/h']})
%     title('Consumption curve for specific velocity')
%     xlabel('Velocity [km/h]')
%     ylabel('Consumption [Wh/km]')
%     subplot(2,1,2)
%     plot(speed,(55000./consumption))
%     title('Travel distance for specific velocity')
%     xlabel('Velocity [km/h]')
%     ylabel('Distance [km]')
%     fprintf('\nThe consumption for %0.2f km/h is %0.2f Wh/km.\nBattery is expected to last %0.1f km given same velocity. \n\n', VELOCITY, CONSUMPTION, MILEAGE); 
% end

end
