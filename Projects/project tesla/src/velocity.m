function VELOCITY = velocity(DISTANCE,ROUTE)
%VELOCITY Compute the velocity at a given distance
%  VELOCITY(X,ROUTE) returns a value by evaluating a piecewise polynomial
%  at the desired distance obtained using cubic splines on an 
%  imported data set inside the working directory.
%
%  The function will also accept matrices which enables for computation of
%  several consumption values of their corresponding velocity
%
%  Example:
%        velocity(10,'speed_anna') returns velocity at 10km in route 'speed_anna'
%
%        92.3228
%
%  Example:
%        velocity(24.5,'speed_elsa') returns velocity at 24.5km for route 'speed_elsa'
%
%        94.8067
%
%  See also CONSUMPTION, TIME_TO_DESTINATION, TOTAL_CONSUMPTION

%  Hassan Al N 5-9-20

%checks for allowed input and correct amount of input
if nargin ~= 2 || ~isnumeric(DISTANCE) || nargin > 2
    error('Usage: velocity = velocity(distance,route)')
end

%tries to load ROUTE and returns if it cant find anything
try load (num2str(ROUTE))
catch
    warning('specified route data not found')
    return
end

%checks if the velocity is bounded to the data
if  ~ismatrix(DISTANCE)
    if DISTANCE < min(distance_km) || DISTANCE > max(distance_km)
        warning('distance is outside of defined interpolated boundary');
        return
    end
end

%checks if the velocity (vector) is bounded to the data
if ismatrix(DISTANCE)
    if sum(DISTANCE < min(distance_km)) >= 1 || sum(DISTANCE > max(distance_km)) >= 1
        warning('distance vector is outside defined boundary for interpolated function');
        return
    end
end

% create Piece-wise Polynomial using imported data and evaluate at DISTANCE
pp_speed = spline(distance_km, speed_kmph);
VELOCITY = ppval(pp_speed,DISTANCE);

%  % ~ plot subroutine ~
% if nargin == 2
%     distance = min(distance_km):0.01:max(distance_km);
%     speed = spline(distance_km, speed_kmph, distance);
%     clf
%     plot(distance,speed)
%     hold on
%     p2 = plot(DISTANCE,VELOCITY, 'xr','LineWidth',5, 'MarkerSize', 10);
%     hold off
%     legend(p2, 'current location', 'Location', 'best')
%     xlabel('Distance [km]')
%     ylabel('Velocity [km/h]')
%     split = strsplit(num2str(ROUTE),'_');
%     title(['Velocity for route' split(2) ])
%     fprintf('\nThe velocity at %0.1f km is %0.2f km/h.\n\n', DISTANCE, VELOCITY);
% end

end