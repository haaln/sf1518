function EXPECTED_DISTANCE = distance(T, ROUTE)
%DISTANCE Calculate the expected distance travelled for a given time
%  DISTANCE(T,ROUTE) returns a value by finding the root of the polynomial
%  constructed from the Time (in hours) subtracted by TIME_TO_DESTINATION using
%  Newtons iteration formula.
%
%  The function will also accept arguments in matrix form and in turn
%  return another matrix as an output.
%
%  The function will also warn you if the Time exceeds the time needed to
%  complete the chosen route.
%
%  Example:
%       distance(0.5, 'speed_anna')
%
%       51.0755
%
%  Example:
%       distance([0.5 1],'speed_elsa')
%
%       33.5967   67.1934
%
%  For some reason using the function with matrices returns different
%  answers, it only occurs if the function extrapolates. It gives the right
%  answer for the extrapolated answer but the wrong for the interpolated
%  answer.
%
%  See also REACH

%check input arguments and load ROUTE
if nargin ~= 2 || ~isnumeric(T) || nargin > 2
    error('Usage: distance distance(T,route)')
end
try load (num2str(ROUTE))
catch
    warning('specified route data not found')
    return
end
%construct polynomial using anonymous functions.
g = @(DISTANCE) time_to_destination(DISTANCE,ROUTE,1e4) - T;
gprime = @(DISTANCE) 1./velocity(DISTANCE,ROUTE);

%calculate a reasonable inital guess
TOT_T = time_to_destination(max(distance_km),ROUTE,1e3);
ratio = T ./ TOT_T;
DISTANCE = max(distance_km) .* ratio;

i = 0;
TOL = 0.5 * 1e-1; %tolerance
ERROR = 1;

%the try statement is to catch warnings if the function begins to extrapolate
try
    warning('off')
    while abs(ERROR) > TOL
        DISTANCE_2 = DISTANCE - g(DISTANCE)./gprime(DISTANCE); %Newton's Method
        ERROR = DISTANCE_2 - DISTANCE; %the difference between old value and new value
        DISTANCE = DISTANCE_2;
        i = i + 1;
    end
catch ME %error handling
    switch ME.identifier
        case 'MATLAB:unassignedOutputs' % if the function has begun to extrapolate
            warning('on')
            if ~ismatrix(DISTANCE) % if DISTANCE is integer
                warning('The Time exceeds the amount needed to complete the route %s, extrapolated results may vary.\nYou are projected to complete the chosen route and be able to go for another %0.1f km', ROUTE, DISTANCE-max(distance_km))
            else % if DISTANCE is matrix
                warning('The Time exceeds the amount needed to complete the route %s, extrapolated results may vary.\nYou are projected to complete the chosen route and go even further.', ROUTE)
            end
    end
    warning('on')
end

EXPECTED_DISTANCE = DISTANCE;

end
