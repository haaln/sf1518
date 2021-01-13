function ETA = time_to_destination(DISTANCE, ROUTE, N)
%TIME_TO_DESTINATION Calculate Estimated Time Arrival
%  TIME_TO_DESTIONATION(DISTANCE,ROUTE,N) returns a value for the time (h) it
%  takes to travel in km using composite-trapezoidal integration with N
%  partitions if N is given, if not, it is set to 1.
%
%  the function spans a vector from zero to x with x/n step length and then
%  computes the velocity by calling the VELOCITY(DISTANCE,ROUTE) function 
%  for every index. The integral is then computed using the
%  composite-trapezoidal formula.
%
%  Example:
%       time_to_destination(20,'speed_anna')
%
%       1.2185e+3
%
%  Example:
%       time_to_destinati   on(20,'speed_anna',20)
%
%       1.7539e+03       
%
%
% See also CONSUMPTION, VELOCITY, TOTAL_CONSUMPTION


%checks allowed input and amount of input
    if nargin == 3
        if ~isnumeric(DISTANCE) || isempty(ROUTE) || ~isnumeric(N)
            error('Invalid input')
        end
    elseif nargin == 2 && ~isempty(DISTANCE) && ~isempty(ROUTE)
        N = 1;
    else
        error('Invalid input parameters, number of arguments must be at least two')
    end
    
    % assign step-length, h, bounds, X, and get function value, Y,
    h = DISTANCE/N;
    X = 0:h:DISTANCE;
    Y = 1./(velocity(X,ROUTE));
    ETA = ((h/2)*(Y(1)+Y(end)+2.*sum(Y(2:(end-1))))); % trapezoidal integration
    
%%%%% using 'integral' for verification
%     load (num2str(ROUTE))
%     FUNCTION = @(X) 1./interp1(distance_km,speed_kmph, X, 'spline');
%     REAL_INTEGRAL = integral(FUNCTION, 0, DISTANCE);
%     
%%%%% How many intervals are needed such that the trapezoid method returns
%     % the correct amount of minutes?
%     TOL = 1e-2;
%     Tn = 1;
%     N = 1;
%     while abs(REAL_INTEGRAL - Tn) > TOL
%         N = N +1;
%         h = DISTANCE/N;
%         X = 0:h:DISTANCE;
%         Y = 1./(velocity(X,ROUTE));
%         Tn = ((h/2)*(Y(1)+Y(end)+2.*sum(Y(2:(end-1)))));
%     end
%     iterations = N;
% 
%     
%     
%%%%% How long time for Anna & Elsa to fully complete their routes?
%     load speed_anna.mat
%     ROUTE = 'speed_anna';
%     %using integral
%     FUNC_TIME_A = @(X) 1./interp1(distance_km,speed_kmph, X, 'spline');
%     TIME_ANNA = integral(FUNC_TIME_A, 0, max(distance_km));
%%%%% using trapezoid and not integral
%     h = max(distance_km)/N;
%     X = 0:h:max(distance_km);
%     Y = 1./(velocity(X,ROUTE));
%     TIME_ANNA_TRAPZ = ((h/2)*(Y(1)+Y(end)+2.*sum(Y(2:(end-1)))));
%%%%% comparing accuracy
%     TOL = 1e-2;
%     Tn = 1e9;
%     N = 1;
%     while abs(TIME_ANNA - Tn) > TOL
%         N = N + 1;
%         h = max(distance_km)/N;
%         X = 0:h:max(distance_km);
%         Y = 1./(velocity(X,ROUTE));
%         Tn = ((h/2)*(Y(1)+Y(end)+2.*sum(Y(2:(end-1)))));
%     end
%%%%% number of iterations for anna
%     ANNA_N = N;


%%%%% ELSA
%     load speed_elsa.mat
%     ROUTE = 'speed_elsa';
%%%%% using integral
%     FUNC_TIME_E = @(X) 1./interp1(distance_km,speed_kmph, X, 'spline');
%     TIME_ELSA = integral(FUNC_TIME_E, 0, max(distance_km));
%     %using trapezoid and not integral
%     h = max(distance_km)/N;
%     X = 0:h:max(distance_km);
%     Y = 1./(velocity(X,ROUTE));
%     TIME_ELSA_TRAPZ = ((h/2)*(Y(1)+Y(end)+2.*sum(Y(2:(end-1)))));
%%%%% comparing accuracy
%     TOL = 1e-2;
%     Tn = 1e6;
%     N = 1;
%     while abs(TIME_ELSA - Tn) > TOL
%         N = N + 1;
%         h = max(distance_km)/N;
%         X = 0:h:max(distance_km);
%         Y = 1./(velocity(X,ROUTE));
%         Tn = ((h/2)*(Y(1)+Y(end)+2.*sum(Y(2:(end-1)))));
%     end
%%%%% number of iterations for ELSA
%     ELSA_N = N;
%     
%     
%     
%         
%     
%     
%    fprintf('\nreal_integral = %0.6f and it takes %0.f iterations using Trapezoid composite method to achieve precision %0.e',REAL_INTEGRAL, iterations, TOL)
% %    fprintf('\nTime for ANNA/ELSA = %0.2f and %0.2f', TIME_ANNA, TIME_ELSA)
end
