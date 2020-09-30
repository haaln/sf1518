function E = total_consumption(DISTANCE, ROUTE, N)
%TOTAL_CONSUMPTION Calculate accumulated power consumption for Tesla Roadster
%  TOTAL_CONSUMPTION(DISTANCE, ROUTE, N) returns the total consumption over a
%  distance by calling velocity.m and inserting the returned value
%  into consumption.m where thereafter it integrates the compounded result
%  using trapezoidal composite method with n intervals. 
%
%  If the third argument N is not specified, it is automatically set to 1.
%  
%  Example:
%       total_consumption(20,'speed_anna')
%
%       3.6338e+03
%
%  Example:
%       total_consumption(20,'speed_anna',20)
%
%       3.6338e+03
%
%  See also CONSUMPTION, VELOCITY, TIME_TO_DESTINATION, 

%  Hassan Al N 5-9-20

%checks inputs and returns if it's not an integer
    if nargin == 3
        if ~isnumeric(DISTANCE) || isempty(ROUTE) || ~isnumeric(N)
            error('Invalid input')
        end
    elseif nargin == 2 && ~isempty(DISTANCE) && ~isempty(ROUTE)
        N = 1;
    else
        error('Invalid input parameters, number of arguments must be at least two')
    end
    
    % assign step-length, h, bounds, X, and get function values VELOCITY and CONSUMPTION,
    h = DISTANCE/N;    
    X = 0:h:DISTANCE;
    VELOCITY = velocity(X,ROUTE);
    CONSUMPTION = consumption(VELOCITY);
    E = (h/2)*(CONSUMPTION(1)+CONSUMPTION(end)+2.*sum(CONSUMPTION(2:(end-1)))); %trapezoidal integration
% 
%%%%% using 'integral' for verification
%     load (num2str(ROUTE)) 
%     pp_speed = spline(distance_km, speed_kmph); 
%     FUNC_VELOCITY = @(x) ppval(pp_speed,x); 
%     load roadster.mat
%     pp_consumption = spline(speed_kmph,consumption_Whpkm);
%     FUNC_CONSUMPTION =@(x) ppval(pp_consumption, FUNC_VELOCITY(x));
%     REAL_INTEGRAL = integral(FUNC_CONSUMPTION,0,DISTANCE);
%     fprintf('The actual consumption value over the distance %0.f is %0.3f Wh.\n',DISTANCE,REAL_INTEGRAL)
%         
%%%%% How many intervals are needed to get the correct amount of Wh?    
%     load (num2str(ROUTE)) 
%     TOL = 1e-3;
%     Tn = 1;
%     N = 1;
%     while abs(REAL_INTEGRAL - Tn) > TOL
%         N = N + 1;
%         h = DISTANCE/N;
%         X = 0:h:DISTANCE;
%         VELOCITY = velocity(X,ROUTE);
%         CONSUMPTION = consumption(VELOCITY);
%         Tn = ((h/2)*(CONSUMPTION(1)+CONSUMPTION(end)+2.*sum(CONSUMPTION(2:(end-1)))));
%     end
%     N_ACC = N;
%     
% 
%%%%% How much energy do Anna & Elsa consume after the completion of their route?
%     load speed_anna
%     ROUTE = 'speed_anna';
%     pp_speed = spline(distance_km, speed_kmph);
%     FUNC_VELOCITY = @(x) ppval(pp_speed,x);
%     FUNC_CONSUMPTION =@(x) ppval(pp_consumption, FUNC_VELOCITY(x));
%     POWER_ANNA = integral(FUNC_CONSUMPTION,0,max(distance_km));
%     fprintf('Anna consumes a total of %0.f Wh at completion of speed_anna.\n',POWER_ANNA)
%     
%%%%% calculates the same thing using trapezoid and not integral
%     h = max(distance_km)/N;
%     X = 0:h:max(distance_km);
%     VELOCITY = velocity(X,ROUTE);
%     CONSUMPTION = consumption(VELOCITY);
%     POWER_ANNA_TRAPZ = (h/2)*(CONSUMPTION(1)+CONSUMPTION(end)+2.*sum(CONSUMPTION(2:(end-1)))); 
% %   fprintf('Anna consumes a total of %0.f Wh at completion of speed_anna for %0.f partitions.\n',POWER_ANNA_TRAPZ, N)
%     %how many iterations needed for accurate answer?
%     TOL = 1e-4;
%     Tn = 1;
%     N = 1;
%     while abs(POWER_ANNA- Tn) > TOL
%         N = N + 1;
%         h = max(distance_km)/N;
%         X = 0:h:max(distance_km);
%         VELOCITY = velocity(X,ROUTE);
%         CONSUMPTION = consumption(VELOCITY);
%         Tn = ((h/2)*(CONSUMPTION(1)+CONSUMPTION(end)+2.*sum(CONSUMPTION(2:(end-1)))));
%     end
%     ANNA_N = N;
%     
%%%%% ELSA
%     load speed_elsa.mat
%     ROUTE = 'speed_elsa';
%     pp_speed = spline(distance_km, speed_kmph);
%     FUNC_VELOCITY = @(x) ppval(pp_speed,x);
%     FUNC_CONSUMPTION =@(x) ppval(pp_consumption, FUNC_VELOCITY(x));
%     POWER_ELSA = integral(FUNC_CONSUMPTION,0,max(distance_km));
%     fprintf('Elsa consumes a total of %0.f Wh at completion of speed_elsa. \n',POWER_ELSA)
%     
%%%%% calculates the same thing using trapezoid and not integral
%     h = max(distance_km)/N;
%     X = 0:h:max(distance_km);
%     VELOCITY = velocity(X,ROUTE);
%     CONSUMPTION = consumption(VELOCITY);
%     POWER_ELSA_TRAPZ = (h/2)*(CONSUMPTION(1)+CONSUMPTION(end)+2.*sum(CONSUMPTION(2:(end-1))));
% %     fprintf('Elsa consumes a total of %0.f Wh at completion of speed_elsa for %0.f partitions. \n',POWER_ELSA_TRAPZ, N)
%%%%% how many iterations needed for accurate answer?
%     TOL = 1e-4;
%     Tn = 1;
%     N = 1;
%     while abs(POWER_ELSA - Tn) > TOL
%         N = N + 1;
%         h = max(distance_km)/N;
%         X = 0:h:max(distance_km);
%         VELOCITY = velocity(X,ROUTE);
%         CONSUMPTION = consumption(VELOCITY);
%         Tn = ((h/2)*(CONSUMPTION(1)+CONSUMPTION(end)+2.*sum(CONSUMPTION(2:(end-1)))));
%     end
%     ELSA_N = N;
%     
%     fprintf('The number of intervals needed are %.f to achieve the precision %0.e\n', N_ACC, TOL)
end
