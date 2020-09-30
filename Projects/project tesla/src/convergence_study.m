clear,clc
format long

% Matlab 'integral' to cross-check error convergence:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROUTE = 'speed_anna';
load (num2str(ROUTE))
DISTANCE = max(distance_km);
xq = min(speed_kmph):0.001:max(speed_kmph);
fun = @(xq) 1./interp1(distance_km,speed_kmph, xq, 'spline');
REAL_INTEGRAL = integral(fun, 0, DISTANCE);




% numerical convergence analysis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 100000;
TOL = 1e-8;
ROUND = 1;
BREAK_CONDITION = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    while(true)
        % assign iteration parameters
        h = DISTANCE/N;
        X = 0:h:DISTANCE;
        Y = 1./(velocity(X,ROUTE));
        
        TRAPEZOID_APPROXIMATION(ROUND) = ((h/2)*(Y(1)+Y(end)+2.*sum(Y(2:(end-1))))); % calculate approximation using trapezoid rule
        REAL_ERROR(ROUND) = abs(REAL_INTEGRAL - TRAPEZOID_APPROXIMATION(ROUND));
        STEP_LENGTH(ROUND) = h;
        ROUND_V(ROUND) = ROUND;
        N = 2*N;
        
        if ROUND > BREAK_CONDITION
            break
        elseif REAL_ERROR(ROUND) < TOL
            break
        end
        
        ROUND = ROUND + 1;
    end
    
    % |e(h)-e(h/2)|
    for i = 1:(ROUND-1)
        APPROXIMATIVE_ERROR(i) = (abs(TRAPEZOID_APPROXIMATION(i) - TRAPEZOID_APPROXIMATION(i+1))); 
        h_APPROXIMATIVE(i) = STEP_LENGTH(i);
    end
    
    % log2(|e(h)-e(h/2)| / |e(h/2)-e(h/4)|)= log(2^p) = p + O(h²)
    for i = 1:(ROUND-2)
        NUMERICAL_CONVERGENCE(i) = log2(APPROXIMATIVE_ERROR(i)/APPROXIMATIVE_ERROR(i+1));
        h_CONVERGENCE(i) = STEP_LENGTH(i);
    end
    
    
    subplot(2,1,1) %plot order of convergence with 'help-lines'
    loglog(h_CONVERGENCE,h_CONVERGENCE.^NUMERICAL_CONVERGENCE,h_CONVERGENCE,h_CONVERGENCE.^2,'-.',h_CONVERGENCE,h_CONVERGENCE.^1,'--')
    legend('Experimental P', 'P=2','P=1','Location','best')
    xlabel('h')
    title('Order of convergence P')
    grid on
    
    subplot(2,1,2) % plot error as a function of step-length
    loglog(h_APPROXIMATIVE,APPROXIMATIVE_ERROR,h_APPROXIMATIVE, h_APPROXIMATIVE.^2,'-.',h_APPROXIMATIVE,h_APPROXIMATIVE,'--')
    grid on
    legend('Error(h)','O(h^2)','O(h)','Location','best')
    xlabel('h')
    title('Error difference vs step-length, h')
    
    ERROR_TABLE = [ROUND_V' STEP_LENGTH' TRAPEZOID_APPROXIMATION' REAL_ERROR']
    NUMERICAL_CONVERGENCE = [NUMERICAL_CONVERGENCE']
    
    
