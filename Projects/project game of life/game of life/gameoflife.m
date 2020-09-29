%% disclaimer about the fidelity of this script vs project requirements.

% you might notice that the function 'antalgrannar' has 'NEIGHBORS' as an
% input and an output parameter, this is to avoid reallocating memory each
% time it enters the nested for-loop, but mainly for matlab to shut up
% about the errors.

% both functions furthermore don't use a specific set of coordinates as an input
% and rely solely on the BOARD and NEIGHBOR variables mainly, this design
% choice is of no other reason than its intuitive ease-of-use.

% If, for some reason you would want to change this you would need to remove
% the nested loops from both of the functions and paste them just before
% the functions and assign the i and j variables as input parameters
% instead of BOARD, that would yield the desired outcome for both
% functions, but solely for that coordinate. Thus you would have to put in those
% values in a matrix outside of the function in the main script to have a
% complete overview of the next game state.

% the antalgrannar function also counts ALL the neighbors, including
% itself, but this issue is easily overcome by changing the levnadsregler.m
% i.e cell dies if >4 rather than >3

%%
clear, clc

try n = input('Enter size of map: ');
    while n < 2 || round(n) ~= n && isinteger == 0
        warning('Invalid value, only natural numbers larger than one');
        n = input('Please choose another value: ');
    end
catch
    warning('Invalid value');
    warning('Exiting script')
    return
end


BOARD_X = zeros(1,n);
BOARD_Y = zeros(1,n);
BOARD = meshgrid(BOARD_X, BOARD_Y);
TEMP_BOARD = meshgrid(BOARD_X, BOARD_Y);
NEIGHBORS = meshgrid(BOARD_X, BOARD_Y);
ROUND = 0;

%%%%% n > 2 %%%%%%%
% BOARD(1,2) = 1;
% BOARD(1,1) = 1;
% BOARD(2,2) = 1;

% %%% n ~ 9  %%%%%%
% BOARD(5,5) = 1;
% BOARD(5,4) = 1;
% BOARD(4,4) = 1;
% BOARD(6,6) = 1;

% %%%% n ~ 15 %%%%%%%
% BOARD(5:9,5) = 1;
% BOARD(5,7:11) = 1;
% BOARD(11,5:9) = 1;
% BOARD(7:11,11) = 1;

%%%%% n == 10 %%%%%%%
% BOARD(1:3,1) = 1;
% BOARD(1,1:3) = 1;
% BOARD(8:10,10) =1;
% BOARD(10,8:10) =1;



spy(BOARD,50,'k')
title(['CURRENT ROUND: ', num2str(ROUND)])
pause(1)

while(ROUND < 50)
    
    NEIGHBORS = antalgrannar(BOARD, NEIGHBORS);
    TEMP_BOARD = levnadsregler(BOARD, NEIGHBORS, TEMP_BOARD);
    
    ROUND = ROUND + 1;
    BOARD = TEMP_BOARD;
    TEMP_BOARD = meshgrid(BOARD_X, BOARD_Y);
    spy(BOARD,50,'k')
    title(['CURRENT ROUND: ', num2str(ROUND)])
    pause(0.1)
    
    if sum(sum(BOARD)) == 0
        disp(['all pieces annihiliated at at round ', num2str(ROUND)])
        break
    end
    
end
if ROUND == 50
    disp('system achieved steady-state equilibrium, probably')
end