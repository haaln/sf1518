% changes that need to be done:
%
% draw map function from gui
% automatic exit if game reaches steady-state solution
% separate everything into functions and objects
% debugging
%

clear, clc
n = input('Enter size of map: ');
BOARD_X = zeros(1,n);
BOARD_Y = zeros(1,n);
BOARD = meshgrid(BOARD_X, BOARD_Y);
TEMP_BOARD = meshgrid(BOARD_X, BOARD_Y);
ROUND = 0;
NEIGHBORS = 1337;
BOARD = TEMP_BOARD;

%%%%%%%%%%%%%%%%%%%%%%%%
% 
% BOARD(1) = 1;
% BOARD(2,1) = 1;
% BOARD(1,2) = 1;
% 
% BOARD(5,5) = 1;
% BOARD(6,5) = 1;
% BOARD(5,6) = 1;
% 
% BOARD(1,13) = 1;
% BOARD(1,14) = 1;
% BOARD(2,14) = 1;
% BOARD(3,14) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%



BOARD(5,5) = 1;
BOARD(5,4) = 1;
BOARD(4,4) = 1;
BOARD(6,6) = 1;


spy(BOARD,70,'k')
title(['CURRENT ROUND: ', num2str(ROUND)])
pause(2.5)
while(true)
    for i = 1:length(BOARD_Y)
        for j = 1:length(BOARD_X)
            if i == 1 && j == 1
                NEIGHBORS = sum(sum(BOARD(i:(i+1),j:(j+1))));
            elseif i > 1 && i < length(BOARD_Y) && j == 1
                NEIGHBORS = sum(sum(BOARD((i-1):(i+1),j:(j+1)))); 
            elseif i == 1 && j > 1 && j < length(BOARD_X)
                NEIGHBORS = sum(sum(BOARD((i):(i+1),(j-1):(j+1)))); 
            elseif i > 1 && i < length(BOARD_Y) && j > 1 && j < length(BOARD_X)
                NEIGHBORS = sum(sum(BOARD((i-1):(i+1),(j-1):(j+1))));
            elseif i == length(BOARD_Y) && j == 1
                NEIGHBORS = sum(sum(BOARD((i-1):i,j:(j+1))));
            elseif i == length(BOARD_Y) && j > 1 && j < length(BOARD_X)
                NEIGHBORS = sum(sum(BOARD((i-1):i,(j-1):(j+1))));
            elseif i == 1 && j == length(BOARD_X)
                NEIGHBORS = sum(sum(BOARD(i,(j-1):j)));
            elseif i > 1 && i < length(BOARD_Y) && j == length(BOARD_X)
                NEIGHBORS = sum(sum(BOARD((i-1):(i+1),(j-1):j))); 
            end        

            % rules of engagement
            if NEIGHBORS < 3 && BOARD(i,j) == 1
                TEMP_BOARD(i,j) = 0;
            end
            if NEIGHBORS == 3 
                TEMP_BOARD(i,j) = 1;
            end

            if NEIGHBORS > 4
                TEMP_BOARD(i,j) = 0;
            end        
            % rules of engagement end

        end
    end
ROUND = ROUND + 1;
BOARD = TEMP_BOARD;
TEMP_BOARD = meshgrid(BOARD_X, BOARD_Y);

% spy(BOARD,70,'k')
title(['CURRENT ROUND: ', num2str(ROUND)])
pause(2.5)
end
