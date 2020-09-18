clear,clf, clc

n = input('Enter size of map: ');
BOARD_X = zeros(1,n);
BOARD_Y = zeros(1,n);
BOARD = meshgrid(BOARD_X, BOARD_Y);
TEMP_BOARD = meshgrid(BOARD_X, BOARD_Y);
NEIGHBORS = meshgrid(BOARD_X, BOARD_Y);
ROUND = 0;
BOARD = TEMP_BOARD;

% %%%%%%%%%%%%%%%%%%
% BOARD(5,5) = 1;
% BOARD(5,4) = 1;
% BOARD(4,4) = 1;
% BOARD(6,6) = 1;
% %%%%%%%%%%%%%%%%%%
% BOARD(5:9,5) = 1;
% BOARD(5,7:11) = 1;
% BOARD(11,5:9) = 1;
% BOARD(7:11,11) = 1;
%%%%%%%%%%%%%%%%%%%%
% BOARD(1:3,1) = 1;
% BOARD(1,1:3) = 1;
% BOARD(8:10,10) =1;
% BOARD(10,8:10) =1;
%%%%%%%%%%%%%%%%%%%


spy(BOARD,70,'k')
title(['CURRENT ROUND: ', num2str(ROUND)])
pause(1)


while(ROUND < 30)
    
    NEIGHBORS = antalgrannar(BOARD, NEIGHBORS);
    TEMP_BOARD = levnadsregler(BOARD, NEIGHBORS, TEMP_BOARD);
    
    ROUND = ROUND + 1;
    BOARD = TEMP_BOARD;
    TEMP_BOARD = meshgrid(BOARD_X, BOARD_Y);
    spy(BOARD,70,'k')
    title(['CURRENT ROUND: ', num2str(ROUND)])
    pause(0.1)
    
    if sum(sum(BOARD)) == 0
        disp(['all pieces anhiliated at at round ', num2str(ROUND)])
        break
    end
    
end
if ROUND == 30
    disp('system achieved steady-state equilibrium, probably')
end

function NEIGHBORS = antalgrannar(BOARD, NEIGHBORS)
    for i = 1:length(BOARD)
        for j = 1:length(BOARD)

            if i == 1 && j == 1
                NEIGHBORS(i,j) = sum(sum(BOARD(i:(i+1),j:(j+1))));

            elseif i > 1 && i < length(BOARD) && j == 1
                NEIGHBORS(i,j) = sum(sum(BOARD((i-1):(i+1),j:(j+1))));

            elseif i == 1 && j > 1 && j < length(BOARD)
                NEIGHBORS(i,j) = sum(sum(BOARD((i):(i+1),(j-1):(j+1))));

            elseif i > 1 && i < length(BOARD) && j > 1 && j < length(BOARD)
                NEIGHBORS(i,j) = sum(sum(BOARD((i-1):(i+1),(j-1):(j+1))));

            elseif i == length(BOARD) && j == 1
                NEIGHBORS(i,j) = sum(sum(BOARD((i-1):i,j:(j+1))));

            elseif i == length(BOARD) && j > 1 && j < length(BOARD)
                NEIGHBORS(i,j) = sum(sum(BOARD((i-1):i,(j-1):(j+1))));

            elseif i == 1 && j == length(BOARD)
                NEIGHBORS(i,j) = sum(sum(BOARD(i:(i+1),(j-1):j)));

            elseif i == length(BOARD) && j == length(BOARD)
                NEIGHBORS(i,j) = sum(sum(BOARD((i-1):i,(j-1):j)));
                
            elseif i > 1 && i < length(BOARD) && j == length(BOARD)
                NEIGHBORS(i,j) = sum(sum(BOARD((i-1):(i+1),(j-1):j)));

            end
        end
    end
end

function TEMP_BOARD = levnadsregler(BOARD, NEIGHBORS, TEMP_BOARD)

    for i = 1:length(BOARD)
        for j = 1:length(BOARD)

            if NEIGHBORS(i,j) < 3 && BOARD(i,j) == 1
                TEMP_BOARD(i,j) = 0;
                
            elseif NEIGHBORS(i,j) == 3
                TEMP_BOARD(i,j) = 1;
                
            elseif NEIGHBORS(i,j) == 4 && BOARD(i,j) == 1
                TEMP_BOARD(i,j) = 1;
                
            elseif NEIGHBORS(i,j) > 4
                TEMP_BOARD(i,j) = 0;
                
            end
        end
    end
end