function TEMP_BOARD = levnadsregler(BOARD, NEIGHBORS, TEMP_BOARD)

% function to decide whether the first matrix input whether to live or die
% using the second matrix input then putting the result in the third matrix
% input.

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