function NEIGHBORS = antalgrannar(BOARD, NEIGHBORS)

% function to calculate number of neighbors around each index within a
% matrix. it uses a primitive algorithm by categorization of indeces within
% a matrix using a bunch of if-statements then proceeding to count the neighbors around it, including itself.

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