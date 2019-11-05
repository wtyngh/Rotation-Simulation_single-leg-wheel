function stair_landscape = craete_stair_landscape(x_partition, stair_level, level_height)
    % testing
%     clc; clear variables;
%     x_range = [0, 10];
%     x_partition_diff = 0.1;
%     x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition
%     stair_level = 3;
%     level_height = 1;
    %
    
    pointer_index = 1;

    inputsize = length(x_partition);
    stair_landscape = zeros(1,inputsize);
    
    partition_grid = floor(inputsize/stair_level);
    
    fillet_region = 10;
    for level_index = 1:stair_level
        for i = 1:partition_grid
            if i < fillet_region % first few points
                if level_index == 1
                    stair_landscape(1,pointer_index) = level_height * (level_index - 1);
                else
                    stair_landscape(1,pointer_index) = level_height * (level_index - 1)...
                                - (fillet_region - i)/fillet_region*level_height*0.5;
                end
            elseif i > (partition_grid - fillet_region)
                stair_landscape(1,pointer_index) = level_height * (level_index - 1)...
                    + (i - (partition_grid - fillet_region))/fillet_region*level_height*0.5;
            else
                
                stair_landscape(1,pointer_index) = level_height * (level_index - 1);
            end
            pointer_index = pointer_index + 1;
        end
    end

    % fill the last few points
    if pointer_index <= inputsize
        for p = pointer_index:inputsize
            stair_landscape(1,pointer_index) = level_height * (level_index - 1);
            pointer_index = pointer_index + 1;
        end
    end
    
    % testing
%     plot(x_partition,stair_landscape);
  

    

% 
end