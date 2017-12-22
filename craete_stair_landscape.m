function stair_landscape = craete_stair_landscape(x_partition, stair_level, level_height)
    
%     clc; clear variables;
%     x_range = [0, 150];
%     x_partition_diff = 0.1;
%     x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition
%     stair_level = 3;
%     level_height = 1;

    pointer_index = 0;

    inputsize = length(x_partition);
    stair_landscape = zeros(1,inputsize);
    
    partition_grid = floor(inputsize/stair_level);
    
    
    for level_index = 1:stair_level

        for i = 1:partition_grid
            pointer_index = pointer_index + 1;
            stair_landscape(1,pointer_index) = level_height * (level_index - 1);
            
            
        end
    end

    if pointer_index < inputsize
        for p = pointer_index:inputsize
            stair_landscape(1,pointer_index) = level_height * (level_index - 1);
            pointer_index = pointer_index + 1;
        end
    end
%     plot(x_partition,stair_landscape);
  

    


end