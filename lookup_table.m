function interpolation_y = lookup_table(data_x,data_y,input_x)

    % data_x = linspace(0,8,21);
    % data_y = sin(data(1,:));
    % x = [4 5 6];
    inputsize = size(input_x);
    interpolation_y = zeros(1,inputsize(2));

    for loop_index = 1:inputsize(2)

        x = input_x(1,loop_index);

        previous_index = find(data_x <= x, 1,'last');

        previous_value_x = data_x(previous_index);
        next_value_x = data_x(previous_index + 1);



        previous_value_y = data_y(previous_index);
        next_value_y = data_y(previous_index + 1);


        % use linear interpolation

        interpolation_y(1,loop_index) = previous_value_y + ...
            (x - previous_value_x)/(next_value_x - previous_value_x) * (next_value_y - previous_value_y);


    %         plot(data(1,:),data(2,:));
    %         hold on;
    % 
    %         plot(previous_value_x,previous_value_y,'marker','*');
    %         plot(next_value_x,next_value_y,'marker','x');
    % 
    %         plot(x,interpolation_y,'marker','o');


    end

end