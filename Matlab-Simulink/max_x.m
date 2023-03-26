function [x1_max,x2_max,x3_max,x4_max] = max_x(uMR_values)

    for i = 1 : length(uMR_values)

        uMR = uMR_values(i);

        x = extract_x(uMR);

        if i == 1
            x1_max = x(1);
            x2_max = x(2);
            x3_max = x(3);
            x4_max = x(4);
        else
            if x(1) > x1_max
                x1_max = x(1);
            end
            if x(2) > x2_max
                x2_max = x(2);
            end
            if x(3) > x3_max
                x3_max = x(3);
            end
            if x(4) > x4_max
                x4_max = x(4);
            end    
        end
    end
    
end

