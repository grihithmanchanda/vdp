% inputs:
%   currents: the list of current measurements in all given directions,
%       specified by the order below
%   voltages: the list of voltage measurements specified in the same order
%   prec: how precise must we be to end the iteration
%       (error from y=1 < 10^-prec)
%   start_val: an initial value to start our iterative process from. Doesn't
%       matter too much as long as it's reasonable and above zero
%   t: thickness of material to determine resistivity, in meters
function res = vdp(currents, voltages, t, prec, start_val)
    res_list = setup(currents, voltages);
    [res, vals, err] = prop(res_list, prec, start_val);
    res = t * res;
    vals = t * vals;
    err = t * err;
    fprintf('The resistivity of the material is %4.3f \x03a9 m\n', res);
    plot(vals);
    hold on;
    plot(err);
end

% caculate all the correct resistances
function res_list = setup(currents, voltages)
    % order of arg list:
    % I_12
    % I_21
    % I_34
    % I_43
    % I_23
    % I_32
    % I_41
    % I_14
    % 
    % Same order goes for voltages. Thus, opposites are adjacent to each other
    
    % get 8 individual resistances, two columns end up to be vertical and
    % horizontal
    resistances = zeros(length(currents) / 2, 2);
    
    for i = 1:length(currents)
        offset = 2;
        if mod(floor((i - 1) / 2), 2) == 1
            offset = -2;
        end
        resistances(i) = voltages(i + offset) / currents(i);
    end

    % get the horizontal and vertical resistances
    res_list = sum(resistances, 1) / (length(currents) / 2);
end

function [res, vals, err] = prop(res_list, prec, start_val)
    % use linear approximation about current point to determine x change to
    % reach the ideal value
    res = start_val;
    vals = [res];
    y_dist = 1 - fn(res_list, res);
    err = [y_dist];
    while abs(y_dist) > (10^(-prec))
        new_res = res + (y_dist / Dfn(res_list, res));
        % if linear approximation goes negative, half original value and try
        % again
        if new_res <= 0
            res = res / 2;
        else
            res = new_res;
        end
        % keep track of all the values and errors (y distances)
        vals = [vals res];
        y_dist = (1 - fn(res_list, res));
        err = [err y_dist];
    end
end

% function to solve
function opt = fn(res_list, R)
    opt = exp(-pi * (res_list(1) / R)) + exp(-pi * (res_list(2) / R));
end

% derivative of function
function opt = Dfn(res_list, R)
    opt = (pi / (R^2)) * ((res_list(1) * exp(-pi * (res_list(1) / R))) ...
        + (res_list(2) * exp(-pi * (res_list(2) / R))));
end
