function [c, ceq] = nonlinear_constraints_mixture(params,x,y)

    % the sum of the weigth must be = 1
    ceq = sum(params(end,:))-1;
    c = [];

end
