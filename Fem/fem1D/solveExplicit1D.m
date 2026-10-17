function [T] = solveExplicit1D(C,K,F,model,dirichletNodes,dirichletValues)
% Resuelve la evolucion con Euler explicito.

    timeMass = model.rho*model.cp*C;
    temperature = model.T0(:);
    temperature(dirichletNodes) = dirichletValues;
    T = temperature;

    for iteration = 1 : model.maxit
        previousTemperature = temperature;
        temperature = timeMass\(timeMass*previousTemperature + ...
            model.dt*(F-K*previousTemperature));
        temperature(dirichletNodes) = dirichletValues;
        T = [T temperature];

        if relativeChange(temperature,previousTemperature) < model.tol
            return;
        end
    end
end

function [value] = relativeChange(currentValue,previousValue)
    denominator = norm(currentValue,2);
    if denominator == 0
        value = norm(currentValue-previousValue,2);
    else
        value = norm(currentValue-previousValue,2)/denominator;
    end
end
