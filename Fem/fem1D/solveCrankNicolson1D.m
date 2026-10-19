function [T] = solveCrankNicolson1D(C,K,F,model,dirichletNodes,dirichletValues)
% Resuelve la evolucion con el esquema de Crank-Nicolson.

    timeMass = model.rho*model.cp*C;
    halfTimeStep = model.dt/2;
    leftMatrix = timeMass + halfTimeStep*K;
    rightMatrix = timeMass - halfTimeStep*K;
    temperature = model.T0(:);
    temperature(dirichletNodes) = dirichletValues;
    T = temperature;

    for iteration = 1 : model.maxit
        previousTemperature = temperature;
        temperature = leftMatrix\(rightMatrix*previousTemperature + ...
            model.dt*F);
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