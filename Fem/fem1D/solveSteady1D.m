function [T] = solveSteady1D(K,F)
% Resuelve el problema estacionario K*T = F.

    T = K\F;
end
