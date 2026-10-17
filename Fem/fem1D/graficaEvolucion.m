function graficaEvolucion(xnode,T,dt)
% Muestra la evolucion temporal de la temperatura 1D.
%
% xnode: coordenadas nodales.
% T: matriz devuelta por elemFinitos.
% dt: paso temporal. Para una solucion estacionaria puede omitirse.

    xnode = xnode(:);
    if nargin < 3
        dt = 1;
    end

    figure('Name','Evolucion de temperatura 1D');
    temperatureMinimum = min(T(:));
    temperatureMaximum = max(T(:));

    for state = 1 : size(T,2)
        plot(xnode,T(:,state),'o-','LineWidth',1.5);
        xlabel('x');
        ylabel('Temperatura');
        title(sprintf('Estado %d - t = %g',state-1,(state-1)*dt));
        grid on;
        xlim([min(xnode), max(xnode)]);
        ylim([temperatureMinimum, temperatureMaximum]);
        drawnow;
        pause(0.03);
    end
end
