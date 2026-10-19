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
    numberOfStates = size(T,2);
    frameStep = max(1,ceil(numberOfStates/120));

    plotHandle = plot(xnode,T(:,1),'o-','LineWidth',1.5);
    xlabel('x');
    ylabel('Temperatura');
    grid on;
    xlim([min(xnode), max(xnode)]);
    ylim([temperatureMinimum, temperatureMaximum]);

    for state = 1 : frameStep : numberOfStates
        set(plotHandle,'YData',T(:,state));
        title(sprintf('Estado %d - t = %g',state-1,(state-1)*dt));
        drawnow;
        pause(0.005);
    end

    if state ~= numberOfStates
        set(plotHandle,'YData',T(:,numberOfStates));
        title(sprintf('Estado %d - t = %g',numberOfStates-1,...
            (numberOfStates-1)*dt));
        drawnow;
    end
end
