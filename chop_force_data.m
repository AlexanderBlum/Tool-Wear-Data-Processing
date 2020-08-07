function chop_force_data(force_data)

Fs = 10000;
dt = 1/Fs;

%% rough bracketing and chopping

% select each subset from left to right:
% click one and click two bracket plunge 1
% click two and click three bracket facing 1
% click three and click four bracket plunge 2
% repeat thru 12 total clicks for six plunges and five facings

% figure
hold on
plot(force_data{1});
plot(force_data{2});
plot(force_data{3});

[x, ~] = ginput(12);

% close(gcf);

x = floor(x);

plunges = cell(6,1);
for ii = 1:6
    % [Fx, Fy, Fz, t]
    plunges{ii} = nan(x(ii*2)   - x(ii*2-1) + 1, 4);
    % Fx, Fy, Fz
    for jj = 1:3
        plunges{ii}(:,jj) = force_data{jj}(x(ii*2-1):x(ii*2));
    end
    plunges{ii}(:,4) = linspace(0, dt*length(plunges{ii}(:,jj)), length(plunges{ii}(:,jj)));    

end
facings = cell(5,1);
for ii = 1:5
    facings{ii} = nan(x(ii*2+1) - x(ii*2) + 1 , 4);
    for jj = 1:3
        facings{ii}(:,jj) = force_data{jj}(x(ii*2):x(ii*2+1));
    end
    facings{ii}(:,4) = linspace(0, dt*length(facings{ii}(:,jj)), length(facings{ii}(:,jj)));    
end

%% fine bracketing and chopping
figure
for ii = 1:6
    hold on
    for jj = 1:3
        plot(plunges{ii}(:, jj))
    end
    hold off
    [x, ~]  = ginput(2);
    x = floor(x);
    for jj = 1:3
        plunges{ii} = plunges{ii}(x(1):x(2), jj);
    end
    plunges{ii}(:,4) = linspace(0, dt*length(plunges{ii}(:,jj)), length(plunges{ii}(:,jj)));     
end

for ii = 1:5
    hold on
    for jj = 1:3
        plot(facings{ii}(:, jj))
    end
    hold off
    [x, ~]  = ginput(2);
    x = floor(x);
    for jj = 1:3
        facings{ii} = facings{ii}(x(1):x(2), jj);
    end
    facings{ii}(:,4) = linspace(0, dt*length(facings{ii}(:,jj)), length(facings{ii}(:,jj)));   
end
end