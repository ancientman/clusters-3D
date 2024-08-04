function spotDist = distTime(c1c2SpotProp)
timePoints = length(c1c2SpotProp);
nNuc = length(c1c2SpotProp{1}.C1C2SpotVolUM);
minDist = cell(1, nNuc);
for i = 1:nNuc
    minDist{i} = zeros(timePoints, 1);
    for t = 1:timePoints    
%     volTemp = c1c2SpotProp{t}.C1C2SpotVolUM{i};
%     [volTemp, volIdx] = sort(volTemp);
    distTemp = c1c2SpotProp{t}.C1C2SpotDistUM{i};
%     distTemp = distTemp(volIdx);
    if ~isempty(distTemp)
        minDist{i}(t) = min(distTemp);
    else
        minDist{i}(t) = NaN;
    end
    end
% figure('Color', 'w');
% p1 = plot (1:timePoints, minDist{i});
% p1.Color = [0.1 0.1 0.1];
% p1.LineWidth = 1.5;
% p1.LineStyle = '-';
% legend('off');
% ax = gca;
% ax.FontSize = 12;
% ax.LineWidth = 1.5;
% grid ('on');
% xlabel('Time (s)');
% ylabel('Dist. from mRNA center {\mu}m');
% title('Distance of closest bicoid from mRNA');
% hold on;
% x = 0:1:timePoints;
% c=mean(minDist{i});
% const = @(x)(c).*x.^(0);
% p2 = plot(x, const(x));
% p2.Color = [0.7 0.2 0.2];
% p2.LineWidth = 1.5;
% p2.LineStyle = '-.';
% x0 = 100;
% y0= 100;
% plotWidth=300;
% plotHeight=300;
% set(gcf,'position',[x0,y0,plotWidth,plotHeight])
end
spotDist = minDist;
end