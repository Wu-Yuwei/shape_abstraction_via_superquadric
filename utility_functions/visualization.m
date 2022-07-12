function visualization(Z, K, theta, point, fig1, fig2, viewAngle)

meshColor = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840];[1.0000, 0.0745, 0.6510]];
pointColor = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];...
    [0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];[0.6350, 0.0780, 0.1840];[1.0000, 0.0745, 0.6510]];

xLim = [min(point(1,:)) max(point(1,:))];
yLim = [min(point(2,:)) max(point(2,:))];
zLim = [min(point(3,:)) max(point(3,:))];

colorSize = size(pointColor,1);

f1 = figure(fig1);
clf(f1)
hold on
for ii = 1:K
    c = mod(ii,colorSize);
    c(c==0) = colorSize;
    showPoints(point(:,Z(:,ii)==1),'MarkerSize',5,'Color',pointColor(c,:))
end
hold off
view(viewAngle)
xlim(xLim*1.5)
ylim(yLim*1.5)
zlim(zLim*1.5)

f2 = figure(fig2);
clf(f2)
hold on
for ii = 1:K
    c = mod(ii,colorSize);
    c(c==0) = colorSize;
    showSuperquadric_mesh(theta(ii,:),'Color', meshColor(c,:),'Taper',true)
end
hold off
view(viewAngle)
xlim(xLim*1.5)
ylim(yLim*1.5)
zlim(zLim*1.5)