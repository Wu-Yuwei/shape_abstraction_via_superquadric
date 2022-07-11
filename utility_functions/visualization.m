function visualization(Z, K, theta, point, pointColor, meshColor, fig1, fig2, viewAngle)

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