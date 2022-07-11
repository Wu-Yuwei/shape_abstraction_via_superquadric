function pcd = downSample(pc,ub,lb)
ga = 0.001;
if length(pc.Location) > lb && length(pc.Location) < ub
    pcd = pc;
else
    while 1
        if length(pc.Location) < lb
            pcd = pc;
            break
        end
        pcd = pcdownsample(pc,'gridAverage',ga);
        if length(pcd.Location) > ub
            ga = ga + 0.001;
        elseif length(pcd.Location) > lb
            break
        else
            ga = ga - 0.001;
            pcd = pcdownsample(pc,'gridAverage',ga);
            break
        end
    end
end

        