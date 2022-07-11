function D = radial_distance_relative(X, para)

alpha_phi = 1e-4;

e1 = para(1);
e2 = para(2);
a1 = para(3);
a2 = para(4);
a3 = para(5);
gamma = para(6);
beta = para(7);
alpha = para(8);
tx = para(9);
ty = para(10);
tz = para(11);

kx = para(12);
ky = para(13);

R = eul2rotm([gamma, beta, alpha]);
t = [tx;ty;tz];

Xc = R \ (X-t);
xc = Xc(1,:);
yc = Xc(2,:);
zc = Xc(3,:);

xc = xc ./ (kx/a3 * zc + 1);
yc = yc ./ (ky/a3 * zc + 1);
Xc = [xc;yc;zc];

HH = H(xc,yc,zc, a1, a2, a3, e1, e2, alpha_phi);

D = abs(HH .^ (-1/2) - 1);


end


function value = phi(x, y, t, alpha)
    u = max(x,y);
    v = min(x,y);
    
    if t > alpha
        value = u .* g(x,y,t);
    else
        value = u .* ((g(x,y,alpha)-1).*t/alpha + 1);
    end
    assert(sum(isnan(value)) == 0) 
end

function value = g(x, y, t)
    idx = (x == 0) & (y == 0);
    u = max(x,y);
    v = min(x,y);
    value = (1 + (v./u).^(1/t)).^t;
    value(1,idx) = 0;
    assert(sum(isnan(value)) == 0)
end

function value = H(x,y,z, a1, a2, a3, e1, e2, alpha)
    x_bar = (x/a1).^2;
    y_bar = (y/a2).^2;
    z_bar = (z/a3).^2;
    
    value = phi(phi(x_bar, y_bar, e2, alpha), z_bar, e1, alpha);
    assert(sum(isnan(value)) == 0)
end