function [D,J] = radial_distance_with_gradient(X, para)

%% Approximation parameters for the auxiliary functions
alpha_phi = 1e-4;
alpha_l = 1e-12;
alpha_psi = 0.1;
beta_psi = 0.5;
alpha_h = 1e-7;

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

xx = X(1,:);
yy = X(2,:);
zz = X(3,:);

R = eul2rotm([gamma, beta, alpha]);
t = [tx;ty;tz];

Xc = R \ (X-t);
xc = Xc(1,:);
yc = Xc(2,:);
zc = Xc(3,:);

xc = xc ./ (kx./a3 .* zc + 1);
yc = yc ./ (ky./a3 .* zc + 1);

Xc = [xc;yc;zc];

phi1 = phi((xc./a1).^2, (yc./a2).^2, e2, alpha_phi);
psi1 = psi(phi1, (zc./a3).^2, e1, alpha_psi, beta_psi);
psi1_ = psi((zc./a3).^2, phi1, e1, alpha_psi, beta_psi);

psi2 = psi((xc./a1).^2, (yc./a2).^2, e2, alpha_psi, beta_psi);
psi2_ = psi((yc./a2).^2, (xc./a1).^2, e2, alpha_psi, beta_psi);

HH = H(xc,yc,zc, a1, a2, a3, e1, e2, alpha_phi,phi1);
r_0 = vecnorm(Xc);
D = r_0 .* (HH .^ (-1./2) - 1);
DH = dH_dX(xc,yc,zc,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH, psi1, psi2,psi1_,psi2_);


J = [dG_de1(xc,yc,zc,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH, psi1,psi2,psi1_,psi2_);
     dG_de2(xc,yc,zc,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH, psi1,psi2,psi1_,psi2_);
     dG_da1(xc,yc,zc,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, psi1,psi2,psi1_,psi2_);
     dG_da2(xc,yc,zc,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, psi1,psi2,psi1_,psi2_);
     dG_da3(xc,yc,zc,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, kx, ky, alpha, beta, gamma, tx,ty,tz, psi1,psi2,psi1_,psi2_);
     dG_dgamma(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH,kx,ky)
     dG_dbeta(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH,kx,ky)
     dG_dalpha(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH,a1,a2,a3, HH,kx,ky)
     dG_dtx(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH,kx,ky);
     dG_dty(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH,kx,ky);
     dG_dtz(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH,kx,ky);
     dG_dkx(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH,kx,ky);
     dG_dky(xc,yc,zc,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH,kx,ky)];

% D = D;
J = J';



end

%% auxiliary functions
function value = phi(x, y, t, alpha)
    u = max(x,y);
    v = min(x,y);
    
    if t > alpha
        value = u .* g(x,y,t);
    else
        value = u .* ((g(x,y,alpha)-1).*t./alpha + 1);
    end
    assert(sum(isnan(value)) == 0) 
end

function value = g(x, y, t)
    idx = (x == 0) & (y == 0);
    u = max(x,y);
    v = min(x,y);
    value = (1 + (v./u).^(1./t)).^t;
    value(1,idx) = 0;
    assert(sum(isnan(value)) == 0)
end

function value = l(t, alpha)
    value = zeros(size(t));
    idx = find(t >= alpha);
    value(1, idx) = t(1, idx) .* log(t(1, idx));
    assert(sum(isnan(value)) == 0)
end

function value = f(r, t, alpha, beta)
    value = 1 ./ (1 + r.^(1./t));
    if t < alpha
        idx = find((1-r) < beta);
        value(1,idx) = 1 ./ (1 + r(1,idx).^(1./alpha));
    end
    assert(sum(isnan(value)) == 0)
end        

function value = psi(x, y, t, alpha, beta)
    value = nan(size(x));
    value(1,(x == 0) & (y == 0)) = 1./2;
    
    idx1 = find((x-y) >= 0 & x ~=0);
    x1 = x(1, idx1);
    y1 = y(1, idx1);
    rr = y1./x1;
    value(1, idx1) = f(rr, t, alpha, beta);
    
    idx2 = find((y-x) > 0 & y ~= 0);
    x2 = x(1, idx2);
    y2 = y(1, idx2);
    rr = x2./y2;
    value(1, idx2) = 1 - f(rr, t, alpha, beta);
    assert(sum(isnan(value)) == 0)
end

function value = H(x,y,z, a1, a2, a3, e1, e2, alpha, phi1)
    x_bar = (x./a1).^2;
    y_bar = (y./a2).^2;
    z_bar = (z./a3).^2;
    
    value = phi(phi1, z_bar, e1, alpha);
    assert(sum(isnan(value)) == 0)
end

% function value = F(x,y,z, a1, a2, a3, e1, e2, alpha)
%     HH = H(x,y,z, a1, a2, a3, e1, e2, alpha);
%     r_0 = vecnorm([x;y;z]);
%     value = r_0 .* abs(HH .^ (-1./2) - 1);
%     assert(sum(isnan(value)) == 0)
% end

function value = dH_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, psi1,psi2,psi1_,psi2_)
%     x_bar = (x./a1).^2;
%     y_bar = (y./a2).^2;
%     z_bar = (z./a3).^2;
    
    value = -2./a1 .* HH .*  ...
        psi1 .* psi2;
    assert(sum(isnan(value)) == 0)
end

function value = dH_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, psi1,psi2,psi1_,psi2_)
%     x_bar = (x./a1).^2;
%     y_bar = (y./a2).^2;
%     z_bar = (z./a3).^2;
    
    value = -2./a2 .* HH .*  ...
        psi1 .* psi2_;
    assert(sum(isnan(value)) == 0)
end

function value = dH_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, psi1,psi2,psi1_,psi2_)
%     x_bar = (x./a1).^2;
%     y_bar = (y./a2).^2;
%     z_bar = (z./a3).^2;
    
    value = -2./a3 .* HH .* ...
        psi2_; 
    assert(sum(isnan(value)) == 0)
end

function value = dH_de1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH, psi1,psi2,psi1_,psi2_)
%     x_bar = (x./a1).^2;
%     y_bar = (y./a2).^2;
%     z_bar = (z./a3).^2;
    
    value = -(HH .* ...
        (l(psi1_, alpha_l) + ...
         l(psi1, alpha_l)));
    assert(sum(isnan(value)) == 0)
end

function value = dH_de2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH, psi1,psi2,psi1_,psi2_)
%     x_bar = (x./a1).^2;
%     y_bar = (y./a2).^2;
%     z_bar = (z./a3).^2;
    
    value = (-HH .* ...
        psi1 .* ...
        (l(psi2,alpha_l) + l(psi2_, alpha_l)));
    assert(sum(isnan(value)) == 0)
end

function value = dH_dx(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h,HH, psi1,psi2,psi1_,psi2_)
    dHda1 = dH_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi,HH, psi1,psi2,psi1_,psi2_);
    value = -x ./ (alpha_h.^2.*a1) .* dHda1;
    idx = abs(x./a1) > alpha_h;
    value(1,idx) = -a1./x(1,idx) .* dHda1(1,idx);
    assert(sum(isnan(value)) == 0)
end

function value = dH_dy(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h,HH, psi1,psi2,psi1_,psi2_)
    dHda2 = dH_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi,HH, psi1,psi2,psi1_,psi2_);
    value = -y ./ (alpha_h.^2.*a2) .* dHda2;
    idx = abs(y./a2) > alpha_h;
    value(1,idx) = -a2./y(1,idx) .* dHda2(1,idx);
    assert(sum(isnan(value)) == 0)
end

function value = dH_dz(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h,HH, psi1,psi2,psi1_,psi2_)
    dHda3 = dH_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi,HH, psi1,psi2,psi1_,psi2_);
    value = -z ./ (alpha_h.^2.*a3) .* dHda3;
    idx = abs(z./a3) > alpha_h;
    value(1,idx) = -a3./z(1,idx) .* dHda3(1,idx);
    assert(sum(isnan(value)) == 0)
end

function value = dT_dgamma(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = (cos(beta).*cos(gamma).*(ty - y) - cos(beta).*sin(gamma).*(tx - x))./((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (kx.*((cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(tx - x) + (sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(ty - y)).*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)))./(a3.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(2,:) = - ((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(tx - x) + (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(ty - y))./((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (ky.*((cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(tx - x) + (sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(ty - y)).*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)))./(a3.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(3,:) = - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(tx - x) - (sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(ty - y);

    assert(sum(isnan(value),'all') == 0)
end

function value = dT_dbeta(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = - (cos(beta).*(tz - z) + cos(gamma).*sin(beta).*(tx - x) + sin(beta).*sin(gamma).*(ty - y))./((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (kx.*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)).*(cos(alpha).*cos(beta).*cos(gamma).*(tx - x) - cos(alpha).*sin(beta).*(tz - z) + cos(alpha).*cos(beta).*sin(gamma).*(ty - y)))./(a3.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(2,:) = (cos(beta).*cos(gamma).*sin(alpha).*(tx - x) - sin(alpha).*sin(beta).*(tz - z) + cos(beta).*sin(alpha).*sin(gamma).*(ty - y))./((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (ky.*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)).*(cos(alpha).*cos(beta).*cos(gamma).*(tx - x) - cos(alpha).*sin(beta).*(tz - z) + cos(alpha).*cos(beta).*sin(gamma).*(ty - y)))./(a3.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(3,:) = cos(alpha).*sin(beta).*(tz - z) - cos(alpha).*cos(beta).*cos(gamma).*(tx - x) - cos(alpha).*cos(beta).*sin(gamma).*(ty - y);

    assert(sum(isnan(value),'all') == 0)
end

function value = dT_dalpha(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = (kx.*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)).*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)))./(a3.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(2,:) = ((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z))./((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) + (ky.*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)).^2)./(a3.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(3,:) = (cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z);

    assert(sum(isnan(value),'all') == 0)
end

function value = dT_dtx(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = (cos(beta).*cos(gamma))./((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (kx.*(sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)))./(a3.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(2,:) = - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta))./((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (ky.*(sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)))./(a3.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(3,:) = (- sin(alpha).*sin(gamma) - cos(alpha).*cos(gamma).*sin(beta)) .* ones(size(x));
    assert(sum(isnan(value),'all') == 0)
end

function value = dT_dty(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = (cos(beta).*sin(gamma))./((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) + (kx.*(cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)))./(a3.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(2,:) = (cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma))./((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) + (ky.*(cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)))./(a3.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(3,:) = cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma) .* ones(size(x));
    assert(sum(isnan(value),'all') == 0)
end

function value = dT_dtz(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = - sin(beta)./((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (kx.*cos(alpha).*cos(beta).*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)))./(a3.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(2,:) = (cos(beta).*sin(alpha))./((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1) - (ky.*cos(alpha).*cos(beta).*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)))./(a3.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(3,:) = (-cos(alpha).*cos(beta)) .* ones(size(x));
    assert(sum(isnan(value),'all') == 0)
end

function value = dT_dkx(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = -(((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)).*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)))./(a3.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(2,:) = 0 .* ones(size(x));
    value(3,:) = 0 .* ones(size(x));
    assert(sum(isnan(value),'all') == 0)
end

function value = dT_dky(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3)
    value(1,:) = 0 .* ones(size(x));
    value(2,:) = -(((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)).*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)))./(a3.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    value(3,:) = 0 .* ones(size(x));
    assert(sum(isnan(value),'all') == 0)
end

function value = dH_dX(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH, psi1,psi2,psi1_,psi2_)
    value = [dH_dx(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH, psi1,psi2,psi1_,psi2_);
             dH_dy(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH, psi1,psi2,psi1_,psi2_);
             dH_dz(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_h, HH, psi1,psi2,psi1_,psi2_)];
    assert(sum(isnan(value),'all') == 0)
end

function [value,DT] = dH_dgamma(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)
    DT = dT_dgamma(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function [value,DT] = dH_dbeta(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)     
    DT = dT_dbeta(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function [value,DT] = dH_dalpha(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)
    DT = dT_dalpha(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function [value,DT] = dH_dtx(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)
    DT = dT_dtx(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function [value,DT] = dH_dty(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)
    DT = dT_dty(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function [value,DT] = dH_dtz(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)
    DT = dT_dtz(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function [value,DT] = dH_dkx(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)
    DT = dT_dkx(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function [value,DT] = dH_dky(gamma,beta,alpha,tx,ty,tz, kx,ky,x,y,z, DH,a1,a2,a3)
    DT = dT_dky(gamma,beta,alpha,tx,ty,tz,kx,ky,x,y,z,a1,a2,a3);
    value = dot(DH, DT);
    assert(sum(isnan(value)) == 0)
end

function value = dG_de1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH, psi1,psi2,psi1_,psi2_)
    value = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2) .* dH_de1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l,HH, psi1,psi2,psi1_,psi2_));
    assert(sum(isnan(value)) == 0)
end

function value = dG_de2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l, HH, psi1,psi2,psi1_,psi2_)
    value = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* dH_de2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, alpha_l,HH, psi1,psi2,psi1_,psi2_);
    assert(sum(isnan(value)) == 0)
end

function value = dG_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, psi1,psi2,psi1_,psi2_)
    value = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* dH_da1(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi,HH, psi1,psi2,psi1_,psi2_);
    assert(sum(isnan(value)) == 0)
end

function value = dG_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, psi1,psi2,psi1_,psi2_)
    value = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* dH_da2(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi,HH, psi1,psi2,psi1_,psi2_);
    assert(sum(isnan(value)) == 0)
end

function value = dG_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi, HH, kx, ky, alpha, beta, gamma, tx,ty,tz, psi1,psi2,psi1_,psi2_)
    dxc_da3 = (kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)).*(cos(beta).*cos(gamma).*(tx - x) - sin(beta).*(tz - z) + cos(beta).*sin(gamma).*(ty - y)))./(a3.^2.*((kx.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    dyc_da3 = (ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)).*((cos(alpha).*cos(gamma) + sin(alpha).*sin(beta).*sin(gamma)).*(ty - y) - (cos(alpha).*sin(gamma) - cos(gamma).*sin(alpha).*sin(beta)).*(tx - x) + cos(beta).*sin(alpha).*(tz - z)))./(a3.^2.*((ky.*((sin(alpha).*sin(gamma) + cos(alpha).*cos(gamma).*sin(beta)).*(tx - x) - (cos(gamma).*sin(alpha) - cos(alpha).*sin(beta).*sin(gamma)).*(ty - y) + cos(alpha).*cos(beta).*(tz - z)))./a3 - 1).^2);
    dzc_da3 = 0;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*dxc_da3 + y.*dyc_da3 + z.*dzc_da3).* (HH.^(-1./2)-1);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* dH_da3(x,y,z,a1,a2,a3,e1,e2, alpha_phi, alpha_psi, beta_psi,HH, psi1,psi2,psi1_,psi2_);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dgamma(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH,a1,a2,a3, HH, kx,ky)
    [value2,DT] = dH_dgamma(gamma,beta,alpha,tx,ty,tz, kx,ky, xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)).* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dbeta(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH, kx,ky)
    [value2,DT] = dH_dbeta(gamma,beta,alpha,tx,ty,tz, kx,ky, xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)) .* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dalpha(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH, kx,ky)
    [value2,DT] = dH_dalpha(gamma,beta,alpha,tx,ty,tz, kx,ky,xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)) .* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dtx(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH, kx,ky)
    [value2,DT] = dH_dtx(gamma,beta,alpha,tx,ty,tz, kx,ky, xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)) .* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dty(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH, kx,ky)
    [value2,DT] = dH_dty(gamma,beta,alpha,tx,ty,tz, kx,ky,xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)) .* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dtz(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3,HH, kx,ky)
    [value2,DT] = dH_dtz(gamma,beta,alpha,tx,ty,tz, kx,ky, xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)) .* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dkx(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH, kx,ky)
    [value2,DT] = dH_dkx(gamma,beta,alpha,tx,ty,tz, kx,ky, xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)) .* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end

function value = dG_dky(x,y,z,gamma,beta,alpha,tx,ty,tz, xx,yy,zz, DH, a1,a2,a3, HH, kx,ky)
    [value2,DT] = dH_dky(gamma,beta,alpha,tx,ty,tz,kx,ky,xx,yy,zz, DH,a1,a2,a3);
    value2 = vecnorm([x;y;z]) .* (-1./2.*HH.^(-3./2)) .* value2;
    value1 = (x.^2+y.^2+z.^2).^(-1./2) .* (x.*DT(1) + y.*DT(2) + z.*DT(3)) .* (HH.^(-1./2)-1);
    value =value1 + value2;
    assert(sum(isnan(value)) == 0)
end



































