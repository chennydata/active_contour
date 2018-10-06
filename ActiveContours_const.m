function phi = ActiveContours_const(u0, nu,lambda1,lambda2,TimeSteps)
% EXAMPLE:
% u0 is the gray scale image
% This program is written by Chiu-Yen Kao
% Reference: T.F. Chan and L.A. Vese, Active Contours without edges, IEEE
% transactions on Image Processing, 2001, 10(2): 266-277.


[M,N] = size(u0);    % possibly Color picture

% plot figure every plot_iter
plot_iter = 1;
max_iter_phi = 1;

% create an initial guess for phi
x = 1:M; y = 1:N;
dx = 1; dy = 1; dt = (10e7) * 0.9 / ((2*nu/(dx^2)) + (2*nu/(dy^2)));
[X,Y] = meshgrid(x,y);
x0 = (M+1)/2;
y0 = (N+1)/2;
r0 = (min(M,N))/3;
phi = sqrt((X-x0).^2 + (Y-y0).^2) - r0;
phi = phi';


for(n=1:TimeSteps)
    disp(['The number of iteration is ' num2str(n)])
    if mod(n-1,plot_iter) ==0
        figure(1);clf
        imagesc(u0, [0 255]); colormap(gray);
        axis image;
        hold on
        axis image;
        contour(Y,X,phi',[0 0],'m');drawnow
    end

    c1 = C1(u0,phi);
    c2 = C2(u0,phi);
    
    % Gauss-Jacobi iteration
    iter = 0;
    error_phi = 1;
    while (error_phi > 1.e-1) & (iter<max_iter_phi)
        iter = iter+1;        
        phi_xp = [phi(2:M,:); phi(M,:)];      % vertical concatenation
        phi_xm = [phi(1,:);   phi(1:M-1,:)];        % (since x values are rows)
        phi_yp = [phi(:,2:N) phi(:,N)];       % horizontal concatenation
        phi_ym = [phi(:,1)   phi(:,1:N-1)];         % (since y values are columns)
        Dx_m = (phi - phi_xm)/dx;           % first derivatives
        Dx_p = (phi_xp - phi)/dx;
        Dxx = (Dx_p - Dx_m)/dx;     % second derivative
        Dy_m = (phi - phi_ym)/dy;
        Dy_p = (phi_yp - phi)/dy;
        Dyy = (Dy_p - Dy_m)/dy;
        Dx_0 = (phi_xp - phi_xm) / (2*dx);
        Dy_0 = (phi_yp - phi_ym) / (2*dy);
        C1x = 1 ./ sqrt(Dx_p.^2  + Dy_0.^2    + (10e-7)^2);
        C2x = 1 ./ sqrt(Dx_m.^2  + Dy_0.^2    + (10e-7)^2);
        C3y = 1 ./ sqrt(Dx_0.^2  + Dy_p.^2    + (10e-7)^2);
        C4y = 1 ./ sqrt(Dx_0.^2  + Dy_m.^2    + (10e-7)^2);
        Grad = DiracDelta( phi, 1);      
        m = (dt/(dx*dy)) * Grad .* nu;
        C = 1 + m.*(C1x + C2x + C3y + C4y);
        phi_temp = (1 ./ C) .* ( phi + m.*(C1x.*phi_xp + C2x.*phi_xm+C3y.*phi_yp + C4y.*phi_ym) + (dt*Grad).*( -lambda1*(u0 - c2).^2 + lambda2*(u0 - c1).^2) );
        error_phi = norm(phi_temp-phi);
        phi = phi_temp;
    end
end
    
figure(1)
imagesc(u0, [0 255]); colormap(gray);
axis image;
hold on
axis image;
contour(Y,X,phi',[0 0],'m');
