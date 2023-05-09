function dydt = non_linear_dynamics(~,x,u,m,I_z,l_r,l_f,g,B,C,mu)
    F_z = m * g/2; % Normal force on each tire
    D = mu * F_z; 
    
    % Define the slip angles:
    alpha_R = atan((x(5) - l_r *x(6))/x(4));
    alpha_F = atan((x(5) + l_f *x(6))/x(4)) - x(7);
    
    F_Ry = D * sin(C*atan(B*alpha_R));
    F_Fy = D * sin(C*atan(B*alpha_F));


    dydt = [x(4)*cos(x(3)) - x(5)*sin(x(3))
       x(4)*sin(x(3)) + x(5)*cos(x(3))
       x(6)
       u(1)/m
       (x(4)*u(2) + x(7)*u(1)/m)*(l_r/(l_r+l_f))
       (x(4)*u(2) + x(7)*u(1)/m)*(1/(l_r+l_f))
       u(2)
        ];

end