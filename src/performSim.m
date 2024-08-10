% Function that performs simulation
function [t, Hs, Vs] = performSim(C_d, C_v, A_x, Q, m_0, T_f, n)
   
    % Time set-up
    inter = [0,T_f];
    dt = (inter(2)-inter(1))/n;
    t(1) = inter(1);

    % Initial Conditions
    Vs(1) = 0;
    Hs(1) = 0;

    % Main Euler Method Loop
    for i = 1:n
        t(i+1) = t(i) + dt;
        [dv, dh] = rocket_step(t(i),Vs(i),Hs(i),m_0, C_d, C_v, A_x, Q, dt);
        
        Vs(i+1) = Vs(i) + dv;
        Hs(i+1) = Hs(i) + dh;
    end
end