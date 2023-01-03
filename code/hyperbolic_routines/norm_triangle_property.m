clear;
GEOM = 'HYP';
sigma = 1;
num = 5000;
mu_0 = [0 0 1];
a = -10; b=10;
for i=1:num 

    if ~strcmp(GEOM,'HYP')
        ux = a + (b-a)*rand();
        uy = a + (b-a)*rand();
        uz = a + (b-a)*rand();
        if uz < 0 
            continue
        end
        u = [ux, uy, uz];

        vx = a + (b-a)*rand();
        vy = a + (b-a)*rand();
        vz = a + (b-a)*rand();
        if vz < 0 
            continue
        end

        v = [vx, vy, vz];
        
        if ux*ux + uy*uy - uz*uz > 0 || ux*ux + uy*uy - uz*uz < -1
            continue
        end

        if vx*vx + vy*vy - vz*vz > 0 || vx*vx + vy*vy - vz*vz < -1
            continue
        end
        
    else
        [~,u] = h_pseudo_hyp_Gauss(mu_0,sigma);
        [~,v] = h_pseudo_hyp_Gauss(mu_0,sigma);
    end
    
%     if h_squared_Lorentz_distance(u,v) <= h_squared_Lorentz_distance(u,mu_0) + h_squared_Lorentz_distance(v,mu_0)
%         'triangle'
%     else
%         'not triangle'
%     end

        h_cosh_angle(u,v)
end