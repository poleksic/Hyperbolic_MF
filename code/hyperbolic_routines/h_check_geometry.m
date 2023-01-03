clear all;
% a point on hyperboloid (vertex)
mu = [0;0;1];
% a point y in the tangent plane at mu
y = [1;2;0];

% a point mu1 on hyperboloid
mu1 = [4;4;sqrt(33)];
% a point y1 in the tangent plane at mu1
y1 = [4;4;32/sqrt(33)];
% an arbitrary point from the ambient space
z = [20;21;22];

% CHECK EXPONENTIAL MAP
normy = h_norm(y1);

hpoint = h_exp_map(mu1,y1);
hdist_mu_y = h_distance(mu1,hpoint);

fprintf('normy=%f hdist_mu_y=%f\n',normy, hdist_mu_y);

% CHECK PROJECTION

proj = h_projection(mu1,z);
% check orthogonality of projecion
a = h_inner(proj,proj-z);
b = h_inner(proj,proj);
fprintf('a=%f b=%f\n',a,b);

% CHECK IF POINT IN HYPERBOLIC SPACE (HYPERBOLOID)

% h_vector_check(hpoint)

% number of samples
m = 10000;
% mu = [0 0 1];
mu = [1 1 sqrt(3)]
% Gaussian
G = zeros(m,length(mu));
% Pseudo-hyperbolic-Gaussian
PHG = zeros(m,length(mu));

for i=1:m
    [xnorm,xhyper] = h_pseudo_hyp_Gauss(mu,0.5);
    G(i,:) = xnorm;
    PHG(i,:) = xhyper;
end 

scatter3(G(:,1),G(:,2),G(:,3),1);

hold on 

scatter3(PHG(:,1),PHG(:,2),PHG(:,3),1);