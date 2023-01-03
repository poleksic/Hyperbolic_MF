function [P] = uploadParams(METHOD, CV, TARG) 

    BASE_DIR = '/home/aleksandar';
    addpath(strcat(BASE_DIR,'/HYPERBOLIC_MF/data/yam08/'));

    P.eta = {0.05, 0.25};
    P.iter = {200, 800};
    P.J = {1, 2, 4, 8};
    P.c = {1, 4, 16};
    P.theta = {2, 3, 5, 8, 12};
    P.rnk = {5, 10, 20, 40};

    P.lM = {0.0625, 0.125, 0.25, 0.5, 1, 2};
    P.lN = {0.0625, 0.125, 0.25, 0.5, 1, 2};
    P.alphaU = {0.5};
    P.alphaV = {0.5};
    P.pseudo_std = {0.1};
    P.gamma = {0.8, 0.9};

    if strcmp(METHOD,'HDNILMF')
            P.alpha = {0.5, 1.0};
            P.beta = {0.25, 0.5, 1.0};
            P.gamma = {0.25, 0.5, 1.0};
    end

    fn = strcat(TARG,'_R08');
    load ([fn '.dat']);
    P.R = full(spconvert(eval(fn)));
    if strcmp(CV,'LOO')
        P.R = P.R';
    end
end
