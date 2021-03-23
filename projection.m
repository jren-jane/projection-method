
%%%%%%%%%%%%%%%%%%%%%% code by Jing Ren 01/896410 %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% parameters and known variable values %%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% For the ease of the viewer, I am keeping the %%%%%%%%%%%%%%%
%%%%%%%%%%%%% the state spaces and the value functions for %%%%%%%%%%%%%%%
%%%%%%%%%%%%% the good and bad states separately, although %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% they are in fact, always stacked together %%%%%%%%%%%%%%%%

  alpha = 0.36;
  beta = 0.99;
  delta = 0.025;

  pr_gg = 0.977;
  pr_gb = 0.074;
  pr_bg = 0.023;
  pr_bb = 0.926;
  
%  P = [pr_gg , pr_bg ; pr_gb , pr_bb]; % Set up the transition matrix
  
  zg = 1.05;      % parameter in the production function for the good state
  zb = 0.842;     % parameter in the production function for the bad state
  eg = 1.088;     % future expectation of z when the current state is good
  eb = 0.804;     % future expectation of z when the current state is bad
% The above values are obtained from the previous question in problem set 2
  
%%%%%%%%%%%%%%%% Initial guess for the value functions %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% The current values are stored in columns %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% The future values are stored in rows %%%%%%%%%%%%%%%%5%%%

  kg = (0.05 : 0.05 : 20)';     % state space for capital in the good state
  kb = (0.05 : 0.05 : 20)';     % state space for capital in the bad state
  % Caution: Here 0 cannot be included, otherwise a "-inf" will show up in 
  % the value function
  N = size(kg , 1);             % number of states for capital
%  k = [kg ; kb];                % Stack the two vectors together
  
  kgcs = (0.01 : 0.01 : 20);    % "cs" is to denote points in cubic spline
  kbcs = (0.01 : 0.01 : 20);    % "cs" is to denote points in cubic spline
  M = size(kgcs , 2);           % number of points in cubic spline
%  kcs = [kgcs , kbcs];          % Stack the two vectors together
  
  vg = ones(N , 1);          % Each element corresponds to a current state
  vb = ones(N , 1);
  v = [vg ; vb];             % stack the two vectors together
  
%  vgnew = ones(N , 1);       % the updated value function
%  vbnew = ones(N , 1);       % the updated value function
%  vnew = [vgnew ; vbnew];    % Stack the two vectors together
  
  vgcs = ones(1 , M);        
  vbcs = ones(1 , M);
%  vcs = [vgcs , vbcs];
  
%%%%%%%%%%%%%%%%%% Set up the current return matrix %%%%%%%%%%%%%%%%%%%%%%

%   R = ones(2 * N , 2 * M);   % Initialize the current return matrix.
%                              % This part does not need to be updated
%                              % through iterations and therefore is written
%                              % separately.
%                              
  C = [zg * kg.^alpha ; zb * kb.^alpha] * ones(1 , 2 * M) ...
      - ones(2 * N , 1) * [kgcs , kbcs] ...
      + (1 - delta) * [kg ; kb] * ones(1 , 2 * M); % the consumption matrix
  
  R = log(C .* (C > 0)) .* (C > 0) - 999 * (C <= 0);
  
  % If utility or consumption is negative, utility = -999

%%%%%%%%%%%%%%%%%%%%%%%%%%% Start iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  MxIt = 1000;             % maximum number of iterations
  eps = 0.01;              % convergence criterion

for l = 1 : MxIt
    
    vgcs = csapi(kg , vg , kgcs); % values on the points in the cubic spline
    vbcs = csapi(kb , vb , kbcs); % values on the points in the cubic spline
    
    ezg = pr_gg * vgcs + pr_bg * vbcs; % expectation for z in the good state
    ezb = pr_gb * vgcs + pr_bb * vbcs; % expectation for z in the bad state
    
    W = R + beta * ones(2 * N , 1) * [ezg , ezb]; % RHS inside the maximum
    
    [W , index] = max(W');% Because the future states is in the rows,
                          % we need to transpose the matrix to maximize
                          % over the columns
    
    vnew = W';            % Transpose back so that it is a column vector
    Pol = index';         % Transpose the index in row vector back to 
                          % a column vector for the policy function
    
%%%%%%%%%%%%%%%%%%%%%%%% convergence criterion %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if max(abs(vnew - v)) < eps
        fprintf(1 , 'The value function has converged in %2.0f iterations.\n', l);
        break;
    end

    v = vnew;            % If criteria not met, update the value function
    
    vg = v(1 : N);       % Pass on the values to each value function
    vb = v(N + 1 : 2 * N);% Pass on the values to each value function
    
end

%%%%%%%%%%%%%%%%%%%%%%% Plot the value functions %%%%%%%%%%%%%%%%%%%%%%%%%

  subplot(2,2,1)
  plot(kg , vg)
  xlabel('value function for the good state')

  subplot(2,2,2)
  plot(kb , vb)
  xlabel('value function for the bad state')
  
%%%%%%%%%%%%%%%%%%%%%% Plot the policy functions %%%%%%%%%%%%%%%%%%%%%%%%%

  subplot(2,2,3)
  plot(kg , Pol(1:N))
  xlabel('policy function for the good state')

  subplot(2,2,4)
  plot(kb , Pol(N+1:2*N))
  xlabel('policy function for the bad state')