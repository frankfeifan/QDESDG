addpath ../

if ~exist('do_not_clear')
  do_not_clear = false;
end

if ~do_not_clear
  clear
end

% {{{ Init git stuff
[~, git_hash] = system('git rev-parse --short HEAD');
git_hash = strtrim(git_hash);

[~,git_clean] = system('git status --untracked-files=no --porcelain');
if numel(git_clean) > 0
  git_hash = [git_hash,'_unclean'];
  warning('git is dirty')
end
% }}}

% {{{ Parameter);
if ~exist('mesh_version')
  mesh_version = 1;
end
if ~exist('N')
  N = 4;
end
if ~exist('R')
  R = 0;
end
if ~exist('V_p')
  V_p = 1e-9;
end
if ~exist('P_beta')
  P_beta = 10;
end

sim_years = 1200;

P.plt = usejava('jvm');

data_base_name = ['data/BP2', ...
                  '_N_', num2str(N), ...
                  '_R_', num2str(R), ...
                  '_P_beta_', num2str(P_beta), ...
                  '_mesh_version_', num2str(mesh_version), ...
                  '_', git_hash];
disp(data_base_name)
% }}}

% {{{ elasticity parameters
P.rho = 2.670;
P.cs = 3.464;
P.mu = P.cs^2 * P.rho;
P.lam = P.cs^2 * P.rho;
P.eta = P.mu./(2*P.cs); %radiation damping parameter
% }}}

% {{{ rate-and-state parameters etc.
P.V_p = 1e-9;
P.H1 = 15;
P.H2 = 18;
P.a0 = 0.01;
P.amax = 0.025;
P.b0 = 0.015;
P.s_NORM = 50;
P.D_c = 0.008;
P.slip_init = 0;
P.f0 = 0.6;
P.v0 = 1e-6;
P.vinit = 1e-9;
P.beta = P_beta * N^2;
% }}}

% {{{ Create the dip mesh:
mesh_file = ['meshes/BP2_V', num2str(mesh_version),'.inp'];
[xg, yg, EToV, EToE, EToF, EToO, EToB, EToT] = inp_read_2d(mesh_file);
[O, M, G] = init_mesh(N, [], R, xg, yg, EToV, EToE, EToF, EToO, EToB, EToT);
if P.plt
  plot_field([], M, G);
  axis image
end
% }}}

% {{{ Create the operator and factor
A = assemble_matrix(O, M, P, G);
Ntotal = M.Np * M.num_elm;
Az = A(2 * Ntotal + (1:Ntotal), 2 * Ntotal + (1:Ntotal));
disp('factoring...')
[fact.L, fact.U, fact.P, fact.Q, fact.R] = lu(Az);
%{
[fact.P, g, fact.L] = chol(Az, 'lower');
if g ~= 0
  error ('Matrix must be positive definite.') ;
end
%}
% }}}

% {{{ Plotting stuff
% map from fault data to which fault
fault_map{1} = reshape(1:numel(M.mapB{7}), size(M.mapB{7}));
M.fault_map = fault_map;
mapF = [M.mapB{7}];
M.mapF = mapF;

yF = M.face_interp * M.y;
yfault = yF(mapF(:));

xF = M.face_interp * M.x;
xfault = xF(mapF(:));

plot_fault{1} = reshape(find(M.nxM(mapF(fault_map{1}))>0), ...
                size(fault_map{1}, 1), size(fault_map{1}, 2) / 2);

[~,I] = sort(yfault(plot_fault{1}(1, :)));
plot_fault{1} = plot_fault{1}(:, I);

M.plot_fault = plot_fault;
% hold on
% plot(xfault, yfault, '*-')
% hold off
% }}}

% {{{ Create the fault data

ndof_fault = numel(xfault);

P.a = P.a0 - (P.a0 - P.amax) * min(1, max(0, (P.H1 + yfault)/(P.H1 - P.H2)));
P.b = P.b0 * ones(ndof_fault, 1);

P.tauz0 = P.s_NORM * P.amax * ...
          asinh(P.vinit/(2*P.v0) *...
                exp((P.f0 + P.b0 * log(P.v0 / P.vinit))/P.amax)) + ...
          P.eta * P.vinit;

slip = zeros(ndof_fault, 1);
Q = (P.D_c / P.v0) *...
    exp((P.a ./ P.b) .* log((2*P.v0 ./ P.vinit) * ...
        sinh((P.tauz0 - P.eta * P.vinit) ./ (P.a * P.s_NORM))) - P.f0 ./ P.b);
psi = P.f0 + P.b .* log(P.v0 * Q ./ P.D_c);

P.tauz0 = P.tauz0 * sign(M.nxM(M.mapF));
P.s_NORM =  P.s_NORM * ones(ndof_fault, 1);
P.eta =  P.eta * ones(ndof_fault, 1);
P.D_c =  P.D_c * ones(ndof_fault, 1);

%{
figure(1)
plot(psi(plot_fault{1}(:)), yfault(plot_fault{1}(:)))
title('psi0')
axis tight
xlim(xlim + diff(xlim)*0.1*[-1 1])
print -dpng psi.png

figure(2)
plot(Q(plot_fault{1}(:)), yfault(plot_fault{1}(:)))
title('Q0')
axis tight
xlim(xlim + diff(xlim)*0.1*[-1 1])
print -dpng theta.png

figure(3)
plot(P.tauz0(plot_fault{1}(:)), yfault(plot_fault{1}(:)))
title('tau0')
axis tight
xlim(xlim + diff(xlim)*0.1*[-1 1])
print -dpng tau.png

figure(4)
plot(P.a(plot_fault{1}(:)), yfault(plot_fault{1}(:)))
title('a')
axis tight
xlim(xlim + diff(xlim)*0.1*[-1 1])
print -dpng a.png

figure(5)
plot(P.b(plot_fault{1}(:)), yfault(plot_fault{1}(:)))
title('b')
axis tight
xlim(xlim + diff(xlim)*0.1*[-1 1])
print -dpng b.png

%}

% }}}

% {{{ Globals
global initial_guess
Vnew = zeros(ndof_fault, 1);
tauz = zeros(ndof_fault, 1);
initial_guess = Vnew;

global ig
ig.t_n = 0;
ig.ctr = 0;
ig.save_stride_fields = 5;

global ssd
% }}}

save([data_base_name,'_data.mat'], 'fault_map', 'plot_fault', ...
      'xfault', 'yfault', 'N', 'R', 'P', 'V_p', 'mesh_version');

% {{{ Run the ODE solver
options = odeset('OutputFcn', @my_output_fcn_A, ...
                 'reltol', 1e-7, ...
                 'abstol', 1e-7, ...
                 'InitialStep', 1e-6);
disp('starting simulation')
tic
t = 0;
tf = t + sim_years*31556926;
w = [slip;psi];
ssd = initialize_fields(w, t, Vnew, tauz(:), 0, 0, [data_base_name, '_']);
ode45(@rhs_A, [t tf], w, options, P, O, M, G, fact, data_base_name);
toc
% }}}
