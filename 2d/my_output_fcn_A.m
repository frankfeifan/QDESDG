function stop = my_output_fcn_A(t, w, done, P, O, M, G, fact, data_base_name)

stop = false;

global initial_guess
global ssd
global ig

if length(t) ~= 4
   ig.ctr = 0;
   ig.tn = 0;
   tic
   return
else
    t = t(end);
    w = w(:,end);
end

N = length(w);

slip = w(1:N/2,end);
psi = w(N/2+1:end);

res = create_res_fault_antiplane(M, t, slip, P);


[~, ~, bz] = elasticAx_antiplane([], [], [], O, M, P, res);
bz = -bz;
ux = zeros(M.Np, M.num_elm);
uy = zeros(M.Np, M.num_elm);

if isfield(fact, 'Q')
  uz = fact.Q * (fact.U \ (fact.L \ (fact.P * (fact.R \ bz(:)))));
else
  uz = fact.P * (fact.L' \ (fact.L \ (fact.P' * bz(:))));
end
uz = reshape(uz, M.Np, M.num_elm);

[~, ~, ~, sxz, syz] = elasticAx_antiplane(ux, uy, uz, O, M, P, res);

if ~isfield(M, 'mapF')
  M.mapF = M.mapB{7};
end

syzF = M.face_interp*syz;
syzM = syzF(M.mapM(M.mapF));
syzP = syzF(M.mapP(M.mapF));
syzF = (syzM + syzP)./2;


sxzF = M.face_interp*sxz;
sxzM = sxzF(M.mapM(M.mapF));
sxzP = sxzF(M.mapP(M.mapF));
sxzF = (sxzM + sxzP)./2;

tauz = M.nxM(M.mapF).*sxzF + M.nyM(M.mapF).*syzF + P.tauz0;

V0 = initial_guess;

%V0 = zeros(N/2,1); %initial guess = change to previous V?
ffn = @(V,m) rateStateFriction(V,psi,P,tauz(:),m);
VL = -abs(tauz(:))./P.eta(:);
VR = abs(tauz(:))./P.eta(:);
opt = struct('PassMask',true);
[V,f,iter,error_code] = fnewtbndv(ffn,VL,VR,V0,opt);

initial_guess = V;


tauz = tauz(:);

dt = t-ig.t_n;

ig.t_n = t;

fprintf('Avergage time per step = %e sec\n', toc / (ig.ctr+1));
fprintf('time    = %e years :: dt         = %e years\n', ...
        t / 31556926, dt / 31556926);
fprintf('max |V| = %e m / s :: max |tauz| = %e MPa\n\n', ...
        max(abs(V)), max(abs(tauz)));

if (mod(ig.ctr, ig.save_stride_fields) == 0) && P.plt
  yF = M.face_interp * M.y; yF = yF(M.mapF); yF = yF(:);
  subplot(3, 2, 1);
  plot_field(uz, M, G, false), colorbar
  title(sprintf('uz: year = %e, ', t./31556926));
  caxis(P.V_p * t/2 * [-1 1]);
  subplot(3, 2, 2);
  plot_field(sxz, M, G, false), colorbar
  title(sprintf('sxz: year = %e, ', t./31556926));
  if ~isfield(M, 'plot_fault')
    M.plot_fault{1} = 1:N/4;
  end
  for k = 1:length(M.plot_fault)
    PF = M.plot_fault{k};
    subplot(3, 2, 3)
    plot(yF(PF), slip(PF))
    axis tight
    title(sprintf('slip: year = %e, ', t./31556926));
    subplot(3, 2, 4)
    plot(yF(PF), initial_guess(PF))
    axis tight
    title(sprintf('V: year = %e, ', t./31556926));
    subplot(3, 2, 5)
    plot(yF(PF), tauz(PF))
    axis tight
    title(sprintf('tau: year = %e, ', t./31556926));
    subplot(3, 2, 6)
    plot(yF(PF), psi(PF))
    axis tight
    title(sprintf('psi: year = %e, ', t./31556926));
    drawnow
  end
end

%% SAVE STUFF %%%
if (mod(ig.ctr, ig.save_stride_fields) == 0)
  if isempty(t) == 0
    fprintf('---------\nSaving fields: (%s)\n---------\n', data_base_name);
    ssd.w       = SaveStreamData('Write', ssd.w,       w);
    ssd.time    = SaveStreamData('Write', ssd.time,    t);
    ssd.vel     = SaveStreamData('Write', ssd.vel,     V);
    ssd.tau     = SaveStreamData('Write', ssd.tau,     tauz);
    ssd.time_step = SaveStreamData('Write',ssd.time_step, dt);
  end
% plot_field(sxz, M, G),colorbar, %caxis([7 12])
% axis image
% drawnow
end

ig.ctr = ig.ctr + 1;

end



