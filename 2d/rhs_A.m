function dw = rhs_A(t, w, P, O, M, G, fact, data_base_name)


global initial_guess

dw = zeros(length(w),1);

N = length(dw);

slip = w(1:N/2,end);
psi = w(N/2+1:end);

% Create the residual
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

syzF = M.face_interp*syz;
if ~isfield(M, 'mapF')
  M.mapF = M.mapB{7};
end
syzM = syzF(M.mapM(M.mapF));
syzP = syzF(M.mapP(M.mapF));
syz = (syzM + syzP)./2;


sxzF = M.face_interp*sxz;
sxzM = sxzF(M.mapM(M.mapF));
sxzP = sxzF(M.mapP(M.mapF));
sxz = (sxzM + sxzP)./2;

% tauz = M.nxM(M.mapF).*sxz + M.nyM(M.mapF).*syz;
% I = find(abs(tauz) > 1e-5);
% disp([tauz(I), P.tauz0(I)])
tauz = M.nxM(M.mapF).*sxz + M.nyM(M.mapF).*syz + P.tauz0;


V0 = initial_guess;

%V0 = zeros(N/2,1); %initial guess = change to previous V?
ffn = @(V,m) rateStateFriction(V,psi,P,tauz(:),m);
VL = -abs(tauz(:))./P.eta(:);
VR = abs(tauz(:))./P.eta(:);
opt = struct('PassMask',true);
[V,f,iter,error_code] = fnewtbndv(ffn,VL,VR,V0,opt);

initial_guess = V;
dw(1:N/2,1) = V;


for i = 1:N/2
  if P.b(i) ~= 0
    dw(N/2+i) = (P.b(i).*P.v0./P.D_c(i)).*(exp((P.f0-psi(i))./P.b(i)) - abs(V(i))./P.v0);
  else
    dw(N/2+i) = 0;
  end
end


if any(~isfinite(dw)) == 1
  vals = find(~isfinite(dw))
  y = M.face_interp*M.y;
  y = y(M.mapF);
  disp('y');
  y(mod(vals, length(y)))
  disp('V');
  V(mod(vals, length(y)))
  disp('psi');
  psi(mod(vals, length(y)))
  disp('slip');
  psi(mod(vals, length(y)))
  keyboard
end

