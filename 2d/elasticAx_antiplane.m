function [Ax, Ay, Az, vol_sxz, vol_syz] = elasticAx_antiplane(ux, uy, uz, O, M, P, res)

Ax = zeros(M.Np, M.num_elm);
Ay = zeros(M.Np, M.num_elm);

lam = P.lam;
mu = P.mu;

do_res = nargin > 6;

% {{{ Compute derivatives of ux, uy, uz
if isempty(uz)
  uz = zeros(M.Np, M.num_elm);
end
uz_r = O.Dr * uz;
uz_s = O.Ds * uz;
uz_x = M.drdx .* uz_r + M.dsdx .* uz_s;
uz_y = M.drdy .* uz_r + M.dsdy .* uz_s;
% }}}

% {{{ Compute the stress
%% Volume terms
vol_sxz = mu .* uz_x;
vol_syz = mu .* uz_y;

%% Get plus and minus side face values
uzF = M.face_interp * uz;

uzM = uzF(M.mapM);

uzP = uzF(M.mapP);

%% surface terms
uzS = (uzM + uzP) / 2;
for b = 1:length(M.mapB)
  if M.bc(b) == 1 % Displacement
    uzS(M.mapB{b}) = 0;
    if do_res
      uzS(M.mapB{b}) = res.bc{b}.uz;
    end
  elseif M.bc(b) == 2 % Traction
    uzS(M.mapB{b}) = uzM(M.mapB{b});
  elseif M.bc(b) == 7 || M.bc(b) == 8 || M.bc(b) == 3% Displacement jump
    if do_res
      uzS(M.mapB{b}) = uzS(M.mapB{b}) + res.bc{b}.duz/2;
    end
  else
    if ~isempty(M.mapB{b})
      error('invalid bc')
    end
  end
end

gxz = (uzS - uzM) .* M.nxM / 2;
gyz = (uzS - uzM) .* M.nyM / 2;

sxzM  = O.Mf * (M.sJ .* (2 * mu .* gxz));
syzM  = O.Mf * (M.sJ .* (2 * mu .* gyz));

%% Lift surface terms back to the volume
sxz = vol_sxz + O.MI * (M.JI .* (M.face_interp' * sxzM));
syz = vol_syz + O.MI * (M.JI .* (M.face_interp' * syzM));

% }}}

% {{{ elastic matrix

%% volume terms
Az = 0;
Az = Az - O.Dr' * ( O.M  * ( M.drdx .* M.J .* sxz));
Az = Az - O.Ds' * ( O.M  * ( M.dsdx .* M.J .* sxz));
Az = Az - O.Dr' * ( O.M  * ( M.drdy .* M.J .* syz));
Az = Az - O.Ds' * ( O.M  * ( M.dsdy .* M.J .* syz));

if do_res && isfield(res, 'fz')
  Az = Az - O.M  * ( M.J .* res.fz);
end

%% surface flux term
% Get the plus and minus side stresses
% SIPG Flux
sxzF = M.face_interp * vol_sxz;
syzF = M.face_interp * vol_syz;

sxzM = sxzF(M.mapM);
syzM = syzF(M.mapM);

sxzP = sxzF(M.mapP);
syzP = syzF(M.mapP);

% Compute average traction across the face
TzM = M.nxM .* sxzM + M.nyM .* syzM;

TzP = M.nxM .* sxzP + M.nyM .* syzP;

TzS = (TzM + TzP) / 2;
for b = 1:length(M.mapB)
  if M.bc(b) == 1 % Displacement
    TzS(M.mapB{b}) = TzM(M.mapB{b});
  elseif M.bc(b) == 2 % Traction
    TzS(M.mapB{b}) = 0;
    if do_res
      TzS(M.mapB{b}) = res.bc{b}.Tz;
    end
  elseif M.bc(b) == 7 || M.bc(b) == 8 || M.bc(b) == 3 % Displacement jump
  else
    if ~isempty(M.mapB{b})
      error('invalid bc')
    end
  end
end

Az = Az + M.face_interp' * (O.Mf * (M.sJ .* TzS));

%% surface penalty term
uzD = uzM - uzP;
for b = 1:length(M.mapB)
  if M.bc(b) == 1 % Displacement
    uzD(M.mapB{b}) = uzM(M.mapB{b});
    if do_res
      uzD(M.mapB{b}) = uzD(M.mapB{b}) - res.bc{b}.uz;
    end
  elseif M.bc(b) == 2 % Traction
    uzD(M.mapB{b}) = 0;
  elseif M.bc(b) == 7 || M.bc(b) == 8 || M.bc(b) == 3 % Displacement jump
    if do_res
      uzD(M.mapB{b}) = uzD(M.mapB{b}) - res.bc{b}.duz;
    end
  else
    if ~isempty(M.mapB{b})
      error('invalid bc')
    end
  end
end

gxz = uzD .* M.nxM / 2;
gyz = uzD .* M.nyM / 2;

fxzD  = (2 * mu .* gxz);
fyzD  = (2 * mu .* gyz);

Az = Az - P.beta * M.face_interp' * ( O.Mf * ( (M.nxM .* fxzD + ...
                                              M.nyM .* fyzD)));
% }}}

Az = -Az;
