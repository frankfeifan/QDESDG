function res = create_res_fault_antiplane(M, t, slip, p)

xF = M.face_interp * M.x;
for b = 1:length(M.mapB)
  if b == 1 % Displacement
    res.bc{b}.ux = zeros(size(M.mapB{b}));
    res.bc{b}.uy = zeros(size(M.mapB{b}));
    if isfield(p, 'tau_inf')
      res.bc{b}.uz = sign(xF(M.mapB{b})) * p.V_p * t / 2 + (p.tau_inf./p.mu) .* xF(M.mapB{b});
    else
      res.bc{b}.uz = sign(xF(M.mapB{b})) * p.V_p * t / 2;
    end
  elseif b == 2 % Traction
    res.bc{b}.Tx = zeros(size(M.mapB{b}));
    res.bc{b}.Ty = zeros(size(M.mapB{b}));
    res.bc{b}.Tz = zeros(size(M.mapB{b}));
  elseif b == 3 % Vp on the fault
    res.bc{b}.dux = zeros(size(M.mapB{b}));
    res.bc{b}.duy = zeros(size(M.mapB{b}));
    res.bc{b}.duz = -sign(M.nxM(M.mapB{b})) * p.V_p * t;
  elseif b == 7 || b == 8 % Displacement jumps
    res.bc{b}.dux = zeros(size(M.mapB{b}));
    res.bc{b}.duy = zeros(size(M.mapB{b}));
    if isfield(M, 'fault_map')
      res.bc{b}.duz = -slip(M.fault_map{b-6});
    else
      res.bc{b}.duz = reshape(-slip, size(M.mapB{b}));
    end
  else
    if ~isempty(M.mapB{b})
      error('invalid bc')
    end
  end
end
