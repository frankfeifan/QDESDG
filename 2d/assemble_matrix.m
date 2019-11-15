function A = assemble_matrix(O, M, P, G)

lam = P.lam;
mu = P.mu;

% Storage for the matrix elements
NNZ = 5 * M.num_elm * M.Np^2;
[Ive, Jve] = meshgrid(1:M.Np);
Qve = 1:(M.Np)^2;

% Add in penalty term storage
NNZ = NNZ + 4 * 5 * M.num_elm * M.Np^2;

k = 0;

Ntotal = M.num_elm * M.Np;
FP{1} = 0 * M.Nq + (1:M.Nq);
FP{2} = 1 * M.Nq + (1:M.Nq);
FP{3} = 2 * M.Nq + (1:M.Nq);
FP{4} = 3 * M.Nq + (1:M.Nq);

MI{1} = M.face_interp(FP{1}, :);
MI{2} = M.face_interp(FP{2}, :);
MI{3} = M.face_interp(FP{3}, :);
MI{4} = M.face_interp(FP{4}, :);

PI{1,1} = M.face_interp(       FP{1} , :);
PI{2,1} = M.face_interp(       FP{2} , :);
PI{3,1} = M.face_interp(       FP{3} , :);
PI{4,1} = M.face_interp(       FP{4} , :);
PI{1,2} = M.face_interp(fliplr(FP{1}), :);
PI{2,2} = M.face_interp(fliplr(FP{2}), :);
PI{3,2} = M.face_interp(fliplr(FP{3}), :);
PI{4,2} = M.face_interp(fliplr(FP{4}), :);

V = zeros(NNZ, 1);
I = zeros(NNZ, 1);
J = zeros(NNZ, 1);
tic
for e = 1:M.num_elm
  if mod(e, 1000) == 0
    time = toc;
    time = ((M.num_elm - e) * (time / e));
    fprintf('\n');
    fprintf('Assembling element %d of %d\n',e, M.num_elm);
    fprintf('Estimated time remaining = %.2e sec (%.2e min)\n', time, time / 60);
  end
  Dx = diag(M.drdx(:,e)) * O.Dr + diag(M.dsdx(:,e)) * O.Ds;
  Dy = diag(M.drdy(:,e)) * O.Dr + diag(M.dsdy(:,e)) * O.Ds;

  MJ = O.M * diag(M.J(:,e));
  Axx = Dx' * MJ * (lam + 2 * mu) * Dx + Dy' * MJ * mu * Dy;
  Axy = Dy' * MJ * mu * Dx + Dx' * MJ * lam * Dy;
  Ayy = Dx' * MJ * mu * Dx + Dy' * MJ * (lam + 2 * mu) * Dy;
  Ayx = Dx' * MJ * mu * Dy + Dy' * MJ * lam * Dx;
  Azz = Dx' * MJ * mu * Dx + Dy' * MJ * mu * Dy;

  for f = 1:4
    ep = G.EToE(f,e);
    fp = G.EToF(f,e);
    op = G.EToO(f,e)+1;

    DxP = diag(M.drdx(:,ep)) * O.Dr + diag(M.dsdx(:,ep)) * O.Ds;
    DyP = diag(M.drdy(:,ep)) * O.Dr + diag(M.dsdy(:,ep)) * O.Ds;

    Bxx = zeros(M.Np, M.Np);
    Byx = zeros(M.Np, M.Np);
    Byy = zeros(M.Np, M.Np);
    Bxy = zeros(M.Np, M.Np);
    Bzz = zeros(M.Np, M.Np);

    % (u^{*} - u)
    if G.EToB(f,e) == 0 || G.EToB(f,e) == 7 || G.EToB(f,e) == 8 || G.EToB(f,e) == 3
      Mx = -MI{f}' * O.Mf(FP{f},FP{f}) * ...
            diag(M.sJ(FP{f},e) .* M.nxM(FP{f},e)) * MI{f} / 2;
      My = -MI{f}' * O.Mf(FP{f},FP{f}) * ...
            diag(M.sJ(FP{f},e) .* M.nyM(FP{f},e)) * MI{f} / 2;

      Px = -MI{f}' * O.Mf(FP{f},FP{f}) * ...
            diag(M.sJ(FP{f},e) .* M.nxM(FP{f},e)) * PI{fp, op} / 2;
      Py = -MI{f}' * O.Mf(FP{f},FP{f}) * ...
            diag(M.sJ(FP{f},e) .* M.nyM(FP{f},e)) * PI{fp, op} / 2;
    elseif M.bc(G.EToB(f,e)) == 1
      Mx = -MI{f}' * O.Mf(FP{f},FP{f}) * ...
            diag(M.sJ(FP{f},e) .* M.nxM(FP{f},e)) * MI{f};
      My = -MI{f}' * O.Mf(FP{f},FP{f}) * ...
            diag(M.sJ(FP{f},e) .* M.nyM(FP{f},e)) * MI{f};
      Px = zeros(M.Np, M.Np);
      Py = zeros(M.Np, M.Np);
    elseif M.bc(G.EToB(f,e)) == 2
      Mx = zeros(M.Np, M.Np);
      My = zeros(M.Np, M.Np);

      Px = zeros(M.Np, M.Np);
      Py = zeros(M.Np, M.Np);
    else
      error('invalid bc')
    end
    Axx = Axx + Dx' * (lam + 2 * mu) * Mx + Dy' * mu * My;
    Ayx = Ayx + Dx' * mu * My + Dy' * lam * Mx;
    Ayy = Ayy + Dx' * mu * Mx + Dy' * (lam + 2 * mu) * My;
    Axy = Axy + Dy' * mu * Mx + Dx' * lam * My;
    Azz = Azz + Dx' * mu * Mx + Dy' * mu * My;

    Bxx = Bxx - (Dx' * (lam + 2 * mu) * Px + Dy' * mu * Py);
    Byx = Byx - (Dx' * mu * Py + Dy' * lam * Px);
    Byy = Byy - (Dx' * mu * Px + Dy' * (lam + 2 * mu) * Py);
    Bxy = Bxy - (Dy' * mu * Px + Dx' * lam * Py);
    Bzz = Bzz - (Dx' * mu * Px + Dy' * mu * Py);

    % T^{*}
    % Since we know the operator is symmetric we can add the transpose on
    Axx = Axx + (Dx' * (lam + 2 * mu) * Mx + Dy' * mu * My)';
    Ayx = Ayx + (Dy' * mu * Mx + Dx' * lam * My)';
    Ayy = Ayy + (Dx' * mu * Mx + Dy' * (lam + 2 * mu) * My)';
    Axy = Axy + (Dx' * mu * My + Dy' * lam * Mx)';
    Azz = Azz + (Dx' * mu * Mx + Dy' * mu * My)';

    Bxx = Bxx + (DxP' * (lam + 2 * mu) * Px' + DyP' * mu * Py')';
    Byx = Byx + (DyP' * mu * Px' + DxP' * lam * Py')';
    Byy = Byy + (DxP' * mu * Px' + DyP' * (lam + 2 * mu) * Py')';
    Bxy = Bxy + (DxP' * mu * Py' + DyP' * lam * Px')';
    Bzz = Bzz + (DxP' * mu * Px' + DyP' * mu * Py')';

    % penalty term
    if G.EToB(f,e) == 0 || G.EToB(f,e) == 7 || G.EToB(f,e) == 8 || G.EToB(f,e) == 3
      Mxx = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nxM(FP{f},e).^2) * MI{f};
      Myy = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nyM(FP{f},e).^2) * MI{f};
      Mzz = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nzM(FP{f},e).^2) * MI{f};
      Myx = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nxM(FP{f},e) .* ...
                                              M.nyM(FP{f},e)) * MI{f};

      Pxx = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nxM(FP{f},e).^2) * PI{fp, op};
      Pyy = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nyM(FP{f},e).^2) * PI{fp, op};
      Pzz = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nzM(FP{f},e).^2) * PI{fp, op};
      Pyx = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nxM(FP{f},e) .* ...
                                              M.nyM(FP{f},e)) * PI{fp, op};
    elseif M.bc(G.EToB(f,e)) == 1
      Mxx = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nxM(FP{f},e).^2) * MI{f};
      Myy = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nyM(FP{f},e).^2) * MI{f};
      Mzz = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nzM(FP{f},e).^2) * MI{f};
      Myx = MI{f}' * O.Mf(FP{f},FP{f}) * diag(M.nxM(FP{f},e) .* ...
                                              M.nyM(FP{f},e)) * MI{f};
      Pxx = zeros(M.Np, M.Np);
      Pyx = zeros(M.Np, M.Np);
      Pyy = zeros(M.Np, M.Np);
      Pzz = zeros(M.Np, M.Np);
    elseif M.bc(G.EToB(f,e)) == 2
      Mxx = zeros(M.Np, M.Np);
      Myx = zeros(M.Np, M.Np);
      Myy = zeros(M.Np, M.Np);
      Mzz = zeros(M.Np, M.Np);

      Pxx = zeros(M.Np, M.Np);
      Pyx = zeros(M.Np, M.Np);
      Pyy = zeros(M.Np, M.Np);
      Pzz = zeros(M.Np, M.Np);
    else
      error('invalid bc')
    end
    Axx = Axx + P.beta * (lam + 2 * mu) * Mxx + P.beta * mu * Myy;
    Axy = Axy + P.beta * mu * Myx + P.beta * lam * Myx;
    Ayy = Ayy + P.beta * mu * Mxx + P.beta * (lam + 2 * mu) * Myy;
    Ayx = Ayx + P.beta * mu * Myx + P.beta * lam * Myx;
    Azz = Azz + P.beta * mu * Mxx + P.beta * mu * Myy;

    Bxx = Bxx - P.beta * (lam + 2 * mu) * Pxx - P.beta * mu * Pyy;
    Bxy = Bxy - P.beta * mu * Pyx - P.beta * lam * Pyx;
    Byy = Byy - P.beta * mu * Pxx - P.beta * (lam + 2 * mu) * Pyy;
    Byx = Byx - P.beta * mu * Pyx - P.beta * lam * Pyx;
    Bzz = Bzz - P.beta * mu * Pxx - P.beta * mu * Pyy;

    V(k + Qve, 1) = Bxx(:);
    I(k + Qve, 1) = Ive(:) + (ep - 1) * M.Np;
    J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np;
    k = k + Qve(end);

    V(k + Qve, 1) = Byy(:);
    I(k + Qve, 1) = Ive(:) + (ep - 1) * M.Np + Ntotal;
    J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np + Ntotal;
    k = k + Qve(end);

    V(k + Qve, 1) = Bzz(:);
    I(k + Qve, 1) = Ive(:) + (ep - 1) * M.Np + 2 * Ntotal;
    J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np + 2 * Ntotal;
    k = k + Qve(end);

    V(k + Qve, 1) = Byx(:);
    I(k + Qve, 1) = Ive(:) + (ep - 1) * M.Np;
    J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np + Ntotal;
    k = k + Qve(end);

    V(k + Qve, 1) = Bxy(:);
    I(k + Qve, 1) = Ive(:) + (ep - 1) * M.Np + Ntotal;
    J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np;
    k = k + Qve(end);
  end

  V(k + Qve, 1) = Axx(:);
  I(k + Qve, 1) = Ive(:) + (e - 1) * M.Np;
  J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np;
  k = k + Qve(end);

  V(k + Qve, 1) = Ayy(:);
  I(k + Qve, 1) = Ive(:) + (e - 1) * M.Np + Ntotal;
  J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np + Ntotal;
  k = k + Qve(end);

  V(k + Qve, 1) = Azz(:);
  I(k + Qve, 1) = Ive(:) + (e - 1) * M.Np + 2 * Ntotal;
  J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np + 2 * Ntotal;
  k = k + Qve(end);

  V(k + Qve, 1) = Ayx(:);
  I(k + Qve, 1) = Ive(:) + (e - 1) * M.Np;
  J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np + Ntotal;
  k = k + Qve(end);

  V(k + Qve, 1) = Axy(:);
  I(k + Qve, 1) = Ive(:) + (e - 1) * M.Np + Ntotal;
  J(k + Qve, 1) = Jve(:) + (e - 1) * M.Np;
  k = k + Qve(end);
end

assert(k == NNZ)

fprintf('\n');
disp('Building matrix...')
A = sparse(I, J, V);

