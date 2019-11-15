function [O, M, G] = init_mesh(N, filename, num_ref, ...
  x_grid, y_grid, EToV, EToE, EToF, EToO, EToB, EToT)


if(N > 22)
  error('cannot handle N > 22 (need more Chebfun routines)')
end

% Number of quadrature points
Nq = N + 1;
Np = Nq^2;

% Read in the mesh
if nargin < 4
  [x_grid, y_grid, EToV, EToE, EToF, EToO, EToB, EToT] = inp_read_2d(filename);
end
for r = 1:num_ref
  [x_grid, y_grid, EToV, EToE, EToF, EToO, EToB, EToT] = ...
    refine(x_grid, y_grid, EToV, EToE, EToF, EToO, EToB, EToT);
end

num_elm = size(EToV, 2);

% Determine the quadrature rule
[r_lgl, w_lgl, v_lgl] = lobpts(Nq);

% fgure out the mesh
r = kron(ones(Nq, 1), r_lgl);
s = kron(r_lgl, ones(Nq, 1));

x = (((1 - r) .* (1 - s)) * x_grid(EToV(1, :))' + ...
     ((1 + r) .* (1 - s)) * x_grid(EToV(2, :))' + ...
     ((1 - r) .* (1 + s)) * x_grid(EToV(3, :))' + ...
     ((1 + r) .* (1 + s)) * x_grid(EToV(4, :))') / 4;

y = (((1 - r) .* (1 - s)) * y_grid(EToV(1, :))' + ...
     ((1 + r) .* (1 - s)) * y_grid(EToV(2, :))' + ...
     ((1 - r) .* (1 + s)) * y_grid(EToV(3, :))' + ...
     ((1 + r) .* (1 + s)) * y_grid(EToV(4, :))') / 4;

%{
% Plot the mesh
plot(x, y, '*')
hold on
plot(x_grid(EToV([1, 2, 4, 3, 1], :)), y_grid(EToV([1, 2, 4, 3, 1], :)), 'k')
hold off
axis image
%}

% Compute face maps
fmap = zeros(Nq, 4);
fmap(:, 1) = (1:Nq:Np);
fmap(:, 2) = (Nq:Nq:Np);
fmap(:, 3) = 1:Nq;
fmap(:, 4) = N * Nq + (1:Nq);

face_interp = [kron(speye(Nq), sparse(1,  1, 1, 1, Nq));...
               kron(speye(Nq), sparse(1, Nq, 1, 1, Nq));...
               kron(sparse(1,  1, 1, 1, Nq), speye(Nq));...
               kron(sparse(1, Nq, 1, 1, Nq), speye(Nq))];

Nfaces = 4;

mapM = reshape(1:Nfaces * Nq * num_elm, Nfaces * Nq, num_elm);

mapP = mapM;
for eM = 1:num_elm
  for fM = 1:Nfaces
    eP = EToE(fM, eM);
    fP = EToF(fM, eM);
    oP = EToO(fM, eM);
    if oP
      mapP((fM - 1) * Nq + (1:Nq), eM) = mapM((fP - 1) * Nq + (Nq:-1:1), eP);
    else
      mapP((fM - 1) * Nq + (1:Nq), eM) = mapM((fP - 1) * Nq + (1:+1:Nq), eP);
    end
  end
end
xf = face_interp * x;
yf = face_interp * y;
assert(max(max(abs(xf(mapM) - xf(mapP)))) == 0)
assert(max(max(abs(yf(mapM) - yf(mapP)))) == 0)

% Compute the boundary maps
bc = 1:max(max(unique(EToB)));
for b = bc(:)'
  [f, e] = find(EToB == b);
  numbc = length(f);
  mapB{b} = zeros(Nq, numbc);
  for k = 1:numbc
    mapB{b}(:, k) = (1:Nq) + Nq * ((f(k) - 1) + Nfaces * (e(k) - 1));
  end
  if numbc > 0
    assert(max(max(abs(xf(mapM(mapB{b})) - xf(mapP(mapB{b}))))) == 0)
    assert(max(max(abs(yf(mapM(mapB{b})) - yf(mapP(mapB{b}))))) == 0)
  end
end

% 1-D derivative operator
D_lgl = spectral_derivative(r_lgl, v_lgl);

Dr = kron(speye(Nq), D_lgl);
Ds = kron(D_lgl, speye(Nq));

% Compute metric terms
dxdr = Dr * x;
dydr = Dr * y;
dxds = Ds * x;
dyds = Ds * y;
J = dxdr .* dyds - dxds .* dydr;
Jsign = face_interp * (2 * (J > 0) - 1);
J = abs(J);
drdx =  dyds ./ J;
drdy = -dxds ./ J;
dsdx = -dydr ./ J;
dsdy =  dxdr ./ J;

% Compute the face normals
nx = Jsign .* (kron(speye(Nfaces), D_lgl) * yf);
ny = Jsign .* (kron(speye(Nfaces), D_lgl) * xf);
nz = Jsign .* (zeros(size(nx)));
sJ = sqrt(nx.^2 + ny.^2);
nx = nx ./ sJ;
ny = ny ./ sJ;
nx(0 * Nq + (1:Nq),:) = -nx(0 * Nq + (1:Nq),:);
ny(1 * Nq + (1:Nq),:) = -ny(1 * Nq + (1:Nq),:);
ny(2 * Nq + (1:Nq),:) = -ny(2 * Nq + (1:Nq),:);
nx(3 * Nq + (1:Nq),:) = -nx(3 * Nq + (1:Nq),:);

%{
for k = 1:num_elm
  plot(x_grid(EToV([1, 2, 4, 3, 1], :)), y_grid(EToV([1, 2, 4, 3, 1], :)), 'k')
  hold on
  quiver(xf(:,k), yf(:,k), nx(:,k), ny(:,k))
  hold off
  axis image
  pause
end
plot(x_grid(EToV([1, 2, 4, 3, 1], :)), y_grid(EToV([1, 2, 4, 3, 1], :)), 'k')
hold on
quiver(xf(:,:), yf(:,:), nx(:,:), ny(:,:))
hold off
axis image
%}

O.Dr = Dr;
O.Ds = Ds;
O.M  = diag(sparse(kron(w_lgl     , w_lgl     )));
O.MI = diag(sparse(kron(1 ./ w_lgl, 1 ./ w_lgl)));
O.Mf = diag(sparse(kron(ones(1, 4), w_lgl     )));
O.r_lgl = r_lgl;
O.v_lgl = v_lgl;
O.w_lgl = w_lgl;

M.num_elm = num_elm;
M.N = N;
M.Nq = Nq;
M.Np = Np;
M.drdx = drdx;
M.dsdx = dsdx;
M.drdy = drdy;
M.dsdy = dsdy;
M.J = J;
M.JI = 1./J;
M.face_interp = face_interp;
M.mapM = mapM;
M.mapP = mapP;
M.mapB = mapB;
M.nxM = nx;
M.nyM = ny;
M.nzM = nz;
M.nxP = -nx;
M.nyP = -ny;
M.nzP = -nz;
M.sJ = sJ;
M.x = x;
M.y = y;
M.bc = bc;

G.x    = x_grid;
G.y    = y_grid;
G.EToV = EToV;
G.EToE = EToE;
G.EToF = EToF;
G.EToO = EToO;
G.EToB = EToB;
G.EToT = EToT;


% Plotting stuff
T1 = [0, 1, Nq]';
Ta = kron(kron(ones(N,1),(1:N)') + kron(Nq*(0:N-1)',ones(N,1)), ones(3,1)) ...
   + kron(ones(N^2,1),T1);
T2 = [1, Nq, Nq+1]';
Tb = kron(kron(ones(N,1),(1:N)') + kron(Nq*(0:N-1)',ones(N,1)), ones(3,1)) ...
   + kron(ones(N^2,1),T2);
Tc = [Ta; Tb];
T = kron(ones(num_elm, 1), Tc) + kron(Np*((1:num_elm)'-1), ones(size(Tc)));
M.T = reshape(T, 3, length(T)/3)';
