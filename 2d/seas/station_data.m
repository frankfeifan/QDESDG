function [] = station_data(yfault, plot_fault, pre, ystation)
disp(pre)

z = yfault;
J = find(z(plot_fault(1,:)) > z(plot_fault(end,:)));
plot_fault(:,J) = plot_fault(end:-1:1,J);
[~,J] = sort(z(plot_fault(1,:)));
plot_fault = plot_fault(:,J);


%%%%%%% READ DATA %%%%%%%%%%%%%%%%%%%%%%%

W = SaveStreamData('Read', [pre,'w.dat']);
V = SaveStreamData('Read', [pre,'Vel.dat']);
S = SaveStreamData('Read', [pre,'tau.dat']);
T = SaveStreamData('Read', [pre,'Time.dat']);
disp(T(end) / 31556926)

N = min([length(T), size(W, 2), size(V, 2), size(S, 2)]);
W = W(:,1:N);
V = V(:,1:N);
rho = 2.670;
cs = 3.464;
mu = cs^2 * rho;
lam = cs^2 * rho;
eta = mu./(2*cs); %radiation damping parameter
S = S(:,1:N) - eta * V;
T = T(1:N);

b0 = 0.015;
Dc = 0.008;
v0 = 1e-6;
f0 = 0.6;

U = W(1:end/2,:);
psi = W(end/2+1:end,:);
Q = Dc * exp((psi - f0) / b0) / v0;

ind = find_indices(pre, V);  %find indices where max(V) > 1e-3 m/s

interval = [0.1 * 31556926 0.1];

TI = [];
for i = 1:2:length(ind)-1
  if(length(ind) < i+1)
    break
  end
  TI = [TI,T(ind(i)):interval(1):T(ind(i+1))];

  if(length(ind) < i+2)
    break
  end
  TI = [TI,T(ind(i+1)+1):interval(2):T(ind(i+2)-1)];
end

for s = 1:length(ystation)
  e = min(find(z(plot_fault(1,:)) <= ystation(s) & ...
          z(plot_fault(end,:)) >= ystation(s)));

  J = plot_fault(:,e);
  I = lagrange_interpolation_matrix(z(J), ystation(s));

  U_s = interp1(T, I * U(J,:), TI);
  V_s = interp1(T, I * V(J,:), TI);
  S_s = interp1(T, I * S(J,:), TI);
  Q_s = interp1(T, I * Q(J,:), TI);
  fid = fopen(sprintf('%sz_fltst_dp%03.0f.station', pre, 10*abs(ystation(s))), 'w+');

  strs = strsplit(pre, '/');
  fprintf(fid, '# Problem:  %s\n', strs{end}(1:end-1));
  fprintf(fid, '# Author:   Jeremy Kozdon\n');
  fprintf(fid, '# Method:   Symmetric Interior Penalty DG\n');
  fprintf(fid, '# Station:  %e\n', ystation(s));
  fprintf(fid, '# aseismic step = %e years\n', interval(1)/31556926);
  fprintf(fid, '# seismic  step = %e seconds\n', interval(2));
  fprintf(fid, 't slip slip_rate shear_stress state\n');

  for n = 2:length(TI)
    fprintf(fid, '%.16e %e %e %e %e\n', ...
            TI(n), U_s(n), log10(V_s(n)), S_s(n), log10(Q_s(n)));
  end

  fclose(fid);
end
