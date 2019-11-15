function [] = station_data(yfault, plot_fault, pre, ystation)

z = yfault;
J = find(z(plot_fault(1,:)) > z(plot_fault(end,:)));
plot_fault(:,J) = plot_fault(end:-1:1,J);
[~,J] = sort(z(plot_fault(1,:)));
plot_fault = plot_fault(:,J);

ind = find_indices(pre);  %find indices where max(V) > 1e-3 m/s

%%%%%%% READ DATA %%%%%%%%%%%%%%%%%%%%%%%

W = SaveStreamData('Read', [pre,'w.dat']);
V = SaveStreamData('Read', [pre,'Vel.dat']);
T = SaveStreamData('Read', [pre,'Time.dat']);

N = min([length(T), size(W, 2), size(V, 2)]);
W = W(:,1:N);
V = V(:,1:N);
U = W(1:end/2,:);

interval = [5 * 31556926 1];

TI = [];
for i = 1:2:length(ind)-1
  if(length(ind) < i+1)
    break
  end
  TI = [TI,T(ind(i)):interval(1):T(ind(i+1))];

  if(length(ind) < i+2)
    break
  end
  TI = [TI,T(ind(i+1)):interval(2):T(ind(i+2))];
end

for s = 1:length(ystation)
  e = min(find(z(plot_fault(1,:)) <= ystation(s) & ...
          z(plot_fault(end,:)) >= ystation(s)));

  J = plot_fault(:,e);
  I = lagrange_interpolation_matrix(z(J), ystation(s));

  U_s{s} = interp1(T, I * U(J,:), TI);
end
Vm = interp1(T, max(abs(V)), TI);

fid = fopen(sprintf('%s.slip', pre(1:end-1)), 'w+');

strs = strsplit(pre, '/');
fprintf(fid, '# Problem:  %s\n', strs{end}(1:end-1));
fprintf(fid, '# Author:   Jeremy Kozdon\n');
fprintf(fid, '# Method:   Symmetric Interior Penalty DG\n');
fprintf(fid, '# slip data\n');
fprintf(fid, '# aseismic step = %e years\n', interval(1)/31556926);
fprintf(fid, '# seismic  step = %e seconds\n', interval(2));
fprintf(fid, 't max_slip_rate slip\n');

fprintf(fid, '%.16e %.16e', 0, 0);
for s = 1:length(ystation)
  fprintf(fid, ' %.16e', 1000*abs(ystation(s)));
end
fprintf(fid, '\n');


for n = 1:length(TI)
  fprintf(fid, '%.16e %.16e', TI(n), log10(Vm(n)));
  for s = 1:length(ystation)
    fprintf(fid, ' %.16e', U_s{s}(n));
  end
  fprintf(fid, '\n');
end
fclose(fid);
