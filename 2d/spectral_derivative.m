function D = spectral_derivative(x, v)
  n = length(x);
  D = zeros(n,n);

  if nargin < 2
    w = lagrange_weights(x);
  end

  for k = 1:n
    for j = 1:n
      if k == j
        for l = 1:n
          if l ~= k
            D(j, k) = D(j, k) + 1 / (x(k) - x(l));
          end
        end
      else
        D(j, k) = (v(k) / v(j)) / (x(j) - x(k));
      end
    end
  end
end
