function [x, y, EToV, EToE, EToF, EToO, EToB, EToT] = inp_read_2d(filename)

  fid = fopen(filename,'r');

  default_bc = 1;

  if(fid < 0)
    error(['Cannot open: ', filename]);
  end

  % {{{ Read the nodes
  node_start = seek_to_substring(fid, 'NSET=ALLNODES', true);

  % First count the number of nodes
  num_nodes = 0;
  while true
    line = fgetl(fid);
    if isempty(sscanf(line, '%d, %f, %f, %f'))
      break
    end
    num_nodes = num_nodes + 1;
  end

  % Read the nodes
  fseek(fid, node_start, 'bof');
  x = nan(num_nodes, 1);
  y = nan(num_nodes, 1);
  z = nan(num_nodes, 1);
  n = 0;
  while true
    line = fgetl(fid);
    A = sscanf(line, '%d, %f, %f, %f');
    if isempty(A)
      break
    end
    n = n + 1;
    x(n) = A(2);
    y(n) = A(3);
    z(n) = A(4);
  end
  assert(n == num_nodes);
  % }}}

  % go back to beginning of the file
  frewind(fid);

  % {{{ Read the elements
  elm_start = seek_to_substring(fid, '*ELEMENT', true);

  % First count the number of elements
  num_elm = 0;
  while true
    line = fgetl(fid);
    if isempty(sscanf(line, '%d, %d, %d, %d, %d'))
      break
    end
    num_elm = num_elm + 1;
  end

  % Read the nodes
  fseek(fid, elm_start, 'bof');
  EToV = zeros(4, num_elm);
  n = 0;
  while true
    line = fgetl(fid);
    A = sscanf(line, '%d, %d, %d, %d, %d');
    if isempty(A)
      break
    end
    n = n + 1;
    EToV([1 2 4 3], n) = A(2:5);
  end
  assert(n == num_elm);
  % }}}

  % {{{ Determine the connectivity
  EToE = zeros(4, num_elm);
  EToF = zeros(4, num_elm);
  EToO = zeros(4, num_elm);
  EToB = zeros(4, num_elm);
  EToT = zeros(1, num_elm);
  for k = 1:num_elm
    [I1, J1] = find(EToV(1,k) == EToV);
    [I2, J2] = find(EToV(2,k) == EToV);
    [I3, J3] = find(EToV(3,k) == EToV);
    [I4, J4] = find(EToV(4,k) == EToV);

    [J1, v] = setdiff(J1, k); I1 = I1(v);
    [J2, v] = setdiff(J2, k); I2 = I2(v);
    [J3, v] = setdiff(J3, k); I3 = I3(v);
    [J4, v] = setdiff(J4, k); I4 = I4(v);

    EToE(:, k) = k;
    EToF(:, k) = 1:4;

    if ~isempty(intersect(J1, J3))
      [EToE(1, k), v, w] = intersect(J1, J3);
      [EToF(1, k), EToO(1,k)] = vert_to_face(I1(v), I3(w));
    else
      EToB(1, k) = default_bc;
    end
    if ~isempty(intersect(J2, J4))
      [EToE(2, k), v, w] = intersect(J2, J4);
      [EToF(2, k), EToO(2, k)] = vert_to_face(I2(v), I4(w));
    else
      EToB(2, k) = default_bc;
    end
    if ~isempty(intersect(J1, J2))
      [EToE(3, k), v, w] = intersect(J1, J2);
      [EToF(3, k), EToO(3, k)] = vert_to_face(I1(v), I2(w));
    else
      EToB(3, k) = default_bc;
    end
    if ~isempty(intersect(J3, J4))
      [EToE(4, k), v, w] = intersect(J3, J4);
      [EToF(4, k), EToO(4, k)] = vert_to_face(I3(v), I4(w));
    else
      EToB(4, k) = default_bc;
    end
  end
  % }}}

  % go back to beginning of the file
  frewind(fid);

  % {{{ Read side set info
  [start, line] = seek_to_substring(fid, '*ELSET', false);
  inp_to_f = [3,  2, 4, 1];
  while ischar(line)
    SSE = str2double( regexp(line,'.*(\d+)_E(\d+)','tokens','once') );
    if ~isempty(SSE)
      f = inp_to_f(SSE(2));
      while true
        line = fgetl(fid);
        if ~isempty(findstr(line, '*'))
          break
        end

        for k = str2num(line)
          EToB(f, k) = SSE(1);
        end
      end
    else
      SSE = str2double( regexp(line,'.*SS(\d+).*','tokens','once') );
      if ~isempty(SSE)
        while true
          line = fgetl(fid);
          if ~isempty(findstr(line, '*'))
            break
          end
          for k = str2num(line)
            EToT(k) = SSE(1);
          end
        end
      end
    end

    if isempty(findstr(line, '*ELSET'))
      [start, line] = seek_to_substring(fid, '*ELSET', false);
    end
  end
  % }}}

  fclose(fid);

end

function [loc, line] = seek_to_substring(fid, key, err)
  line = fgetl(fid);
  while isempty(findstr(line, key))
    line = fgetl(fid);
    if ~ischar(line) 
      if err
        error([''' Did not find ''', key, ''''])
      end
      break
    end
  end

  loc = ftell(fid);
end

function [f, o] = vert_to_face(vin, win)
  f = -1;
  o = (vin > win);
  v = min(vin, win);
  w = max(vin, win);
  if v == 1 && w == 3
    f = 1;
  elseif v == 2 && w == 4
    f = 2;
  elseif v == 1 && w == 2
    f = 3;
  elseif v == 3 && w == 4
    f = 4;
  end

  if i < 0
    error(sprintf('invalid vertex numbers: %d:%d', v, w))
  end
end

