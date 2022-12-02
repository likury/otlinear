## Copyright (C) 2022 lucas
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} simplex (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: lucas <lucas@DESKTOP-GPL3NTT>
## Created: 2022-12-02

function retval = simplex (A, m, b, n)

  ##problema auxiliar

  aux_A = horzcat(A, eye(m))
  bind_aux = [n+1:n + m]

  c_aux = zeros(n + m, 1)
  c_aux(bind_aux) = 1
  c_aux_B = c_aux(bind_aux)
  aux_bfs = zeros(n + m, 1)
  aux_bfs(bind_aux,1) = b
  inv_B = eye(m)

  p = c_aux_B'*inv_B
  reduced_c = c_aux' - p*aux_A

  negative_index = find(reduced_c < 0)

  j = negative_index(1)

  u = inv_B*A(:, j)
  positive_index_u = find(u > 0)
  [theta, theta_index] = min(bsxfun(@rdivide, aux_bfs(bind_aux(positive_index_u)), u(positive_index_u)))

  l = theta_index
  aux_matrix = horzcat(inv_B, u)

  pivot_element = u(l)

  for (i = 1:m)
    if (i != l)
      multiple = - u(i)/u(l);
      new_row = aux_matrix(i, :) + multiple*aux_matrix(l, :)
      aux_matrix(i, :) = new_row
    endif
  endfor

  multiple = 1/u(l)
  aux_matrix(l, :) = multiple*aux_matrix(l, :)

  for (i = 1:m)
    aux_bfs(bind_aux(i)) = aux_bfs(bind_aux(i)) - theta*u(i)
  endfor
  aux_bfs(j) = theta
  bind_aux(j) = l






endfunction
