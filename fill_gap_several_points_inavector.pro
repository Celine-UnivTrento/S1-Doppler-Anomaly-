function fill_gap_several_points_inavector, vector

vector_out = vector
val_def = 0
p = where(vector ne 0)
q = abs(shift(p,-1)-p)
p_disc = where(q gt 1)

n_disc = n_elements(p_disc)-1
n_vector = n_elements(vector)

for i_disc =0, n_disc -1 do begin 
  ;print, vector_out
  ;print, p[p_disc[i_disc]+1],p[p_disc[i_disc]]
  ;print, vector[p[p_disc[i_disc]+1]], vector[p[p_disc[i_disc]]]
  
  if q[p_disc[i_disc]] lt n_vector/4 then begin 
  
 ;   print, 'intervalle pour fit 1D : ', p[p_disc[i_disc]+1],p[p_disc[i_disc]], abs( p[p_disc[i_disc]+1]- p[p_disc[i_disc]])
    
;    if abs( p[p_disc[i_disc]+1]- p[p_disc[i_disc]]) gt 100 then stop
    
    poly_deg = poly_fit([p[p_disc[i_disc]+1],p[p_disc[i_disc]] ], [vector[p[p_disc[i_disc]+1]], vector[p[p_disc[i_disc]]]], 1)
    long_disc = p[p_disc[i_disc]+1]-p[p_disc[i_disc]]+1
    vector_out[p[p_disc[i_disc]]:p[p_disc[i_disc]+1]] = poly_deg[1]*(findgen(long_disc)+p[p_disc[i_disc]])+poly_deg[0]
  endif; else print, 'intervalle de (sur total) : ', q[p_disc[i_disc]], n_vector
endfor

return, vector_out


end