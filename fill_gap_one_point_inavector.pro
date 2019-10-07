function fill_gap_one_point_inavector, vector

vector_out = vector

val_def = 0.

p = where(vector ne val_def)
num_hold = -1 
trans = p[0]

for i = 0L, n_elements(p)-1 do begin
  if trans ne p[i] then num_hold = p[i]
  if num_hold ne -1 then begin 
    diff = num_hold -trans
    if diff eq 1 then vector_out[p[i]-1] = (vector[p[i]]+vector[p[i]-2])/2.
    num_hold = -1
  endif
  trans = p[i] +1  
endfor

return, vector_out

end