pro pointille, x, y, col=col, th=th

IF n_elements(col) EQ 0 THEN col=-1
IF n_elements(th) EQ 0 THEN th=1

min_y = min(y, /nan);round(min(y))
min_x = min(x, /nan);round(min(x))
max_y = max(y, /nan);round(max(y))
max_x = max(x, /nan);round(max(x))

nbre_lig = 6.

dx = ((max_x-min_x)/nbre_lig) 

nbre_col = 6.
dy = ((max_y-min_y)/nbre_col) 

FOR i = 0, nbre_lig do oplot,[-1,1]*0+min_x+i*dx,[min_y, max_y],col=col, lin=3, thick=th
FOR j = 0, nbre_col do oplot,[min_x,max_x],[-1,1]*0+min_y+j*dy,col=col, lin=3, thick=th

END
