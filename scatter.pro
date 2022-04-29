pro scatter, x,y,z ,equal=equal, xrange = xrange, yrange = yrange, zrange=zrange, $
    cltable=cltable, xtitle=xtitle, ytitle=ytitle, title=title, $
    cbar=cbar, cbtitle = cbtitle, cbhoriz = cbhoriz, cbformat = cbformat, $
    psym = psym, symsize = symsize, $
    map=map, world=world, $;inv_ctable=inv_ctable, $
    background=background, colorline=colorline,_EXTRA=e,$
    contour_fill=contour_fill, is_zero=is_zero, isotropic=isotropic
    
    
  ;; dec. 2011
  ;;
  ;;   equivalent of matlab scatter function
  ;;
  ;;  world keyword display data on a world map (x=lon  y=lat)
    
  p = where(finite(z))
  if n_elements(p) lt 2 then return
  p = -1
  
  IF n_elements(xtitle) eq 0 THEN xtitle = ''
  IF n_elements(ytitle) eq 0 THEN ytitle = ''
  IF n_elements(title) eq 0 THEN title = ''
  
  IF n_elements(cltable) EQ 0 THEN cltable = 13 else cltable = cltable
  IF cltable EQ -1 THEN cltable = 0
  if n_elements(background) eq 0 then background=0L
  if n_elements(colorline) eq 0 then colorline=!d.n_colors-1 - background ;else print, colorline
  IF n_elements(cbar) eq 0 THEN cbar = 0
  IF n_elements(cbtitle) eq 0 THEN cbtitle = ''
  IF N_ELEMENTS(cbhoriz) EQ 0 THEN cbhoriz = 0
  IF N_ELEMENTS(cbformat) EQ 0 THEN cbformat = '(F6.2)'
 ; IF N_ELEMENTS(inv_ctable) EQ 0 THEN inv_ctable = 0
  
  mxx = max(x, /nan)
  mix = min(x, /nan)
  mxy = max(y, /nan)
  miy = min(y, /nan)
  
  ;ratio = (mxy-miy)/(mxx-mix)
  ;ysize = 1000.
  ;xsize = fix(float(ysize)/ratio)
  
  ;print, xsize, ysize
  ;window, /free, xs=xsize, ys=ysize
  
  if not(keyword_set(zrange)) then begin
    mxz = max(z, /nan)
    miz = min(z, /nan)
  endif else begin
    mxz = zrange[1]
    miz = zrange[0]
  endelse
  
  if miz eq mxz then begin
    miz = miz-0.10*miz-1.0e-15
    mxz = mxz+0.10*mxz+1.0e-15
  endif
  
  ;  STOP
  
  IF (!D.Flags AND 1) NE 0 THEN scalablePixels = 1 ELSE scalablePixels = 0
  
  if not(keyword_set(xrange)) then xrange = [mix,mxx]
  if not(keyword_set(yrange)) then yrange = [miy,mxy]
  
  
  ; define plotting symbols
  if not(keyword_set(psym)) then psym = 3
  if not(keyword_set(symsize)) then symsize = 1
  
  ; set background color
  !p.background = background
  
  IF cbar NE 0 THEN BEGIN
    position_old = !P.position
    ;    if cbhoriz then !P.position = [.1, .2, 0.9, 0.9] $
    ;    else !P.position = [0.1, 0.1, .85, 0.9]
    if cbhoriz then !P.position = [.1, .2, .95, .95] $
    else !P.position = [.1, .1, .85, .95]
  ENDIF
  
  if (not(keyword_set(map)) and not(keyword_set(world))) then $
    plot, [mix,mxx],[miy,mxy], xrange = xrange, yrange = yrange ,psym = psym, $
    xstyle=1, ystyle=1,color=colorline,$
    xtitle=xtitle, ytitle=ytitle, title=title, _EXTRA=e, isotropic=isotropic, /nodata
    
  if keyword_set(world) then $
    MAP_SET, 0, 0, /CYLINDRICAL, /ISOTROPIC,limit=[-80,-180,80,180], $
    title=title, color=colorline,_EXTRA=e
    
  ;load color table
  tvlct, r_orig, g_orig, b_orig, /get
  
  loadct, cltable, /silent
  if not(scalablePixels) then device, decomposed = 0
  ;  loadct, cltable, /silent
  ncolors = !d.table_size
  
  if keyword_set(equal) then begin
    c=hist_equal(z)
 ;   if inv_ctable then c = rotate(c, 2)
  endif   else begin
  
    c =  ncolors*(z-miz)/(mxz-miz)
    ind_tmp = where(z le miz,cnt)
    if cnt gt 0 then c[ind_tmp]=0l
    ind_tmp=-1
    ind_tmp = where(z ge mxz,cnt)
    if cnt gt 0 then c[ind_tmp]= ncolors-1
    ind_tmp = -1
    
 ;   if inv_ctable then c = rotate(c, 2)
  endelse
  
  ; contour
  ; règlage palette de couleur
  ; 21 == nombre de couleurs de la palette
  IF (keyword_set(contour_fill)) THEN begin
  
    IF (contour_fill EQ 1) THEN step=(mxz-miz)/21. $
    ELSE step=contour_fill
    
    lbeg = floor(miz/step)*step
    lend = ceil(mxz/step)*step
    nlevels = (lend-lbeg)/step+1
    lev = findgen(nlevels)*step+lbeg
    ; ncolors = !d.table_size
    c_col = fix((findgen(nlevels+2))*(ncolors-1)/(nlevels+1))
    
    ;; set background color
    !p.background = background
    
    lev_min = min(lev, /nan) - (lev[1]-lev[0])
    
    contour, z, x, y, $
      /overplot, $
      xstyle = 1, ystyle = 1, xminor = -1, yminor = -1,$
      /cell_fill,background=0L, $
      c_color=[0L, c_col], /irregular,$
      levels=[lev_min,lev] , /closed ;
      
    ; fait un contour autour de la valeur zero
    IF (keyword_set(is_zero)) THEN $
      contour, z,x,y, $
      /overplot, levels=[0], c_labels=[0],$
      /irregular,background=0L,color=-255
      
  ;; Rajoute les points par dessus pour ce qui sont sur la zone de représentation
  ;p = where(z lt mxz and z gt miz )
  ;plots, x[p], y[p], psym = 4, color = c[p] , $
  ;       thick=2,noclip=0
  ;;    plots, x[p], y[p], psym = 4, color = c[p] , $
  ;;          symsize = symsize,_EXTRA=e,noclip=0, thick=2
      
  ENDIF else begin
    ; plot data
    ;p = where(z eq !VALUES.F_NAN)
    ;p = where(finite(z, /nan))
    ;IF p[0] eq -1 THEN $
    ;   plots, x, y, psym = psym, color = c , $
    ;          symsize = symsize,_EXTRA=e,noclip=0
    ;p =-1
    p = where(finite(z) and z ge miz and z le mxz)
    IF p[0] ne -1 THEN $
      plots, x[p], y[p], psym = psym, color = c[p] , $
      symsize = symsize,_EXTRA=e,noclip=0 else print, 'Que des valeurs non finies !!!!!!!!!!!!!!!!!!!!'
    p=-1
  ENDELSE
  
  ;load color table
  if not(scalablePixels) then device, decomposed =1
  tvlct, r_orig, g_orig, b_orig
  ;!p.background = background_old
  
  if keyword_set(world) then begin
    MAP_CONTINENTS,/fill,color=128
    MAP_grid,/box_axes,color=colorline
  endif
  
  if keyword_set(cbar) then begin
    if cbhoriz then begin
      cbarpos = [!X.window(0), .05, !X.window(1),.08]
      loadct, cltable
      colorbar_obs, /bottom, range=[miz,mxz], $
        color=colorline,  format=cbformat, $
        position = cbarpos, title =cbtitle;, cltable = cltable;, _EXTRA=e
    endif else begin
      cbarpos = [.87, !Y.window(0), .90, !Y.window(1)]
      loadct, cltable
      colorbar_obs, /vert, /right, range=[miz,mxz], $
        color=colorline,  format=cbformat,$
        position = cbarpos, title = cbtitle;, cltable = cltable;, _EXTRA=e,cltable = cltable,
    endelse
    ;    colorbar, position = cbpos, range = [miz,mxz],color=colorline, cltable=cltable, $
    ;      DIVISIONS=10,title=cbtit, EXTRA=e
    !P.position = position_old
  endif
  
END
