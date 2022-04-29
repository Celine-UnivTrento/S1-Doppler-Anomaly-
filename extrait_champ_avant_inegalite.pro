FUNCTION extrait_champ_avant_inegalite, chaine

pos_2pts = strpos(chaine, ':')
pos_ine =  strpos(chaine, '>')

return,  strmid(chaine, pos_2pts+1, pos_ine-pos_2pts-1)

end
