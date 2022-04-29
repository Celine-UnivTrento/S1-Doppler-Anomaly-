FUNCTION extrait_entre_inegalite, chaine

str =  stregex(chaine, '>[^*]+<', /FOLD_CASE, /extract)

return,  strmid(str, 1, strlen(str)-2)

end
