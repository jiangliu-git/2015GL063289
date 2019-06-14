pro main_test
;;; test stuff
;nx=20
;ny=10
;posx=randomu(s,500)
;posy=randomu(s,500)
;value=posx^2+posy^2
;;field=interp_cic(value,posx*nx,nx,posy*ny,ny,/average)
;field=interp_ngp(value,posx*nx,nx,posy*ny,ny,/average)
;;field=interp_tsc(value,posx*nx,nx,posy*ny,ny,/average)
;surface,field,/lego

folder_cdf = 'ctcs_values'
filename_p = folder_cdf+'/'+'pressure_isotropic.cdf'
cdf_var_show, filename_p, /ncdf

;cdfid = ncdf_open('ct_currents.cdf')
;ncdf_varget, cdfid, 'jCrossEq', jCrossEq
;ncdf_close, cdfid
;print, jCrossEq
stop
end
