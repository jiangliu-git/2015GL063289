pro main_dy_plot
;;;;;;; plot dy of the satellite orbit
;; this is for single DFBs, for multi DFBs, go to main_multi_plot

@folders

;list_suf = '_fgs'
;list_suf = ''
list_suf = '_earthward_df_binxbout'

;;; limit the separation
dx_cri = 3.
dy_cri = 100.
dz_cri = 2.

dy_range = [0., 2.]

dxyz_all = datain_simple(save_folder+'/dxyz'+list_suf+'.dat', dim=3, type='double')
xyz_all = datain_simple(save_folder+'/xyz_ave'+list_suf+'.dat', dim=3, type='double')
if strmatch(list_suf, '*_df_*') then n_all = datain_simple(save_folder+'/n_sideall'+list_suf+'.dat', dim=3, type='double')

good_fin = finite(dxyz_all(1,*))
good_dx = abs(dxyz_all(0,*)) le dx_cri
good_dy = abs(dxyz_all(1,*)) le dy_cri
good_dz = abs(dxyz_all(2,*)) le dz_cri

base_req = good_fin and good_dx and good_dy and good_dz

;;;;;;;;; dy as a function of y
;;yrange = [-12, 12]
;yrange = [-6, 6]
;ybin_size = 2.5
;
;i_bin = where(base_req)
;stat_plot, transpose(xyz_all(1,i_bin)), abs(transpose(dxyz_all(1,i_bin))), k_c = 3, bin_range = yrange, binsize = ybin_size, qtt_2_range = dy_range, qtt_range = yrange, qtt_2_title = '|'+cap_delta_letter+'Y!n| [R!dE!n]', qtt_title = 'mean Y [R!dE!n]', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick

;;;;;;;; dy as a function of ny, dawn dusk separate
if strmatch(list_suf, '*_df_*') then begin
	dawn_loc = 0.
	dusk_loc = 2.
	ny_crit = 0.0

	bin_range = [-1, 1.]
	ny_binsize = 0.2

	dusk = xyz_all(1,*) gt dusk_loc
	dawn = xyz_all(1,*) lt dawn_loc
	morning = n_all(1,*) lt -ny_crit
	evening = n_all(1,*) gt ny_crit
	
	i_dusk = where(base_req and dusk, n_dusk)
	i_dawn = where(base_req and dawn, n_dawn)
	i_dusk_morning = where(base_req and dusk and morning, n_dusk_morning)
	i_dusk_evening = where(base_req and dusk and evening, n_dusk_evening)
	i_dawn_morning = where(base_req and dawn and morning, n_dawn_morning)
	i_dawn_evening = where(base_req and dawn and evening, n_dawn_evening)

	stat_plot, transpose(n_all(1,i_dusk)), abs(transpose(dxyz_all(1,i_dusk))), k_c = 3, bin_range = bin_range, binsize = ny_binsize, qtt_2_range = dy_range, qtt_range = bin_range, qtt_2_title = '|'+cap_delta_letter+'Y!n| [R!dE!n]', qtt_title = 'n!dy', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick, title = 'Dusk'
	makepng, pic_folder+'/dy_dusk'
	stat_plot, transpose(n_all(1,i_dawn)), abs(transpose(dxyz_all(1,i_dawn))), k_c = 3, bin_range = bin_range, binsize = ny_binsize, qtt_2_range = dy_range, qtt_range = bin_range, qtt_2_title = '|'+cap_delta_letter+'Y!n| [R!dE!n]', qtt_title = 'n!dy', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick, title = 'Dawn'
	makepng, pic_folder+'/dy_dawn'
endif

stop
end
