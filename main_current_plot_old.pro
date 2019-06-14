pro main_current_plot_old
;;; plot the linear density of different current components
;; plot the strength/width of the dfcs/cross-tail cs
thm_init

computer = 'I:'
@folders

;;; choose which magnetic field to use
b_type = 'bin'
;b_type = 'bave'
;b_type = 'bdirave'

;;; choose whether to use only good normal direction events
goodn_only = 'yes'
;goodn_only = 'no'
c_angle = 30 ;; 30 degs as the criterion for good normal
k_c = 5 ;; bins less than 5 points will be marked NaN

;fold = 0 ;; no fold
fold = 2 ;; fold north and south
;fold = 4 ;; fold all 4 places together

;;; conditions
case b_type of
'bin': list_suf = '_earthward_df_binxl_nocare'
'bave': list_suf = '_earthward_df_bavexl_nocare'
'bdirave': list_suf = '_earthward_df_bdiravexl_nocare'
endcase
case b_type of
'bin': b_str = 'B!din!n'
'bave': b_str = jiao_l+'B'+jiao_r
'bdirave': b_str = jiao_l+'b'+jiao_r
endcase

;;;; load data
pos = datain_simple(save_folder+'/pos'+list_suf+'.dat', dim=3, type='double')
h = datain_simple(save_folder+'/h'+list_suf+'.dat', dim=1, type='double')
b_in = datain_simple(save_folder+'/B_in'+list_suf+'.dat', dim=3, type='double')
b_out = datain_simple(save_folder+'/B_out'+list_suf+'.dat', dim=3, type='double')
l = datain_simple(save_folder+'/l'+list_suf+'.dat', dim=3, type='double')
if strcmp(b_type, 'bin') then b_use = b_in
if strcmp(b_type, 'bave') then b_use = datain_simple(save_folder+'/B_ave'+list_suf+'.dat', dim=3, type = 'double')
if strcmp(b_type, 'bdirave') then b_use = datain_simple(save_folder+'/B_dirave'+list_suf+'.dat', dim=3, type = 'double')

;;;;;;;; the angle between B_use and L
b_use_strength = transpose(sqrt(total(b_use^2, 1)))
b_use_dir = b_use/[b_use_strength, b_use_strength, b_use_strength]
dotp_arr = b_use_dir(0,*)*l(0,*)+b_use_dir(1,*)*l(1,*)+b_use_dir(2,*)*l(2,*)
crossp_arr = [b_use_dir(1,*)*l(2,*)-b_use_dir(2,*)*l(1,*),b_use_dir(2,*)*l(0,*)-b_use_dir(0,*)*l(2,*),b_use_dir(0,*)*l(1,*)-b_use_dir(1,*)*l(0,*)]
angle_abs = acos(dotp_arr)*180./!pi
angle_sign = crossp_arr(0,*)/abs(crossp_arr(0,*))
angle_L = angle_abs*angle_sign

;;;;;;;; get normal direction on the fly
n = crossp_long(b_use, l, /normalize, length = crossp_length)
for i = 0, n_elements(n(0,*))-1 do begin
  if n(0,i) lt 0 then n(*,i) = -n(*,i)
endfor

if strcmp(goodn_only, 'yes') then begin
	good_suf = '_goodn'
	i_goodn = where(angle_abs gt c_angle, j_goodn)
	if j_goodn gt 0 then begin
		pos = pos(*, i_goodn)
		h = h(*, i_goodn)
		b_in = b_in(*, i_goodn)
		b_out = b_out(*, i_goodn)
		b_use = b_use(*, i_goodn)
		b_use_dir = b_use_dir(*, i_goodn)
		l = l(*, i_goodn)
		n = n(*, i_goodn)
		angle_abs = angle_abs(*, i_goodn)
		angle_sign = angle_sign(*, i_goodn)
		angle_L = angle_L(*, i_goodn)
	endif
endif else good_suf = ''

x = pos(0,*)
y = pos(1,*)
z = pos(2,*)
;MLT = 
ny = n(1,*)
phi = asin(ny)*180./!pi

dbl_strength = transpose(total((b_in-b_out)*l,1)) 
dbl = rebin(dbl_strength, 3, n_elements(dbl_strength))*l

dbl_par_strength = dotp_long(dbl, b_use_dir)
dbl_par = rebin(dbl_par_strength, 3, n_elements(dbl_par_strength))*b_use_dir
dbl_per = dbl-dbl_par
dbl_per_strength = transpose(sqrt(total(dbl_per^2, 1)))
I_full = abs(dbl_strength)/mu0
I_per_signed = dbl_par_strength/mu0
I_par_signed = dbl_per_strength/mu0*angle_sign
I_per = abs(I_per_signed)
I_par = abs(I_par_signed)
I_par_projecth = I_par/abs(h)*1. ;; linear relationship
I_par_signed_projecth = I_par_signed/abs(h)*1. ;; linear relationship

;;;;;;;;; the angle between Bin and Z, BinXZ to Z
;z = dblarr(3, n_elements(b_in_dir(0,*)))
;z(2,*)=1.
;dotp_arr_z = b_in_dir(0,*)*z(0,*)+b_in_dir(1,*)*z(1,*)+b_in_dir(2,*)*z(2,*)
;crossp_arr_z = [b_in_dir(1,*)*z(2,*)-b_in_dir(2,*)*z(1,*),b_in_dir(2,*)*z(0,*)-b_in_dir(0,*)*z(2,*),b_in_dir(0,*)*z(1,*)-b_in_dir(1,*)*z(0,*)]
;angle_abs_z = acos(dotp_arr_z)*180./!pi
;angle_sign_z = crossp_arr_z(0,*)/abs(crossp_arr_z(0,*))
;angle_Bin2z = angle_abs_z*angle_sign_z
;angle_bxz2z = atan2(b_in_dir(0,*), b_in_dir(2,*))*180./!pi 

;;;;;;;;;;;;;;;;;;; make 2d surface plot to examine the current distribution ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;n_bin_phi = 6
;n_bin_h = 7
;
;if fold eq 0 then h_range_plot = [-1., 1.] else begin
;	h = abs(h)
;	h_range_plot = [0., 1.]
;endelse
;if fold eq 4 then begin
;	phi = abs(phi)
;	phi_range_plot = [0, 90.]
;endif else phi_range_plot = [-90., 90.]
;
;binsize_phi = (phi_range_plot(1)-phi_range_plot(0))/(n_bin_phi-1)
;phi_range_bin = phi_range_plot+[-0.5*binsize_phi, 0.5*binsize_phi+0.01]
;binsize_h = (h_range_plot(1)-h_range_plot(0))/(n_bin_h-1)
;h_range_bin = h_range_plot+[-0.5*binsize_h, 0.5*binsize_h+0.01]
;
;titles = ['i', 'i_perp', 'i_par']
;ztitles = titles+' [nA/m]'
;qtts2bin = [I_full, I_per, I_par]
;
;for i = 0, n_elements(qtts2bin(*,0))-1 do begin
;	qtt2bin = qtts2bin(i, *)
;	;; bin the data
;	bin2d, transpose(phi), transpose(h), transpose(qtt2bin), xrange = phi_range_bin, yrange = h_range_bin, binsize = [binsize_phi, binsize_h], binhistogram = counts, xcenters = phi_cntrs, ycenters = h_cntrs, averages = qtt_avrg, medians = qtt_med, stdevs = qtt_std
;	i_few = where(counts lt k_c, j_few)
;	if j_few gt 0 then begin
;		qtt_avrg(i_few) = !values.f_nan
;		qtt_med(i_few) = !values.f_nan
;		qtt_std(i_few) = !values.f_nan
;	endif
;	;; plot the data
;	plotxyz, phi_cntrs, h_cntrs, qtt_med, xrange = phi_range_plot, yrange = h_range_plot, xstyle = 1, ystyle = 1, title=titles(i), xtitle = 'phi=asin(ny) [deg]', ytitle = 'Bqx/Blobe,q', ztitle = ztitles(i), /noisotropic
;	makepng, pic_folder+'/'+titles(i)+'_'+b_type
;endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;; make stat plot for different MLT events ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
y_dusk = 2.
y_dawn = 2.

qtt_bin = phi
qtt_name = 'phi'
qtt_range = [-90., 90.]
qtt_title = 'phi [deg]'
;bin_bounds = [-90., -60., 60., 90.]
;bin_bounds = [-90., -75., -60., 60., 75., 90.]
bin_bounds = [-90., -60., -30., 0., 30., 60., 90.]
;bin_bounds = [-90., -30., 0., 30., 90.]
;bin_bounds = [-90., 0., 90.]
;bin_bounds = [-90., -10., 10., 90.]

;qtt_name = 'h'
;qtt_bin = h
;qtt_range = [-1., 1.]
;qtt_title = 'h'
;bin_bounds = [-1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.]

;qtt_name = 'h_fold'
;qtt_bin = abs(h)
;qtt_range = [0., 1.]
;qtt_title = '|h|'
;bin_bounds = [0., 0.25, 0.5, 0.75, 1.]
;vertical = 1

;qtt_2_name = 'i_par'
;qtt_2_bin = I_par
;qtt_2_title = 'i_par [nA/m]'
;qtt_2_range = [0., 2.5e7]

;qtt_2_name = 'i_par_signed_fold2n'
;qtt_2_bin = I_par_signed*sign(h)
;qtt_2_title = 'i_par regarding h, folded to north [nA/m]'
;qtt_2_range = [-2e7, 2e7]

;qtt_2_name = 'i_parproh'
;qtt_2_bin = I_par_projecth
;qtt_2_title = 'i_par projected regarding h [nA/m]'
;qtt_2_range = [0., 5e7]

qtt_2_name = 'i_parproh_signed_fold2n'
qtt_2_bin = I_par_signed_projecth*sign(h)
qtt_2_title = 'i_par projected regarding h, folded to north [nA/m]'
qtt_2_range = [-4e7, 4e7]

dusk = y gt y_dusk
dawn = y lt y_dawn
morning = phi lt 0.
evening = phi gt 0.

i_dusk_sec = where(dusk)
i_dusk_morning = where(dusk and morning)
i_dusk_evening = where(dusk and evening)
i_dawn_sec = where(dawn)
i_dawn_morning = where(dawn and morning)
i_dawn_evening = where(dawn and evening)

store_data, 'dusk_sec', data = {i:where(dusk), title:'Dusk sector', bin_boundaries:bin_bounds}
store_data, 'dusk_morning', data = {i:where(dusk and morning), title:'Dusk sector, Morning side', bin_boundaries:bin_bounds}
store_data, 'dusk_evening', data = {i:where(dusk and evening), title:'Dusk sector, Evening side', bin_boundaries:bin_bounds}
store_data, 'dawn_sec', data = {i:where(dawn), title:'Dawn sector', bin_boundaries:bin_bounds}
store_data, 'dawn_morning', data = {i:where(dawn and morning), title:'Dawn sector, Morning side', bin_boundaries:bin_bounds}
store_data, 'dawn_evening', data = {i:where(dawn and evening), title:'Dawn sector, Evening side', bin_boundaries:bin_bounds}

vars_exam = ['dusk_sec', 'dawn_sec']
;vars_exam = ['dusk_morning', 'dusk_evening', 'dawn_morning', 'dawn_evening']

for i = 0, n_elements(vars_exam)-1 do begin
	get_data, vars_exam(i), data = this
	i_this = this.i
	stat_plot, qtt_bin(i_this), qtt_2_bin(i_this), k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = qtt_2_range, qtt_range = qtt_range, qtt_2_title = qtt_2_title, qtt_title = qtt_title, kinbin = kinbin_all, bincntrs_out = hcntrs, vertical=vertical, bin_boundaries = this.bin_boundaries, title = this.title, avrg = avrg, std = std, med = med, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write_pm
	oplot, [0,0], !y.crange 
	oplot, !x.crange, [0,0]
	print, med
	makepng, pic_folder+'/'+qtt_2_name+'_vs_'+qtt_name+'_'+vars_exam(i)+'_'+b_type+good_suf
endfor

stop
end
