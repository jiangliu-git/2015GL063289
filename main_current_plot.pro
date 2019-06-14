pro main_current_plot
;;; plot the linear density of different current components
;; plot the strength/width of the dfcs/cross-tail cs
;; do the plot for main_cs_thick.pro, main_dfcs_thick.pro, main_current_strength.pro
thm_init

computer = 'I:'
@folders

;;; choose which magnetic field to use
b_type = 'bin'
;b_type = 'bave'
;b_type = 'bdirave'

;;; choose which dB representation to use
;cu_type = 'l'
cu_type = 'db'

;;; choose whether to use only good normal direction events
goodn_only = 'yes'
;goodn_only = 'no'
c_angle = 30 ;; 30 degs as the criterion for good normal
k_c = 3 ;; bins less than 5 points will be marked NaN

;fold = 0 ;; no fold
fold = 2 ;; fold north and south
;fold = 4 ;; fold all 4 places together

;;; conditions
case b_type of
'bin': b_str = 'B!din!n'
'bave': b_str = jiao_l+'B'+jiao_r
'bdirave': b_str = jiao_l+'b'+jiao_r
endcase

;;;; load data
pos = datain_simple(save_folder+'/pos.dat', dim=3, type='double')
h = datain_simple(save_folder+'/h.dat', dim=1, type='double')
thick_vec = datain_simple(save_folder+'/df_thick_vec.dat', dim=3, type='double')
b_in = datain_simple(save_folder+'/B_in.dat', dim=3, type='double')
b_out = datain_simple(save_folder+'/B_out.dat', dim=3, type='double')
b_ave = datain_simple(save_folder+'/B_ave.dat', dim=3, type = 'double')
l = datain_simple(save_folder+'/l.dat', dim=3, type='double')
lambda = datain_simple(save_folder+'/lambda.dat', dim=3, type='double')
db = b_in - b_out
if strcmp(b_type, 'bin') then b_use = b_in
if strcmp(b_type, 'bave') then b_use = b_ave
if strcmp(b_type, 'bdirave') then b_use = datain_simple(save_folder+'/B_dirave.dat', dim=3, type = 'double')
if strcmp(cu_type, 'l') then cu = l
if strcmp(cu_type, 'db') then cu = db

b_out_strength = transpose(sqrt(total(b_out^2,1))) 
db_strength = transpose(sqrt(total(db^2,1))) 
dbl_strength = transpose(total(db*l,1)) ;; by definition is always positive
dbl = rebin(dbl_strength, 3, n_elements(dbl_strength))*l

;;;;;;;; the angle between B_use and CU
b_use_strength = transpose(sqrt(total(b_use^2, 1)))
b_use_dir = b_use/[b_use_strength, b_use_strength, b_use_strength]
if strcmp(cu_type, 'l') then cu_dir = l
if strcmp(cu_type, 'db') then begin
	db_dir = db/[db_strength, db_strength, db_strength]
	cu_dir = db_dir
endif
dotp_arr = b_use_dir(0,*)*cu_dir(0,*)+b_use_dir(1,*)*cu_dir(1,*)+b_use_dir(2,*)*cu_dir(2,*)
crossp_arr = [b_use_dir(1,*)*cu_dir(2,*)-b_use_dir(2,*)*cu_dir(1,*),b_use_dir(2,*)*cu_dir(0,*)-b_use_dir(0,*)*cu_dir(2,*),b_use_dir(0,*)*cu_dir(1,*)-b_use_dir(1,*)*cu_dir(0,*)]
angle_abs = acos(dotp_arr)*180./!pi
angle_sign = crossp_arr(0,*)/abs(crossp_arr(0,*))
angle_L = angle_abs*angle_sign

;;;;;;;; the angle between B_use and B_out
b_out_dir = b_out/[b_out_strength, b_out_strength, b_out_strength]
dotp_arr = dotp_long(b_use_dir, b_out_dir)
angle_binbout_abs = acos(dotp_arr)*180./!pi

;;;;;;;; get normal direction on the fly
n = crossp_long(b_use, cu, /normalize, length = crossp_length)
for i = 0, n_elements(n(0,*))-1 do begin
  if n(0,i) lt 0 then n(*,i) = -n(*,i)
endfor

;;;;;;;; calculate DF thickness and re-adjust n
thick = dotp_long(thick_vec, n)
earthward = thick gt 0
tailward = thick lt 0
i_tail = where(tailward, n_tail)
if n_tail gt 0 then begin
	n(*,i_tail) = -n(*,i_tail)
	thick = abs(thick)
endif

;; the "good angle" conditions
good_l = lambda(0,*)/lambda(1,*) gt 3. ;; make sure L is correct
good_db = db_strength gt 2. ;; in fact all of them are gt 2
good_b_out = b_out_strength gt 2. ;; 1296 of them are gt 2
if strcmp(cu_type, 'l') then begin
	good_cu = good_l 
	additional_condition = 0
endif else begin
	good_cu = good_db
	;additional_condition = 0 ;; only dBxBuse
	additional_condition = (angle_binbout_abs gt c_angle) and (angle_binbout_abs lt 180.-c_angle) and good_b_out ;; when dB bad angle, if Bout good also do
endelse
if strcmp(b_type, 'bdirave') then b_strength = transpose(sqrt(total(b_ave^2, 1))) else b_strength = transpose(sqrt(total(b_use^2, 1)))
good_b = b_strength gt 2

good_earth = earthward and good_b and good_cu

;; good normal conditions
if strcmp(goodn_only, 'yes') then begin
	good_suf = '_goodn'
	good_angle = ((angle_abs gt c_angle) and (angle_abs lt 180.-c_angle)) or additional_condition
	;good_angle = (angle_abs gt c_angle)
	good = good_angle and good_earth
endif else begin
	good_suf = ''
	good = good_earth
endelse

i_good = where(good, j_good)
if j_good gt 0 then begin
	pos = pos(*, i_good)
	h = h(*, i_good)
	b_in = b_in(*, i_good)
	b_out = b_out(*, i_good)
	b_use = b_use(*, i_good)
	b_use_dir = b_use_dir(*, i_good)
	l = l(*, i_good)
	n = n(*, i_good)
	thick = thick(*, i_good)
	thick_vec = thick_vec(*, i_good)
	angle_abs = angle_abs(*, i_good)
	angle_sign = angle_sign(*, i_good)
	angle_L = angle_L(*, i_good)
endif

;; number of earthward-moving DFBs
no_use = where(thick_vec(0,*) ge 0, n_earthmoving)
print, 'Percentage of earthward moving:'
print, double(n_earthmoving)/double(j_good)

x = pos(0,*)
y = pos(1,*)
z = pos(2,*)
rho = sqrt(x^2+y^2)
MLT = asin(-y/rho)/(15./180*!pi)
nx = n(0,*)
ny = n(1,*)
phi = asin(ny)*180./!pi
phi_xy = asin(ny/sqrt(nx^2+ny^2))*180./!pi

if strcmp(cu_type, 'l') then begin
	dbcu_strength = dbl_strength(*, i_good) 
	dbcu = dbl(*, i_good) 
endif
if strcmp(cu_type, 'db') then begin
	dbcu_strength = db_strength(*, i_good) 
	dbcu = db(*, i_good) 
endif

dbcu_par_strength = dotp_long(dbcu, b_use_dir)
dbcu_par = rebin(dbcu_par_strength, 3, n_elements(dbcu_par_strength))*b_use_dir
dbcu_per = dbcu-dbcu_par
dbcu_per_strength = transpose(sqrt(total(dbcu_per^2, 1)))

I_full = abs(dbcu_strength)/mu0

;;; signed values
I_per_signed = dbcu_par_strength/mu0
I_par_signed = dbcu_per_strength/mu0*angle_sign
I_par_signed_projecth = pro2h(I_par_signed, h, sign(h), cu_type = cu_type, b_type = b_type, quantity = 'i', method = 'linear_med_abs')

j_par_signed = I_par_signed/(thick*1000.)

;;;;;;;;; the angle between Bin and Z, BinXZ to Z
;z = dblarr(3, n_elements(b_in_dir(0,*)))
;z(2,*)=1.
;dotp_arr_z = b_in_dir(0,*)*z(0,*)+b_in_dir(1,*)*z(1,*)+b_in_dir(2,*)*z(2,*)
;crossp_arr_z = [b_in_dir(1,*)*z(2,*)-b_in_dir(2,*)*z(1,*),b_in_dir(2,*)*z(0,*)-b_in_dir(0,*)*z(2,*),b_in_dir(0,*)*z(1,*)-b_in_dir(1,*)*z(0,*)]
;angle_abs_z = acos(dotp_arr_z)*180./!pi
;angle_sign_z = crossp_arr_z(0,*)/abs(crossp_arr_z(0,*))
;angle_Bin2z = angle_abs_z*angle_sign_z
;angle_bxz2z = atan2(b_in_dir(0,*), b_in_dir(2,*))*180./!pi 

;;;;;;;;; the angle between B_use and the XZ/YZ plane
angle_buse_XZ = asin(b_use_dir[1,*])*180./!pi
angle_buse_YZ = asin(b_use_dir[0,*])*180./!pi
angle_buse_XY = asin(b_use_dir[2,*])*180./!pi
;; the angle from the XY projection to B_use direction
z_dir = [0., 0., 1.]
surf_dir = crossp_long(n, rebin(z_dir, 3, n_elements(n[0,*])), /normalize)
angle_b_surf = (!pi/2-acos(dotp_long(surf_dir, b_use_dir)))*sign(ny)*180./!pi

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
;; choose whether use MLT or Y to seperate events.
;sep = 'Y'
sep = 'MLT'

x_sep_near = -15 ;; separation between near-Earth and far
x_sep_far = -22 ;; separation between far and too far

y_dusk = 2.
y_dawn = 2.
mlt_dusk = -1.
mlt_dawn = -1.

;qtt_name = 'ny'
;qtt_bin = ny
;qtt_range = [-1., 1.]
;qtt_title = 'ny'
;bin_bounds = [-1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.]

qtt_name = 'phi'
qtt_bin = phi
qtt_range = [90., -90.]
qtt_title = phi_letter+' [deg]'
area = 1 ;; make shaded area
;;; same-size bins
n_bins = 7 ;; for paper
;n_bins = 5 ;; for auxilliary
;n_bins = 2 ;; for values
bin_bounds = (findgen(n_bins+1)/n_bins-0.5)*180
;;; custom bins
;bounds_oneside = [6.92308, 20.7692, 34.6154, 51., 67, 90.] ;; very good for dusk side events
;bin_bounds = [-reverse(bounds_oneside), bounds_oneside]

;qtt_name = 'phixy'
;qtt_bin = phi_xy
;qtt_range = [90., -90.]
;qtt_title = 'phi (XY plane) [deg]'
;;; same-size bins
;n_bins = 7
;bin_bounds = (findgen(n_bins+1)/n_bins-0.5)*180
;;;; custom bins
;;bounds_oneside = [6.92308, 20.7692, 34.6154, 51., 67, 90.] ;; very good for dusk side events
;;bin_bounds = [-reverse(bounds_oneside), bounds_oneside]

;qtt_name = 'h'
;qtt_bin = h
;qtt_range = [-1., 1.]
;qtt_title = 'h'
;bin_bounds = [-1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.]
;vertical = 1

;qtt_name = 'h_fold'
;qtt_bin = abs(h)
;qtt_range = [0., 1.]
;;qtt_title = '|B!dqx!n/B!dlobe,q!n|'
;qtt_title = '|h|'
;;bin_bounds = [0., 0.333, 0.666, 1.]
;;bin_bounds = [0., 0.25, 0.5, 0.75, 1.] ;; fo
;bin_bounds = [0., 0.2, 0.4, 0.6, 0.8, 1.] ;; for Bin+dB
;vertical = 1

;;;;;;; i ;;;;;;;;;;;;;;;;
;qtt_2_name = 'i_par'
;qtt_2_bin = abs(I_par_signed)*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE
;qtt_2_title = '|i!d//!n| [MA/R!dE!n]'
;qtt_2_range = [0., 1.8e7]*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE
;;;;; original units, for pro2h values
;;qtt_2_name = 'i_par'
;;qtt_2_bin = abs(I_par_signed) ;; factor: transform to MA/RE
;;qtt_2_title = '|i!d//!n| [nA/m]'
;;qtt_2_range = [0., 1.8e7] ;; factor: transform to MA/RE

;qtt_2_name = 'i_par_signed'
;qtt_2_bin = I_par_signed*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE
;qtt_2_title = 'i!d//!n [MA/R!dE!n]'
;qtt_2_range = [-1.5e7, 1.5e7]*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE

qtt_2_name = 'i_par_signed_fold2n'
qtt_2_bin = I_par_signed*sign(h)*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE
qtt_2_title = 'i!d//!n (folded)!n [MA/R!dE!n]'
qtt_2_range = [-1.5e7, 1.5e7]*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE

;qtt_2_name = 'i_par_signed_fold2n_projectXY'
;qtt_2_bin = I_par_signed*sign(h)*abs(cos(angle_b_surf))*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE
;qtt_2_title = 'i!d//!n to Earth (XY) [MA/R!dE,c!n]'
;qtt_2_range = [-1e7, 1e7]*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE

;qtt_2_name = 'i_par_signed_fold2evening'
;qtt_2_bin = I_par_signed*sign(ny)
;qtt_2_title = 'i_par, folded to evening [nA/m]'
;qtt_2_range = [-2e7, 2e7]
;fit = 1

;qtt_2_name = 'i_par_signed_fold2n2evening'
;qtt_2_bin = I_par_signed*sign(h)*sign(ny)
;qtt_2_title = 'i_par, folded to north and evening [nA/m]'
;qtt_2_range = [-2e7, 2e7]

;qtt_2_name = 'i_parproh'
;qtt_2_bin = abs(I_par_signed_projecth)
;qtt_2_title = 'i_par projected regarding h [nA/m]'
;qtt_2_range = [0., 5e7]

;qtt_2_name = 'i_parproh_signed_fold2n'
;qtt_2_bin = I_par_signed_projecth*sign(h)*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE
;qtt_2_title = 'i!d//!n to Earth [MA/R!dE!n]'
;qtt_2_range = [-2e7, 2e7]*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE

;qtt_2_name = 'i_parproh_signed_fold2n_projectXY' ;; this is not physical: cannot project XYsurf after project regarding h; because the angles at the two locations are different.
;qtt_2_bin = I_par_signed_projecth*sign(h)*abs(cos(angle_b_surf))*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE
;qtt_2_title = 'i!d//!n to Earth (XY) [MA/R!dE,c!n]'
;qtt_2_range = [-2e7, 2e7]*1e-9*6371.*1000.*1e-6 ;; factor: transform to MA/RE

;;;;;; j ;;;;;;;;;;;;;;;
;qtt_2_name = 'j_par'
;qtt_2_bin = abs(j_par_signed)
;qtt_2_title = '|j_par| [nA/m!u2!n]'
;qtt_2_range = [0., 20]

;qtt_2_name = 'j_par_signed_fold2n'
;qtt_2_bin = I_par_signed*sign(h)/(thick*1000.)
;qtt_2_title = 'j_par regarding h, folded to north [nA/m!u2!n]'
;qtt_2_range = [-10, 10]

;qtt_2_name = 'j_par_signed_fold2evening'
;qtt_2_bin = j_par_signed*sign(ny)
;qtt_2_title = 'j_par, folded to evening [nA/m]'
;qtt_2_range = [-20, 20]
;fit = 1

;qtt_2_name = 'j_parproh_signed_fold2n'
;qtt_2_bin = I_par_signed_projecth*sign(h)/(thick*1000.)
;qtt_2_title = 'j_par projected regarding h, folded to north [nA/m!u2!n]'
;qtt_2_range = [-20, 20]

;;;;; other ;;;;;;;;;;;
;qtt_2_name = b_type+'_XZ'
;qtt_2_bin = angle_buse_XZ 
;qtt_2_title = 'the angle between '+b_type+' and XZ plane [deg]'
;qtt_2_range = [-90, 90]

;qtt_2_name = b_type+'_XZ_abs'
;qtt_2_bin = abs(angle_buse_XZ)
;qtt_2_title = 'the |angle| between '+b_type+' and XZ plane [deg]'
;qtt_2_range = [0, 90]

;qtt_2_name = b_type+'_XY'
;qtt_2_bin = angle_buse_XY 
;qtt_2_title = 'the angle between '+b_type+' and XY plane [deg]'
;qtt_2_range = [-90, 90]

;qtt_2_name = b_type+'_XY_abs'
;qtt_2_bin = abs(angle_buse_XY)
;qtt_2_title = 'the |angle| between '+b_type+' and XY plane [deg]'
;qtt_2_range = [0, 90]

;qtt_2_name = b_type+'_surf'
;qtt_2_bin = angle_b_surf
;qtt_2_title = 'the angle between '+b_type+' and surface proj [deg]'
;qtt_2_range = [-90, 90]

;qtt_2_name = b_type+'_surf_abs'
;qtt_2_bin = abs(angle_b_surf)
;qtt_2_title = 'the |angle| between '+b_type+' and surface proj [deg]'
;qtt_2_range = [0, 90]


;;; finite
finite_q = finite(qtt_bin) and finite(qtt_2_bin)

;;; dawn or dusk
if strcmp(sep, 'Y') then begin
	dusk = y gt y_dusk
	dawn = y lt y_dawn
endif

if strcmp(sep, 'MLT') then begin
	dusk = MLT lt mlt_dusk
	dawn = MLT gt mlt_dawn
endif

;;; morning or evening
morning = phi lt 0.
evening = phi gt 0.

;;; x location
far = (x lt x_sep_near) and (x gt x_sep_far)

;;; north or south
north = h gt 0.
south = h lt 0.

store_data, 'all', data = {i:where(finite_q), title:'All', bin_boundaries:bin_bounds}
store_data, 'dusk_sec', data = {i:where(dusk), title:'Dusk sector', bin_boundaries:bin_bounds}
store_data, 'dusk_far', data = {i:where(dusk and far), title:'Dusk sector, far', bin_boundaries:bin_bounds}
store_data, 'dusk_morning', data = {i:where(dusk and morning), title:'Dusk sector, Morning side', bin_boundaries:bin_bounds}
store_data, 'dusk_evening', data = {i:where(dusk and evening), title:'Dusk sector, Evening side', bin_boundaries:bin_bounds}
store_data, 'dusk_north', data = {i:where(dusk and north), title:'Dusk sector, North', bin_boundaries:bin_bounds}
store_data, 'dusk_south', data = {i:where(dusk and south), title:'Dusk sector, South', bin_boundaries:bin_bounds}
store_data, 'dawn_sec', data = {i:where(dawn), title:'Dawn sector', bin_boundaries:bin_bounds}
store_data, 'dawn_far', data = {i:where(dawn and far), title:'Dawn sector, far', bin_boundaries:bin_bounds}
store_data, 'dawn_morning', data = {i:where(dawn and morning), title:'Dawn sector, Morning side', bin_boundaries:bin_bounds}
store_data, 'dawn_evening', data = {i:where(dawn and evening), title:'Dawn sector, Evening side', bin_boundaries:bin_bounds}
store_data, 'dawn_north', data = {i:where(dawn and north), title:'Dawn sector, North', bin_boundaries:bin_bounds}
store_data, 'dawn_south', data = {i:where(dawn and south), title:'Dawn sector, South', bin_boundaries:bin_bounds}

;vars_exam = 'all'
vars_exam = ['dusk_sec', 'dawn_sec']
;vars_exam = ['dusk_far', 'dawn_far']
;vars_exam = ['dusk_morning', 'dusk_evening', 'dawn_morning', 'dawn_evening']
;vars_exam = ['dusk_north', 'dusk_south', 'dawn_north', 'dawn_south']

n_panels = n_elements(vars_exam)

left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
space_horiz = 0.007
space_vert = 0.007
if n_panels eq 2 then begin
	n_p_horiz = n_panels
	n_p_vert = 1
	xsize = 6.
	ysize = 2.6
	abc = ['(a)', '(b)']
endif
if n_panels eq 1 then begin
	n_p_horiz = n_panels
	n_p_vert = 1
	xsize = 3.3
	ysize = 2.8
	abc = ['']
endif
if n_panels eq 4 then begin
	n_p_horiz = 2
	n_p_vert = 2
	bot_margin = 0.1
	right_margin = 0.05
	xsize = 6
	ysize = 5
	abc = ['(a)', '(c)', '(b)', '(d)']
endif

positions = panel_positions([n_p_horiz, n_p_vert], lr_margins = [left_margin, right_margin], bt_margins = [bot_margin, top_margin], space = [space_horiz, space_vert])

popen, pic_folder+'/'+qtt_2_name+'_vs_'+qtt_name+'_'+cu_type+'_vs_'+b_type+'_'+sep+good_suf
print_options,xsize=xsize,ysize=ysize
for i = 0, n_panels-1 do begin
	get_data, vars_exam[i], data = this
	print, 'Group: '+vars_exam[i]
	;; title settings
	if n_panels ge 4 then begin
		if (i eq 0) or (i eq 1) then begin
		    qtt_2_ticknames = ''
			qtt_2_title_show = qtt_2_title
		endif else begin
		    qtt_2_ticknames = replicate(' ', 59)
			qtt_2_title_show = ''
		endelse
		if (i eq 0) or (i eq 2) then begin
		    title = strmid(this.title,0,11)
		    qtt_ticknames = replicate(' ', 59)
			qtt_title_show = ''
		endif else begin
			title = ''
		    qtt_ticknames = ''
			qtt_title_show = qtt_title
		endelse
	endif else begin
		title = this.title
		qtt_title_show = qtt_title
		if i eq 0 then begin
		    qtt_2_ticknames = ''
			qtt_2_title_show = qtt_2_title
		endif else begin
		    qtt_2_ticknames = replicate(' ', 59)
			qtt_2_title_show = ''
		endelse
	endelse

	i_this = this.i
	qtt_this = qtt_bin(i_this)
	qtt2_this = qtt_2_bin(i_this)
	;;; compute average values at morning and evening
	if strmatch(qtt_name, '*phi*') then begin
		i_evening = where(qtt_this gt 0, n_evening)
		if n_evening gt 0 then begin
			print, '- Average at evening: '+strcompress(string(mean(qtt2_this[i_evening], /nan)))
			print, '- Median at evening: '+strcompress(string(median(qtt2_this[i_evening])))
		endif
		i_morning = where(qtt_this lt 0, n_morning)
		if n_morning gt 0 then begin
			print, '- Average values at morning: '+strcompress(string(mean(qtt2_this[i_morning], /nan)))
			print, '- median values at morning: '+strcompress(string(median(qtt2_this[i_morning])))
		endif
	endif
	stat_plot, qtt_this, qtt2_this, k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = qtt_2_range, qtt_range = qtt_range, qtt_2_title = qtt_2_title_show, qtt_title = qtt_title_show, kinbin = kinbin_all, bincntrs_out = hcntrs, vertical=vertical, qtt_tickname = qtt_ticknames, qtt_2_tickname = qtt_2_ticknames, bin_boundaries = this.bin_boundaries, title = title, avrg = avrg, std = std, med = med, /no_mean, color_med = 0, color_quar = 0, type_med = 'square', type_quar = 'ebar', /no_write_pm, area = area, position = positions[*,i], /noerase
	oplot, [0,0], !y.crange 
	oplot, !x.crange, [0,0]
	;;; get the linear dependence of i to h
	if keyword_set(fit) then begin
		n_pts = 300 ;; must be even
		;; first, make the med points a trend with 300 points
		qtt_all = dindgen(n_pts)/(n_pts-1.)*(qtt_range(1)-qtt_range(0))+qtt_range(0) 
		qtt2_all = interpol(med, hcntrs, qtt_all, /nan)
		;; take abs values because log does not alow abs.
		qtt_all_fit = abs(qtt_all)
		qtt2_all_fit = abs(qtt2_all)
		;; fit log
		fit = poly_fit(alog10(qtt_all_fit), alog10(qtt2_all_fit), 1)
		print, 'fit:'
		print, fit(1)
		;; plot fit to plot
		qtt_plot_half = dindgen(n_pts/2)/(n_pts/2-1.)*(qtt_range(1)-qtt_range(0))/2.
		qtt2_plot_half = -10.^fit(0)*qtt_plot_half^fit(1)
		qtt_plot = [-reverse(qtt_plot_half), qtt_plot_half]
		qtt2_plot = [-reverse(qtt2_plot_half), qtt2_plot_half]
		oplot, qtt2_plot, qtt_plot, color = 4, thick = 2
	endif

	;;; get the net area and positive, negative parts
	med_all = [med[0],med,med[-1]]
	if keyword_set(vertical) then begin
		ctr_all = [!y.crange[0],hcntrs,!y.crange[1]]
	endif else begin
		ctr_all = [!x.crange[0],hcntrs,!x.crange[1]]
	endelse
	ctr_all = ctr_all[sort(ctr_all)]
	area_net = total((ctr_all[1:*]-ctr_all[0:-2])*(med_all[0:-2]+med_all[1:*])*0.5)
	if strmatch(qtt_name, '*phi*') then area_net = area_net*!pi/180. ;; deg to rad
	;xyouts, !x.crange[0]+0.5*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.8*(!y.crange[1]-!y.crange[0]), 'Net area:'+strcompress(string(area_net)), /data, align = 0.5
	print, 'Net area:'+strcompress(string(area_net))
	;; get the point that is zero
	i_trans = where(med_all[1:*]*med_all[0:-2] lt 0, n_trans)
	if n_trans eq 1 then begin
		ctr_mid = interpol([ctr_all[i_trans], ctr_all[i_trans+1]], [med_all[i_trans], med_all[i_trans+1]], 0.)
		i_low = where(ctr_all lt ctr_mid)
		i_high = where(ctr_all gt ctr_mid)
		ctr_low = [ctr_all[i_low], ctr_mid]
		ctr_high = [ctr_mid, ctr_all[i_high]]
		med_low = [med_all[i_low], 0.]
		med_high = [0., med_all[i_high]]
		area_low = total((ctr_low[1:*]-ctr_low[0:-2])*(med_low[0:-2]+med_low[1:*])*0.5)
		area_high = total((ctr_high[1:*]-ctr_high[0:-2])*(med_high[0:-2]+med_high[1:*])*0.5)
		;;; times DFB radius
		area_xr_low = total((ctr_low[1:*]-ctr_low[0:-2])*(med_low[0:-2]+med_low[1:*])*0.5*findr(0.5*(ctr_low[1:*]+ctr_low[0:-2]), sector = vars_exam(i)))
		area_xr_high = total((ctr_high[1:*]-ctr_high[0:-2])*(med_high[0:-2]+med_high[1:*])*0.5*findr(0.5*(ctr_high[1:*]+ctr_high[0:-2]), sector = vars_exam(i)))
		;;; deg to rad
		if strmatch(qtt_name, '*phi*') then begin
			area_low = area_low*!pi/180. 
			area_high = area_high*!pi/180. 
			area_xr_low = area_xr_low*!pi/180. 
			area_xr_high = area_xr_high*!pi/180. 
		endif
		area_sm = min(abs([area_low, area_high]))
		area_xr_sm = min(abs([area_xr_low, area_xr_high]))
		area_xr_net = area_xr_high+area_xr_low
		rate = abs(area_net)/area_sm
		rate_xr = abs(area_xr_net)/area_xr_sm
		;;; write down the areas of the shadows
		;xyouts, !x.crange[0]+0.7*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.58*(!y.crange[1]-!y.crange[0]), strcompress(string(area_low, format = '(g9.2)'))+'!cMA/[R!dE,r!n]', /data, charsize = 0.7
		;xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.43*(!y.crange[1]-!y.crange[0]), strcompress(string(area_high, format = '(g9.2)'))+'!cMA/[R!dE,r!n]', /data, charsize = 0.7
		;xyouts, !x.crange[0]+0.3*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.3*(!y.crange[1]-!y.crange[0]), 'net is'+strcompress(string(rate))+'of smaller', /data
		print, 'NOT CONSIDERING R:'
		print, 'low area:'+strcompress(string(area_low))
		print, 'high area:'+strcompress(string(area_high))
		print, 'net is'+strcompress(string(rate))+'of smaller'
		print, 'CONSIDERING R:'
		print, 'low current:'+strcompress(string(area_xr_low))+' MA'
		print, 'high area:'+strcompress(string(area_xr_high))+' MA'
		print, 'net is'+strcompress(string(area_xr_net))+' MA'
		print, 'net is'+strcompress(string(rate_xr))+' of smaller'
	endif
	
	;;; draw labels
	;; panel number
	xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.89*(!y.crange[1]-!y.crange[0]), abc[i], /data
	;; evening and morning
	if strmatch(qtt_name, '*phi*') and ((n_panels eq 2) or ((n_panels eq 4) and ((i eq 0) or (i eq 2)))) then begin
		xyouts, !x.crange[0], !y.crange[1]+0.01*(!y.crange[1]-!y.crange[0]), 'Evening', /data, align = 0., charsize = 0.7
		xyouts, !x.crange[1], !y.crange[1]+0.01*(!y.crange[1]-!y.crange[0]), 'Morning', /data, align = 1., charsize = 0.7
	endif
	;; north and south
	if n_p_horiz eq 2 then begin
		str = ''
		if i eq 2 then str = 'North'
		if i eq 3 then str = 'South'
		xyouts, !x.crange[1]+0.03*(!x.crange[1]-!x.crange[0]), !y.crange[0]+0.63*(!y.crange[1]-!y.crange[0]), str, /data, align = 0., charsize = 1.3, orientation = -90
	endif
	;;; save figure
	;makepng, pic_folder+'/'+qtt_2_name+'_vs_'+qtt_name+'_'+vars_exam(i)+'_'+cu_type+'_vs_'+b_type+'_'+sep+good_suf
endfor
pclose

stop
end
