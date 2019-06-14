pro main_dy
;; find out the dy between satellites during events or non events
thm_init

computer = 'I:'
;computer = '/home/jliu'

@folders

;list_suf = '_fgs'
;list_suf = ''
list_suf = '_earthward_df_binxbout'

v_suffix = ''

xrange = [-12., -6.]*RE
rho_req = 12.*RE

if strmatch(list_suf, '*_df_*') then n = datain_simple(save_folder+'/n'+list_suf+'.dat', dim=3, type='double')
;;; load events
listname = 'dfb_list_lead_tail'+list_suf
events = load_list(listname+'.txt', folder = listfolder)
;;; diagnose
;events = events(*, 0:3)

;;;;;;;;;;;;;;;;;;; get the satellite distances for all events ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
probes_more = ['a', 'b', 'c', 'd', 'e']
probes_less = ['a', 'd', 'e']
all_dxyz = [0., 0., 0.]
all_xyz = [0., 0., 0.] ;; this is the average XYZ between the two probes.
if strmatch(list_suf, '*_df_*') then all_n = [0., 0., 0.] ;; this is the average XYZ between the two probes.
for i = 0, n_elements(events(0,*))-1 do begin
	time = time_double(events(0,i))
	if time le time_double('2010 1 1') then probes = probes_more else probes = probes_less
	sc = events(3,i)
	sc_num = probes_num_string(sc, all = probes)
	other_num = fw_array_uid(probes_num_string(probes, all = probes), sc_num, /diff)
	del_data, '*'
	load_bin_data, probe = probes, trange = [time-90., time+90.], datatype='pos', /tclip, datafolder = pos_folder
	pos = dblarr(3, n_elements(probes))
	;;; get satellite propagation velocity and X*Y* coordinates
	if ~strcmp(v_suffix, '') then begin
		events_do = [time_string(time), time_string(time), 'm', sc]
		no_use = df_thick(events_do, normal_method = 'binxbout', v_use = 'efs', seconds_check = 15., b_type = 'fgl', fgs_folder = fgs_folder, fgl_folder = fgl_folder, e_folder = efs_folder, normal_dir = n_this, c_angle = 30., maxgap = 10, vperp = vperp, vxy = vxy)
		case v_suffix of
		'vperp': vp = vperp
		'vxy': vp = vxy
		endcase
	endif else begin
		vp = [100., 0., 0.]
	endelse
	vp_strgh_xy = sqrt(total(vp(0:1)^2))
	vdir_xy = vp(0:1)/vp_strgh_xy
	front_dir_xy = [-vdir_xy(1), vdir_xy(0)]
	;;; get probe locations
	for j = 0, n_elements(probes)-1 do begin
		get_data, 'th'+probes(j)+'_state_pos_tclip', data = data
		pos(0,j) = mean(data.y(*,0), /nan)
		pos(1,j) = mean(data.y(*,1), /nan)
		pos(2,j) = mean(data.y(*,2), /nan)
	endfor
	dxyz_this = dblarr(3, n_elements(other_num))+!values.f_nan
	xyz_this = dxyz_this
	if strmatch(list_suf, '*_df_*') then n_this = dxyz_this
	for j = 0, n_elements(other_num)-1 do begin
		;;; limit the positions of the probes to be within requirement
		x_other = pos(0, other_num(j))
		y_other = pos(1, other_num(j))
		z_other = pos(2, other_num(j))
		rho_other = sqrt(y_other^2+z_other^2)
		if (x_other gt xrange(0)) and (x_other lt xrange(1)) and (rho_other lt rho_req) then begin
			;;; dxyz
			dxyz_temp = pos(*, other_num(j))-pos(*, sc_num)
			dx = dxyz_temp(0)*vdir_xy(0)+dxyz_temp(1)*vdir_xy(1)
			dy = dxyz_temp(0)*front_dir_xy(0)+dxyz_temp(1)*front_dir_xy(1)
			dxyz_this(*,j) = [dx, dy, dxyz_temp(2)]
			;;; average xyz
			xyz_this(*,j) = 0.5*(pos(*, other_num(j))+pos(*, sc_num))
			;;; normal direction
			if strmatch(list_suf, '*_df_*') then n_this(*,j) = n(*,i)
		endif
	endfor
	all_dxyz = [[all_dxyz], [dxyz_this]]
	all_xyz = [[all_xyz], [xyz_this]]
	all_n = [[all_n], [n_this]]
endfor
if n_elements(all_dxyz(0,*)) gt 1 then begin
	all_dxyz = all_dxyz(*, 1:*)
	all_xyz = all_xyz(*, 1:*)
	all_n = all_n(*, 1:*)
endif

all_dxyz = all_dxyz/RE
all_xyz = all_xyz/RE

dataout_simple, save_folder+'/dxyz'+list_suf, all_dxyz
dataout_simple, save_folder+'/xyz_ave'+list_suf, all_xyz
dataout_simple, save_folder+'/n_sideall'+list_suf, all_n
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
stop

end
