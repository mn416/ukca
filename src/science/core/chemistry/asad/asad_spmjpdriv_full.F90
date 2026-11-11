module asad_spmjpdriv_full_mod
  implicit none
  public

  contains
  subroutine asad_spmjpdriv(n_points)
    use asad_mod, only : asad_state_type, base_tracer, ctype, ffrac, fracterms, frpx, ftol, h_in_ss, ilcf, ilct, ilft, ilftr, ilss, ilstmin, jpcf, jpco, jpcspf, jpfm, jpif, jpmsp, jpna, jpnr, jpoo, jpro2, jpsp, jpspec, jsro2, madvtr, majors, moffam, n_in_ss, ndepd, ndepw, negterms, nf, nfracterms, ngrp, nitfg, nitnr, njcoth, nldepx, nlf, nlfro2, nlmajmin, nlstst, nmsjac, nnegterms, nnfrp, nodd, nonzero_map, nonzero_map_unordered, nposterms, nprdx1, nprdx2, nprdx3, nsjac1, nsph, nspho2, nspi, nspn, nspno, nspo1d, nspo3, nspo3p, nspoh, nss_h, nss_n, nss_o1d, nss_o3p, nssi, nsspi, nsspt, nssri, nssrt, nssrx, nsst, nstst, ntabfp, ntro3, nuni, o3p_in_ss, peps, posterms, ptol, reorder, s=>asad_state, speci, spfjsize_max, total
    use ukca_config_specification_mod, only : ukca_config
    use parkind1, only : jpim, jprb
    use yomhook, only : dr_hook, lhook
    use ereport_mod, only : ereport
    use asad_diffun_mod, only : asad_diffun
    use asad_spimpmjp_mod, only : asad_spimpmjp
    use asad_ftoy_mod, only : asad_ftoy, ift, istmin
    use asad_sparse_vars, only : spfuljac, splinslv2, spresolv2
    use ukca_um_legacy_mod, only : mype
    use asad_fuljac_mod, only : asad_fuljac
    use asad_steady_mod, only : asad_steady
    use asad_fyfixr_mod, only : asad_fyfixr
    use asad_fyinit_mod, only : asad_fyinit, icf, ict, iftr, iss
    use asad_prls_mod, only : asad_prls
    use asad_fyself_mod, only : asad_fyself
    use ukca_environment_fields_mod, only : atmospheric_ch4, atmospheric_co2, atmospheric_h2, atmospheric_n2, atmospheric_o2
    use ukca_constants, only : c_ch4, c_co2, c_h2, c_ho2, c_n2, c_no3, c_o2, c_o3, c_oh
    use asad_inicnt_mod, only : asad_inicnt
    integer, parameter :: max_redo = 128
    logical, parameter :: not_first_call = .false.
    integer, parameter :: maxneg = 2000
    real, parameter :: f_max = 1.0e30
    real, parameter :: max_val = -1.0
    logical, parameter :: not_first_call_1 = .false.
    integer, parameter :: ltrig_iter = 51
    integer, parameter :: incr_limiter = 7
    integer, parameter :: n_o3 = 1
    integer, parameter :: n_oh = 2
    integer, parameter :: n_ho2 = 3
    integer, parameter :: n_no = 4
    integer, parameter :: n_o3_1 = 1
    integer, parameter :: n_oh_1 = 2
    integer, parameter :: n_ho2_1 = 3
    integer, parameter :: n_no_1 = 4
    integer, parameter :: n_o3_2 = 1
    integer, parameter :: n_oh_2 = 2
    integer, parameter :: n_ho2_2 = 3
    integer, parameter :: n_no_2 = 4
    real, parameter :: fco2_default = 350.0e-6
    real, parameter :: fh2_default = 5.0e-7
    real, parameter :: fn2_default = 0.78084
    real, parameter :: fo2_default = 0.20945
    real, parameter :: fch4_default = 1.76e-6
    real, parameter :: fco2_default_1 = 350.0e-6
    real, parameter :: fh2_default_1 = 5.0e-7
    real, parameter :: fn2_default_1 = 0.78084
    real, parameter :: fo2_default_1 = 0.20945
    real, parameter :: fch4_default_1 = 1.76e-6
    real, parameter :: fco2_default_2 = 350.0e-6
    real, parameter :: fh2_default_2 = 5.0e-7
    real, parameter :: fn2_default_2 = 0.78084
    real, parameter :: fo2_default_2 = 0.20945
    real, parameter :: fch4_default_2 = 1.76e-6
    integer, intent(in) :: n_points
    integer :: ix
    integer :: jy
    integer :: nlev
    integer :: exit_code
    integer :: ncsteps_initial
    integer :: iredo
    integer :: js
    integer :: iter
    integer :: location
    integer :: solver_iter
    integer :: num_iter
    integer :: jit
    integer :: errcode
    real :: cdt_initial
    real, dimension(n_points,jpcspf) :: f_initial
    integer :: iter_1
    integer :: jl
    integer :: jtr
    integer :: ifi
    integer :: i
    integer :: itr
    integer :: count_negatives
    real :: ztmp
    real :: f_min
    real :: reltol_residual_error
    real :: reltol_error
    real :: rafmin
    real :: rafmax
    real :: deltt
    real :: damp1
    real :: error_norm
    real :: residual_error
    real, dimension(n_points,jpcspf) :: f_initial_1
    real, dimension(n_points,jpcspf) :: f_incr
    real, dimension(n_points,jpcspf) :: g_f
    real, dimension(n_points,jpcspf) :: g_f_old
    real, dimension(n_points,jpcspf) :: g_ftmp
    real, dimension(n_points,jpcspf) :: delta_g
    real :: coeff
    real, dimension(jpcspf) :: zsum
    integer :: widx2
    integer :: widx1
    integer :: jr
    integer :: ix_1
    integer :: i_1
    integer :: j
    integer :: x
    real, dimension(n_points) :: ssnum
    real, dimension(n_points) :: ssden
    real, dimension(n_points,4) :: dssnum
    real, dimension(n_points,4) :: dssden
    real :: deltt_1
    character(len=2) :: ityped
    integer :: p
    integer :: irj
    integer :: ifamd
    integer :: itrd
    integer :: i_2
    integer :: jc
    integer :: j3
    integer :: jn
    integer :: js_16
    integer :: jl_1
    integer :: ij
    integer :: jl_2
    integer :: jtr_1
    integer :: ip
    real, dimension(n_points) :: fro2
    real, dimension(n_points) :: fro2_old
    real :: damp_ro2
    integer :: j_1
    integer :: jit_2
    integer :: jl_3
    integer :: js_1
    integer :: ifam
    integer :: imaj
    integer :: istart
    integer :: iend
    integer :: iodd
    integer :: itr_1
    integer :: icode
    integer :: iro2
    integer :: errcode_6
    real :: zthresh
    real :: sl
    real, dimension(n_points,jpspec) :: zy
    real, dimension(n_points) :: zb
    real, dimension(n_points) :: zc
    real, dimension(n_points) :: zd
    logical :: gconv
    integer :: jr_1
    integer :: ix_2
    integer :: i_3
    integer :: j_2
    integer :: x_1
    real, dimension(n_points) :: ssnum_1
    real, dimension(n_points) :: ssden_1
    real, dimension(n_points,4) :: dssnum_1
    real, dimension(n_points,4) :: dssden_1
    integer :: ifam_1
    integer :: itr_2
    integer :: j_3
    integer :: jl_4
    integer :: jtr_2
    integer :: js_2
    logical :: gfam
    logical :: gtr
    real, dimension(n_points,jpcspf) :: tmprc
    integer :: jl_5
    integer :: jtr_3
    integer :: j_4
    real :: deltt_2
    character(len=2) :: ityped_1
    integer :: p_1
    integer :: irj_1
    integer :: ifamd_1
    integer :: itrd_1
    integer :: i_4
    integer :: jc_1
    integer :: j3_1
    integer :: jn_1
    integer :: js_3
    integer :: jl_6
    integer :: ij_1
    integer :: total1
    integer :: kr
    integer :: jl_7
    integer :: i_5
    integer :: j_5
    integer :: ikr
    integer :: krj
    integer :: ij_2
    real, dimension(n_points,jpcspf) :: bb1
    real, dimension(n_points,jpcspf) :: xx1
    real, dimension(n_points) :: pivot
    real, dimension(n_points) :: kfact
    integer :: jl_8
    integer :: jtr_4
    real, dimension(n_points) :: fro2_1
    real, dimension(n_points) :: fro2_old_1
    real :: damp_ro2_1
    integer :: j_6
    integer :: jit_1
    integer :: jl_9
    integer :: js_4
    integer :: ifam_2
    integer :: imaj_1
    integer :: istart_1
    integer :: iend_1
    integer :: iodd_1
    integer :: itr_3
    integer :: icode_1
    integer :: iro2_1
    integer :: errcode_1
    real :: zthresh_1
    real :: sl_1
    real, dimension(n_points,jpspec) :: zy_1
    real, dimension(n_points) :: zb_1
    real, dimension(n_points) :: zc_1
    real, dimension(n_points) :: zd_1
    logical :: gconv_1
    integer :: jr_2
    integer :: ix_3
    integer :: i_6
    integer :: j_7
    integer :: x_2
    real, dimension(n_points) :: ssnum_2
    real, dimension(n_points) :: ssden_2
    real, dimension(n_points,4) :: dssnum_2
    real, dimension(n_points,4) :: dssden_2
    integer :: ifam_3
    integer :: itr_4
    integer :: j_8
    integer :: jl_10
    integer :: jtr_5
    integer :: js_5
    logical :: gfam_1
    logical :: gtr_1
    integer :: kr_1
    integer :: i_7
    integer :: j_9
    integer :: krj_1
    integer :: ikr_1
    real, dimension(n_points,jpcspf) :: bb1_1
    real, dimension(n_points,jpcspf) :: xx1_1
    real, dimension(n_points) :: pivot_1
    integer :: j_10
    integer :: jl_11
    integer :: js_6
    integer :: istart_2
    integer :: iend_2
    integer :: ifam_4
    integer :: itr_5
    integer :: errcode_2
    real :: peps10
    real :: fch4
    real :: fco2
    real :: fh2
    real :: fn2
    real :: fo2
    real :: zcnst
    character(len=10) :: species_tmp
    integer :: j_11
    integer :: jl_12
    integer :: jr_3
    integer :: js_7
    integer :: isp
    integer :: istart_3
    integer :: iend_3
    integer :: j_12
    integer :: jl_13
    integer :: js_8
    integer :: ifam_5
    integer :: imaj_2
    integer :: itr_6
    integer :: j_13
    integer :: j1
    integer :: j2
    integer :: j3_2
    integer :: jl_14
    integer :: jr_4
    integer :: js_9
    integer :: jskip
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: ir1
    integer :: ir2
    integer :: ip1
    integer :: ip2
    integer :: ip3
    integer :: iunip1
    integer :: iss_1
    integer :: ir
    integer :: ix_4
    integer :: istart_4
    integer :: iend_4
    real :: fr
    integer :: j_14
    integer :: j1_1
    integer :: j2_1
    integer :: j3_3
    integer :: jl_15
    integer :: jr_5
    integer :: js_10
    integer :: jskip_1
    integer :: i1_1
    integer :: i2_1
    integer :: i3_1
    integer :: ir1_1
    integer :: ir2_1
    integer :: ip1_1
    integer :: ip2_1
    integer :: ip3_1
    integer :: iunip1_1
    integer :: iss_2
    integer :: ir_1
    integer :: ix_5
    integer :: istart_5
    integer :: iend_5
    real :: fr_1
    integer :: kr_2
    integer :: i_8
    integer :: j_15
    integer :: krj_2
    integer :: ikr_2
    real, dimension(n_points,jpcspf) :: bb1_2
    real, dimension(n_points,jpcspf) :: xx1_2
    real, dimension(n_points) :: pivot_2
    integer :: j_16
    integer :: jl_16
    integer :: js_11
    integer :: istart_6
    integer :: iend_6
    integer :: ifam_6
    integer :: itr_7
    integer :: errcode_3
    real :: peps10_1
    real :: fch4_1
    real :: fco2_1
    real :: fh2_1
    real :: fn2_1
    real :: fo2_1
    real :: zcnst_1
    character(len=10) :: species_tmp_1
    integer :: j_17
    integer :: jl_17
    integer :: jr_6
    integer :: js_12
    integer :: isp_1
    integer :: istart_7
    integer :: iend_7
    integer :: j_18
    integer :: jl_18
    integer :: js_13
    integer :: ifam_7
    integer :: imaj_3
    integer :: itr_8
    integer :: j_19
    integer :: j1_2
    integer :: j2_2
    integer :: j3_4
    integer :: jl_19
    integer :: jr_7
    integer :: js_14
    integer :: jskip_2
    integer :: i1_2
    integer :: i2_2
    integer :: i3_2
    integer :: ir1_2
    integer :: ir2_2
    integer :: ip1_2
    integer :: ip2_2
    integer :: ip3_2
    integer :: iunip1_2
    integer :: iss_3
    integer :: ir_2
    integer :: ix_6
    integer :: istart_8
    integer :: iend_8
    real :: fr_2
    integer :: j_20
    integer :: j1_3
    integer :: j2_3
    integer :: j3_5
    integer :: jl_20
    integer :: jr_8
    integer :: js_15
    integer :: jskip_3
    integer :: i1_3
    integer :: i2_3
    integer :: i3_3
    integer :: ir1_3
    integer :: ir2_3
    integer :: ip1_3
    integer :: ip2_3
    integer :: ip3_3
    integer :: iunip1_3
    integer :: iss_4
    integer :: ir_3
    integer :: ix_7
    integer :: istart_9
    integer :: iend_9
    real :: fr_3
    integer :: errcode_4
    integer :: row_length
    integer :: rows
    integer :: iro2_2
    integer :: j_21
    real, dimension(n_points) :: fro2_2
    integer :: errcode_5
    integer :: row_length_1
    integer :: rows_1
    integer :: iro2_3
    integer :: j_22
    real, dimension(n_points) :: fro2_3
    real, dimension(n_points) :: fro2_4
    real, dimension(n_points) :: fro2_old_2
    real :: damp_ro2_2
    integer :: j_23
    integer :: jit_3
    integer :: jl_21
    integer :: js_17
    integer :: ifam_8
    integer :: imaj_4
    integer :: istart_10
    integer :: iend_10
    integer :: iodd_2
    integer :: itr_9
    integer :: icode_2
    integer :: iro2_4
    integer :: errcode_7
    real :: zthresh_2
    real :: sl_2
    real, dimension(n_points,jpspec) :: zy_2
    real, dimension(n_points) :: zb_2
    real, dimension(n_points) :: zc_2
    real, dimension(n_points) :: zd_2
    logical :: gconv_2
    integer :: ifam_9
    integer :: itr_10
    integer :: j_24
    integer :: jl_22
    integer :: jtr_6
    integer :: js_18
    logical :: gfam_2
    logical :: gtr_2
    integer :: j_25
    integer :: jl_23
    integer :: js_19
    integer :: istart_11
    integer :: iend_11
    integer :: ifam_10
    integer :: itr_11
    integer :: errcode_8
    real :: peps10_2
    real :: fch4_2
    real :: fco2_2
    real :: fh2_2
    real :: fn2_2
    real :: fo2_2
    real :: zcnst_2
    character(len=10) :: species_tmp_2
    integer :: j_26
    integer :: jl_24
    integer :: jr_9
    integer :: js_20
    integer :: isp_2
    integer :: istart_12
    integer :: iend_12
    integer :: j_27
    integer :: jl_25
    integer :: js_21
    integer :: ifam_11
    integer :: imaj_5
    integer :: itr_12
    integer :: j_28
    integer :: j1_4
    integer :: j2_4
    integer :: j3_6
    integer :: jl_26
    integer :: jr_10
    integer :: js_22
    integer :: jskip_4
    integer :: i1_4
    integer :: i2_4
    integer :: i3_4
    integer :: ir1_4
    integer :: ir2_4
    integer :: ip1_4
    integer :: ip2_4
    integer :: ip3_4
    integer :: iunip1_4
    integer :: iss_5
    integer :: ir_4
    integer :: ix_8
    integer :: istart_13
    integer :: iend_13
    real :: fr_4
    integer :: j_29
    integer :: j1_5
    integer :: j2_5
    integer :: j3_7
    integer :: jl_27
    integer :: jr_11
    integer :: js_23
    integer :: jskip_5
    integer :: i1_5
    integer :: i2_5
    integer :: i3_5
    integer :: ir1_5
    integer :: ir2_5
    integer :: ip1_5
    integer :: ip2_5
    integer :: ip3_5
    integer :: iunip1_5
    integer :: iss_6
    integer :: ir_5
    integer :: ix_9
    integer :: istart_14
    integer :: iend_14
    real :: fr_5
    integer :: errcode_9
    integer :: row_length_2
    integer :: rows_2
    integer :: iro2_5
    integer :: j_30
    real, dimension(n_points) :: fro2_5
    integer :: t
    integer :: n_states
    integer :: idx
    integer :: idx_1
    integer :: idx_2
    integer :: idx_3
    integer :: idx_4
    integer :: idx_5
    integer :: idx_6
    integer :: idx_7
    integer :: idx_8
    integer :: idx_9
    integer :: idx_10
    integer :: idx_11
    integer :: idx_12
    integer :: idx_13
    integer :: idx_14
    integer :: idx_15
    integer :: idx_16
    integer :: idx_17
    integer :: idx_18
    integer :: idx_19
    integer :: idx_20
    integer :: idx_21
    integer :: idx_22
    integer :: idx_23
    integer :: idx_24
    integer :: idx_25
    integer :: idx_26
    integer :: idx_27
    integer :: idx_28
    integer :: idx_29
    integer :: idx_30
    integer :: idx_31
    integer :: idx_32
    integer :: idx_33
    integer :: idx_34
    integer :: idx_35
    integer :: idx_36
    integer :: idx_37
    integer :: idx_38
    integer :: idx_39
    integer :: idx_40
    integer :: idx_41
    integer :: idx_42
    integer :: idx_43
    integer :: idx_44
    integer :: idx_45
    integer :: idx_46
    integer :: idx_47
    integer :: idx_48
    integer :: idx_49
    integer :: idx_50
    integer :: idx_51
    integer :: idx_52
    integer :: idx_53
    integer :: idx_54
    integer :: idx_55
    integer :: idx_56
    integer :: idx_57
    integer :: idx_58
    integer :: idx_59
    integer :: idx_60
    integer :: idx_61
    integer :: idx_62
    integer :: idx_63
    integer :: idx_64
    integer :: idx_65
    integer :: idx_66
    integer :: idx_67
    integer :: idx_68
    integer :: idx_69
    integer :: idx_70
    integer :: idx_71
    integer :: idx_72
    integer :: idx_73
    integer :: idx_74
    integer :: idx_75
    integer :: idx_76
    integer :: idx_77
    integer :: idx_78
    integer :: idx_79
    integer :: idx_80
    integer :: idx_81
    integer :: idx_82
    integer :: idx_83
    integer :: idx_84
    integer :: idx_85
    integer :: idx_86
    integer :: idx_87
    integer :: idx_88
    integer :: idx_89
    integer :: idx_90
    integer :: idx_91
    integer :: idx_92
    integer :: idx_93
    integer :: idx_94
    integer :: idx_95
    integer :: idx_96
    integer :: idx_97
    integer :: idx_98
    integer :: idx_99
    integer :: idx_100
    integer :: idx_101
    integer :: idx_102
    integer :: idx_103
    integer :: idx_104
    integer :: idx_105
    integer :: idx_106
    integer :: idx_107
    integer :: idx_108
    integer :: idx_109
    integer :: idx_110
    integer :: idx_111
    integer :: idx_112
    integer :: idx_113
    integer :: idx_114
    integer :: idx_115
    integer :: idx_116
    integer :: idx_117
    integer :: idx_118
    integer :: idx_119
    integer :: idx_120
    integer :: idx_121
    integer :: idx_122
    integer :: idx_123
    integer :: idx_124
    integer :: idx_125
    integer :: idx_126
    integer :: idx_127
    integer :: idx_128
    integer :: idx_129
    integer :: idx_130
    integer :: idx_131
    integer :: idx_132
    integer :: idx_133
    integer :: idx_134
    integer :: idx_135
    integer :: idx_136
    integer :: idx_137
    integer :: idx_138
    integer :: idx_139
    integer :: idx_140
    integer :: idx_141
    integer :: idx_142
    integer :: idx_143
    integer :: idx_144
    integer :: idx_145
    integer :: idx_146
    integer :: idx_147
    integer :: idx_148
    integer :: idx_149
    integer :: idx_150
    integer :: idx_151
    integer :: idx_152
    integer :: idx_153
    integer :: idx_154
    integer :: idx_155
    integer :: idx_156
    integer :: idx_157
    integer :: idx_158
    integer :: idx_159
    integer :: idx_160
    integer :: idx_161
    integer :: idx_162
    integer :: idx_163
    integer :: idx_164
    integer :: idx_165
    integer :: idx_166
    integer :: idx_167
    integer :: idx_168
    integer :: idx_169
    integer :: idx_170
    integer :: idx_171
    integer :: idx_172
    integer :: idx_173
    integer :: idx_174
    integer :: idx_175
    integer :: idx_176
    integer :: idx_177
    integer :: idx_178
    integer :: idx_179
    integer :: idx_180
    integer :: idx_181
    integer :: idx_182
    integer :: idx_183
    integer :: idx_184
    integer :: idx_185
    integer :: idx_186
    integer :: idx_187
    integer :: idx_188
    integer :: idx_189
    integer :: idx_190
    integer :: idx_191
    integer :: idx_192
    integer :: idx_193
    integer :: idx_194
    integer :: idx_195
    integer :: idx_196
    integer :: idx_197
    integer :: idx_198
    integer :: idx_199
    integer :: idx_200
    integer :: idx_201
    integer :: idx_202
    integer :: idx_203
    integer :: idx_204
    integer :: idx_205
    integer :: idx_206
    integer :: idx_207
    integer :: idx_208
    integer :: idx_209
    integer :: idx_210
    integer :: idx_211
    integer :: idx_212
    integer :: idx_213
    integer :: idx_214
    integer :: idx_215
    integer :: idx_216
    integer :: idx_217
    integer :: idx_218
    integer :: idx_219
    integer :: idx_220
    integer :: idx_221
    integer :: idx_222
    integer :: idx_223
    integer :: idx_224
    integer :: idx_225
    integer :: idx_226
    integer :: idx_227
    integer :: idx_228

    ix = 0
    jy = 0
    nlev = 0
    n_states = SIZE(s)
  !$ACC  PARALLEL LOOP GANG &
  !$ACC& PRIVATE(exit_code, ncsteps_initial, iredo, js, iter, location, &
  !$ACC& solver_iter, num_iter, jit, errcode, cdt_initial, f_initial, &
  !$ACC& iter_1, jl, jtr, ifi, i, itr, count_negatives, &
  !$ACC& ztmp, f_min, reltol_residual_error, reltol_error, rafmin, &
  !$ACC& rafmax, deltt, damp1, error_norm, residual_error, &
  !$ACC& f_initial_1, f_incr, g_f, g_f_old, g_ftmp, delta_g, coeff, &
  !$ACC& zsum,  widx2, widx1, jr, ix_1, i_1, j, x, &
  !$ACC& ssnum, ssden, dssnum, dssden,  deltt_1, &
  !$ACC& ityped, p, irj, ifamd, itrd, i_2, jc, j3, jn, js_16, jl_1, &
  !$ACC& ij,  jl_2, jtr_1, ip,  fro2, &
  !$ACC& fro2_old, damp_ro2, j_1, jit_2, jl_3, js_1, ifam, imaj, &
  !$ACC& istart, iend, iodd, itr_1, icode, iro2, errcode_6, zthresh, &
  !$ACC& sl, zy, zb, zc, zd, gconv,  jr_1, ix_2, i_3, &
  !$ACC& j_2, x_1, ssnum_1, ssden_1, dssnum_1, dssden_1, &
  !$ACC&  ifam_1, itr_2, j_3, jl_4, jtr_2, js_2, gfam, &
  !$ACC& gtr,  tmprc, jl_5, jtr_3, j_4, &
  !$ACC&  deltt_2, ityped_1, p_1, irj_1, ifamd_1, &
  !$ACC& itrd_1, i_4, jc_1, j3_1, jn_1, js_3, jl_6, ij_1, &
  !$ACC&  total1, kr, jl_7, i_5, j_5, ikr, krj, ij_2, &
  !$ACC& bb1, xx1, pivot, kfact,  jl_8, jtr_4, &
  !$ACC&  fro2_1, fro2_old_1, damp_ro2_1, j_6, jit_1, &
  !$ACC& jl_9, js_4, ifam_2, imaj_1, istart_1, iend_1, iodd_1, itr_3, &
  !$ACC& icode_1, iro2_1, errcode_1, zthresh_1, sl_1, zy_1, zb_1, &
  !$ACC& zc_1, zd_1, gconv_1,  jr_2, ix_3, i_6, j_7, &
  !$ACC& x_2, ssnum_2, ssden_2, dssnum_2, dssden_2,  &
  !$ACC& ifam_3, itr_4, j_8, jl_10, jtr_5, js_5, gfam_1, gtr_1, &
  !$ACC&  kr_1, i_7, j_9, krj_1, ikr_1, bb1_1, xx1_1, &
  !$ACC& pivot_1,  j_10, jl_11, js_6, istart_2, &
  !$ACC& iend_2, ifam_4, itr_5, errcode_2, peps10, fch4, fco2, fh2, &
  !$ACC& fn2, fo2, zcnst,  species_tmp, j_11, jl_12, &
  !$ACC& jr_3, js_7, isp,  istart_3, iend_3, j_12, &
  !$ACC& jl_13, js_8, ifam_5, imaj_2, itr_6,  j_13, &
  !$ACC& j1, j2, j3_2, jl_14, jr_4, js_9, jskip, i1, i2, i3, ir1, &
  !$ACC& ir2, ip1, ip2, ip3, iunip1, iss_1, ir, ix_4, istart_4, &
  !$ACC& iend_4, fr,  j_14, j1_1, j2_1, j3_3, jl_15, &
  !$ACC& jr_5, js_10, jskip_1, i1_1, i2_1, i3_1, ir1_1, ir2_1, ip1_1, &
  !$ACC& ip2_1, ip3_1, iunip1_1, iss_2, ir_1, ix_5, istart_5, iend_5, &
  !$ACC& fr_1,  kr_2, i_8, j_15, krj_2, ikr_2, bb1_2, &
  !$ACC& xx1_2, pivot_2,  j_16, jl_16, js_11, &
  !$ACC& istart_6, iend_6, ifam_6, itr_7, errcode_3, peps10_1, &
  !$ACC& fch4_1, fco2_1, fh2_1, fn2_1, fo2_1, zcnst_1, &
  !$ACC&  species_tmp_1, j_17, jl_17, jr_6, js_12, &
  !$ACC& isp_1,  istart_7, iend_7, j_18, jl_18, &
  !$ACC& js_13, ifam_7, imaj_3, itr_8,  j_19, j1_2, &
  !$ACC& j2_2, j3_4, jl_19, jr_7, js_14, jskip_2, i1_2, i2_2, i3_2, &
  !$ACC& ir1_2, ir2_2, ip1_2, ip2_2, ip3_2, iunip1_2, iss_3, ir_2, &
  !$ACC& ix_6, istart_8, iend_8, fr_2,  j_20, j1_3, &
  !$ACC& j2_3, j3_5, jl_20, jr_8, js_15, jskip_3, i1_3, i2_3, i3_3, &
  !$ACC& ir1_3, ir2_3, ip1_3, ip2_3, ip3_3, iunip1_3, iss_4, ir_3, &
  !$ACC& ix_7, istart_9, iend_9, fr_3,  errcode_4, &
  !$ACC& row_length, rows, iro2_2, j_21,  fro2_2, &
  !$ACC& errcode_5, row_length_1, rows_1, iro2_3, j_22, &
  !$ACC&  fro2_3, fro2_4, fro2_old_2, damp_ro2_2, &
  !$ACC& j_23, jit_3, jl_21, js_17, ifam_8, imaj_4, istart_10, &
  !$ACC& iend_10, iodd_2, itr_9, icode_2, iro2_4, errcode_7, &
  !$ACC& zthresh_2, sl_2, zy_2, zb_2, zc_2, zd_2, gconv_2, &
  !$ACC&  ifam_9, itr_10, j_24, jl_22, jtr_6, js_18, &
  !$ACC& gfam_2, gtr_2,  j_25, jl_23, js_19, &
  !$ACC& istart_11, iend_11, ifam_10, itr_11, errcode_8, peps10_2, &
  !$ACC& fch4_2, fco2_2, fh2_2, fn2_2, fo2_2, zcnst_2, &
  !$ACC&  species_tmp_2, j_26, jl_24, jr_9, js_20, &
  !$ACC& isp_2,  istart_12, iend_12, j_27, jl_25, &
  !$ACC& js_21, ifam_11, imaj_5, itr_12,  j_28, j1_4, &
  !$ACC& j2_4, j3_6, jl_26, jr_10, js_22, jskip_4, i1_4, i2_4, i3_4, &
  !$ACC& ir1_4, ir2_4, ip1_4, ip2_4, ip3_4, iunip1_4, iss_5, ir_4, &
  !$ACC& ix_8, istart_13, iend_13, fr_4,  j_29, j1_5, &
  !$ACC& j2_5, j3_7, jl_27, jr_11, js_23, jskip_5, i1_5, i2_5, i3_5, &
  !$ACC& ir1_5, ir2_5, ip1_5, ip2_5, ip3_5, iunip1_5, iss_6, ir_5, &
  !$ACC& ix_9, istart_14, iend_14, fr_5,  errcode_9, &
  !$ACC& row_length_2, rows_2, iro2_5, j_30,  fro2_5, &
  !$ACC& t)
  do t = 1, n_states, 1
      ncsteps_initial = s(t)%ncsteps
      cdt_initial = s(t)%cdt
      iredo = s(t)%ncsteps_factor
      do idx = 1, UBOUND(f_initial, dim=2), 1
        do idx_1 = 1, UBOUND(f_initial, dim=1), 1
          f_initial(idx_1,idx) = s(t)%f(idx_1,idx )
        enddo
      enddo
      s(t)%ltrig = .false.
      if (ukca_config%l_ukca_asad_full) then
        location = 1
      elseif (ukca_config%l_ukca_asad_columns) then
        location = ix + (jy - 1) * ukca_config%row_length
      else
        location = nlev
      end if
      num_iter = 0
      iter = 1
      do while (iter <= iredo)
        f_min = SQRT(peps)
        reltol_residual_error = 1.0e-10
        reltol_error = 10 * ptol
        rafmin = 1.0e-01
        rafmax = 1.0e+04
        exit_code = 0
        error_norm = 1.0
        deltt = 1.0 / s(t)%cdt
        damp1 = 0.5
        do idx_2 = 1, jpcspf, 1
          do idx_3 = 1, n_points, 1
            f_initial_1(idx_3,idx_2) = s(t)%f(idx_3 ,idx_2 )
          enddo
        enddo
        do widx2 = 1, jpcspf, 1
          do widx1 = 1, n_points, 1
            if (s(t)%f(widx1,widx2) < f_min) then
              s(t)%f(widx1,widx2) = f_min
            end if
          enddo
        enddo
        if (nstst /= 0) then
          do idx_4 = 1, UBOUND(s(t)%deriv, dim=3), 1
            do idx_5 = 1, UBOUND(s(t)%deriv, dim=2), 1
              do idx_6 = 1, UBOUND(s(t)%deriv, dim=1), 1
                s(t)%deriv(idx_6,idx_5,idx_4) = 1.0
              enddo
            enddo
          enddo
          do ix_1 = 1, nsst, 1
            do idx_7 = 1, n_points, 1
              ssnum(idx_7) = 0.0
            enddo
            do idx_8 = 1, n_points, 1
              ssden(idx_8) = 0.0
            enddo
            do idx_9 = 1, 4, 1
              do idx_10 = 1, n_points, 1
                dssden(idx_10,idx_9) = 0.0
              enddo
            enddo
            do idx_11 = 1, 4, 1
              do idx_12 = 1, n_points, 1
                dssnum(idx_12,idx_11) = 0.0
              enddo
            enddo
            do jr = 1, nsspt(ix_1), 1
              i_1 = nsspi(ix_1,jr)
              if (i_1 <= nuni) then
                do idx_13 = 1, n_points, 1
                  ssnum(idx_13) = ssnum(idx_13) + s(t)%rk(idx_13,i_1) * s(t)%y(idx_13,nspi(i_1,1))
                enddo
                if (ix_1 < 5 .AND. nspo3 == nspi(i_1,1)) then
                  do idx_14 = 1, n_points, 1
                    dssnum(idx_14,n_o3) = dssnum(idx_14,n_o3) + s(t)%rk(idx_14,i_1)
                  enddo
                end if
                if (ix_1 < 5 .AND. nspno == nspi(i_1,1)) then
                  do idx_15 = 1, n_points, 1
                    dssnum(idx_15,n_no) = dssnum(idx_15,n_no) + s(t)%rk(idx_15,i_1)
                  enddo
                end if
              else
                do idx_16 = 1, n_points, 1
                  ssnum(idx_16) = ssnum(idx_16) + s(t)%rk(idx_16,i_1) * s(t)%y(idx_16,nspi(i_1,1)) * s(t)%y(idx_16,nspi(i_1,2))
                enddo
                if (ix_1 < 5) then
                  if (nspo1d == nspi(i_1,1)) then
                    do idx_17 = 1, n_points, 1
                      dssnum(idx_17,n_o3) = dssnum(idx_17,n_o3) + s(t)%rk(idx_17,i_1) * s(t)%y(idx_17,nspi(i_1,2)) * s(t)%deriv(idx_17,nss_o1d,n_o3)
                    enddo
                  end if
                  if (nspo1d == nspi(i_1,2)) then
                    do idx_18 = 1, n_points, 1
                      dssnum(idx_18,n_o3) = dssnum(idx_18,n_o3) + s(t)%rk(idx_18,i_1) * s(t)%y(idx_18,nspi(i_1,1)) * s(t)%deriv(idx_18,nss_o1d,n_o3)
                    enddo
                  end if
                  if (nspo3p == nspi(i_1,1)) then
                    do idx_19 = 1, n_points, 1
                      dssnum(idx_19,n_o3) = dssnum(idx_19,n_o3) + s(t)%rk(idx_19,i_1) * s(t)%y(idx_19,nspi(i_1,2)) * s(t)%deriv(idx_19,nss_o3p,n_o3)
                    enddo
                  end if
                  if (nspo3p == nspi(i_1,2)) then
                    do idx_20 = 1, n_points, 1
                      dssnum(idx_20,n_o3) = dssnum(idx_20,n_o3) + s(t)%rk(idx_20,i_1) * s(t)%y(idx_20,nspi(i_1,1)) * s(t)%deriv(idx_20,nss_o3p,n_o3)
                    enddo
                  end if
                  if (nspo3 == nspi(i_1,1)) then
                    do idx_21 = 1, n_points, 1
                      dssnum(idx_21,n_o3) = dssnum(idx_21,n_o3) + s(t)%rk(idx_21,i_1) * s(t)%y(idx_21,nspi(i_1,2))
                    enddo
                  end if
                  if (nspo3 == nspi(i_1,2)) then
                    do idx_22 = 1, n_points, 1
                      dssnum(idx_22,n_o3) = dssnum(idx_22,n_o3) + s(t)%rk(idx_22,i_1) * s(t)%y(idx_22,nspi(i_1,1))
                    enddo
                  end if
                  if (nspo3p == nspi(i_1,1)) then
                    do idx_23 = 1, n_points, 1
                      dssnum(idx_23,n_oh) = dssnum(idx_23,n_oh) + s(t)%rk(idx_23,i_1) * s(t)%y(idx_23,nspi(i_1,2)) * s(t)%deriv(idx_23,nss_o3p,n_oh)
                    enddo
                  end if
                  if (nspo3p == nspi(i_1,2)) then
                    do idx_24 = 1, n_points, 1
                      dssnum(idx_24,n_oh) = dssnum(idx_24,n_oh) + s(t)%rk(idx_24,i_1) * s(t)%y(idx_24,nspi(i_1,1)) * s(t)%deriv(idx_24,nss_o3p,n_oh)
                    enddo
                  end if
                  if (nspoh == nspi(i_1,1)) then
                    do idx_25 = 1, n_points, 1
                      dssnum(idx_25,n_oh) = dssnum(idx_25,n_oh) + s(t)%rk(idx_25,i_1) * s(t)%y(idx_25,nspi(i_1,2))
                    enddo
                  end if
                  if (nspoh == nspi(i_1,2)) then
                    do idx_26 = 1, n_points, 1
                      dssnum(idx_26,n_oh) = dssnum(idx_26,n_oh) + s(t)%rk(idx_26,i_1) * s(t)%y(idx_26,nspi(i_1,1))
                    enddo
                  end if
                  if (nspo3p == nspi(i_1,1)) then
                    do idx_27 = 1, n_points, 1
                      dssnum(idx_27,n_ho2) = dssnum(idx_27,n_ho2) + s(t)%rk(idx_27,i_1) * s(t)%y(idx_27,nspi(i_1,2)) * s(t)%deriv(idx_27,nss_o3p,n_ho2)
                    enddo
                  end if
                  if (nspo3p == nspi(i_1,2)) then
                    do idx_28 = 1, n_points, 1
                      dssnum(idx_28,n_ho2) = dssnum(idx_28,n_ho2) + s(t)%rk(idx_28,i_1) * s(t)%y(idx_28,nspi(i_1,1)) * s(t)%deriv(idx_28,nss_o3p,n_ho2)
                    enddo
                  end if
                  if (nspho2 == nspi(i_1,1)) then
                    do idx_29 = 1, n_points, 1
                      dssnum(idx_29,n_ho2) = dssnum(idx_29,n_ho2) + s(t)%rk(idx_29,i_1) * s(t)%y(idx_29,nspi(i_1,2))
                    enddo
                  end if
                  if (nspho2 == nspi(i_1,2)) then
                    do idx_30 = 1, n_points, 1
                      dssnum(idx_30,n_ho2) = dssnum(idx_30,n_ho2) + s(t)%rk(idx_30,i_1) * s(t)%y(idx_30,nspi(i_1,1))
                    enddo
                  end if
                  if (nspno == nspi(i_1,1)) then
                    do idx_31 = 1, n_points, 1
                      dssnum(idx_31,n_no) = dssnum(idx_31,n_no) + s(t)%rk(idx_31,i_1) * s(t)%y(idx_31,nspi(i_1,2))
                    enddo
                  end if
                  if (nspno == nspi(i_1,2)) then
                    do idx_32 = 1, n_points, 1
                      dssnum(idx_32,n_no) = dssnum(idx_32,n_no) + s(t)%rk(idx_32,i_1) * s(t)%y(idx_32,nspi(i_1,1))
                    enddo
                  end if
                end if
              end if
            enddo
            do jr = 1, nssrt(ix_1), 1
              i_1 = nssri(ix_1,jr)
              j = nssrx(ix_1,jr)
              if (i_1 <= nuni) then
                do idx_33 = 1, n_points, 1
                  ssden(idx_33) = ssden(idx_33) + s(t)%rk(idx_33,i_1)
                enddo
              else
                do idx_34 = 1, n_points, 1
                  ssden(idx_34) = ssden(idx_34) + s(t)%rk(idx_34,i_1) * s(t)%y(idx_34,nspi(i_1,j))
                enddo
                if (ix_1 < 5) then
                  if (nspo3 == nspi(i_1,j)) then
                    do idx_35 = 1, n_points, 1
                      dssden(idx_35,n_o3) = dssden(idx_35,n_o3) + s(t)%rk(idx_35,i_1)
                    enddo
                  end if
                  if (nspoh == nspi(i_1,j)) then
                    do idx_36 = 1, n_points, 1
                      dssden(idx_36,n_oh) = dssden(idx_36,n_oh) + s(t)%rk(idx_36,i_1)
                    enddo
                  end if
                  if (nspho2 == nspi(i_1,j)) then
                    do idx_37 = 1, n_points, 1
                      dssden(idx_37,n_ho2) = dssden(idx_37,n_ho2) + s(t)%rk(idx_37,i_1)
                    enddo
                  end if
                  if (nspno == nspi(i_1,j)) then
                    do idx_38 = 1, n_points, 1
                      dssden(idx_38,n_no) = dssden(idx_38,n_no) + s(t)%rk(idx_38,i_1)
                    enddo
                  end if
                end if
              end if
            enddo
            do idx_39 = 1, n_points, 1
              s(t)%y(idx_39,nssi(ix_1)) = ssnum(idx_39) / ssden(idx_39)
            enddo
            if (ix_1 < 5) then
              do jr = 1, 4, 1
                do idx_40 = 1, n_points, 1
                  s(t)%deriv(idx_40,ix_1,jr) = (ssden(idx_40) * dssnum(idx_40,jr) - ssnum(idx_40) * dssden(idx_40,jr)) / (ssden(idx_40) * ssden(idx_40))
                enddo
              enddo
            end if
          enddo
          do x = 1, n_points, 1
            if (peps < s(t)%y(x,nspo1d)) then
              s(t)%deriv(x,nss_o1d,n_o3) = s(t)%deriv(x,nss_o1d,n_o3) * s(t)%y(x,nspo3) / s(t)%y(x,nspo1d)
              s(t)%deriv(x,nss_o1d,n_oh) = s(t)%deriv(x,nss_o1d,n_oh) * s(t)%y(x,nspoh) / s(t)%y(x,nspo1d)
              s(t)%deriv(x,nss_o1d,n_ho2) = s(t)%deriv(x,nss_o1d,n_ho2) * s(t)%y(x,nspho2) / s(t)%y(x,nspo1d)
              s(t)%deriv(x,nss_o1d,n_no) = s(t)%deriv(x,nss_o1d,n_no) * s(t)%y(x,nspno) / s(t)%y(x,nspo1d)
            else
              s(t)%deriv(x,nss_o1d,n_o3) = 1.0
              s(t)%deriv(x,nss_o1d,n_oh) = 1.0
              s(t)%deriv(x,nss_o1d,n_ho2) = 1.0
              s(t)%deriv(x,nss_o1d,n_no) = 1.0
            end if
          enddo
          if (o3p_in_ss) then
            do x = 1, n_points, 1
              if (peps < s(t)%y(x,nspo3p)) then
                s(t)%deriv(x,nss_o3p,n_o3) = s(t)%deriv(x,nss_o3p,n_o3) * s(t)%y(x,nspo3) / s(t)%y(x,nspo3p)
                s(t)%deriv(x,nss_o3p,n_oh) = s(t)%deriv(x,nss_o3p,n_oh) * s(t)%y(x,nspoh) / s(t)%y(x,nspo3p)
                s(t)%deriv(x,nss_o3p,n_ho2) = s(t)%deriv(x,nss_o3p,n_ho2) * s(t)%y(x,nspho2) / s(t)%y(x,nspo3p)
                s(t)%deriv(x,nss_o3p,n_no) = s(t)%deriv(x,nss_o3p,n_no) * s(t)%y(x,nspno) / s(t)%y(x,nspo3p)
              else
                s(t)%deriv(x,nss_o3p,n_o3) = 1.0
                s(t)%deriv(x,nss_o3p,n_oh) = 1.0
                s(t)%deriv(x,nss_o3p,n_ho2) = 1.0
                s(t)%deriv(x,nss_o3p,n_no) = 1.0
              end if
            enddo
          end if
          if (n_in_ss) then
            do x = 1, n_points, 1
              if (peps < s(t)%y(x,nspn)) then
                s(t)%deriv(x,nss_n,n_o3) = s(t)%deriv(x,nss_n,n_o3) * s(t)%y(x,nspo3) / s(t)%y(x,nspn)
                s(t)%deriv(x,nss_n,n_oh) = s(t)%deriv(x,nss_n,n_oh) * s(t)%y(x,nspoh) / s(t)%y(x,nspn)
                s(t)%deriv(x,nss_n,n_ho2) = s(t)%deriv(x,nss_n,n_ho2) * s(t)%y(x,nspho2) / s(t)%y(x,nspn)
                s(t)%deriv(x,nss_n,n_no) = s(t)%deriv(x,nss_n,n_no) * s(t)%y(x,nspno) / s(t)%y(x,nspn)
              else
                s(t)%deriv(x,nss_n,n_o3) = 1.0
                s(t)%deriv(x,nss_n,n_oh) = 1.0
                s(t)%deriv(x,nss_n,n_ho2) = 1.0
                s(t)%deriv(x,nss_n,n_no) = 1.0
              end if
            enddo
          end if
          if (h_in_ss) then
            do x = 1, n_points, 1
              if (peps < s(t)%y(x,nsph)) then
                s(t)%deriv(x,nss_h,n_o3) = s(t)%deriv(x,nss_h,n_o3) * s(t)%y(x,nspo3) / s(t)%y(x,nsph)
                s(t)%deriv(x,nss_h,n_oh) = s(t)%deriv(x,nss_h,n_oh) * s(t)%y(x,nspoh) / s(t)%y(x,nsph)
                s(t)%deriv(x,nss_h,n_ho2) = s(t)%deriv(x,nss_h,n_ho2) * s(t)%y(x,nspho2) / s(t)%y(x,nsph)
                s(t)%deriv(x,nss_h,n_no) = s(t)%deriv(x,nss_h,n_no) * s(t)%y(x,nspno) / s(t)%y(x,nsph)
              else
                s(t)%deriv(x,nss_h,n_o3) = 1.0
                s(t)%deriv(x,nss_h,n_oh) = 1.0
                s(t)%deriv(x,nss_h,n_ho2) = 1.0
                s(t)%deriv(x,nss_h,n_no) = 1.0
              end if
            enddo
          end if
        end if
        do i_2 = 1, jpcspf, 1
          do jl_1 = 1, n_points, 1
            if (f_min > s(t)%f(jl_1,i_2)) then
              s(t)%f(jl_1,i_2) = f_min
            end if
          enddo
        enddo
        deltt_1 = 1.0 / s(t)%cdt
        do idx_41 = 1, UBOUND(s(t)%spfj, dim=2), 1
          do idx_42 = 1, UBOUND(s(t)%spfj, dim=1), 1
            s(t)%spfj(idx_42,idx_41) = 0.0
          enddo
        enddo
        do i_2 = 1, jpcspf, 1
          do idx_43 = 1, UBOUND(s(t)%spfj, dim=1), 1
            s(t)%spfj(idx_43,nonzero_map(i_2,i_2)) = -deltt_1 * s(t)%f(idx_43 ,i_2)
          enddo
        enddo
        do p = 1, total, 1
          do i_2 = nposterms(p), 1, -1
            do idx_44 = 1, UBOUND(s(t)%spfj, dim=1), 1
              s(t)%spfj(idx_44,p) = s(t)%spfj(idx_44,p) + s(t)%prk(idx_44 ,posterms(p,i_2))
            enddo
          enddo
          do i_2 = nnegterms(p), 1, -1
            do idx_45 = 1, UBOUND(s(t)%spfj, dim=1), 1
              s(t)%spfj(idx_45,p) = s(t)%spfj(idx_45,p) - s(t)%prk(idx_45 ,negterms(p,i_2))
            enddo
          enddo
          do i_2 = nfracterms(p), 1, -1
            do idx_46 = 1, UBOUND(s(t)%spfj, dim=1), 1
              s(t)%spfj(idx_46,p) = s(t)%spfj(idx_46,p) + ffrac(p,i_2) * s(t)%prk(idx_46 ,fracterms(p,i_2))
            enddo
          enddo
        enddo
        do jc = 1, nstst, 1
          do j3 = 1, nmsjac(jc), 1
            irj = nsjac1(j3,jc)
            do jn = 1, jpmsp, 1
              ij = njcoth(irj,jn)
              if (ij /= 0) then
                p = nonzero_map(ij,ntro3)
                if (jn < 3) then
                  do idx_47 = 1, UBOUND(s(t)%spfj, dim=1), 1
                    s(t)%spfj(idx_47,p) = s(t)%spfj(idx_47,p) - s(t)%prk(idx_47 ,irj) * s(t)%deriv(idx_47 ,jc,1)
                  enddo
                else
                  do idx_48 = 1, UBOUND(s(t)%spfj, dim=1), 1
                    s(t)%spfj(idx_48,p) = s(t)%spfj(idx_48,p) + s(t)%prk(idx_48 ,irj) * s(t)%deriv(idx_48 ,jc,1)
                  enddo
                end if
              end if
            enddo
          enddo
        enddo
        if (ndepw /= 0 .OR. ndepd /= 0) then
          do js_16 = 1, jpspec, 1
            ifamd = moffam(js_16)
            itrd = madvtr(js_16)
            ityped = ctype(js_16)
            if (ifamd /= 0) then
              p = nonzero_map(ifamd,ifamd)
              do jl_1 = 1, n_points, 1
                if (ityped == jpfm .OR. ityped == jpif .AND. s(t)%linfam(jl_1,itrd)) then
                  s(t)%spfj(jl_1,p) = s(t)%spfj(jl_1,p) - nodd(js_16) * (s(t)%dpd(jl_1,js_16) + s(t)%dpw(jl_1,js_16)) * s(t)%y(jl_1,js_16)
                end if
              enddo
            end if
            if (itrd /= 0) then
              p = nonzero_map(itrd,itrd)
              do idx_49 = 1, UBOUND(s(t)%spfj, dim=1), 1
                s(t)%spfj(idx_49,p) = s(t)%spfj(idx_49,p) - (s(t)%dpd(idx_49 ,js_16) + s(t)%dpw(idx_49 ,js_16)) * s(t)%y(idx_49 ,js_16)
              enddo
            end if
          enddo
        end if
        do p = 1, total, 1
          do idx_50 = 1, UBOUND(s(t)%spfj, dim=1), 1
            s(t)%spfj(idx_50,p) = s(t)%spfj(idx_50,p) / s(t)%f(idx_50 ,base_tracer(p))
          enddo
        enddo
        do jtr_1 = 1, jpcspf, 1
          ip = nonzero_map(jtr_1,jtr_1)
          do jl_2 = 1, n_points, 1
            if (s(t)%spfj(jl_2,ip) > 0.0) then
              s(t)%f(jl_2,jtr_1) = f_initial_1(jl_2,jtr_1) + s(t)%cdt * s(t)%fdot(jl_2,jtr_1)
            else
              s(t)%f(jl_2,jtr_1) = f_initial_1(jl_2,jtr_1) + (s(t)%cdt * s(t)%fdot(jl_2,jtr_1)) / (1.0 - s(t)%cdt * s(t)%spfj(jl_2,ip))
            end if
            if (f_min > s(t)%f(jl_2,jtr_1)) then
              s(t)%f(jl_2,jtr_1) = f_min
            end if
          enddo
        enddo
        do iter_1 = 1, ukca_config%nrsteps, 1
          ifi = 0
          if (iter_1 == 1) then
            ifi = nitnr
          end if
          ifam = 0
          imaj = 0
          sl = 0.0
          zthresh = 2.0 / s(t)%cdt
          do idx_51 = 1, jpspec, 1
            do idx_52 = 1, n_points, 1
              zy(idx_52,idx_51) = 0.0
            enddo
          enddo
          do idx_53 = 1, n_points, 1
            zb(idx_53) = 0.0
          enddo
          do idx_54 = 1, n_points, 1
            zc(idx_54) = 0.0
          enddo
          do idx_55 = 1, n_points, 1
            zd(idx_55) = 0.0
          enddo
          if (not_first_call_1) then
            peps10 = 10.0 * peps
            istart_2 = nlmajmin(1)
            iend_2 = nlmajmin(2)
            do j_10 = istart_2, iend_2, 1
              js_6 = nlmajmin(j_10)
              ifam_4 = moffam(js_6)
              do jl_11 = 1, n_points, 1
                s(t)%y(jl_11,js_6) = s(t)%f(jl_11,ifam_4)
              enddo
            enddo
            istart_2 = nlmajmin(3)
            iend_2 = nlmajmin(4)
            do j_10 = istart_2, iend_2, 1
              js_6 = nlmajmin(j_10)
              do jl_11 = 1, n_points, 1
                s(t)%y(jl_11,js_6) = peps10
              enddo
            enddo
            do j_10 = 1, iss, 1
              js_6 = ilss(j_10)
              do jl_11 = 1, n_points, 1
                s(t)%y(jl_11,js_6) = peps10
              enddo
            enddo
            fco2 = fco2_default
            fch4 = fch4_default
            fh2 = fh2_default
            fo2 = fo2_default
            fn2 = fn2_default
            if (ukca_config%l_chem_environ_co2_scalar) then
              fco2 = atmospheric_co2 / c_co2
              if (ukca_config%l_chem_environ_co2_fld) then
                fco2 = 0.0
              end if
            end if
            if (ukca_config%l_chem_environ_h2_scalar) then
              fh2 = atmospheric_h2 / c_h2
            end if
            if (ukca_config%l_chem_environ_n2_scalar) then
              fn2 = atmospheric_n2 / c_n2
            end if
            if (ukca_config%l_chem_environ_o2_scalar) then
              fo2 = atmospheric_o2 / c_o2
            end if
            if (ukca_config%l_chem_environ_ch4_scalar) then
              fch4 = atmospheric_ch4 / c_ch4
            end if
            do j_10 = 1, ict, 1
              js_6 = ilct(j_10)
              if (speci(js_6) == 'N2        ') then
                zcnst = fn2
              elseif (speci(js_6) == 'O2        ') then
                zcnst = fo2
              elseif (speci(js_6) == 'CO2       ') then
                zcnst = fco2
              elseif (speci(js_6) == 'H2        ') then
                zcnst = fh2
              elseif (speci(js_6) == 'CH4       ') then
                zcnst = fch4
              else
                errcode_2 = -1
                zcnst = 0.0
              end if
              do jl_11 = 1, n_points, 1
                s(t)%y(jl_11,js_6) = zcnst * s(t)%tnd(jl_11)
              enddo
            enddo
          end if
          do j_10 = 1, iftr, 1
            js_6 = ilftr(j_10)
            itr_5 = madvtr(js_6)
            do jl_11 = 1, n_points, 1
              s(t)%y(jl_11,js_6) = s(t)%f(jl_11,itr_5)
            enddo
          enddo
          if (not_first_call_1) then
            do j_10 = 1, icf, 1
              js_6 = ilcf(j_10)
              species_tmp = speci(js_6)
              row_length = ukca_config%row_length
              rows = ukca_config%rows
              if (species_tmp(1:4) == 'CO2 ' .AND. ukca_config%l_chem_environ_co2_fld) then
                do idx_56 = 1, n_points, 1
                  s(t)%y(idx_56,js_6) = s(t)%co2(idx_56 )
                enddo
              elseif (species_tmp(1:4) == 'H2O ') then
                do idx_57 = 1, n_points, 1
                  s(t)%y(idx_57,js_6) = s(t)%wp(idx_57 )
                enddo
              elseif (species_tmp(1:4) == 'RO2 ') then
                if (ukca_config%l_ukca_ro2_perm) then
                  do idx_58 = 1, UBOUND(fro2_2, dim=1), 1
                    fro2_2(idx_58) = 0.0
                  enddo
                  do j_21 = 1, jpro2, 1
                    iro2_2 = nlfro2(j_21)
                    do idx_59 = 1, UBOUND(fro2_2, dim=1), 1
                      fro2_2(idx_59) = fro2_2(idx_59) + s(t)%f(idx_59 ,iro2_2)
                    enddo
                  enddo
                  do idx_60 = 1, n_points, 1
                    s(t)%y(idx_60,js_6) = fro2_2(idx_60) / s(t)%tnd(idx_60 )
                  enddo
                else
                  errcode_4 = 126
                end if
              else
                errcode_4 = 124
              end if
              do idx_61 = 1, n_points, 1
                s(t)%y(idx_61,js_6) = s(t)%y(idx_61,js_6) * s(t)%tnd(idx_61 )
              enddo
            enddo
          end if
          if (ukca_config%l_ukca_ro2_perm) then
            do idx_62 = 1, n_points, 1
              fro2_old(idx_62) = s(t)%y(idx_62,jsro2)
            enddo
            do idx_63 = 1, n_points, 1
              fro2(idx_63) = 0.0
            enddo
            do j_1 = 1, jpro2, 1
              iro2 = nlfro2(j_1)
              do idx_64 = 1, n_points, 1
                fro2(idx_64) = fro2(idx_64) + s(t)%f(idx_64,iro2)
              enddo
            enddo
            if (iter_1 <= 0) then
              damp_ro2 = 1.0
            elseif (iter_1 == 1) then
              damp_ro2 = 0.5
            elseif (iter_1 == 2) then
              damp_ro2 = 0.9
            else
              damp_ro2 = 0.99
            end if
            do idx_65 = 1, n_points, 1
              s(t)%y(idx_65,jsro2) = fro2_old(idx_65) + damp_ro2 * (fro2(idx_65) - fro2_old(idx_65))
            enddo
          end if
          if (nstst == 0) then
          end if
          if (not_first_call_1 .AND. ifi < 5) then
            ifi = 5
            icode = -1
          end if
          if (not_first_call_1) then
            do j_11 = 1, nstst, 1
              js_7 = nlstst(j_11)
              do jl_12 = 1, n_points, 1
                s(t)%qa(jl_12,js_7) = 0.0
              enddo
            enddo
            do jr_3 = 1, jpnr, 1
              isp = nspi(jr_3,1)
              if (isp == nspi(jr_3,2)) then
                do jl_12 = 1, n_points, 1
                  s(t)%qa(jl_12,isp) = s(t)%qa(jl_12,isp) + 2.0 * s(t)%rk(jl_12,jr_3)
                enddo
              end if
            enddo
          end if
          if (ifi == 0) then
            istart_3 = nlmajmin(1)
            iend_3 = nlmajmin(2)
            do j_12 = istart_3, iend_3, 1
              js_8 = nlmajmin(j_12)
              ifam_5 = moffam(js_8)
              do jl_13 = 1, n_points, 1
                s(t)%y(jl_13,js_8) = s(t)%f(jl_13,ifam_5) * s(t)%ratio(jl_13,js_8)
              enddo
            enddo
            istart_3 = nlmajmin(3)
            iend_3 = nlmajmin(4)
            do j_12 = istart_3, iend_3, 1
              js_8 = nlmajmin(j_12)
              ifam_5 = moffam(js_8)
              imaj_2 = majors(ifam_5)
              if (ctype(js_8) /= jpif) then
                do jl_13 = 1, n_points, 1
                  s(t)%y(jl_13,js_8) = s(t)%y(jl_13,imaj_2) * s(t)%ratio(jl_13,js_8)
                enddo
              else
                itr_6 = madvtr(js_8)
                do jl_13 = 1, n_points, 1
                  if (s(t)%linfam(jl_13,itr_6)) then
                    s(t)%y(jl_13,js_8) = s(t)%y(jl_13,imaj_2) * s(t)%ratio(jl_13,js_8)
                  end if
                enddo
              end if
            enddo
          end if
          do jit_2 = 1, ifi, 1
            if (jit_2 > 1) then
              do j_1 = 1, nstst, 1
                js_1 = nlstst(j_1)
                do jl_3 = 1, n_points, 1
                  if (peps < s(t)%pd(jl_3,js_1) .AND. peps > ABS(s(t)%y(jl_3,js_1))) then
                    s(t)%y(jl_3,js_1) = peps
                  end if
                enddo
              enddo
            end if
            do js_9 = 1, jpspec, 1
              do jl_14 = 1, n_points, 1
                s(t)%pd(jl_14,js_9) = 0.0
                s(t)%pd(jl_14,js_9 + jpspec) = 0.0
              enddo
            enddo
            do jr_4 = 1, nuni, 1
              ir1 = nspi(jr_4,1)
              do jl_14 = 1, n_points, 1
                s(t)%prk(jl_14,jr_4) = s(t)%rk(jl_14,jr_4) * s(t)%y(jl_14,ir1)
              enddo
            enddo
            iunip1 = nuni + 1
            do jr_4 = iunip1, jpnr, 1
              ir1 = nspi(jr_4,1)
              ir2 = nspi(jr_4,2)
              do jl_14 = 1, n_points, 1
                s(t)%prk(jl_14,jr_4) = s(t)%rk(jl_14,jr_4) * s(t)%y(jl_14,ir1) * s(t)%y(jl_14,ir2)
              enddo
            enddo
            do jskip = 0, jpspec, jpspec
              do j_13 = 1, istmin, 1
                js_9 = ilstmin(j_13) + jskip
                ip3 = ngrp(js_9,3)
                ip2 = ngrp(js_9,2)
                ip1 = ngrp(js_9,1)
                do j3_2 = 1, ip3, 1
                  i1 = nprdx3(1,j3_2,js_9)
                  i2 = nprdx3(2,j3_2,js_9)
                  i3 = nprdx3(3,j3_2,js_9)
                  do jl_14 = 1, n_points, 1
                    s(t)%pd(jl_14,js_9) = s(t)%pd(jl_14,js_9) + s(t)%prk(jl_14,i1) + s(t)%prk(jl_14,i2) + s(t)%prk(jl_14,i3)
                  enddo
                enddo
                do j2 = 1, ip2, 1
                  i1 = nprdx2(1,js_9)
                  i2 = nprdx2(2,js_9)
                  do jl_14 = 1, n_points, 1
                    s(t)%pd(jl_14,js_9) = s(t)%pd(jl_14,js_9) + s(t)%prk(jl_14,i1) + s(t)%prk(jl_14,i2)
                  enddo
                enddo
                do j1 = 1, ip1, 1
                  i1 = nprdx1(js_9)
                  do jl_14 = 1, n_points, 1
                    s(t)%pd(jl_14,js_9) = s(t)%pd(jl_14,js_9) + s(t)%prk(jl_14,i1)
                  enddo
                enddo
              enddo
            enddo
            do jr_4 = 1, nnfrp, 1
              iss_1 = ntabfp(jr_4,1)
              ir = ntabfp(jr_4,2)
              ix_4 = ntabfp(jr_4,3)
              fr = frpx(ix_4)
              do jl_14 = 1, n_points, 1
                s(t)%pd(jl_14,iss_1) = s(t)%pd(jl_14,iss_1) + fr * s(t)%prk(jl_14,ir)
              enddo
            enddo
            if (.false.) then
              istart_4 = nldepx(1)
              iend_4 = nldepx(2)
              do j_13 = istart_4, iend_4, 1
                js_9 = nldepx(j_13)
                do jl_14 = 1, n_points, 1
                  s(t)%pd(jl_14,js_9 + jpspec) = s(t)%pd(jl_14,js_9 + jpspec) + (s(t)%dpd(jl_14,js_9) + s(t)%dpw(jl_14,js_9)) * s(t)%y(jl_14,js_9)
                enddo
              enddo
              istart_4 = nldepx(3)
              iend_4 = nldepx(4)
              do j_13 = istart_4, iend_4, 1
                js_9 = nldepx(j_13)
                do jl_14 = 1, n_points, 1
                  s(t)%pd(jl_14,js_9 + jpspec) = s(t)%pd(jl_14,js_9 + jpspec) + s(t)%dpd(jl_14,js_9) * s(t)%y(jl_14,js_9)
                enddo
              enddo
              istart_4 = nldepx(5)
              iend_4 = nldepx(6)
              do j_13 = istart_4, iend_4, 1
                js_9 = nldepx(j_13)
                do jl_14 = 1, n_points, 1
                  s(t)%pd(jl_14,js_9 + jpspec) = s(t)%pd(jl_14,js_9 + jpspec) + s(t)%dpw(jl_14,js_9) * s(t)%y(jl_14,js_9)
                enddo
              enddo
            end if
            do j_1 = 1, nstst, 1
              js_1 = nlstst(j_1)
              do jl_3 = 1, n_points, 1
                s(t)%ratio(jl_3,js_1) = 0.0
              enddo
            enddo
            do j_1 = 1, istmin, 1
              js_1 = ilstmin(j_1)
              do jl_3 = 1, n_points, 1
                if (peps > s(t)%y(jl_3,js_1)) then
                  sl = 0.0
                else
                  sl = s(t)%pd(jl_3,js_1 + jpspec) / s(t)%y(jl_3,js_1)
                end if
                if (peps < s(t)%qa(jl_3,js_1)) then
                  zb(jl_3) = sl - s(t)%qa(jl_3,js_1) * s(t)%y(jl_3,js_1)
                  zc(jl_3) = s(t)%pd(jl_3,js_1)
                  zd(jl_3) = zb(jl_3) * zb(jl_3) + 4.0 * s(t)%qa(jl_3,js_1) * zc(jl_3)
                  if (zd(jl_3) > 0.0) then
                    s(t)%y(jl_3,js_1) = (zb(jl_3) - SQRT(zd(jl_3))) / (-2.0 * s(t)%qa(jl_3,js_1))
                    if (s(t)%y(jl_3,js_1) < 0.0) then
                      s(t)%y(jl_3,js_1) = 10.0 * peps
                    end if
                  else
                    s(t)%y(jl_3,js_1) = zb(jl_3) / (-2.0 * s(t)%qa(jl_3,js_1))
                  end if
                elseif (peps < sl .AND. peps >= s(t)%qa(jl_3,js_1)) then
                  s(t)%y(jl_3,js_1) = s(t)%pd(jl_3,js_1) / sl
                else
                  s(t)%y(jl_3,js_1) = 0.0
                end if
              enddo
            enddo
            istart = nlmajmin(3)
            iend = nlmajmin(4)
            do j_1 = nlmajmin(3), nlmajmin(4), 1
              js_1 = nlmajmin(j_1)
              iodd = nodd(js_1)
              ifam = moffam(js_1)
              imaj = 0
              if (ifam /= 0) then
                imaj = majors(ifam)
              end if
              do jl_3 = 1, n_points, 1
                if (peps < s(t)%y(jl_3,imaj)) then
                  s(t)%ratio(jl_3,js_1) = s(t)%y(jl_3,js_1) / s(t)%y(jl_3,imaj)
                  s(t)%ratio(jl_3,imaj) = s(t)%ratio(jl_3,imaj) + iodd * s(t)%ratio(jl_3,js_1)
                else
                  s(t)%ratio(jl_3,js_1) = 0.0
                end if
              enddo
            enddo
            do j_1 = 1, ift, 1
              js_1 = ilft(j_1)
              iodd = nodd(js_1)
              ifam = moffam(js_1)
              itr_1 = madvtr(js_1)
              imaj = majors(ifam)
              if (not_first_call_1 .AND. jit_2 == 1) then
                do jl_3 = 1, n_points, 1
                  if (peps < s(t)%y(jl_3,js_1)) then
                    sl = s(t)%pd(jl_3,js_1 + jpspec) / s(t)%y(jl_3,js_1)
                  else
                    sl = 0.0
                  end if
                  s(t)%linfam(jl_3,itr_1) = sl > zthresh
                  if (s(t)%linfam(jl_3,itr_1)) then
                    s(t)%f(jl_3,ifam) = s(t)%f(jl_3,ifam) + iodd * s(t)%f(jl_3,itr_1)
                  end if
                enddo
              end if
              do jl_3 = 1, n_points, 1
                if (s(t)%linfam(jl_3,itr_1)) then
                  if (peps < s(t)%y(jl_3,imaj)) then
                    s(t)%ratio(jl_3,js_1) = s(t)%y(jl_3,js_1) / s(t)%y(jl_3,imaj)
                    s(t)%ratio(jl_3,imaj) = s(t)%ratio(jl_3,imaj) + iodd * s(t)%ratio(jl_3,js_1)
                  else
                    s(t)%ratio(jl_3,js_1) = 0.0
                  end if
                else
                  s(t)%y(jl_3,js_1) = s(t)%f(jl_3,itr_1)
                end if
              enddo
            enddo
            istart = nlmajmin(1)
            iend = nlmajmin(2)
            do j_1 = istart, iend, 1
              js_1 = nlmajmin(j_1)
              iodd = nodd(js_1)
              ifam = moffam(js_1)
              do jl_3 = 1, n_points, 1
                s(t)%ratio(jl_3,js_1) = 1.0 / (iodd + s(t)%ratio(jl_3,js_1))
                s(t)%y(jl_3,js_1) = s(t)%f(jl_3,ifam) * s(t)%ratio(jl_3,js_1)
              enddo
            enddo
            istart = nlmajmin(3)
            iend = nlmajmin(4)
            do j_1 = istart, iend, 1
              js_1 = nlmajmin(j_1)
              ifam = moffam(js_1)
              imaj = majors(ifam)
              if (ctype(js_1) /= jpif) then
                do jl_3 = 1, n_points, 1
                  s(t)%y(jl_3,js_1) = s(t)%y(jl_3,imaj) * s(t)%ratio(jl_3,js_1)
                enddo
              else
                itr_1 = madvtr(js_1)
                do jl_3 = 1, n_points, 1
                  if (s(t)%linfam(jl_3,itr_1)) then
                    s(t)%y(jl_3,js_1) = s(t)%y(jl_3,imaj) * s(t)%ratio(jl_3,js_1)
                  end if
                enddo
              end if
            enddo
            gconv = .true.
            do j_1 = 1, nstst, 1
              js_1 = nlstst(j_1)
              do jl_3 = 1, n_points, 1
                if (s(t)%pmintnd(jl_3) < s(t)%y(jl_3,js_1) .AND. ftol * s(t)%y(jl_3,js_1) < ABS(s(t)%y(jl_3,js_1) - zy(jl_3,js_1))) then
                  gconv = .false.
                end if
              enddo
              if (.NOT.gconv) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Exit_Stmt
                EXIT
              end if
            enddo
            if (gconv) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Exit_Stmt
              EXIT
            end if
            do j_1 = 1, nstst, 1
              js_1 = nlstst(j_1)
              do jl_3 = 1, n_points, 1
                zy(jl_3,js_1) = s(t)%y(jl_3,js_1)
              enddo
            enddo
          enddo
          if (nstst /= 0 .AND. ifi == 0) then
            do idx_66 = 1, UBOUND(s(t)%deriv, dim=3), 1
              do idx_67 = 1, UBOUND(s(t)%deriv, dim=2), 1
                do idx_68 = 1, UBOUND(s(t)%deriv, dim=1), 1
                  s(t)%deriv(idx_68,idx_67,idx_66) = 1.0
                enddo
              enddo
            enddo
            do ix_2 = 1, nsst, 1
              do idx_69 = 1, n_points, 1
                ssnum_1(idx_69) = 0.0
              enddo
              do idx_70 = 1, n_points, 1
                ssden_1(idx_70) = 0.0
              enddo
              do idx_71 = 1, 4, 1
                do idx_72 = 1, n_points, 1
                  dssden_1(idx_72,idx_71) = 0.0
                enddo
              enddo
              do idx_73 = 1, 4, 1
                do idx_74 = 1, n_points, 1
                  dssnum_1(idx_74,idx_73) = 0.0
                enddo
              enddo
              do jr_1 = 1, nsspt(ix_2), 1
                i_3 = nsspi(ix_2,jr_1)
                if (i_3 <= nuni) then
                  do idx_75 = 1, n_points, 1
                    ssnum_1(idx_75) = ssnum_1(idx_75) + s(t)%rk(idx_75,i_3) * s(t)%y(idx_75,nspi(i_3,1))
                  enddo
                  if (ix_2 < 5 .AND. nspo3 == nspi(i_3,1)) then
                    do idx_76 = 1, n_points, 1
                      dssnum_1(idx_76,n_o3_1) = dssnum_1(idx_76,n_o3_1) + s(t)%rk(idx_76,i_3)
                    enddo
                  end if
                  if (ix_2 < 5 .AND. nspno == nspi(i_3,1)) then
                    do idx_77 = 1, n_points, 1
                      dssnum_1(idx_77,n_no_1) = dssnum_1(idx_77,n_no_1) + s(t)%rk(idx_77,i_3)
                    enddo
                  end if
                else
                  do idx_78 = 1, n_points, 1
                    ssnum_1(idx_78) = ssnum_1(idx_78) + s(t)%rk(idx_78,i_3) * s(t)%y(idx_78,nspi(i_3,1)) * s(t)%y(idx_78,nspi(i_3,2))
                  enddo
                  if (ix_2 < 5) then
                    if (nspo1d == nspi(i_3,1)) then
                      do idx_79 = 1, n_points, 1
                        dssnum_1(idx_79,n_o3_1) = dssnum_1(idx_79,n_o3_1) + s(t)%rk(idx_79,i_3) * s(t)%y(idx_79,nspi(i_3,2)) * s(t)%deriv(idx_79,nss_o1d,n_o3_1)
                      enddo
                    end if
                    if (nspo1d == nspi(i_3,2)) then
                      do idx_80 = 1, n_points, 1
                        dssnum_1(idx_80,n_o3_1) = dssnum_1(idx_80,n_o3_1) + s(t)%rk(idx_80,i_3) * s(t)%y(idx_80,nspi(i_3,1)) * s(t)%deriv(idx_80,nss_o1d,n_o3_1)
                      enddo
                    end if
                    if (nspo3p == nspi(i_3,1)) then
                      do idx_81 = 1, n_points, 1
                        dssnum_1(idx_81,n_o3_1) = dssnum_1(idx_81,n_o3_1) + s(t)%rk(idx_81,i_3) * s(t)%y(idx_81,nspi(i_3,2)) * s(t)%deriv(idx_81,nss_o3p,n_o3_1)
                      enddo
                    end if
                    if (nspo3p == nspi(i_3,2)) then
                      do idx_82 = 1, n_points, 1
                        dssnum_1(idx_82,n_o3_1) = dssnum_1(idx_82,n_o3_1) + s(t)%rk(idx_82,i_3) * s(t)%y(idx_82,nspi(i_3,1)) * s(t)%deriv(idx_82,nss_o3p,n_o3_1)
                      enddo
                    end if
                    if (nspo3 == nspi(i_3,1)) then
                      do idx_83 = 1, n_points, 1
                        dssnum_1(idx_83,n_o3_1) = dssnum_1(idx_83,n_o3_1) + s(t)%rk(idx_83,i_3) * s(t)%y(idx_83,nspi(i_3,2))
                      enddo
                    end if
                    if (nspo3 == nspi(i_3,2)) then
                      do idx_84 = 1, n_points, 1
                        dssnum_1(idx_84,n_o3_1) = dssnum_1(idx_84,n_o3_1) + s(t)%rk(idx_84,i_3) * s(t)%y(idx_84,nspi(i_3,1))
                      enddo
                    end if
                    if (nspo3p == nspi(i_3,1)) then
                      do idx_85 = 1, n_points, 1
                        dssnum_1(idx_85,n_oh_1) = dssnum_1(idx_85,n_oh_1) + s(t)%rk(idx_85,i_3) * s(t)%y(idx_85,nspi(i_3,2)) * s(t)%deriv(idx_85,nss_o3p,n_oh_1)
                      enddo
                    end if
                    if (nspo3p == nspi(i_3,2)) then
                      do idx_86 = 1, n_points, 1
                        dssnum_1(idx_86,n_oh_1) = dssnum_1(idx_86,n_oh_1) + s(t)%rk(idx_86,i_3) * s(t)%y(idx_86,nspi(i_3,1)) * s(t)%deriv(idx_86,nss_o3p,n_oh_1)
                      enddo
                    end if
                    if (nspoh == nspi(i_3,1)) then
                      do idx_87 = 1, n_points, 1
                        dssnum_1(idx_87,n_oh_1) = dssnum_1(idx_87,n_oh_1) + s(t)%rk(idx_87,i_3) * s(t)%y(idx_87,nspi(i_3,2))
                      enddo
                    end if
                    if (nspoh == nspi(i_3,2)) then
                      do idx_88 = 1, n_points, 1
                        dssnum_1(idx_88,n_oh_1) = dssnum_1(idx_88,n_oh_1) + s(t)%rk(idx_88,i_3) * s(t)%y(idx_88,nspi(i_3,1))
                      enddo
                    end if
                    if (nspo3p == nspi(i_3,1)) then
                      do idx_89 = 1, n_points, 1
                        dssnum_1(idx_89,n_ho2_1) = dssnum_1(idx_89,n_ho2_1) + s(t)%rk(idx_89,i_3) * s(t)%y(idx_89,nspi(i_3,2)) * s(t)%deriv(idx_89,nss_o3p,n_ho2_1)
                      enddo
                    end if
                    if (nspo3p == nspi(i_3,2)) then
                      do idx_90 = 1, n_points, 1
                        dssnum_1(idx_90,n_ho2_1) = dssnum_1(idx_90,n_ho2_1) + s(t)%rk(idx_90,i_3) * s(t)%y(idx_90,nspi(i_3,1)) * s(t)%deriv(idx_90,nss_o3p,n_ho2_1)
                      enddo
                    end if
                    if (nspho2 == nspi(i_3,1)) then
                      do idx_91 = 1, n_points, 1
                        dssnum_1(idx_91,n_ho2_1) = dssnum_1(idx_91,n_ho2_1) + s(t)%rk(idx_91,i_3) * s(t)%y(idx_91,nspi(i_3,2))
                      enddo
                    end if
                    if (nspho2 == nspi(i_3,2)) then
                      do idx_92 = 1, n_points, 1
                        dssnum_1(idx_92,n_ho2_1) = dssnum_1(idx_92,n_ho2_1) + s(t)%rk(idx_92,i_3) * s(t)%y(idx_92,nspi(i_3,1))
                      enddo
                    end if
                    if (nspno == nspi(i_3,1)) then
                      do idx_93 = 1, n_points, 1
                        dssnum_1(idx_93,n_no_1) = dssnum_1(idx_93,n_no_1) + s(t)%rk(idx_93,i_3) * s(t)%y(idx_93,nspi(i_3,2))
                      enddo
                    end if
                    if (nspno == nspi(i_3,2)) then
                      do idx_94 = 1, n_points, 1
                        dssnum_1(idx_94,n_no_1) = dssnum_1(idx_94,n_no_1) + s(t)%rk(idx_94,i_3) * s(t)%y(idx_94,nspi(i_3,1))
                      enddo
                    end if
                  end if
                end if
              enddo
              do jr_1 = 1, nssrt(ix_2), 1
                i_3 = nssri(ix_2,jr_1)
                j_2 = nssrx(ix_2,jr_1)
                if (i_3 <= nuni) then
                  do idx_95 = 1, n_points, 1
                    ssden_1(idx_95) = ssden_1(idx_95) + s(t)%rk(idx_95,i_3)
                  enddo
                else
                  do idx_96 = 1, n_points, 1
                    ssden_1(idx_96) = ssden_1(idx_96) + s(t)%rk(idx_96,i_3) * s(t)%y(idx_96,nspi(i_3,j_2))
                  enddo
                  if (ix_2 < 5) then
                    if (nspo3 == nspi(i_3,j_2)) then
                      do idx_97 = 1, n_points, 1
                        dssden_1(idx_97,n_o3_1) = dssden_1(idx_97,n_o3_1) + s(t)%rk(idx_97,i_3)
                      enddo
                    end if
                    if (nspoh == nspi(i_3,j_2)) then
                      do idx_98 = 1, n_points, 1
                        dssden_1(idx_98,n_oh_1) = dssden_1(idx_98,n_oh_1) + s(t)%rk(idx_98,i_3)
                      enddo
                    end if
                    if (nspho2 == nspi(i_3,j_2)) then
                      do idx_99 = 1, n_points, 1
                        dssden_1(idx_99,n_ho2_1) = dssden_1(idx_99,n_ho2_1) + s(t)%rk(idx_99,i_3)
                      enddo
                    end if
                    if (nspno == nspi(i_3,j_2)) then
                      do idx_100 = 1, n_points, 1
                        dssden_1(idx_100,n_no_1) = dssden_1(idx_100,n_no_1) + s(t)%rk(idx_100,i_3)
                      enddo
                    end if
                  end if
                end if
              enddo
              do idx_101 = 1, n_points, 1
                s(t)%y(idx_101,nssi(ix_2)) = ssnum_1(idx_101) / ssden_1(idx_101)
              enddo
              if (ix_2 < 5) then
                do jr_1 = 1, 4, 1
                  do idx_102 = 1, n_points, 1
                    s(t)%deriv(idx_102,ix_2,jr_1) = (ssden_1(idx_102) * dssnum_1(idx_102,jr_1) - ssnum_1(idx_102) * dssden_1(idx_102,jr_1)) / (ssden_1(idx_102) * ssden_1(idx_102))
                  enddo
                enddo
              end if
            enddo
            do x_1 = 1, n_points, 1
              if (peps < s(t)%y(x_1,nspo1d)) then
                s(t)%deriv(x_1,nss_o1d,n_o3_1) = s(t)%deriv(x_1,nss_o1d,n_o3_1) * s(t)%y(x_1,nspo3) / s(t)%y(x_1,nspo1d)
                s(t)%deriv(x_1,nss_o1d,n_oh_1) = s(t)%deriv(x_1,nss_o1d,n_oh_1) * s(t)%y(x_1,nspoh) / s(t)%y(x_1,nspo1d)
                s(t)%deriv(x_1,nss_o1d,n_ho2_1) = s(t)%deriv(x_1,nss_o1d,n_ho2_1) * s(t)%y(x_1,nspho2) / s(t)%y(x_1,nspo1d)
                s(t)%deriv(x_1,nss_o1d,n_no_1) = s(t)%deriv(x_1,nss_o1d,n_no_1) * s(t)%y(x_1,nspno) / s(t)%y(x_1,nspo1d)
              else
                s(t)%deriv(x_1,nss_o1d,n_o3_1) = 1.0
                s(t)%deriv(x_1,nss_o1d,n_oh_1) = 1.0
                s(t)%deriv(x_1,nss_o1d,n_ho2_1) = 1.0
                s(t)%deriv(x_1,nss_o1d,n_no_1) = 1.0
              end if
            enddo
            if (o3p_in_ss) then
              do x_1 = 1, n_points, 1
                if (peps < s(t)%y(x_1,nspo3p)) then
                  s(t)%deriv(x_1,nss_o3p,n_o3_1) = s(t)%deriv(x_1,nss_o3p,n_o3_1) * s(t)%y(x_1,nspo3) / s(t)%y(x_1,nspo3p)
                  s(t)%deriv(x_1,nss_o3p,n_oh_1) = s(t)%deriv(x_1,nss_o3p,n_oh_1) * s(t)%y(x_1,nspoh) / s(t)%y(x_1,nspo3p)
                  s(t)%deriv(x_1,nss_o3p,n_ho2_1) = s(t)%deriv(x_1,nss_o3p,n_ho2_1) * s(t)%y(x_1,nspho2) / s(t)%y(x_1,nspo3p)
                  s(t)%deriv(x_1,nss_o3p,n_no_1) = s(t)%deriv(x_1,nss_o3p,n_no_1) * s(t)%y(x_1,nspno) / s(t)%y(x_1,nspo3p)
                else
                  s(t)%deriv(x_1,nss_o3p,n_o3_1) = 1.0
                  s(t)%deriv(x_1,nss_o3p,n_oh_1) = 1.0
                  s(t)%deriv(x_1,nss_o3p,n_ho2_1) = 1.0
                  s(t)%deriv(x_1,nss_o3p,n_no_1) = 1.0
                end if
              enddo
            end if
            if (n_in_ss) then
              do x_1 = 1, n_points, 1
                if (peps < s(t)%y(x_1,nspn)) then
                  s(t)%deriv(x_1,nss_n,n_o3_1) = s(t)%deriv(x_1,nss_n,n_o3_1) * s(t)%y(x_1,nspo3) / s(t)%y(x_1,nspn)
                  s(t)%deriv(x_1,nss_n,n_oh_1) = s(t)%deriv(x_1,nss_n,n_oh_1) * s(t)%y(x_1,nspoh) / s(t)%y(x_1,nspn)
                  s(t)%deriv(x_1,nss_n,n_ho2_1) = s(t)%deriv(x_1,nss_n,n_ho2_1) * s(t)%y(x_1,nspho2) / s(t)%y(x_1,nspn)
                  s(t)%deriv(x_1,nss_n,n_no_1) = s(t)%deriv(x_1,nss_n,n_no_1) * s(t)%y(x_1,nspno) / s(t)%y(x_1,nspn)
                else
                  s(t)%deriv(x_1,nss_n,n_o3_1) = 1.0
                  s(t)%deriv(x_1,nss_n,n_oh_1) = 1.0
                  s(t)%deriv(x_1,nss_n,n_ho2_1) = 1.0
                  s(t)%deriv(x_1,nss_n,n_no_1) = 1.0
                end if
              enddo
            end if
            if (h_in_ss) then
              do x_1 = 1, n_points, 1
                if (peps < s(t)%y(x_1,nsph)) then
                  s(t)%deriv(x_1,nss_h,n_o3_1) = s(t)%deriv(x_1,nss_h,n_o3_1) * s(t)%y(x_1,nspo3) / s(t)%y(x_1,nsph)
                  s(t)%deriv(x_1,nss_h,n_oh_1) = s(t)%deriv(x_1,nss_h,n_oh_1) * s(t)%y(x_1,nspoh) / s(t)%y(x_1,nsph)
                  s(t)%deriv(x_1,nss_h,n_ho2_1) = s(t)%deriv(x_1,nss_h,n_ho2_1) * s(t)%y(x_1,nspho2) / s(t)%y(x_1,nsph)
                  s(t)%deriv(x_1,nss_h,n_no_1) = s(t)%deriv(x_1,nss_h,n_no_1) * s(t)%y(x_1,nspno) / s(t)%y(x_1,nsph)
                else
                  s(t)%deriv(x_1,nss_h,n_o3_1) = 1.0
                  s(t)%deriv(x_1,nss_h,n_oh_1) = 1.0
                  s(t)%deriv(x_1,nss_h,n_ho2_1) = 1.0
                  s(t)%deriv(x_1,nss_h,n_no_1) = 1.0
                end if
              enddo
            end if
          end if
          do jtr_2 = 1, jpcspf, 1
            do jl_4 = 1, n_points, 1
              s(t)%fdot(jl_4,jtr_2) = 0.0
            enddo
          enddo
          do js_10 = 1, jpspec, 1
            do jl_15 = 1, n_points, 1
              s(t)%pd(jl_15,js_10) = 0.0
              s(t)%pd(jl_15,js_10 + jpspec) = 0.0
            enddo
          enddo
          do jr_5 = 1, nuni, 1
            ir1_1 = nspi(jr_5,1)
            do jl_15 = 1, n_points, 1
              s(t)%prk(jl_15,jr_5) = s(t)%rk(jl_15,jr_5) * s(t)%y(jl_15,ir1_1)
            enddo
          enddo
          iunip1_1 = nuni + 1
          do jr_5 = iunip1_1, jpnr, 1
            ir1_1 = nspi(jr_5,1)
            ir2_1 = nspi(jr_5,2)
            do jl_15 = 1, n_points, 1
              s(t)%prk(jl_15,jr_5) = s(t)%rk(jl_15,jr_5) * s(t)%y(jl_15,ir1_1) * s(t)%y(jl_15,ir2_1)
            enddo
          enddo
          do jskip_1 = 0, jpspec, jpspec
            do j_14 = 1, nf, 1
              js_10 = nlf(j_14) + jskip_1
              ip3_1 = ngrp(js_10,3)
              ip2_1 = ngrp(js_10,2)
              ip1_1 = ngrp(js_10,1)
              do j3_3 = 1, ip3_1, 1
                i1_1 = nprdx3(1,j3_3,js_10)
                i2_1 = nprdx3(2,j3_3,js_10)
                i3_1 = nprdx3(3,j3_3,js_10)
                do jl_15 = 1, n_points, 1
                  s(t)%pd(jl_15,js_10) = s(t)%pd(jl_15,js_10) + s(t)%prk(jl_15,i1_1) + s(t)%prk(jl_15,i2_1) + s(t)%prk(jl_15,i3_1)
                enddo
              enddo
              do j2_1 = 1, ip2_1, 1
                i1_1 = nprdx2(1,js_10)
                i2_1 = nprdx2(2,js_10)
                do jl_15 = 1, n_points, 1
                  s(t)%pd(jl_15,js_10) = s(t)%pd(jl_15,js_10) + s(t)%prk(jl_15,i1_1) + s(t)%prk(jl_15,i2_1)
                enddo
              enddo
              do j1_1 = 1, ip1_1, 1
                i1_1 = nprdx1(js_10)
                do jl_15 = 1, n_points, 1
                  s(t)%pd(jl_15,js_10) = s(t)%pd(jl_15,js_10) + s(t)%prk(jl_15,i1_1)
                enddo
              enddo
            enddo
          enddo
          do jr_5 = 1, nnfrp, 1
            iss_2 = ntabfp(jr_5,1)
            ir_1 = ntabfp(jr_5,2)
            ix_5 = ntabfp(jr_5,3)
            fr_1 = frpx(ix_5)
            do jl_15 = 1, n_points, 1
              s(t)%pd(jl_15,iss_2) = s(t)%pd(jl_15,iss_2) + fr_1 * s(t)%prk(jl_15,ir_1)
            enddo
          enddo
          if (.true.) then
            istart_5 = nldepx(1)
            iend_5 = nldepx(2)
            do j_14 = istart_5, iend_5, 1
              js_10 = nldepx(j_14)
              do jl_15 = 1, n_points, 1
                s(t)%pd(jl_15,js_10 + jpspec) = s(t)%pd(jl_15,js_10 + jpspec) + (s(t)%dpd(jl_15,js_10) + s(t)%dpw(jl_15,js_10)) * s(t)%y(jl_15,js_10)
              enddo
            enddo
            istart_5 = nldepx(3)
            iend_5 = nldepx(4)
            do j_14 = istart_5, iend_5, 1
              js_10 = nldepx(j_14)
              do jl_15 = 1, n_points, 1
                s(t)%pd(jl_15,js_10 + jpspec) = s(t)%pd(jl_15,js_10 + jpspec) + s(t)%dpd(jl_15,js_10) * s(t)%y(jl_15,js_10)
              enddo
            enddo
            istart_5 = nldepx(5)
            iend_5 = nldepx(6)
            do j_14 = istart_5, iend_5, 1
              js_10 = nldepx(j_14)
              do jl_15 = 1, n_points, 1
                s(t)%pd(jl_15,js_10 + jpspec) = s(t)%pd(jl_15,js_10 + jpspec) + s(t)%dpw(jl_15,js_10) * s(t)%y(jl_15,js_10)
              enddo
            enddo
          end if
          do j_3 = 1, nf, 1
            js_2 = nlf(j_3)
            ifam_1 = moffam(js_2)
            itr_2 = madvtr(js_2)
            gfam = ifam_1 /= 0
            gtr = itr_2 /= 0
            do jl_4 = 1, n_points, 1
              s(t)%ydot(jl_4,js_2) = s(t)%pd(jl_4,js_2) - s(t)%pd(jl_4,js_2 + jpspec)
              if (gfam .AND. (.NOT.gtr .OR. gtr .AND. s(t)%linfam(jl_4,itr_2))) then
                s(t)%fdot(jl_4,ifam_1) = s(t)%fdot(jl_4,ifam_1) + nodd(js_2) * s(t)%ydot(jl_4,js_2)
              end if
              if (gtr) then
                s(t)%fdot(jl_4,itr_2) = s(t)%ydot(jl_4,js_2)
              end if
            enddo
          enddo
          if (error_norm < reltol_error) then
            exit_code = 0
            solver_iter = iter_1 - 1

            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
          if (iter_1 == ukca_config%nrsteps) then
            exit_code = 2

            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
          do idx_103 = 1, jpcspf, 1
            do idx_104 = 1, n_points, 1
              g_f(idx_104,idx_103) = (s(t)%f(idx_104 ,idx_103 ) - f_initial_1(idx_104,idx_103)) * deltt - s(t)%fdot(idx_104 ,idx_103 )
            enddo
          enddo
          residual_error = 0.0
          do jtr_3 = 1, jpcspf, 1
            j_4 = nlf(jtr_3)
            do idx_105 = 1, n_points, 1
              tmprc(idx_105,jtr_3) = s(t)%pd(idx_105,j_4) + s(t)%pd(idx_105,j_4 + jpspec)
            enddo
          enddo
          do jtr_3 = 1, jpcspf, 1
            do jl_5 = 1, n_points, 1
              if (f_min < ABS(tmprc(jl_5,jtr_3))) then
                residual_error = MAX(residual_error, ABS(g_f(jl_5,jtr_3) / tmprc(jl_5,jtr_3)))
              end if
            enddo
          enddo
          if (residual_error < reltol_residual_error) then
            exit_code = 0
            solver_iter = iter_1

            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
          do jl = 1, n_points, 1
            do jtr = 1, jpcspf, 1
              if (s(t)%f(jl,jtr) > f_max) then
                exit_code = 3
              end if
            enddo
          enddo
          if (exit_code == 3) then
            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
          do i_4 = 1, jpcspf, 1
            do jl_6 = 1, n_points, 1
              if (f_min > s(t)%f(jl_6,i_4)) then
                s(t)%f(jl_6,i_4) = f_min
              end if
            enddo
          enddo
          deltt_2 = 1.0 / s(t)%cdt
          do idx_106 = 1, UBOUND(s(t)%spfj, dim=2), 1
            do idx_107 = 1, UBOUND(s(t)%spfj, dim=1), 1
              s(t)%spfj(idx_107,idx_106) = 0.0
            enddo
          enddo
          do i_4 = 1, jpcspf, 1
            do idx_108 = 1, UBOUND(s(t)%spfj, dim=1), 1
              s(t)%spfj(idx_108,nonzero_map(i_4,i_4)) = -deltt_2 * s(t)%f(idx_108 ,i_4)
            enddo
          enddo
          do p_1 = 1, total, 1
            do i_4 = nposterms(p_1), 1, -1
              do idx_109 = 1, UBOUND(s(t)%spfj, dim=1), 1
                s(t)%spfj(idx_109,p_1) = s(t)%spfj(idx_109,p_1) + s(t)%prk(idx_109 ,posterms(p_1,i_4))
              enddo
            enddo
            do i_4 = nnegterms(p_1), 1, -1
              do idx_110 = 1, UBOUND(s(t)%spfj, dim=1), 1
                s(t)%spfj(idx_110,p_1) = s(t)%spfj(idx_110,p_1) - s(t)%prk(idx_110 ,negterms(p_1,i_4))
              enddo
            enddo
            do i_4 = nfracterms(p_1), 1, -1
              do idx_111 = 1, UBOUND(s(t)%spfj, dim=1), 1
                s(t)%spfj(idx_111,p_1) = s(t)%spfj(idx_111,p_1) + ffrac(p_1,i_4) * s(t)%prk(idx_111 ,fracterms(p_1,i_4))
              enddo
            enddo
          enddo
          do jc_1 = 1, nstst, 1
            do j3_1 = 1, nmsjac(jc_1), 1
              irj_1 = nsjac1(j3_1,jc_1)
              do jn_1 = 1, jpmsp, 1
                ij_1 = njcoth(irj_1,jn_1)
                if (ij_1 /= 0) then
                  p_1 = nonzero_map(ij_1,ntro3)
                  if (jn_1 < 3) then
                    do idx_112 = 1, UBOUND(s(t)%spfj, dim=1), 1
                      s(t)%spfj(idx_112,p_1) = s(t)%spfj(idx_112,p_1) - s(t)%prk(idx_112 ,irj_1) * s(t)%deriv(idx_112 ,jc_1,1)
                    enddo
                  else
                    do idx_113 = 1, UBOUND(s(t)%spfj, dim=1), 1
                      s(t)%spfj(idx_113,p_1) = s(t)%spfj(idx_113,p_1) + s(t)%prk(idx_113 ,irj_1) * s(t)%deriv(idx_113 ,jc_1,1)
                    enddo
                  end if
                end if
              enddo
            enddo
          enddo
          if (ndepw /= 0 .OR. ndepd /= 0) then
            do js_3 = 1, jpspec, 1
              ifamd_1 = moffam(js_3)
              itrd_1 = madvtr(js_3)
              ityped_1 = ctype(js_3)
              if (ifamd_1 /= 0) then
                p_1 = nonzero_map(ifamd_1,ifamd_1)
                do jl_6 = 1, n_points, 1
                  if (ityped_1 == jpfm .OR. ityped_1 == jpif .AND. s(t)%linfam(jl_6,itrd_1)) then
                    s(t)%spfj(jl_6,p_1) = s(t)%spfj(jl_6,p_1) - nodd(js_3) * (s(t)%dpd(jl_6,js_3) + s(t)%dpw(jl_6,js_3)) * s(t)%y(jl_6,js_3)
                  end if
                enddo
              end if
              if (itrd_1 /= 0) then
                p_1 = nonzero_map(itrd_1,itrd_1)
                do idx_114 = 1, UBOUND(s(t)%spfj, dim=1), 1
                  s(t)%spfj(idx_114,p_1) = s(t)%spfj(idx_114,p_1) - (s(t)%dpd(idx_114 ,js_3) + s(t)%dpw(idx_114 ,js_3)) * s(t)%y(idx_114 ,js_3)
                enddo
              end if
            enddo
          end if
          do p_1 = 1, total, 1
            do idx_115 = 1, UBOUND(s(t)%spfj, dim=1), 1
              s(t)%spfj(idx_115,p_1) = s(t)%spfj(idx_115,p_1) / s(t)%f(idx_115 ,base_tracer(p_1))
            enddo
          enddo
          do j_5 = 1, total, 1
            do jl_7 = 1, n_points, 1
              s(t)%spfj(jl_7,j_5) = MIN(MAX(s(t)%spfj(jl_7,j_5), -f_max), f_max)
            enddo
          enddo
          total1 = total
          do idx_116 = 1, UBOUND(s(t)%modified_map, dim=2), 1
            do idx_117 = 1, UBOUND(s(t)%modified_map, dim=1), 1
              s(t)%modified_map(idx_117,idx_116) = nonzero_map_unordered(idx_117 ,idx_116 )
            enddo
          enddo
          do kr = 1, jpcspf, 1
            do idx_118 = 1, UBOUND(pivot, dim=1), 1
              pivot(idx_118) = s(t)%spfj(idx_118 ,s(t)%modified_map(kr,kr))
            enddo
            do i_5 = 1, n_points, 1
              if (f_min < ABS(pivot(i_5))) then
                pivot(i_5) = 1.0 / pivot(i_5)
              else
                pivot(i_5) = f_max
              end if
            enddo
            do i_5 = kr + 1, jpcspf, 1
              ikr = s(t)%modified_map(i_5,kr)
              if (ikr > 0) then
                do idx_119 = 1, n_points, 1
                  kfact(idx_119) = s(t)%spfj(idx_119 ,ikr) * pivot(idx_119)
                enddo
                do idx_120 = 1, UBOUND(s(t)%spfj, dim=1), 1
                  s(t)%spfj(idx_120,ikr) = kfact(idx_120 )
                enddo
                do j_5 = kr + 1, jpcspf, 1
                  krj = s(t)%modified_map(kr,j_5)
                  if (krj > 0) then
                    ij_2 = s(t)%modified_map(i_5,j_5)
                    if (ij_2 > 0) then
                      do idx_121 = 1, UBOUND(s(t)%spfj, dim=1), 1
                        s(t)%spfj(idx_121,ij_2) = s(t)%spfj(idx_121,ij_2) - kfact(idx_121 ) * s(t)%spfj(idx_121,krj)
                      enddo
                    else
                      total1 = total1 + 1
                      s(t)%modified_map(i_5,j_5) = total1
                      do idx_122 = 1, UBOUND(s(t)%spfj, dim=1), 1
                        s(t)%spfj(idx_122,total1) = -kfact(idx_122 ) * s(t)%spfj(idx_122,krj)
                      enddo
                    end if
                  end if
                enddo
              end if
            enddo
          enddo
          do j_5 = 1, total1, 1
            do jl_7 = 1, n_points, 1
              s(t)%spfj(jl_7,j_5) = MIN(MAX(s(t)%spfj(jl_7,j_5), -f_max), f_max)
            enddo
          enddo
          do i_8 = 1, jpcspf, 1
            do idx_123 = 1, UBOUND(bb1_2, dim=1), 1
              bb1_2(idx_123,i_8) = g_f(idx_123,reorder(i_8))
            enddo
          enddo
          do kr_2 = 1, jpcspf - 1, 1
            do i_8 = kr_2 + 1, jpcspf, 1
              ikr_2 = s(t)%modified_map(i_8,kr_2)
              if (ikr_2 > 0) then
                do idx_124 = 1, UBOUND(bb1_2, dim=1), 1
                  bb1_2(idx_124,i_8) = bb1_2(idx_124,i_8) - s(t)%spfj(idx_124 ,ikr_2) * bb1_2(idx_124,kr_2)
                enddo
              end if
            enddo
          enddo
          do kr_2 = jpcspf, 1, -1
            do idx_125 = 1, UBOUND(xx1_2, dim=1), 1
              xx1_2(idx_125,kr_2) = bb1_2(idx_125,kr_2)
            enddo
            do j_15 = kr_2 + 1, jpcspf, 1
              krj_2 = s(t)%modified_map(kr_2,j_15)
              if (krj_2 > 0) then
                do idx_126 = 1, UBOUND(xx1_2, dim=1), 1
                  xx1_2(idx_126,kr_2) = xx1_2(idx_126,kr_2) - s(t)%spfj(idx_126 ,krj_2) * xx1_2(idx_126,j_15)
                enddo
              end if
            enddo
            do idx_127 = 1, UBOUND(pivot_2, dim=1), 1
              pivot_2(idx_127) = s(t)%spfj(idx_127 ,s(t)%modified_map(kr_2,kr_2))
            enddo
            do i_8 = 1, n_points, 1
              if (f_min > ABS(pivot_2(i_8))) then
                pivot_2(i_8) = f_min
              end if
            enddo
            if (f_max > 0) then
              do idx_128 = 1, UBOUND(xx1_2, dim=1), 1
                xx1_2(idx_128,kr_2) = MIN(MAX(xx1_2(idx_128,kr_2), -f_max), f_max) / pivot_2(idx_128)
              enddo
            else
              do idx_129 = 1, UBOUND(xx1_2, dim=1), 1
                xx1_2(idx_129,kr_2) = xx1_2(idx_129,kr_2) / pivot_2(idx_129)
              enddo
            end if
          enddo
          do idx_130 = 1, UBOUND(f_incr, dim=2), 1
            do idx_131 = 1, UBOUND(f_incr, dim=1), 1
              f_incr(idx_131,idx_130) = 0.0
            enddo
          enddo
          do i_8 = 1, jpcspf, 1
            do idx_132 = 1, UBOUND(f_incr, dim=1), 1
              f_incr(idx_132,reorder(i_8)) = xx1_2(idx_132,i_8)
            enddo
          enddo
          if (iter_1 == 1) then
            do idx_133 = 1, jpcspf, 1
              do idx_134 = 1, n_points, 1
                f_incr(idx_134,idx_133) = damp1 * f_incr(idx_134,idx_133)
              enddo
            enddo
          end if
          do idx_135 = 1, jpcspf, 1
            do idx_136 = 1, n_points, 1
              f_incr(idx_136,idx_135) = MIN(MAX(f_incr(idx_136,idx_135), -f_max), f_max)
            enddo
          enddo
          error_norm = 0.0
          do jtr_4 = 1, jpcspf, 1
            do jl_8 = 1, n_points, 1
              if (ABS(f_incr(jl_8,jtr_4)) > 1.0e-16) then
                error_norm = MAX(error_norm, ABS(f_incr(jl_8,jtr_4) / MAX(s(t)%f(jl_8,jtr_4), f_min)))
              end if
            enddo
          enddo
          count_negatives = 0
          do jtr = 1, jpcspf, 1
            do jl = 1, n_points, 1
              ztmp = s(t)%f(jl,jtr) + f_incr(jl,jtr)
              if (iter_1 < incr_limiter) then
                ztmp = MAX(rafmin * s(t)%f(jl,jtr), MIN(rafmax * s(t)%f(jl,jtr), ztmp))
              end if
              if (ztmp == 0.0) then
                ztmp = f_min
              end if
              if (ztmp < 0.0) then
                ztmp = f_min
                count_negatives = count_negatives + 1
              end if
              s(t)%f(jl,jtr) = ztmp
            enddo
          enddo
          if (count_negatives > maxneg) then
            exit_code = 2

            ! PSyclone CodeBlock (unsupported code) reason:
            !  - Unsupported statement: Goto_Stmt
            GO TO 9999
          end if
          if (ukca_config%l_ukca_quasinewton .AND. error_norm >= reltol_error) then
            if (iter_1 >= ukca_config%i_ukca_quasinewton_start .AND. iter_1 <= ukca_config%i_ukca_quasinewton_end) then
              ifam_2 = 0
              imaj_1 = 0
              sl_1 = 0.0
              zthresh_1 = 2.0 / s(t)%cdt
              do idx_137 = 1, jpspec, 1
                do idx_138 = 1, n_points, 1
                  zy_1(idx_138,idx_137) = 0.0
                enddo
              enddo
              do idx_139 = 1, n_points, 1
                zb_1(idx_139) = 0.0
              enddo
              do idx_140 = 1, n_points, 1
                zc_1(idx_140) = 0.0
              enddo
              do idx_141 = 1, n_points, 1
                zd_1(idx_141) = 0.0
              enddo
              if (not_first_call_1) then
                peps10_1 = 10.0 * peps
                istart_6 = nlmajmin(1)
                iend_6 = nlmajmin(2)
                do j_16 = istart_6, iend_6, 1
                  js_11 = nlmajmin(j_16)
                  ifam_6 = moffam(js_11)
                  do jl_16 = 1, n_points, 1
                    s(t)%y(jl_16,js_11) = s(t)%f(jl_16,ifam_6)
                  enddo
                enddo
                istart_6 = nlmajmin(3)
                iend_6 = nlmajmin(4)
                do j_16 = istart_6, iend_6, 1
                  js_11 = nlmajmin(j_16)
                  do jl_16 = 1, n_points, 1
                    s(t)%y(jl_16,js_11) = peps10_1
                  enddo
                enddo
                do j_16 = 1, iss, 1
                  js_11 = ilss(j_16)
                  do jl_16 = 1, n_points, 1
                    s(t)%y(jl_16,js_11) = peps10_1
                  enddo
                enddo
                fco2_1 = fco2_default_1
                fch4_1 = fch4_default_1
                fh2_1 = fh2_default_1
                fo2_1 = fo2_default_1
                fn2_1 = fn2_default_1
                if (ukca_config%l_chem_environ_co2_scalar) then
                  fco2_1 = atmospheric_co2 / c_co2
                  if (ukca_config%l_chem_environ_co2_fld) then
                    fco2_1 = 0.0
                  end if
                end if
                if (ukca_config%l_chem_environ_h2_scalar) then
                  fh2_1 = atmospheric_h2 / c_h2
                end if
                if (ukca_config%l_chem_environ_n2_scalar) then
                  fn2_1 = atmospheric_n2 / c_n2
                end if
                if (ukca_config%l_chem_environ_o2_scalar) then
                  fo2_1 = atmospheric_o2 / c_o2
                end if
                if (ukca_config%l_chem_environ_ch4_scalar) then
                  fch4_1 = atmospheric_ch4 / c_ch4
                end if
                do j_16 = 1, ict, 1
                  js_11 = ilct(j_16)
                  if (speci(js_11) == 'N2        ') then
                    zcnst_1 = fn2_1
                  elseif (speci(js_11) == 'O2        ') then
                    zcnst_1 = fo2_1
                  elseif (speci(js_11) == 'CO2       ') then
                    zcnst_1 = fco2_1
                  elseif (speci(js_11) == 'H2        ') then
                    zcnst_1 = fh2_1
                  elseif (speci(js_11) == 'CH4       ') then
                    zcnst_1 = fch4_1
                  else
                    errcode_3 = -1
                    zcnst_1 = 0.0
                  end if
                  do jl_16 = 1, n_points, 1
                    s(t)%y(jl_16,js_11) = zcnst_1 * s(t)%tnd(jl_16)
                  enddo
                enddo
              end if
              do j_16 = 1, iftr, 1
                js_11 = ilftr(j_16)
                itr_7 = madvtr(js_11)
                do jl_16 = 1, n_points, 1
                  s(t)%y(jl_16,js_11) = s(t)%f(jl_16,itr_7)
                enddo
              enddo
              if (not_first_call_1) then
                do j_16 = 1, icf, 1
                  js_11 = ilcf(j_16)
                  species_tmp_1 = speci(js_11)
                  row_length_1 = ukca_config%row_length
                  rows_1 = ukca_config%rows
                  if (species_tmp_1(1:4) == 'CO2 ' .AND. ukca_config%l_chem_environ_co2_fld) then
                    do idx_142 = 1, n_points, 1
                      s(t)%y(idx_142,js_11) = s(t)%co2(idx_142 )
                    enddo
                  elseif (species_tmp_1(1:4) == 'H2O ') then
                    do idx_143 = 1, n_points, 1
                      s(t)%y(idx_143,js_11) = s(t)%wp(idx_143 )
                    enddo
                  elseif (species_tmp_1(1:4) == 'RO2 ') then
                    if (ukca_config%l_ukca_ro2_perm) then
                      do idx_144 = 1, UBOUND(fro2_3, dim=1), 1
                        fro2_3(idx_144) = 0.0
                      enddo
                      do j_22 = 1, jpro2, 1
                        iro2_3 = nlfro2(j_22)
                        do idx_145 = 1, UBOUND(fro2_3, dim=1), 1
                          fro2_3(idx_145) = fro2_3(idx_145) + s(t)%f(idx_145 ,iro2_3)
                        enddo
                      enddo
                      do idx_146 = 1, n_points, 1
                        s(t)%y(idx_146,js_11) = fro2_3(idx_146) / s(t)%tnd(idx_146 )
                      enddo
                    else
                      errcode_5 = 126
                    end if
                  else
                    errcode_5 = 124
                  end if
                  do idx_147 = 1, n_points, 1
                    s(t)%y(idx_147,js_11) = s(t)%y(idx_147,js_11) * s(t)%tnd(idx_147 )
                  enddo
                enddo
              end if
              if (ukca_config%l_ukca_ro2_perm) then
                do idx_148 = 1, n_points, 1
                  fro2_old_1(idx_148) = s(t)%y(idx_148,jsro2)
                enddo
                do idx_149 = 1, n_points, 1
                  fro2_1(idx_149) = 0.0
                enddo
                do j_6 = 1, jpro2, 1
                  iro2_1 = nlfro2(j_6)
                  do idx_150 = 1, n_points, 1
                    fro2_1(idx_150) = fro2_1(idx_150) + s(t)%f(idx_150,iro2_1)
                  enddo
                enddo
                if (iter_1 <= 0) then
                  damp_ro2_1 = 1.0
                elseif (iter_1 == 1) then
                  damp_ro2_1 = 0.5
                elseif (iter_1 == 2) then
                  damp_ro2_1 = 0.9
                else
                  damp_ro2_1 = 0.99
                end if
                do idx_151 = 1, n_points, 1
                  s(t)%y(idx_151,jsro2) = fro2_old_1(idx_151) + damp_ro2_1 * (fro2_1(idx_151) - fro2_old_1(idx_151))
                enddo
              end if
              if (nstst == 0) then
              end if
              if (not_first_call_1 .AND. ifi < 5) then
                ifi = 5
                icode_1 = -1
              end if
              if (not_first_call_1) then
                do j_17 = 1, nstst, 1
                  js_12 = nlstst(j_17)
                  do jl_17 = 1, n_points, 1
                    s(t)%qa(jl_17,js_12) = 0.0
                  enddo
                enddo
                do jr_6 = 1, jpnr, 1
                  isp_1 = nspi(jr_6,1)
                  if (isp_1 == nspi(jr_6,2)) then
                    do jl_17 = 1, n_points, 1
                      s(t)%qa(jl_17,isp_1) = s(t)%qa(jl_17,isp_1) + 2.0 * s(t)%rk(jl_17,jr_6)
                    enddo
                  end if
                enddo
              end if
              if (ifi == 0) then
                istart_7 = nlmajmin(1)
                iend_7 = nlmajmin(2)
                do j_18 = istart_7, iend_7, 1
                  js_13 = nlmajmin(j_18)
                  ifam_7 = moffam(js_13)
                  do jl_18 = 1, n_points, 1
                    s(t)%y(jl_18,js_13) = s(t)%f(jl_18,ifam_7) * s(t)%ratio(jl_18,js_13)
                  enddo
                enddo
                istart_7 = nlmajmin(3)
                iend_7 = nlmajmin(4)
                do j_18 = istart_7, iend_7, 1
                  js_13 = nlmajmin(j_18)
                  ifam_7 = moffam(js_13)
                  imaj_3 = majors(ifam_7)
                  if (ctype(js_13) /= jpif) then
                    do jl_18 = 1, n_points, 1
                      s(t)%y(jl_18,js_13) = s(t)%y(jl_18,imaj_3) * s(t)%ratio(jl_18,js_13)
                    enddo
                  else
                    itr_8 = madvtr(js_13)
                    do jl_18 = 1, n_points, 1
                      if (s(t)%linfam(jl_18,itr_8)) then
                        s(t)%y(jl_18,js_13) = s(t)%y(jl_18,imaj_3) * s(t)%ratio(jl_18,js_13)
                      end if
                    enddo
                  end if
                enddo
              end if
              do jit_1 = 1, ifi, 1
                if (jit_1 > 1) then
                  do j_6 = 1, nstst, 1
                    js_4 = nlstst(j_6)
                    do jl_9 = 1, n_points, 1
                      if (peps < s(t)%pd(jl_9,js_4) .AND. peps > ABS(s(t)%y(jl_9,js_4))) then
                        s(t)%y(jl_9,js_4) = peps
                      end if
                    enddo
                  enddo
                end if
                do js_14 = 1, jpspec, 1
                  do jl_19 = 1, n_points, 1
                    s(t)%pd(jl_19,js_14) = 0.0
                    s(t)%pd(jl_19,js_14 + jpspec) = 0.0
                  enddo
                enddo
                do jr_7 = 1, nuni, 1
                  ir1_2 = nspi(jr_7,1)
                  do jl_19 = 1, n_points, 1
                    s(t)%prk(jl_19,jr_7) = s(t)%rk(jl_19,jr_7) * s(t)%y(jl_19,ir1_2)
                  enddo
                enddo
                iunip1_2 = nuni + 1
                do jr_7 = iunip1_2, jpnr, 1
                  ir1_2 = nspi(jr_7,1)
                  ir2_2 = nspi(jr_7,2)
                  do jl_19 = 1, n_points, 1
                    s(t)%prk(jl_19,jr_7) = s(t)%rk(jl_19,jr_7) * s(t)%y(jl_19,ir1_2) * s(t)%y(jl_19,ir2_2)
                  enddo
                enddo
                do jskip_2 = 0, jpspec, jpspec
                  do j_19 = 1, istmin, 1
                    js_14 = ilstmin(j_19) + jskip_2
                    ip3_2 = ngrp(js_14,3)
                    ip2_2 = ngrp(js_14,2)
                    ip1_2 = ngrp(js_14,1)
                    do j3_4 = 1, ip3_2, 1
                      i1_2 = nprdx3(1,j3_4,js_14)
                      i2_2 = nprdx3(2,j3_4,js_14)
                      i3_2 = nprdx3(3,j3_4,js_14)
                      do jl_19 = 1, n_points, 1
                        s(t)%pd(jl_19,js_14) = s(t)%pd(jl_19,js_14) + s(t)%prk(jl_19,i1_2) + s(t)%prk(jl_19,i2_2) + s(t)%prk(jl_19,i3_2)
                      enddo
                    enddo
                    do j2_2 = 1, ip2_2, 1
                      i1_2 = nprdx2(1,js_14)
                      i2_2 = nprdx2(2,js_14)
                      do jl_19 = 1, n_points, 1
                        s(t)%pd(jl_19,js_14) = s(t)%pd(jl_19,js_14) + s(t)%prk(jl_19,i1_2) + s(t)%prk(jl_19,i2_2)
                      enddo
                    enddo
                    do j1_2 = 1, ip1_2, 1
                      i1_2 = nprdx1(js_14)
                      do jl_19 = 1, n_points, 1
                        s(t)%pd(jl_19,js_14) = s(t)%pd(jl_19,js_14) + s(t)%prk(jl_19,i1_2)
                      enddo
                    enddo
                  enddo
                enddo
                do jr_7 = 1, nnfrp, 1
                  iss_3 = ntabfp(jr_7,1)
                  ir_2 = ntabfp(jr_7,2)
                  ix_6 = ntabfp(jr_7,3)
                  fr_2 = frpx(ix_6)
                  do jl_19 = 1, n_points, 1
                    s(t)%pd(jl_19,iss_3) = s(t)%pd(jl_19,iss_3) + fr_2 * s(t)%prk(jl_19,ir_2)
                  enddo
                enddo
                if (.false.) then
                  istart_8 = nldepx(1)
                  iend_8 = nldepx(2)
                  do j_19 = istart_8, iend_8, 1
                    js_14 = nldepx(j_19)
                    do jl_19 = 1, n_points, 1
                      s(t)%pd(jl_19,js_14 + jpspec) = s(t)%pd(jl_19,js_14 + jpspec) + (s(t)%dpd(jl_19,js_14) + s(t)%dpw(jl_19,js_14)) * s(t)%y(jl_19,js_14)
                    enddo
                  enddo
                  istart_8 = nldepx(3)
                  iend_8 = nldepx(4)
                  do j_19 = istart_8, iend_8, 1
                    js_14 = nldepx(j_19)
                    do jl_19 = 1, n_points, 1
                      s(t)%pd(jl_19,js_14 + jpspec) = s(t)%pd(jl_19,js_14 + jpspec) + s(t)%dpd(jl_19,js_14) * s(t)%y(jl_19,js_14)
                    enddo
                  enddo
                  istart_8 = nldepx(5)
                  iend_8 = nldepx(6)
                  do j_19 = istart_8, iend_8, 1
                    js_14 = nldepx(j_19)
                    do jl_19 = 1, n_points, 1
                      s(t)%pd(jl_19,js_14 + jpspec) = s(t)%pd(jl_19,js_14 + jpspec) + s(t)%dpw(jl_19,js_14) * s(t)%y(jl_19,js_14)
                    enddo
                  enddo
                end if
                do j_6 = 1, nstst, 1
                  js_4 = nlstst(j_6)
                  do jl_9 = 1, n_points, 1
                    s(t)%ratio(jl_9,js_4) = 0.0
                  enddo
                enddo
                do j_6 = 1, istmin, 1
                  js_4 = ilstmin(j_6)
                  do jl_9 = 1, n_points, 1
                    if (peps > s(t)%y(jl_9,js_4)) then
                      sl_1 = 0.0
                    else
                      sl_1 = s(t)%pd(jl_9,js_4 + jpspec) / s(t)%y(jl_9,js_4)
                    end if
                    if (peps < s(t)%qa(jl_9,js_4)) then
                      zb_1(jl_9) = sl_1 - s(t)%qa(jl_9,js_4) * s(t)%y(jl_9,js_4)
                      zc_1(jl_9) = s(t)%pd(jl_9,js_4)
                      zd_1(jl_9) = zb_1(jl_9) * zb_1(jl_9) + 4.0 * s(t)%qa(jl_9,js_4) * zc_1(jl_9)
                      if (zd_1(jl_9) > 0.0) then
                        s(t)%y(jl_9,js_4) = (zb_1(jl_9) - SQRT(zd_1(jl_9))) / (-2.0 * s(t)%qa(jl_9,js_4))
                        if (s(t)%y(jl_9,js_4) < 0.0) then
                          s(t)%y(jl_9,js_4) = 10.0 * peps
                        end if
                      else
                        s(t)%y(jl_9,js_4) = zb_1(jl_9) / (-2.0 * s(t)%qa(jl_9,js_4))
                      end if
                    elseif (peps < sl_1 .AND. peps >= s(t)%qa(jl_9,js_4)) then
                      s(t)%y(jl_9,js_4) = s(t)%pd(jl_9,js_4) / sl_1
                    else
                      s(t)%y(jl_9,js_4) = 0.0
                    end if
                  enddo
                enddo
                istart_1 = nlmajmin(3)
                iend_1 = nlmajmin(4)
                do j_6 = nlmajmin(3), nlmajmin(4), 1
                  js_4 = nlmajmin(j_6)
                  iodd_1 = nodd(js_4)
                  ifam_2 = moffam(js_4)
                  imaj_1 = 0
                  if (ifam_2 /= 0) then
                    imaj_1 = majors(ifam_2)
                  end if
                  do jl_9 = 1, n_points, 1
                    if (peps < s(t)%y(jl_9,imaj_1)) then
                      s(t)%ratio(jl_9,js_4) = s(t)%y(jl_9,js_4) / s(t)%y(jl_9,imaj_1)
                      s(t)%ratio(jl_9,imaj_1) = s(t)%ratio(jl_9,imaj_1) + iodd_1 * s(t)%ratio(jl_9,js_4)
                    else
                      s(t)%ratio(jl_9,js_4) = 0.0
                    end if
                  enddo
                enddo
                do j_6 = 1, ift, 1
                  js_4 = ilft(j_6)
                  iodd_1 = nodd(js_4)
                  ifam_2 = moffam(js_4)
                  itr_3 = madvtr(js_4)
                  imaj_1 = majors(ifam_2)
                  if (not_first_call_1 .AND. jit_1 == 1) then
                    do jl_9 = 1, n_points, 1
                      if (peps < s(t)%y(jl_9,js_4)) then
                        sl_1 = s(t)%pd(jl_9,js_4 + jpspec) / s(t)%y(jl_9,js_4)
                      else
                        sl_1 = 0.0
                      end if
                      s(t)%linfam(jl_9,itr_3) = sl_1 > zthresh_1
                      if (s(t)%linfam(jl_9,itr_3)) then
                        s(t)%f(jl_9,ifam_2) = s(t)%f(jl_9,ifam_2) + iodd_1 * s(t)%f(jl_9,itr_3)
                      end if
                    enddo
                  end if
                  do jl_9 = 1, n_points, 1
                    if (s(t)%linfam(jl_9,itr_3)) then
                      if (peps < s(t)%y(jl_9,imaj_1)) then
                        s(t)%ratio(jl_9,js_4) = s(t)%y(jl_9,js_4) / s(t)%y(jl_9,imaj_1)
                        s(t)%ratio(jl_9,imaj_1) = s(t)%ratio(jl_9,imaj_1) + iodd_1 * s(t)%ratio(jl_9,js_4)
                      else
                        s(t)%ratio(jl_9,js_4) = 0.0
                      end if
                    else
                      s(t)%y(jl_9,js_4) = s(t)%f(jl_9,itr_3)
                    end if
                  enddo
                enddo
                istart_1 = nlmajmin(1)
                iend_1 = nlmajmin(2)
                do j_6 = istart_1, iend_1, 1
                  js_4 = nlmajmin(j_6)
                  iodd_1 = nodd(js_4)
                  ifam_2 = moffam(js_4)
                  do jl_9 = 1, n_points, 1
                    s(t)%ratio(jl_9,js_4) = 1.0 / (iodd_1 + s(t)%ratio(jl_9,js_4))
                    s(t)%y(jl_9,js_4) = s(t)%f(jl_9,ifam_2) * s(t)%ratio(jl_9,js_4)
                  enddo
                enddo
                istart_1 = nlmajmin(3)
                iend_1 = nlmajmin(4)
                do j_6 = istart_1, iend_1, 1
                  js_4 = nlmajmin(j_6)
                  ifam_2 = moffam(js_4)
                  imaj_1 = majors(ifam_2)
                  if (ctype(js_4) /= jpif) then
                    do jl_9 = 1, n_points, 1
                      s(t)%y(jl_9,js_4) = s(t)%y(jl_9,imaj_1) * s(t)%ratio(jl_9,js_4)
                    enddo
                  else
                    itr_3 = madvtr(js_4)
                    do jl_9 = 1, n_points, 1
                      if (s(t)%linfam(jl_9,itr_3)) then
                        s(t)%y(jl_9,js_4) = s(t)%y(jl_9,imaj_1) * s(t)%ratio(jl_9,js_4)
                      end if
                    enddo
                  end if
                enddo
                gconv_1 = .true.
                do j_6 = 1, nstst, 1
                  js_4 = nlstst(j_6)
                  do jl_9 = 1, n_points, 1
                    if (s(t)%pmintnd(jl_9) < s(t)%y(jl_9,js_4) .AND. ftol * s(t)%y(jl_9,js_4) < ABS(s(t)%y(jl_9,js_4) - zy_1(jl_9,js_4))) then
                      gconv_1 = .false.
                    end if
                  enddo
                  if (.NOT.gconv_1) then
                    ! PSyclone CodeBlock (unsupported code) reason:
                    !  - Unsupported statement: Exit_Stmt
                    EXIT
                  end if
                enddo
                if (gconv_1) then
                  ! PSyclone CodeBlock (unsupported code) reason:
                  !  - Unsupported statement: Exit_Stmt
                  EXIT
                end if
                do j_6 = 1, nstst, 1
                  js_4 = nlstst(j_6)
                  do jl_9 = 1, n_points, 1
                    zy_1(jl_9,js_4) = s(t)%y(jl_9,js_4)
                  enddo
                enddo
              enddo
              do idx_152 = 1, UBOUND(s(t)%deriv, dim=3), 1
                do idx_153 = 1, UBOUND(s(t)%deriv, dim=2), 1
                  do idx_154 = 1, UBOUND(s(t)%deriv, dim=1), 1
                    s(t)%deriv(idx_154,idx_153,idx_152) = 1.0
                  enddo
                enddo
              enddo
              do ix_3 = 1, nsst, 1
                do idx_155 = 1, n_points, 1
                  ssnum_2(idx_155) = 0.0
                enddo
                do idx_156 = 1, n_points, 1
                  ssden_2(idx_156) = 0.0
                enddo
                do idx_157 = 1, 4, 1
                  do idx_158 = 1, n_points, 1
                    dssden_2(idx_158,idx_157) = 0.0
                  enddo
                enddo
                do idx_159 = 1, 4, 1
                  do idx_160 = 1, n_points, 1
                    dssnum_2(idx_160,idx_159) = 0.0
                  enddo
                enddo
                do jr_2 = 1, nsspt(ix_3), 1
                  i_6 = nsspi(ix_3,jr_2)
                  if (i_6 <= nuni) then
                    do idx_161 = 1, n_points, 1
                      ssnum_2(idx_161) = ssnum_2(idx_161) + s(t)%rk(idx_161,i_6) * s(t)%y(idx_161,nspi(i_6,1))
                    enddo
                    if (ix_3 < 5 .AND. nspo3 == nspi(i_6,1)) then
                      do idx_162 = 1, n_points, 1
                        dssnum_2(idx_162,n_o3_2) = dssnum_2(idx_162,n_o3_2) + s(t)%rk(idx_162,i_6)
                      enddo
                    end if
                    if (ix_3 < 5 .AND. nspno == nspi(i_6,1)) then
                      do idx_163 = 1, n_points, 1
                        dssnum_2(idx_163,n_no_2) = dssnum_2(idx_163,n_no_2) + s(t)%rk(idx_163,i_6)
                      enddo
                    end if
                  else
                    do idx_164 = 1, n_points, 1
                      ssnum_2(idx_164) = ssnum_2(idx_164) + s(t)%rk(idx_164,i_6) * s(t)%y(idx_164,nspi(i_6,1)) * s(t)%y(idx_164,nspi(i_6,2))
                    enddo
                    if (ix_3 < 5) then
                      if (nspo1d == nspi(i_6,1)) then
                        do idx_165 = 1, n_points, 1
                          dssnum_2(idx_165,n_o3_2) = dssnum_2(idx_165,n_o3_2) + s(t)%rk(idx_165,i_6) * s(t)%y(idx_165,nspi(i_6,2)) * s(t)%deriv(idx_165,nss_o1d,n_o3_2)
                        enddo
                      end if
                      if (nspo1d == nspi(i_6,2)) then
                        do idx_166 = 1, n_points, 1
                          dssnum_2(idx_166,n_o3_2) = dssnum_2(idx_166,n_o3_2) + s(t)%rk(idx_166,i_6) * s(t)%y(idx_166,nspi(i_6,1)) * s(t)%deriv(idx_166,nss_o1d,n_o3_2)
                        enddo
                      end if
                      if (nspo3p == nspi(i_6,1)) then
                        do idx_167 = 1, n_points, 1
                          dssnum_2(idx_167,n_o3_2) = dssnum_2(idx_167,n_o3_2) + s(t)%rk(idx_167,i_6) * s(t)%y(idx_167,nspi(i_6,2)) * s(t)%deriv(idx_167,nss_o3p,n_o3_2)
                        enddo
                      end if
                      if (nspo3p == nspi(i_6,2)) then
                        do idx_168 = 1, n_points, 1
                          dssnum_2(idx_168,n_o3_2) = dssnum_2(idx_168,n_o3_2) + s(t)%rk(idx_168,i_6) * s(t)%y(idx_168,nspi(i_6,1)) * s(t)%deriv(idx_168,nss_o3p,n_o3_2)
                        enddo
                      end if
                      if (nspo3 == nspi(i_6,1)) then
                        do idx_169 = 1, n_points, 1
                          dssnum_2(idx_169,n_o3_2) = dssnum_2(idx_169,n_o3_2) + s(t)%rk(idx_169,i_6) * s(t)%y(idx_169,nspi(i_6,2))
                        enddo
                      end if
                      if (nspo3 == nspi(i_6,2)) then
                        do idx_170 = 1, n_points, 1
                          dssnum_2(idx_170,n_o3_2) = dssnum_2(idx_170,n_o3_2) + s(t)%rk(idx_170,i_6) * s(t)%y(idx_170,nspi(i_6,1))
                        enddo
                      end if
                      if (nspo3p == nspi(i_6,1)) then
                        do idx_171 = 1, n_points, 1
                          dssnum_2(idx_171,n_oh_2) = dssnum_2(idx_171,n_oh_2) + s(t)%rk(idx_171,i_6) * s(t)%y(idx_171,nspi(i_6,2)) * s(t)%deriv(idx_171,nss_o3p,n_oh_2)
                        enddo
                      end if
                      if (nspo3p == nspi(i_6,2)) then
                        do idx_172 = 1, n_points, 1
                          dssnum_2(idx_172,n_oh_2) = dssnum_2(idx_172,n_oh_2) + s(t)%rk(idx_172,i_6) * s(t)%y(idx_172,nspi(i_6,1)) * s(t)%deriv(idx_172,nss_o3p,n_oh_2)
                        enddo
                      end if
                      if (nspoh == nspi(i_6,1)) then
                        do idx_173 = 1, n_points, 1
                          dssnum_2(idx_173,n_oh_2) = dssnum_2(idx_173,n_oh_2) + s(t)%rk(idx_173,i_6) * s(t)%y(idx_173,nspi(i_6,2))
                        enddo
                      end if
                      if (nspoh == nspi(i_6,2)) then
                        do idx_174 = 1, n_points, 1
                          dssnum_2(idx_174,n_oh_2) = dssnum_2(idx_174,n_oh_2) + s(t)%rk(idx_174,i_6) * s(t)%y(idx_174,nspi(i_6,1))
                        enddo
                      end if
                      if (nspo3p == nspi(i_6,1)) then
                        do idx_175 = 1, n_points, 1
                          dssnum_2(idx_175,n_ho2_2) = dssnum_2(idx_175,n_ho2_2) + s(t)%rk(idx_175,i_6) * s(t)%y(idx_175,nspi(i_6,2)) * s(t)%deriv(idx_175,nss_o3p,n_ho2_2)
                        enddo
                      end if
                      if (nspo3p == nspi(i_6,2)) then
                        do idx_176 = 1, n_points, 1
                          dssnum_2(idx_176,n_ho2_2) = dssnum_2(idx_176,n_ho2_2) + s(t)%rk(idx_176,i_6) * s(t)%y(idx_176,nspi(i_6,1)) * s(t)%deriv(idx_176,nss_o3p,n_ho2_2)
                        enddo
                      end if
                      if (nspho2 == nspi(i_6,1)) then
                        do idx_177 = 1, n_points, 1
                          dssnum_2(idx_177,n_ho2_2) = dssnum_2(idx_177,n_ho2_2) + s(t)%rk(idx_177,i_6) * s(t)%y(idx_177,nspi(i_6,2))
                        enddo
                      end if
                      if (nspho2 == nspi(i_6,2)) then
                        do idx_178 = 1, n_points, 1
                          dssnum_2(idx_178,n_ho2_2) = dssnum_2(idx_178,n_ho2_2) + s(t)%rk(idx_178,i_6) * s(t)%y(idx_178,nspi(i_6,1))
                        enddo
                      end if
                      if (nspno == nspi(i_6,1)) then
                        do idx_179 = 1, n_points, 1
                          dssnum_2(idx_179,n_no_2) = dssnum_2(idx_179,n_no_2) + s(t)%rk(idx_179,i_6) * s(t)%y(idx_179,nspi(i_6,2))
                        enddo
                      end if
                      if (nspno == nspi(i_6,2)) then
                        do idx_180 = 1, n_points, 1
                          dssnum_2(idx_180,n_no_2) = dssnum_2(idx_180,n_no_2) + s(t)%rk(idx_180,i_6) * s(t)%y(idx_180,nspi(i_6,1))
                        enddo
                      end if
                    end if
                  end if
                enddo
                do jr_2 = 1, nssrt(ix_3), 1
                  i_6 = nssri(ix_3,jr_2)
                  j_7 = nssrx(ix_3,jr_2)
                  if (i_6 <= nuni) then
                    do idx_181 = 1, n_points, 1
                      ssden_2(idx_181) = ssden_2(idx_181) + s(t)%rk(idx_181,i_6)
                    enddo
                  else
                    do idx_182 = 1, n_points, 1
                      ssden_2(idx_182) = ssden_2(idx_182) + s(t)%rk(idx_182,i_6) * s(t)%y(idx_182,nspi(i_6,j_7))
                    enddo
                    if (ix_3 < 5) then
                      if (nspo3 == nspi(i_6,j_7)) then
                        do idx_183 = 1, n_points, 1
                          dssden_2(idx_183,n_o3_2) = dssden_2(idx_183,n_o3_2) + s(t)%rk(idx_183,i_6)
                        enddo
                      end if
                      if (nspoh == nspi(i_6,j_7)) then
                        do idx_184 = 1, n_points, 1
                          dssden_2(idx_184,n_oh_2) = dssden_2(idx_184,n_oh_2) + s(t)%rk(idx_184,i_6)
                        enddo
                      end if
                      if (nspho2 == nspi(i_6,j_7)) then
                        do idx_185 = 1, n_points, 1
                          dssden_2(idx_185,n_ho2_2) = dssden_2(idx_185,n_ho2_2) + s(t)%rk(idx_185,i_6)
                        enddo
                      end if
                      if (nspno == nspi(i_6,j_7)) then
                        do idx_186 = 1, n_points, 1
                          dssden_2(idx_186,n_no_2) = dssden_2(idx_186,n_no_2) + s(t)%rk(idx_186,i_6)
                        enddo
                      end if
                    end if
                  end if
                enddo
                do idx_187 = 1, n_points, 1
                  s(t)%y(idx_187,nssi(ix_3)) = ssnum_2(idx_187) / ssden_2(idx_187)
                enddo
                if (ix_3 < 5) then
                  do jr_2 = 1, 4, 1
                    do idx_188 = 1, n_points, 1
                      s(t)%deriv(idx_188,ix_3,jr_2) = (ssden_2(idx_188) * dssnum_2(idx_188,jr_2) - ssnum_2(idx_188) * dssden_2(idx_188,jr_2)) / (ssden_2(idx_188) * ssden_2(idx_188))
                    enddo
                  enddo
                end if
              enddo
              do x_2 = 1, n_points, 1
                if (peps < s(t)%y(x_2,nspo1d)) then
                  s(t)%deriv(x_2,nss_o1d,n_o3_2) = s(t)%deriv(x_2,nss_o1d,n_o3_2) * s(t)%y(x_2,nspo3) / s(t)%y(x_2,nspo1d)
                  s(t)%deriv(x_2,nss_o1d,n_oh_2) = s(t)%deriv(x_2,nss_o1d,n_oh_2) * s(t)%y(x_2,nspoh) / s(t)%y(x_2,nspo1d)
                  s(t)%deriv(x_2,nss_o1d,n_ho2_2) = s(t)%deriv(x_2,nss_o1d,n_ho2_2) * s(t)%y(x_2,nspho2) / s(t)%y(x_2,nspo1d)
                  s(t)%deriv(x_2,nss_o1d,n_no_2) = s(t)%deriv(x_2,nss_o1d,n_no_2) * s(t)%y(x_2,nspno) / s(t)%y(x_2,nspo1d)
                else
                  s(t)%deriv(x_2,nss_o1d,n_o3_2) = 1.0
                  s(t)%deriv(x_2,nss_o1d,n_oh_2) = 1.0
                  s(t)%deriv(x_2,nss_o1d,n_ho2_2) = 1.0
                  s(t)%deriv(x_2,nss_o1d,n_no_2) = 1.0
                end if
              enddo
              if (o3p_in_ss) then
                do x_2 = 1, n_points, 1
                  if (peps < s(t)%y(x_2,nspo3p)) then
                    s(t)%deriv(x_2,nss_o3p,n_o3_2) = s(t)%deriv(x_2,nss_o3p,n_o3_2) * s(t)%y(x_2,nspo3) / s(t)%y(x_2,nspo3p)
                    s(t)%deriv(x_2,nss_o3p,n_oh_2) = s(t)%deriv(x_2,nss_o3p,n_oh_2) * s(t)%y(x_2,nspoh) / s(t)%y(x_2,nspo3p)
                    s(t)%deriv(x_2,nss_o3p,n_ho2_2) = s(t)%deriv(x_2,nss_o3p,n_ho2_2) * s(t)%y(x_2,nspho2) / s(t)%y(x_2,nspo3p)
                    s(t)%deriv(x_2,nss_o3p,n_no_2) = s(t)%deriv(x_2,nss_o3p,n_no_2) * s(t)%y(x_2,nspno) / s(t)%y(x_2,nspo3p)
                  else
                    s(t)%deriv(x_2,nss_o3p,n_o3_2) = 1.0
                    s(t)%deriv(x_2,nss_o3p,n_oh_2) = 1.0
                    s(t)%deriv(x_2,nss_o3p,n_ho2_2) = 1.0
                    s(t)%deriv(x_2,nss_o3p,n_no_2) = 1.0
                  end if
                enddo
              end if
              if (n_in_ss) then
                do x_2 = 1, n_points, 1
                  if (peps < s(t)%y(x_2,nspn)) then
                    s(t)%deriv(x_2,nss_n,n_o3_2) = s(t)%deriv(x_2,nss_n,n_o3_2) * s(t)%y(x_2,nspo3) / s(t)%y(x_2,nspn)
                    s(t)%deriv(x_2,nss_n,n_oh_2) = s(t)%deriv(x_2,nss_n,n_oh_2) * s(t)%y(x_2,nspoh) / s(t)%y(x_2,nspn)
                    s(t)%deriv(x_2,nss_n,n_ho2_2) = s(t)%deriv(x_2,nss_n,n_ho2_2) * s(t)%y(x_2,nspho2) / s(t)%y(x_2,nspn)
                    s(t)%deriv(x_2,nss_n,n_no_2) = s(t)%deriv(x_2,nss_n,n_no_2) * s(t)%y(x_2,nspno) / s(t)%y(x_2,nspn)
                  else
                    s(t)%deriv(x_2,nss_n,n_o3_2) = 1.0
                    s(t)%deriv(x_2,nss_n,n_oh_2) = 1.0
                    s(t)%deriv(x_2,nss_n,n_ho2_2) = 1.0
                    s(t)%deriv(x_2,nss_n,n_no_2) = 1.0
                  end if
                enddo
              end if
              if (h_in_ss) then
                do x_2 = 1, n_points, 1
                  if (peps < s(t)%y(x_2,nsph)) then
                    s(t)%deriv(x_2,nss_h,n_o3_2) = s(t)%deriv(x_2,nss_h,n_o3_2) * s(t)%y(x_2,nspo3) / s(t)%y(x_2,nsph)
                    s(t)%deriv(x_2,nss_h,n_oh_2) = s(t)%deriv(x_2,nss_h,n_oh_2) * s(t)%y(x_2,nspoh) / s(t)%y(x_2,nsph)
                    s(t)%deriv(x_2,nss_h,n_ho2_2) = s(t)%deriv(x_2,nss_h,n_ho2_2) * s(t)%y(x_2,nspho2) / s(t)%y(x_2,nsph)
                    s(t)%deriv(x_2,nss_h,n_no_2) = s(t)%deriv(x_2,nss_h,n_no_2) * s(t)%y(x_2,nspno) / s(t)%y(x_2,nsph)
                  else
                    s(t)%deriv(x_2,nss_h,n_o3_2) = 1.0
                    s(t)%deriv(x_2,nss_h,n_oh_2) = 1.0
                    s(t)%deriv(x_2,nss_h,n_ho2_2) = 1.0
                    s(t)%deriv(x_2,nss_h,n_no_2) = 1.0
                  end if
                enddo
              end if
              do jtr_5 = 1, jpcspf, 1
                do jl_10 = 1, n_points, 1
                  s(t)%fdot(jl_10,jtr_5) = 0.0
                enddo
              enddo
              do js_15 = 1, jpspec, 1
                do jl_20 = 1, n_points, 1
                  s(t)%pd(jl_20,js_15) = 0.0
                  s(t)%pd(jl_20,js_15 + jpspec) = 0.0
                enddo
              enddo
              do jr_8 = 1, nuni, 1
                ir1_3 = nspi(jr_8,1)
                do jl_20 = 1, n_points, 1
                  s(t)%prk(jl_20,jr_8) = s(t)%rk(jl_20,jr_8) * s(t)%y(jl_20,ir1_3)
                enddo
              enddo
              iunip1_3 = nuni + 1
              do jr_8 = iunip1_3, jpnr, 1
                ir1_3 = nspi(jr_8,1)
                ir2_3 = nspi(jr_8,2)
                do jl_20 = 1, n_points, 1
                  s(t)%prk(jl_20,jr_8) = s(t)%rk(jl_20,jr_8) * s(t)%y(jl_20,ir1_3) * s(t)%y(jl_20,ir2_3)
                enddo
              enddo
              do jskip_3 = 0, jpspec, jpspec
                do j_20 = 1, nf, 1
                  js_15 = nlf(j_20) + jskip_3
                  ip3_3 = ngrp(js_15,3)
                  ip2_3 = ngrp(js_15,2)
                  ip1_3 = ngrp(js_15,1)
                  do j3_5 = 1, ip3_3, 1
                    i1_3 = nprdx3(1,j3_5,js_15)
                    i2_3 = nprdx3(2,j3_5,js_15)
                    i3_3 = nprdx3(3,j3_5,js_15)
                    do jl_20 = 1, n_points, 1
                      s(t)%pd(jl_20,js_15) = s(t)%pd(jl_20,js_15) + s(t)%prk(jl_20,i1_3) + s(t)%prk(jl_20,i2_3) + s(t)%prk(jl_20,i3_3)
                    enddo
                  enddo
                  do j2_3 = 1, ip2_3, 1
                    i1_3 = nprdx2(1,js_15)
                    i2_3 = nprdx2(2,js_15)
                    do jl_20 = 1, n_points, 1
                      s(t)%pd(jl_20,js_15) = s(t)%pd(jl_20,js_15) + s(t)%prk(jl_20,i1_3) + s(t)%prk(jl_20,i2_3)
                    enddo
                  enddo
                  do j1_3 = 1, ip1_3, 1
                    i1_3 = nprdx1(js_15)
                    do jl_20 = 1, n_points, 1
                      s(t)%pd(jl_20,js_15) = s(t)%pd(jl_20,js_15) + s(t)%prk(jl_20,i1_3)
                    enddo
                  enddo
                enddo
              enddo
              do jr_8 = 1, nnfrp, 1
                iss_4 = ntabfp(jr_8,1)
                ir_3 = ntabfp(jr_8,2)
                ix_7 = ntabfp(jr_8,3)
                fr_3 = frpx(ix_7)
                do jl_20 = 1, n_points, 1
                  s(t)%pd(jl_20,iss_4) = s(t)%pd(jl_20,iss_4) + fr_3 * s(t)%prk(jl_20,ir_3)
                enddo
              enddo
              if (.true.) then
                istart_9 = nldepx(1)
                iend_9 = nldepx(2)
                do j_20 = istart_9, iend_9, 1
                  js_15 = nldepx(j_20)
                  do jl_20 = 1, n_points, 1
                    s(t)%pd(jl_20,js_15 + jpspec) = s(t)%pd(jl_20,js_15 + jpspec) + (s(t)%dpd(jl_20,js_15) + s(t)%dpw(jl_20,js_15)) * s(t)%y(jl_20,js_15)
                  enddo
                enddo
                istart_9 = nldepx(3)
                iend_9 = nldepx(4)
                do j_20 = istart_9, iend_9, 1
                  js_15 = nldepx(j_20)
                  do jl_20 = 1, n_points, 1
                    s(t)%pd(jl_20,js_15 + jpspec) = s(t)%pd(jl_20,js_15 + jpspec) + s(t)%dpd(jl_20,js_15) * s(t)%y(jl_20,js_15)
                  enddo
                enddo
                istart_9 = nldepx(5)
                iend_9 = nldepx(6)
                do j_20 = istart_9, iend_9, 1
                  js_15 = nldepx(j_20)
                  do jl_20 = 1, n_points, 1
                    s(t)%pd(jl_20,js_15 + jpspec) = s(t)%pd(jl_20,js_15 + jpspec) + s(t)%dpw(jl_20,js_15) * s(t)%y(jl_20,js_15)
                  enddo
                enddo
              end if
              do j_8 = 1, nf, 1
                js_5 = nlf(j_8)
                ifam_3 = moffam(js_5)
                itr_4 = madvtr(js_5)
                gfam_1 = ifam_3 /= 0
                gtr_1 = itr_4 /= 0
                do jl_10 = 1, n_points, 1
                  s(t)%ydot(jl_10,js_5) = s(t)%pd(jl_10,js_5) - s(t)%pd(jl_10,js_5 + jpspec)
                  if (gfam_1 .AND. (.NOT.gtr_1 .OR. gtr_1 .AND. s(t)%linfam(jl_10,itr_4))) then
                    s(t)%fdot(jl_10,ifam_3) = s(t)%fdot(jl_10,ifam_3) + nodd(js_5) * s(t)%ydot(jl_10,js_5)
                  end if
                  if (gtr_1) then
                    s(t)%fdot(jl_10,itr_4) = s(t)%ydot(jl_10,js_5)
                  end if
                enddo
              enddo
              do idx_189 = 1, jpcspf, 1
                do idx_190 = 1, n_points, 1
                  g_f_old(idx_190,idx_189) = g_f(idx_190,idx_189)
                enddo
              enddo
              do idx_191 = 1, jpcspf, 1
                do idx_192 = 1, n_points, 1
                  g_f(idx_192,idx_191) = (s(t)%f(idx_192 ,idx_191 ) - f_initial_1(idx_192,idx_191)) * deltt - s(t)%fdot(idx_192 ,idx_191 )
                enddo
              enddo
              do idx_193 = 1, jpcspf, 1
                do idx_194 = 1, n_points, 1
                  delta_g(idx_194,idx_193) = g_f(idx_194,idx_193) - g_f_old(idx_194,idx_193)
                enddo
              enddo
              do jl = 1, n_points, 1
                coeff = DOT_PRODUCT(g_f(jl,:), delta_g(jl,:)) / DOT_PRODUCT(delta_g(jl,:), delta_g(jl,:))
                do idx_195 = 1, UBOUND(g_ftmp, dim=2), 1
                  g_ftmp(jl,idx_195) = g_f(jl,idx_195) * (1.0 - coeff)
                enddo
              enddo
              do i_7 = 1, jpcspf, 1
                do idx_196 = 1, UBOUND(bb1_1, dim=1), 1
                  bb1_1(idx_196,i_7) = g_ftmp(idx_196,reorder(i_7))
                enddo
              enddo
              do kr_1 = 1, jpcspf - 1, 1
                do i_7 = kr_1 + 1, jpcspf, 1
                  ikr_1 = s(t)%modified_map(i_7,kr_1)
                  if (ikr_1 > 0) then
                    do idx_197 = 1, UBOUND(bb1_1, dim=1), 1
                      bb1_1(idx_197,i_7) = bb1_1(idx_197,i_7) - s(t)%spfj(idx_197 ,ikr_1) * bb1_1(idx_197,kr_1)
                    enddo
                  end if
                enddo
              enddo
              do kr_1 = jpcspf, 1, -1
                do idx_198 = 1, UBOUND(xx1_1, dim=1), 1
                  xx1_1(idx_198,kr_1) = bb1_1(idx_198,kr_1)
                enddo
                do j_9 = kr_1 + 1, jpcspf, 1
                  krj_1 = s(t)%modified_map(kr_1,j_9)
                  if (krj_1 > 0) then
                    do idx_199 = 1, UBOUND(xx1_1, dim=1), 1
                      xx1_1(idx_199,kr_1) = xx1_1(idx_199,kr_1) - s(t)%spfj(idx_199 ,krj_1) * xx1_1(idx_199,j_9)
                    enddo
                  end if
                enddo
                do idx_200 = 1, UBOUND(pivot_1, dim=1), 1
                  pivot_1(idx_200) = s(t)%spfj(idx_200 ,s(t)%modified_map(kr_1,kr_1))
                enddo
                do i_7 = 1, n_points, 1
                  if (f_min > ABS(pivot_1(i_7))) then
                    pivot_1(i_7) = f_min
                  end if
                enddo
                if (max_val > 0) then
                  do idx_201 = 1, UBOUND(xx1_1, dim=1), 1
                    xx1_1(idx_201,kr_1) = MIN(MAX(xx1_1(idx_201,kr_1), -max_val), max_val) / pivot_1(idx_201)
                  enddo
                else
                  do idx_202 = 1, UBOUND(xx1_1, dim=1), 1
                    xx1_1(idx_202,kr_1) = xx1_1(idx_202,kr_1) / pivot_1(idx_202)
                  enddo
                end if
              enddo
              do idx_203 = 1, UBOUND(f_incr, dim=2), 1
                do idx_204 = 1, UBOUND(f_incr, dim=1), 1
                  f_incr(idx_204,idx_203) = 0.0
                enddo
              enddo
              do i_7 = 1, jpcspf, 1
                do idx_205 = 1, UBOUND(f_incr, dim=1), 1
                  f_incr(idx_205,reorder(i_7)) = xx1_1(idx_205,i_7)
                enddo
              enddo
              do idx_206 = 1, UBOUND(s(t)%f, dim=2), 1
                do idx_207 = 1, UBOUND(s(t)%f, dim=1), 1
                  s(t)%f(idx_207,idx_206) = s(t)%f(idx_207,idx_206) + f_incr(idx_207 ,idx_206 )
                enddo
              enddo
              do idx_208 = 1, UBOUND(s(t)%f, dim=2), 1
                do idx_209 = 1, UBOUND(s(t)%f, dim=1), 1
                  s(t)%f(idx_209,idx_208) = ABS(s(t)%f(idx_209,idx_208))
                enddo
              enddo
            end if
          end if
        enddo

        ! PSyclone CodeBlock (unsupported code) reason:
        !  - Unsupported labelled statement
        9999 CONTINUE
        if (exit_code /= 0) then
          do idx_210 = 1, UBOUND(s(t)%f, dim=2), 1
            do idx_211 = 1, UBOUND(s(t)%f, dim=1), 1
              s(t)%f(idx_211,idx_210) = f_initial_1(idx_211 ,idx_210 )
            enddo
          enddo
          solver_iter = iter_1
          if (iter_1 >= ltrig_iter) then
            exit_code = 4
          end if
        end if
        num_iter = num_iter + solver_iter
        if (exit_code == 0) then
          iter = iter + 1
        else
          do idx_212 = 1, UBOUND(s(t)%f, dim=2), 1
            do idx_213 = 1, n_points, 1
              s(t)%f(idx_213,idx_212) = f_initial(idx_213,idx_212 )
            enddo
          enddo
          if (exit_code == 4) then
            s(t)%ltrig = .true.
          else
            s(t)%ncsteps = s(t)%ncsteps * 2
            s(t)%cdt = s(t)%cdt / 2.0
            iredo = iredo * 2
            if (iredo >= max_redo) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Exit_Stmt
              EXIT
            end if
          end if
          jit = 0
          ifam_8 = 0
          imaj_4 = 0
          sl_2 = 0.0
          zthresh_2 = 2.0 / s(t)%cdt
          do idx_214 = 1, jpspec, 1
            do idx_215 = 1, n_points, 1
              zy_2(idx_215,idx_214) = 0.0
            enddo
          enddo
          do idx_216 = 1, n_points, 1
            zb_2(idx_216) = 0.0
          enddo
          do idx_217 = 1, n_points, 1
            zc_2(idx_217) = 0.0
          enddo
          do idx_218 = 1, n_points, 1
            zd_2(idx_218) = 0.0
          enddo
          if (not_first_call) then
            peps10_2 = 10.0 * peps
            istart_11 = nlmajmin(1)
            iend_11 = nlmajmin(2)
            do j_25 = istart_11, iend_11, 1
              js_19 = nlmajmin(j_25)
              ifam_10 = moffam(js_19)
              do jl_23 = 1, n_points, 1
                s(t)%y(jl_23,js_19) = s(t)%f(jl_23,ifam_10)
              enddo
            enddo
            istart_11 = nlmajmin(3)
            iend_11 = nlmajmin(4)
            do j_25 = istart_11, iend_11, 1
              js_19 = nlmajmin(j_25)
              do jl_23 = 1, n_points, 1
                s(t)%y(jl_23,js_19) = peps10_2
              enddo
            enddo
            do j_25 = 1, iss, 1
              js_19 = ilss(j_25)
              do jl_23 = 1, n_points, 1
                s(t)%y(jl_23,js_19) = peps10_2
              enddo
            enddo
            fco2_2 = fco2_default_2
            fch4_2 = fch4_default_2
            fh2_2 = fh2_default_2
            fo2_2 = fo2_default_2
            fn2_2 = fn2_default_2
            if (ukca_config%l_chem_environ_co2_scalar) then
              fco2_2 = atmospheric_co2 / c_co2
              if (ukca_config%l_chem_environ_co2_fld) then
                fco2_2 = 0.0
              end if
            end if
            if (ukca_config%l_chem_environ_h2_scalar) then
              fh2_2 = atmospheric_h2 / c_h2
            end if
            if (ukca_config%l_chem_environ_n2_scalar) then
              fn2_2 = atmospheric_n2 / c_n2
            end if
            if (ukca_config%l_chem_environ_o2_scalar) then
              fo2_2 = atmospheric_o2 / c_o2
            end if
            if (ukca_config%l_chem_environ_ch4_scalar) then
              fch4_2 = atmospheric_ch4 / c_ch4
            end if
            do j_25 = 1, ict, 1
              js_19 = ilct(j_25)
              if (speci(js_19) == 'N2        ') then
                zcnst_2 = fn2_2
              elseif (speci(js_19) == 'O2        ') then
                zcnst_2 = fo2_2
              elseif (speci(js_19) == 'CO2       ') then
                zcnst_2 = fco2_2
              elseif (speci(js_19) == 'H2        ') then
                zcnst_2 = fh2_2
              elseif (speci(js_19) == 'CH4       ') then
                zcnst_2 = fch4_2
              else
                errcode_8 = -1
                zcnst_2 = 0.0
              end if
              do jl_23 = 1, n_points, 1
                s(t)%y(jl_23,js_19) = zcnst_2 * s(t)%tnd(jl_23)
              enddo
            enddo
          end if
          do j_25 = 1, iftr, 1
            js_19 = ilftr(j_25)
            itr_11 = madvtr(js_19)
            do jl_23 = 1, n_points, 1
              s(t)%y(jl_23,js_19) = s(t)%f(jl_23,itr_11)
            enddo
          enddo
          if (not_first_call) then
            do j_25 = 1, icf, 1
              js_19 = ilcf(j_25)
              species_tmp_2 = speci(js_19)
              row_length_2 = ukca_config%row_length
              rows_2 = ukca_config%rows
              if (species_tmp_2(1:4) == 'CO2 ' .AND. ukca_config%l_chem_environ_co2_fld) then
                do idx_219 = 1, n_points, 1
                  s(t)%y(idx_219,js_19) = s(t)%co2(idx_219 )
                enddo
              elseif (species_tmp_2(1:4) == 'H2O ') then
                do idx_220 = 1, n_points, 1
                  s(t)%y(idx_220,js_19) = s(t)%wp(idx_220 )
                enddo
              elseif (species_tmp_2(1:4) == 'RO2 ') then
                if (ukca_config%l_ukca_ro2_perm) then
                  do idx_221 = 1, UBOUND(fro2_5, dim=1), 1
                    fro2_5(idx_221) = 0.0
                  enddo
                  do j_30 = 1, jpro2, 1
                    iro2_5 = nlfro2(j_30)
                    do idx_222 = 1, UBOUND(fro2_5, dim=1), 1
                      fro2_5(idx_222) = fro2_5(idx_222) + s(t)%f(idx_222 ,iro2_5)
                    enddo
                  enddo
                  do idx_223 = 1, n_points, 1
                    s(t)%y(idx_223,js_19) = fro2_5(idx_223) / s(t)%tnd(idx_223 )
                  enddo
                else
                  errcode_9 = 126
                end if
              else
                errcode_9 = 124
              end if
              do idx_224 = 1, n_points, 1
                s(t)%y(idx_224,js_19) = s(t)%y(idx_224,js_19) * s(t)%tnd(idx_224 )
              enddo
            enddo
          end if
          if (ukca_config%l_ukca_ro2_perm) then
            do idx_225 = 1, n_points, 1
              fro2_old_2(idx_225) = s(t)%y(idx_225,jsro2)
            enddo
            do idx_226 = 1, n_points, 1
              fro2_4(idx_226) = 0.0
            enddo
            do j_23 = 1, jpro2, 1
              iro2_4 = nlfro2(j_23)
              do idx_227 = 1, n_points, 1
                fro2_4(idx_227) = fro2_4(idx_227) + s(t)%f(idx_227,iro2_4)
              enddo
            enddo
            if (jit <= 0) then
              damp_ro2_2 = 1.0
            elseif (jit == 1) then
              damp_ro2_2 = 0.5
            elseif (jit == 2) then
              damp_ro2_2 = 0.9
            else
              damp_ro2_2 = 0.99
            end if
            do idx_228 = 1, n_points, 1
              s(t)%y(idx_228,jsro2) = fro2_old_2(idx_228) + damp_ro2_2 * (fro2_4(idx_228) - fro2_old_2(idx_228))
            enddo
          end if
          if (nstst == 0) then
          end if
          if (not_first_call .AND. nitfg < 5) then
            nitfg = 5
            icode_2 = -1
          end if
          if (not_first_call) then
            do j_26 = 1, nstst, 1
              js_20 = nlstst(j_26)
              do jl_24 = 1, n_points, 1
                s(t)%qa(jl_24,js_20) = 0.0
              enddo
            enddo
            do jr_9 = 1, jpnr, 1
              isp_2 = nspi(jr_9,1)
              if (isp_2 == nspi(jr_9,2)) then
                do jl_24 = 1, n_points, 1
                  s(t)%qa(jl_24,isp_2) = s(t)%qa(jl_24,isp_2) + 2.0 * s(t)%rk(jl_24,jr_9)
                enddo
              end if
            enddo
          end if
          if (nitfg == 0) then
            istart_12 = nlmajmin(1)
            iend_12 = nlmajmin(2)
            do j_27 = istart_12, iend_12, 1
              js_21 = nlmajmin(j_27)
              ifam_11 = moffam(js_21)
              do jl_25 = 1, n_points, 1
                s(t)%y(jl_25,js_21) = s(t)%f(jl_25,ifam_11) * s(t)%ratio(jl_25,js_21)
              enddo
            enddo
            istart_12 = nlmajmin(3)
            iend_12 = nlmajmin(4)
            do j_27 = istart_12, iend_12, 1
              js_21 = nlmajmin(j_27)
              ifam_11 = moffam(js_21)
              imaj_5 = majors(ifam_11)
              if (ctype(js_21) /= jpif) then
                do jl_25 = 1, n_points, 1
                  s(t)%y(jl_25,js_21) = s(t)%y(jl_25,imaj_5) * s(t)%ratio(jl_25,js_21)
                enddo
              else
                itr_12 = madvtr(js_21)
                do jl_25 = 1, n_points, 1
                  if (s(t)%linfam(jl_25,itr_12)) then
                    s(t)%y(jl_25,js_21) = s(t)%y(jl_25,imaj_5) * s(t)%ratio(jl_25,js_21)
                  end if
                enddo
              end if
            enddo
          end if
          do jit_3 = 1, nitfg, 1
            if (jit_3 > 1) then
              do j_23 = 1, nstst, 1
                js_17 = nlstst(j_23)
                do jl_21 = 1, n_points, 1
                  if (ABS(s(t)%y(jl_21,js_17)) < peps .AND. s(t)%pd(jl_21,js_17) > peps) then
                    s(t)%y(jl_21,js_17) = peps
                  end if
                enddo
              enddo
            end if
            do js_22 = 1, jpspec, 1
              do jl_26 = 1, n_points, 1
                s(t)%pd(jl_26,js_22) = 0.0
                s(t)%pd(jl_26,js_22 + jpspec) = 0.0
              enddo
            enddo
            do jr_10 = 1, nuni, 1
              ir1_4 = nspi(jr_10,1)
              do jl_26 = 1, n_points, 1
                s(t)%prk(jl_26,jr_10) = s(t)%rk(jl_26,jr_10) * s(t)%y(jl_26,ir1_4)
              enddo
            enddo
            iunip1_4 = nuni + 1
            do jr_10 = iunip1_4, jpnr, 1
              ir1_4 = nspi(jr_10,1)
              ir2_4 = nspi(jr_10,2)
              do jl_26 = 1, n_points, 1
                s(t)%prk(jl_26,jr_10) = s(t)%rk(jl_26,jr_10) * s(t)%y(jl_26,ir1_4) * s(t)%y(jl_26,ir2_4)
              enddo
            enddo
            do jskip_4 = 0, jpspec, jpspec
              do j_28 = 1, istmin, 1
                js_22 = ilstmin(j_28) + jskip_4
                ip3_4 = ngrp(js_22,3)
                ip2_4 = ngrp(js_22,2)
                ip1_4 = ngrp(js_22,1)
                do j3_6 = 1, ip3_4, 1
                  i1_4 = nprdx3(1,j3_6,js_22)
                  i2_4 = nprdx3(2,j3_6,js_22)
                  i3_4 = nprdx3(3,j3_6,js_22)
                  do jl_26 = 1, n_points, 1
                    s(t)%pd(jl_26,js_22) = s(t)%pd(jl_26,js_22) + s(t)%prk(jl_26,i1_4) + s(t)%prk(jl_26,i2_4) + s(t)%prk(jl_26,i3_4)
                  enddo
                enddo
                do j2_4 = 1, ip2_4, 1
                  i1_4 = nprdx2(1,js_22)
                  i2_4 = nprdx2(2,js_22)
                  do jl_26 = 1, n_points, 1
                    s(t)%pd(jl_26,js_22) = s(t)%pd(jl_26,js_22) + s(t)%prk(jl_26,i1_4) + s(t)%prk(jl_26,i2_4)
                  enddo
                enddo
                do j1_4 = 1, ip1_4, 1
                  i1_4 = nprdx1(js_22)
                  do jl_26 = 1, n_points, 1
                    s(t)%pd(jl_26,js_22) = s(t)%pd(jl_26,js_22) + s(t)%prk(jl_26,i1_4)
                  enddo
                enddo
              enddo
            enddo
            do jr_10 = 1, nnfrp, 1
              iss_5 = ntabfp(jr_10,1)
              ir_4 = ntabfp(jr_10,2)
              ix_8 = ntabfp(jr_10,3)
              fr_4 = frpx(ix_8)
              do jl_26 = 1, n_points, 1
                s(t)%pd(jl_26,iss_5) = s(t)%pd(jl_26,iss_5) + fr_4 * s(t)%prk(jl_26,ir_4)
              enddo
            enddo
            if (.false.) then
              istart_13 = nldepx(1)
              iend_13 = nldepx(2)
              do j_28 = istart_13, iend_13, 1
                js_22 = nldepx(j_28)
                do jl_26 = 1, n_points, 1
                  s(t)%pd(jl_26,js_22 + jpspec) = s(t)%pd(jl_26,js_22 + jpspec) + (s(t)%dpd(jl_26,js_22) + s(t)%dpw(jl_26,js_22)) * s(t)%y(jl_26,js_22)
                enddo
              enddo
              istart_13 = nldepx(3)
              iend_13 = nldepx(4)
              do j_28 = istart_13, iend_13, 1
                js_22 = nldepx(j_28)
                do jl_26 = 1, n_points, 1
                  s(t)%pd(jl_26,js_22 + jpspec) = s(t)%pd(jl_26,js_22 + jpspec) + s(t)%dpd(jl_26,js_22) * s(t)%y(jl_26,js_22)
                enddo
              enddo
              istart_13 = nldepx(5)
              iend_13 = nldepx(6)
              do j_28 = istart_13, iend_13, 1
                js_22 = nldepx(j_28)
                do jl_26 = 1, n_points, 1
                  s(t)%pd(jl_26,js_22 + jpspec) = s(t)%pd(jl_26,js_22 + jpspec) + s(t)%dpw(jl_26,js_22) * s(t)%y(jl_26,js_22)
                enddo
              enddo
            end if
            do j_23 = 1, nstst, 1
              js_17 = nlstst(j_23)
              do jl_21 = 1, n_points, 1
                s(t)%ratio(jl_21,js_17) = 0.0
              enddo
            enddo
            do j_23 = 1, istmin, 1
              js_17 = ilstmin(j_23)
              do jl_21 = 1, n_points, 1
                if (s(t)%y(jl_21,js_17) < peps) then
                  sl_2 = 0.0
                else
                  sl_2 = s(t)%pd(jl_21,js_17 + jpspec) / s(t)%y(jl_21,js_17)
                end if
                if (s(t)%qa(jl_21,js_17) > peps) then
                  zb_2(jl_21) = sl_2 - s(t)%qa(jl_21,js_17) * s(t)%y(jl_21,js_17)
                  zc_2(jl_21) = s(t)%pd(jl_21,js_17)
                  zd_2(jl_21) = zb_2(jl_21) * zb_2(jl_21) + 4.0 * s(t)%qa(jl_21,js_17) * zc_2(jl_21)
                  if (zd_2(jl_21) > 0.0) then
                    s(t)%y(jl_21,js_17) = (zb_2(jl_21) - SQRT(zd_2(jl_21))) / (-2.0 * s(t)%qa(jl_21,js_17))
                    if (s(t)%y(jl_21,js_17) < 0.0) then
                      s(t)%y(jl_21,js_17) = 10.0 * peps
                    end if
                  else
                    s(t)%y(jl_21,js_17) = zb_2(jl_21) / (-2.0 * s(t)%qa(jl_21,js_17))
                  end if
                elseif (s(t)%qa(jl_21,js_17) <= peps .AND. sl_2 > peps) then
                  s(t)%y(jl_21,js_17) = s(t)%pd(jl_21,js_17) / sl_2
                else
                  s(t)%y(jl_21,js_17) = 0.0
                end if
              enddo
            enddo
            istart_10 = nlmajmin(3)
            iend_10 = nlmajmin(4)
            do j_23 = nlmajmin(3), nlmajmin(4), 1
              js_17 = nlmajmin(j_23)
              iodd_2 = nodd(js_17)
              ifam_8 = moffam(js_17)
              imaj_4 = 0
              if (ifam_8 /= 0) then
                imaj_4 = majors(ifam_8)
              end if
              do jl_21 = 1, n_points, 1
                if (s(t)%y(jl_21,imaj_4) > peps) then
                  s(t)%ratio(jl_21,js_17) = s(t)%y(jl_21,js_17) / s(t)%y(jl_21,imaj_4)
                  s(t)%ratio(jl_21,imaj_4) = s(t)%ratio(jl_21,imaj_4) + iodd_2 * s(t)%ratio(jl_21,js_17)
                else
                  s(t)%ratio(jl_21,js_17) = 0.0
                end if
              enddo
            enddo
            do j_23 = 1, ift, 1
              js_17 = ilft(j_23)
              iodd_2 = nodd(js_17)
              ifam_8 = moffam(js_17)
              itr_9 = madvtr(js_17)
              imaj_4 = majors(ifam_8)
              if (not_first_call .AND. jit_3 == 1) then
                do jl_21 = 1, n_points, 1
                  if (s(t)%y(jl_21,js_17) > peps) then
                    sl_2 = s(t)%pd(jl_21,js_17 + jpspec) / s(t)%y(jl_21,js_17)
                  else
                    sl_2 = 0.0
                  end if
                  s(t)%linfam(jl_21,itr_9) = sl_2 > zthresh_2
                  if (s(t)%linfam(jl_21,itr_9)) then
                    s(t)%f(jl_21,ifam_8) = s(t)%f(jl_21,ifam_8) + iodd_2 * s(t)%f(jl_21,itr_9)
                  end if
                enddo
              end if
              do jl_21 = 1, n_points, 1
                if (s(t)%linfam(jl_21,itr_9)) then
                  if (s(t)%y(jl_21,imaj_4) > peps) then
                    s(t)%ratio(jl_21,js_17) = s(t)%y(jl_21,js_17) / s(t)%y(jl_21,imaj_4)
                    s(t)%ratio(jl_21,imaj_4) = s(t)%ratio(jl_21,imaj_4) + iodd_2 * s(t)%ratio(jl_21,js_17)
                  else
                    s(t)%ratio(jl_21,js_17) = 0.0
                  end if
                else
                  s(t)%y(jl_21,js_17) = s(t)%f(jl_21,itr_9)
                end if
              enddo
            enddo
            istart_10 = nlmajmin(1)
            iend_10 = nlmajmin(2)
            do j_23 = istart_10, iend_10, 1
              js_17 = nlmajmin(j_23)
              iodd_2 = nodd(js_17)
              ifam_8 = moffam(js_17)
              do jl_21 = 1, n_points, 1
                s(t)%ratio(jl_21,js_17) = 1.0 / (iodd_2 + s(t)%ratio(jl_21,js_17))
                s(t)%y(jl_21,js_17) = s(t)%f(jl_21,ifam_8) * s(t)%ratio(jl_21,js_17)
              enddo
            enddo
            istart_10 = nlmajmin(3)
            iend_10 = nlmajmin(4)
            do j_23 = istart_10, iend_10, 1
              js_17 = nlmajmin(j_23)
              ifam_8 = moffam(js_17)
              imaj_4 = majors(ifam_8)
              if (ctype(js_17) /= jpif) then
                do jl_21 = 1, n_points, 1
                  s(t)%y(jl_21,js_17) = s(t)%y(jl_21,imaj_4) * s(t)%ratio(jl_21,js_17)
                enddo
              else
                itr_9 = madvtr(js_17)
                do jl_21 = 1, n_points, 1
                  if (s(t)%linfam(jl_21,itr_9)) then
                    s(t)%y(jl_21,js_17) = s(t)%y(jl_21,imaj_4) * s(t)%ratio(jl_21,js_17)
                  end if
                enddo
              end if
            enddo
            gconv_2 = .true.
            do j_23 = 1, nstst, 1
              js_17 = nlstst(j_23)
              do jl_21 = 1, n_points, 1
                if (ABS(s(t)%y(jl_21,js_17) - zy_2(jl_21,js_17)) > ftol * s(t)%y(jl_21,js_17) .AND. s(t)%y(jl_21,js_17) > s(t)%pmintnd(jl_21)) then
                  gconv_2 = .false.
                end if
              enddo
              if (.NOT.gconv_2) then
                ! PSyclone CodeBlock (unsupported code) reason:
                !  - Unsupported statement: Exit_Stmt
                EXIT
              end if
            enddo
            if (gconv_2) then
              ! PSyclone CodeBlock (unsupported code) reason:
              !  - Unsupported statement: Exit_Stmt
              EXIT
            end if
            do j_23 = 1, nstst, 1
              js_17 = nlstst(j_23)
              do jl_21 = 1, n_points, 1
                zy_2(jl_21,js_17) = s(t)%y(jl_21,js_17)
              enddo
            enddo
          enddo
          do jtr_6 = 1, jpcspf, 1
            do jl_22 = 1, n_points, 1
              s(t)%fdot(jl_22,jtr_6) = 0.0
            enddo
          enddo
          do js_23 = 1, jpspec, 1
            do jl_27 = 1, n_points, 1
              s(t)%pd(jl_27,js_23) = 0.0
              s(t)%pd(jl_27,js_23 + jpspec) = 0.0
            enddo
          enddo
          do jr_11 = 1, nuni, 1
            ir1_5 = nspi(jr_11,1)
            do jl_27 = 1, n_points, 1
              s(t)%prk(jl_27,jr_11) = s(t)%rk(jl_27,jr_11) * s(t)%y(jl_27,ir1_5)
            enddo
          enddo
          iunip1_5 = nuni + 1
          do jr_11 = iunip1_5, jpnr, 1
            ir1_5 = nspi(jr_11,1)
            ir2_5 = nspi(jr_11,2)
            do jl_27 = 1, n_points, 1
              s(t)%prk(jl_27,jr_11) = s(t)%rk(jl_27,jr_11) * s(t)%y(jl_27,ir1_5) * s(t)%y(jl_27,ir2_5)
            enddo
          enddo
          do jskip_5 = 0, jpspec, jpspec
            do j_29 = 1, nf, 1
              js_23 = nlf(j_29) + jskip_5
              ip3_5 = ngrp(js_23,3)
              ip2_5 = ngrp(js_23,2)
              ip1_5 = ngrp(js_23,1)
              do j3_7 = 1, ip3_5, 1
                i1_5 = nprdx3(1,j3_7,js_23)
                i2_5 = nprdx3(2,j3_7,js_23)
                i3_5 = nprdx3(3,j3_7,js_23)
                do jl_27 = 1, n_points, 1
                  s(t)%pd(jl_27,js_23) = s(t)%pd(jl_27,js_23) + s(t)%prk(jl_27,i1_5) + s(t)%prk(jl_27,i2_5) + s(t)%prk(jl_27,i3_5)
                enddo
              enddo
              do j2_5 = 1, ip2_5, 1
                i1_5 = nprdx2(1,js_23)
                i2_5 = nprdx2(2,js_23)
                do jl_27 = 1, n_points, 1
                  s(t)%pd(jl_27,js_23) = s(t)%pd(jl_27,js_23) + s(t)%prk(jl_27,i1_5) + s(t)%prk(jl_27,i2_5)
                enddo
              enddo
              do j1_5 = 1, ip1_5, 1
                i1_5 = nprdx1(js_23)
                do jl_27 = 1, n_points, 1
                  s(t)%pd(jl_27,js_23) = s(t)%pd(jl_27,js_23) + s(t)%prk(jl_27,i1_5)
                enddo
              enddo
            enddo
          enddo
          do jr_11 = 1, nnfrp, 1
            iss_6 = ntabfp(jr_11,1)
            ir_5 = ntabfp(jr_11,2)
            ix_9 = ntabfp(jr_11,3)
            fr_5 = frpx(ix_9)
            do jl_27 = 1, n_points, 1
              s(t)%pd(jl_27,iss_6) = s(t)%pd(jl_27,iss_6) + fr_5 * s(t)%prk(jl_27,ir_5)
            enddo
          enddo
          if (.true.) then
            istart_14 = nldepx(1)
            iend_14 = nldepx(2)
            do j_29 = istart_14, iend_14, 1
              js_23 = nldepx(j_29)
              do jl_27 = 1, n_points, 1
                s(t)%pd(jl_27,js_23 + jpspec) = s(t)%pd(jl_27,js_23 + jpspec) + (s(t)%dpd(jl_27,js_23) + s(t)%dpw(jl_27,js_23)) * s(t)%y(jl_27,js_23)
              enddo
            enddo
            istart_14 = nldepx(3)
            iend_14 = nldepx(4)
            do j_29 = istart_14, iend_14, 1
              js_23 = nldepx(j_29)
              do jl_27 = 1, n_points, 1
                s(t)%pd(jl_27,js_23 + jpspec) = s(t)%pd(jl_27,js_23 + jpspec) + s(t)%dpd(jl_27,js_23) * s(t)%y(jl_27,js_23)
              enddo
            enddo
            istart_14 = nldepx(5)
            iend_14 = nldepx(6)
            do j_29 = istart_14, iend_14, 1
              js_23 = nldepx(j_29)
              do jl_27 = 1, n_points, 1
                s(t)%pd(jl_27,js_23 + jpspec) = s(t)%pd(jl_27,js_23 + jpspec) + s(t)%dpw(jl_27,js_23) * s(t)%y(jl_27,js_23)
              enddo
            enddo
          end if
          do j_24 = 1, nf, 1
            js_18 = nlf(j_24)
            ifam_9 = moffam(js_18)
            itr_10 = madvtr(js_18)
            gfam_2 = ifam_9 /= 0
            gtr_2 = itr_10 /= 0
            do jl_22 = 1, n_points, 1
              s(t)%ydot(jl_22,js_18) = s(t)%pd(jl_22,js_18) - s(t)%pd(jl_22,js_18 + jpspec)
              if (gfam_2 .AND. (.NOT.gtr_2 .OR. gtr_2 .AND. s(t)%linfam(jl_22,itr_10))) then
                s(t)%fdot(jl_22,ifam_9) = s(t)%fdot(jl_22,ifam_9) + nodd(js_18) * s(t)%ydot(jl_22,js_18)
              end if
              if (gtr_2) then
                s(t)%fdot(jl_22,itr_10) = s(t)%ydot(jl_22,js_18)
              end if
            enddo
          enddo
          iter = 1
        end if
      end do
      s(t)%ncsteps = ncsteps_initial
      s(t)%cdt = cdt_initial
    enddo
  !$ACC END PARALLEL LOOP

  end subroutine asad_spmjpdriv

end module asad_spmjpdriv_full_mod
