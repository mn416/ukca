! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module defining ASAD arrays, variables, and parameters
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 programming standards.
!
! ----------------------------------------------------------------------

MODULE asad_mod

USE ukca_config_specification_mod, ONLY: ukca_config, int_method_nr

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ukca_missing_data_mod, ONLY: imdi

IMPLICIT NONE
PUBLIC

TYPE :: asad_state_type
  ! chemistry timestep
  REAL :: cdt

  ! chem. timestep for writing diagnostics
  REAL :: cdt_diag

  REAL, ALLOCATABLE :: co3(:)
  REAL, ALLOCATABLE :: deriv(:,:,:)
  REAL, ALLOCATABLE :: dpd(:,:)
  REAL, ALLOCATABLE :: dpw(:,:)
  REAL, ALLOCATABLE :: ej(:,:)
  REAL, ALLOCATABLE :: emr(:,:)

  ! Concentrations of chemically active species (molecule/cc)
  REAL, ALLOCATABLE :: f(:,:)

  ! Tendency of species calculated from rate equations
  REAL, ALLOCATABLE :: fdot(:,:)

  ! Full jacobian
  REAL, ALLOCATABLE :: fj(:,:,:)

  ! 1.0 if type 1 psc's are present, else 0
  REAL, ALLOCATABLE :: fpsc1(:)

  ! 1.0 if type 2 psc's are present, else 0
  REAL, ALLOCATABLE :: fpsc2(:)

  ! lower order solution
  REAL, ALLOCATABLE :: ftilde(:,:)

  ! interval in timesteps between calls to chemistry
  INTEGER :: interval = imdi

  ! Pivot information for solving jacobian
  INTEGER, ALLOCATABLE :: ipa(:,:)

  INTEGER :: jsubs

  ! Latitude
  REAL, ALLOCATABLE :: lati(:)

  LOGICAL, ALLOCATABLE :: linfam(:,:)

  ! ltrig set to debug slow convergence systems shared between
  ! asad_spimpmjp and asad_spmjpdriv
  LOGICAL :: ltrig

  ! modified map (after decomposition)
  INTEGER, ALLOCATABLE :: modified_map(:,:)

  ! No of chemical steps
  INTEGER :: ncsteps

  ! Factor for ncsteps halvings
  INTEGER :: ncsteps_factor

  REAL, ALLOCATABLE :: p(:)

  ! pd[:] = [prod[:], slos[:]]
  REAL, ALLOCATABLE :: pd(:,:)

  REAL, ALLOCATABLE :: pmintnd(:)
  REAL, ALLOCATABLE :: prk(:,:)

  ! Production rate for each species (molecule /cc /s)
  REAL, POINTER     :: prod(:,:) => NULL()

  REAL, ALLOCATABLE :: qa(:,:)
  REAL, ALLOCATABLE :: ratio(:,:)
  REAL, ALLOCATABLE :: rk(:,:)

  ! No. density type 2 psc solid phase h2o
  REAL, ALLOCATABLE :: sh2o(:)

  ! No. density type 1 psc solid phase hno3
  REAL, ALLOCATABLE :: shno3(:)

  ! Loss rate for each species (molecule /cc /s)
  REAL, POINTER     :: slos(:,:) => NULL()

  ! Sparse full Jacobian
  REAL, ALLOCATABLE :: spfj(:,:)

  ! Amount of H2O in solid phase
  REAL, ALLOCATABLE :: sph2o(:)

  ! Amount of HNO3 in solid phase
  REAL, ALLOCATABLE :: sphno3(:)

  REAL, ALLOCATABLE :: t(:)
  REAL, ALLOCATABLE :: t300(:)

  ! total number density (molecule/cc)
  REAL, ALLOCATABLE :: tnd(:)

  ! water vapour field (vmr)
  REAL, ALLOCATABLE :: wp(:)

  ! CO2 field (vmr)
  REAL, ALLOCATABLE :: co2(:)

  ! Concentration array for all species
  REAL, ALLOCATABLE :: y(:,:)

  ! Rate of change of y
  REAL, ALLOCATABLE :: ydot(:,:)

  ! Aerosol surface area
  REAL, ALLOCATABLE :: za(:)
END TYPE

! Per-thread/team state
TYPE(asad_state_type), ALLOCATABLE, TARGET :: asad_state(:)

REAL :: init_cdt       ! Initial chemistry timestep
REAL :: init_cdt_diag  ! Initial chem. timestep for writing diagnostics
! Initial interval in timesteps between calls to chemistry
INTEGER :: init_interval = imdi 
INTEGER :: init_ncsteps  ! Initial number of chemical steps
INTEGER :: init_ncsteps_factor ! Initial factor for ncsteps halvings

INTEGER :: jpctr = 0                 ! No. of transported chemical tracers
INTEGER :: jpspec = 0                ! No. of chemical species
INTEGER :: jpro2 = 0                 ! No. of RO2-type species
INTEGER :: jpcspf = 0                ! No. of active chem. species in f array
INTEGER :: jpbk = 0                  ! No. of bimolecular reactions
INTEGER :: jptk = 0                  ! No. of termolecular reactions
INTEGER :: jppj = 0                  ! No. of photolytic reactions
INTEGER :: jphk = 0                  ! No. of heterogeneous reactions
INTEGER :: jpnr = 0                  ! jpbk + jptk + jppj + jphk
INTEGER :: jpdd = 0                  ! No. of dry deposited species
INTEGER :: jpdw = 0                  ! No. of wet deposited species

REAL, ALLOCATABLE :: depvel(:,:,:)
REAL, ALLOCATABLE :: k298(:)         ! K(298) for Henry law (M/atm)
REAL, ALLOCATABLE :: dhr(:)          ! deltaH/R (K^-1)
REAL, ALLOCATABLE :: kd298(:,:)      ! dissociation constant (M)
REAL, ALLOCATABLE :: ddhr(:,:)       ! deltaH/R
REAL, ALLOCATABLE :: ct_k298(:)      ! As above, but for constant species
REAL, ALLOCATABLE :: ct_dhr(:)       ! Allocated in ukca_chem_offline
REAL, ALLOCATABLE :: ct_kd298(:,:)   !
REAL, ALLOCATABLE :: ct_ddhr(:,:)    !
REAL, ALLOCATABLE :: ab(:,:)
REAL, ALLOCATABLE :: at(:,:)
REAL, ALLOCATABLE :: aj(:,:)
REAL, ALLOCATABLE :: ah(:,:)
REAL, ALLOCATABLE :: ztabpd(:,:)

INTEGER, ALLOCATABLE :: madvtr(:)    ! Family array - major advected tracers
INTEGER, ALLOCATABLE :: majors(:)    ! Family array
INTEGER, ALLOCATABLE :: moffam(:)    ! Family array
INTEGER, ALLOCATABLE :: nodd(:)      ! Number of odd atoms => chch_defs%nodd
INTEGER, ALLOCATABLE :: nltrf(:)
INTEGER, ALLOCATABLE :: nltr3(:)
INTEGER, ALLOCATABLE :: nlnaro2(:)   ! Indices of RO2 species in nadv array
INTEGER, ALLOCATABLE :: nlfro2(:)    ! Indices of RO2 species in f-array
INTEGER, ALLOCATABLE :: ipa2(:)
INTEGER, ALLOCATABLE :: nltrim(:,:)
INTEGER, ALLOCATABLE :: nlpdv(:,:)
INTEGER, ALLOCATABLE :: nfrpx(:)     ! index to fractional product array nfrpx:
! one entry for each reaction. If zero, there are no fractional products.
! If non-zero, contains the array element in frpx for the first coefficient
! for that reaction.
INTEGER, ALLOCATABLE :: ntabfp(:,:)  ! Table used for indexing the fractional
!    products:  ntabfp(i,1) contains the species no.,
!               ntabfp(i,2) contains the reaction no., and
!               ntabfp(i,3) contains the array location in the frpx array.
INTEGER, ALLOCATABLE :: ntabpd(:,:)
INTEGER, ALLOCATABLE :: npdfr(:,:)
INTEGER, ALLOCATABLE :: ngrp(:,:)     ! For each species, holds no. of 3 sum,
                                      ! 2 sum and single sum terms, therefore
                                      ! controls the loop limits (asad_prls)
INTEGER, ALLOCATABLE :: njcgrp(:,:)
INTEGER, ALLOCATABLE :: nprdx3(:,:,:)
INTEGER, ALLOCATABLE :: nprdx2(:,:)
INTEGER, ALLOCATABLE :: nprdx1(:)
INTEGER, ALLOCATABLE :: njacx3(:,:,:)
INTEGER, ALLOCATABLE :: njacx2(:,:)
INTEGER, ALLOCATABLE :: njacx1(:)
INTEGER, ALLOCATABLE :: nmpjac(:)
INTEGER, ALLOCATABLE :: npjac1(:,:)
INTEGER, ALLOCATABLE :: nbrkx(:)
INTEGER, ALLOCATABLE :: ntrkx(:)
INTEGER, ALLOCATABLE :: nprkx(:)
INTEGER, ALLOCATABLE :: nhrkx(:)
INTEGER, ALLOCATABLE :: nlall(:)
INTEGER, ALLOCATABLE :: nlstst(:)   ! Index array of all steady state/family
                                    ! member species
INTEGER, ALLOCATABLE :: nlf(:)      ! Index array of all species in f array
INTEGER, ALLOCATABLE :: nlmajmin(:)
INTEGER, ALLOCATABLE :: nldepd(:)
INTEGER, ALLOCATABLE :: nldepw(:)
INTEGER, ALLOCATABLE :: nldepx(:)
INTEGER, ALLOCATABLE :: njcoth(:,:)
INTEGER, ALLOCATABLE :: nmzjac(:)
INTEGER, ALLOCATABLE :: nzjac1(:,:)
INTEGER, ALLOCATABLE :: njcoss(:,:)
INTEGER, ALLOCATABLE :: nmsjac(:)
INTEGER, ALLOCATABLE :: nsjac1(:,:)
INTEGER, ALLOCATABLE :: nsspt(:)
INTEGER, ALLOCATABLE :: nspi(:,:)
INTEGER, ALLOCATABLE :: nsspi(:,:)
INTEGER, ALLOCATABLE :: nssi(:)
INTEGER, ALLOCATABLE :: nssrt(:)
INTEGER, ALLOCATABLE :: nssri(:,:)
INTEGER, ALLOCATABLE :: nssrx(:,:)
REAL, ALLOCATABLE :: frpb(:)         ! fractional product array (bimol)
REAL, ALLOCATABLE :: frpt(:)         ! fractional product array (trimol)
REAL, ALLOCATABLE :: frpj(:)         ! fractional product array (phot)
REAL, ALLOCATABLE :: frph(:)         ! fractional product array (het)
REAL, ALLOCATABLE :: frpx(:)         ! fractional product array (total)
! sparse algebra
INTEGER, ALLOCATABLE :: nonzero_map_unordered(:,:)  ! Map of nonzero entries
INTEGER, ALLOCATABLE :: nonzero_map(:,:)
                                    ! Map of nonzero entries, before reordering

INTEGER, ALLOCATABLE :: reorder(:)   ! reordering of tracers to minimize fill-in

INTEGER, ALLOCATABLE :: ilcf(:)
INTEGER, ALLOCATABLE :: ilss(:)
INTEGER, ALLOCATABLE :: ilct(:)
INTEGER, ALLOCATABLE :: ilftr(:)
INTEGER, ALLOCATABLE :: ilft(:)
INTEGER, ALLOCATABLE :: ilstmin(:)

LOGICAL, ALLOCATABLE :: ldepd(:)     ! T for dry deposition
LOGICAL, ALLOCATABLE :: ldepw(:)     ! T for wet deposition

CHARACTER(LEN=10), ALLOCATABLE :: advt(:)      ! advected tracers names
CHARACTER(LEN=10), ALLOCATABLE :: nadvt(:)     ! non-advected species names
CHARACTER(LEN=10), ALLOCATABLE :: family(:)    ! family names
CHARACTER(LEN=10), ALLOCATABLE :: speci(:)     ! species names
CHARACTER(LEN=10), ALLOCATABLE :: spro2(:)     ! RO2-type species names
CHARACTER(LEN=10), ALLOCATABLE :: specf(:)     ! Names of species in f array

CHARACTER(LEN=2),  ALLOCATABLE :: ctype(:)     ! species type
CHARACTER(LEN=10), ALLOCATABLE :: spb(:,:)     ! species from bimolecular rates
CHARACTER(LEN=10), ALLOCATABLE :: spt(:,:)     ! species from termolecular rates
CHARACTER(LEN=10), ALLOCATABLE :: spj(:,:)     ! species from photolysis rates
CHARACTER(LEN=10), ALLOCATABLE :: sph(:,:)     ! species from heterogeneous
                                               ! rates

REAL, PARAMETER    :: pmin = 1.0e-20
REAL, PARAMETER    :: ptol = 1.0e-5            ! tolerance for time integration
REAL, PARAMETER    :: ftol = 1.0e-3            ! tolerance in family member
                                               ! iteration

INTEGER, PARAMETER :: kfphot=0
INTEGER, PARAMETER :: jpss = 16
INTEGER, PARAMETER :: jpssr = 51
INTEGER, PARAMETER :: nvar = 17
INTEGER, PARAMETER :: nllv = 17
INTEGER, PARAMETER :: ninv = 23
INTEGER, PARAMETER :: nout=71
!     nout:        Fortran channel for output in subroutine OUTVMR
INTEGER, PARAMETER :: jpem = 9      ! IS THIS A CONSTANT/USED?
INTEGER, PARAMETER :: jpeq = 2      ! dimension for dissociation arrays
INTEGER, PARAMETER :: jddept = 6
!     jddept:      Number of time periods used in dry deposition i.e.
!                  summer(day,night,24h ave), winter(day,night,24h ave)
INTEGER, PARAMETER :: jddepc = 5
!     jddepc:      Number of land use categories used in dry dep.
INTEGER, PARAMETER :: jpdwio = 56
!     jpdwio       Fortran i/o unit to read/write anything to do with
!                  wet/dry deposition
INTEGER, PARAMETER :: jpemio = 57
!     jpemio       Fortran i/o unit to read in emissions
INTEGER, PARAMETER :: jpfrpd=100
INTEGER, PARAMETER :: jpab    = 3
INTEGER, PARAMETER :: jpat    = 7
INTEGER, PARAMETER :: jpaj    = 3
INTEGER, PARAMETER :: jpah    = 3
INTEGER, PARAMETER :: jpspb   = 6
INTEGER, PARAMETER :: jpspt   = 4
INTEGER, PARAMETER :: jpspj   = 6
INTEGER, PARAMETER :: jpsph   = 6
INTEGER, PARAMETER :: jpmsp   = jpspb
INTEGER, PARAMETER :: jppjac  = 10
INTEGER, PARAMETER :: jpkargs = 10
INTEGER, PARAMETER :: jprargs = 10
INTEGER, PARAMETER :: jpcargs = 1

! Number of non-zero elements in sparse Jacobian
INTEGER :: total

! Fractional product parameters -
! these are initialsed to jp values in asad_mod_init
INTEGER :: jpfrpb
INTEGER :: jpfrpt
INTEGER :: jpfrpj
INTEGER :: jpfrph
INTEGER :: jpfrpx

INTEGER, PARAMETER :: jpcio   = 55

! Values for size of arrays set in asad_mod_init in order to
! allow larger size for CRI-Strat mechanism whilst preserving
! speed for other mechanisms
INTEGER, SAVE :: spfjsize_max  ! maximum number of
                               ! nonzero matrix elements
INTEGER, SAVE :: maxterms      ! maximum number of nonzero
                               ! terms for individual species
INTEGER, SAVE :: maxfterms     ! maximum number of terms
                               ! involving fractional  products

!---------------------------------------------------------------------
! Production and loss Arrays ALLOCATABLE to be able to use flexible
! Jacobian size. Must be saveable as only defined once in setup_spfuljac,
! which is called on the first call to spfuljac in asad_sparse_vars.
! These 8 arrays are used to calculate the Jacobian, which is fixed
! throughout the run and is only calculated on the first call.
! These arrays should NOT be deallocated.
INTEGER, ALLOCATABLE, SAVE :: nposterms(:)
INTEGER, ALLOCATABLE, SAVE :: nnegterms(:)
INTEGER, ALLOCATABLE, SAVE :: nfracterms(:)

INTEGER, ALLOCATABLE, SAVE :: posterms(:,:)
INTEGER, ALLOCATABLE, SAVE :: negterms(:,:)
INTEGER, ALLOCATABLE, SAVE :: fracterms(:,:)
INTEGER, ALLOCATABLE, SAVE :: base_tracer(:)

REAL, ALLOCATABLE, SAVE :: ffrac(:,:)
!---------------------------------------------------------------------

! Character codes for different chemical species types

CHARACTER(LEN=2),  PARAMETER :: jpfm = 'FM'  ! Family member
CHARACTER(LEN=2),  PARAMETER :: jpif = 'FT'  ! Family member dpnding in timestep
CHARACTER(LEN=2),  PARAMETER :: jpsp = 'TR'  ! Independent tracer
CHARACTER(LEN=2),  PARAMETER :: jpna = 'SS'  ! Steady-state species
CHARACTER(LEN=2),  PARAMETER :: jpco = 'CT'  ! Constant
CHARACTER(LEN=2),  PARAMETER :: jpcf = 'CF'  ! Constant with spatial field
CHARACTER(LEN=2),  PARAMETER :: jpoo = 'OO'  ! RO2-type species

LOGICAL, PARAMETER :: lvmr=.TRUE.    ! T for volume mixing ratio
LOGICAL      :: o1d_in_ss      ! T for steady state,
LOGICAL      :: o3p_in_ss      ! these are set in routine:
LOGICAL      :: n_in_ss        ! asad_inrats_set_sp_lists
LOGICAL      :: h_in_ss        !
INTEGER, PARAMETER :: nss_o1d=1      ! indicies of deriv array
INTEGER, PARAMETER :: nss_o3p=2      !    "         "
INTEGER, PARAMETER :: nss_n=3        !    "         "
INTEGER, PARAMETER :: nss_h=4        !    "         "

REAL    :: peps                      !
REAL, PARAMETER :: tslimit = 1200.0  ! timestep limit for some solvers

INTEGER, PARAMETER :: nrsteps_max = 200  ! max iterations for N-R solver

INTEGER :: nitnr          ! Iterations in ftoy for IMPACT solver
INTEGER :: nitfg          ! Max no of iterations in ftoy
INTEGER :: ntrf           ! Counter for species in 'f' array
INTEGER :: ntr3
INTEGER :: nnaf           ! Counter for non-advected tracers
INTEGER :: nro2           ! Counter for Non-transported RO2 species
INTEGER :: nuni
INTEGER :: nsst           ! No of steady-state species
INTEGER :: nit0=20        ! ftoy iterations with method=0
INTEGER :: nfphot
INTEGER :: method          ! chemistry integration method
INTEGER :: nnfrp           ! Total number of fractional products
INTEGER :: nstst           ! No of steady state species
INTEGER :: nf
INTEGER :: ndepd           ! No of dry deposited species
INTEGER :: ndepw           ! No of wet deposited species
INTEGER :: ntro3, ntroh, ntrho2, ntrno
INTEGER :: nspo1d, nspo3p, nspo3, nspoh
INTEGER :: nspho2, nspno, nspn, nsph
INTEGER :: ih_o3, ih_h2o2, ih_so2, ih_hno3 ! index for soluble species
INTEGER :: ih_o3_const                     ! index for soluble species
                                           !  as constant species
INTEGER :: ihso3_h2o2                      ! Index for HSO3- + H2O2(aq) reaction
INTEGER :: ihso3_o3                        ! Index for HSO3- + O3(aq) reaction
INTEGER :: iso3_o3                         ! Index for SO3-- + O3(aq) reaction
INTEGER :: ih2so4_hv                       ! Index for H2SO4 + hv reaction
INTEGER :: iso2_oh                         ! Index for SO2 + OH reaction
INTEGER :: ih2o2_oh                        ! Index for H2O2 + OH reaction
INTEGER :: ihno3_oh                        ! Index for HNO3 + OH reaction
INTEGER :: in2o5_h                         ! Index for N2O5 => HONO2 heterog.
                                           ! reaction
INTEGER :: iho2_h                          ! Index for HO2 + HO2 => H2O2
                                           ! heterogeneous reaction

INTEGER :: jsro2                           ! Index for summed RO2 concentration

LOGICAL :: ljacx

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ASAD_MOD'

CONTAINS

! ######################################################################
SUBROUTINE asad_mod_pre_setup_init()

! To allocate ASAD arrays used during the initial setup processing in
! asad_inrats_set_sp_lists

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_MOD_PRE_SETUP_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! All arrays within ASAD which previously were of size jpctr have been
! resized to jpcspf, unless they specifically refer to number of transported
! tracers rather than number of species to solve chemistry for, in order to
! include non-transported RO2 species.

IF (.NOT. ALLOCATED(nodd)) ALLOCATE(nodd(jpspec))
IF (.NOT. ALLOCATED(nltrf)) ALLOCATE(nltrf(jpcspf))
IF (.NOT. ALLOCATED(nltr3)) ALLOCATE(nltr3(jpcspf))
IF (.NOT. ALLOCATED(nlnaro2)) ALLOCATE(nlnaro2(jpro2))
IF (.NOT. ALLOCATED(nlfro2)) ALLOCATE(nlfro2(jpro2))
IF (.NOT. ALLOCATED(advt)) ALLOCATE(advt(jpctr))
IF (.NOT. ALLOCATED(nadvt)) ALLOCATE(nadvt(jpspec-jpctr))
IF (.NOT. ALLOCATED(spro2)) ALLOCATE(spro2(jpro2))
IF (.NOT. ALLOCATED(specf)) ALLOCATE(specf(jpcspf))
IF (.NOT. ALLOCATED(family)) ALLOCATE(family(jpspec))
IF (.NOT. ALLOCATED(speci)) ALLOCATE(speci(jpspec))
IF (.NOT. ALLOCATED(ctype)) ALLOCATE(ctype(jpspec))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mod_pre_setup_init


! ######################################################################
SUBROUTINE asad_mod_init(n_points)

! To allocate and initialise ASAD arrays and variables

!$ USE omp_lib

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points

INTEGER :: n_parallel, i
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_MOD_INIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Fractional product parameters - set total number allowed to be
! total number of reactions * 2 as each fractional product has 4 potentials!!
jpfrpb  = (jpspb-2)*jpbk
jpfrpt  = (jpspt-2)*jptk
jpfrpj  = (jpspj-2)*jppj
jpfrph  = (jpsph-2)*jphk
jpfrpx  = jpfrpb + jpfrpt + jpfrpj + jpfrph  ! Total FPs

! All arrays within ASAD which previously were of size jpctr have been
! resized to jpcspf, unless they specifically refer to number of transported
! tracers rather than number of species to solve chemistry for, in order to
! include non-transported RO2 species.
IF (.NOT. ALLOCATED(madvtr)) ALLOCATE(madvtr(jpspec))
IF (.NOT. ALLOCATED(majors)) ALLOCATE(majors(jpcspf))
IF (.NOT. ALLOCATED(moffam)) ALLOCATE(moffam(jpspec))
IF (.NOT. ALLOCATED(ipa2)) ALLOCATE(ipa2(jpctr))
IF (.NOT. ALLOCATED(nspi)) ALLOCATE(nspi(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nsspt)) ALLOCATE(nsspt(jpss))
IF (.NOT. ALLOCATED(nsspi)) ALLOCATE(nsspi(jpss,jpssr))
IF (.NOT. ALLOCATED(nssi)) ALLOCATE(nssi(jpss))
IF (.NOT. ALLOCATED(nssrt)) ALLOCATE(nssrt(jpss))
IF (.NOT. ALLOCATED(nssri)) ALLOCATE(nssri(jpss,jpssr))
IF (.NOT. ALLOCATED(nssrx)) ALLOCATE(nssrx(jpss,jpssr))
IF (.NOT. ALLOCATED(depvel)) ALLOCATE(depvel(jddept,jddepc,jpdd))
IF (.NOT. ALLOCATED(k298)) ALLOCATE(k298(jpdw))
IF (.NOT. ALLOCATED(dhr)) ALLOCATE(dhr(jpdw))
IF (.NOT. ALLOCATED(kd298)) ALLOCATE(kd298(jpdw,jpeq))
IF (.NOT. ALLOCATED(ddhr)) ALLOCATE(ddhr(jpdw,jpeq))
IF (.NOT. ALLOCATED(ldepd)) ALLOCATE(ldepd(jpspec))
IF (.NOT. ALLOCATED(ldepw)) ALLOCATE(ldepw(jpspec))
IF (.NOT. ALLOCATED(nltrim)) ALLOCATE(nltrim(0:jpcspf,3))
IF (.NOT. ALLOCATED(nlpdv)) ALLOCATE(nlpdv((jpspj-2)*jppj,2))
IF (.NOT. ALLOCATED(ab)) ALLOCATE(ab(jpbk,jpab))
IF (.NOT. ALLOCATED(at)) ALLOCATE(at(jptk,jpat))
IF (.NOT. ALLOCATED(aj)) ALLOCATE(aj(jppj,jpaj))
IF (.NOT. ALLOCATED(ah)) ALLOCATE(ah(jphk,jpah))
IF (.NOT. ALLOCATED(spb)) ALLOCATE(spb(jpbk+1,jpspb))
IF (.NOT. ALLOCATED(spt)) ALLOCATE(spt(jptk+1,jpspt))
IF (.NOT. ALLOCATED(spj)) ALLOCATE(spj(jppj+1,jpspj))
IF (.NOT. ALLOCATED(sph)) ALLOCATE(sph(jphk+1,jpsph))
IF (.NOT. ALLOCATED(frpb)) ALLOCATE(frpb(jpfrpb))
IF (.NOT. ALLOCATED(frpt)) ALLOCATE(frpt(jpfrpt))
IF (.NOT. ALLOCATED(frpj)) ALLOCATE(frpj(jpfrpj))
IF (.NOT. ALLOCATED(frph)) ALLOCATE(frph(jpfrph))
IF (.NOT. ALLOCATED(frpx)) ALLOCATE(frpx(jpfrpx))
IF (.NOT. ALLOCATED(ztabpd)) ALLOCATE(ztabpd(jpfrpd,2))
IF (.NOT. ALLOCATED(nfrpx)) ALLOCATE(nfrpx(jpnr))
IF (.NOT. ALLOCATED(ntabfp)) ALLOCATE(ntabfp(jpfrpx,3))
IF (.NOT. ALLOCATED(ntabpd)) ALLOCATE(ntabpd(jpfrpd,3))
IF (.NOT. ALLOCATED(npdfr)) ALLOCATE(npdfr(jpnr,2))
IF (.NOT. ALLOCATED(ngrp)) ALLOCATE(ngrp(2*jpspec,3))
IF (.NOT. ALLOCATED(njcgrp)) ALLOCATE(njcgrp(jpcspf,3))
IF (.NOT. ALLOCATED(nprdx3)) ALLOCATE(nprdx3(3,(jpnr/(3*3))+3*3,2*jpspec))
IF (.NOT. ALLOCATED(nprdx2)) ALLOCATE(nprdx2(2,2*jpspec))
IF (.NOT. ALLOCATED(nprdx1)) ALLOCATE(nprdx1(2*jpspec))
IF (.NOT. ALLOCATED(njacx3)) ALLOCATE(njacx3(3,(jpnr/(3*3))+3*3,jpcspf))
IF (.NOT. ALLOCATED(njacx2)) ALLOCATE(njacx2(2,jpcspf))
IF (.NOT. ALLOCATED(njacx1)) ALLOCATE(njacx1(jpcspf))
IF (.NOT. ALLOCATED(nmpjac)) ALLOCATE(nmpjac(jpcspf))
IF (.NOT. ALLOCATED(npjac1)) ALLOCATE(npjac1(jppjac,jpcspf))

IF (.NOT. ALLOCATED(nbrkx)) ALLOCATE(nbrkx(jpbk+1))
IF (.NOT. ALLOCATED(ntrkx)) ALLOCATE(ntrkx(jptk+1))
IF (.NOT. ALLOCATED(nprkx)) ALLOCATE(nprkx(jppj+1))
IF (.NOT. ALLOCATED(nhrkx)) ALLOCATE(nhrkx(jphk+1))
IF (.NOT. ALLOCATED(nlall)) ALLOCATE(nlall(jpspec))
IF (.NOT. ALLOCATED(nlstst)) ALLOCATE(nlstst(jpspec))
IF (.NOT. ALLOCATED(nlf)) ALLOCATE(nlf(jpspec))
IF (.NOT. ALLOCATED(nlmajmin)) ALLOCATE(nlmajmin(jpspec))
IF (.NOT. ALLOCATED(nldepd)) ALLOCATE(nldepd(jpspec))
IF (.NOT. ALLOCATED(nldepw)) ALLOCATE(nldepw(jpspec))
IF (.NOT. ALLOCATED(nldepx)) ALLOCATE(nldepx(jpspec+6))
! nldepx(1:2) = start, end indices (in nldepx) of dry+wet deposited species
! nldepx(3:4) = start, end indices of species undergoing dry deposition only
! nldepx(5:6) = start, end indices of wet deposited species
IF (.NOT. ALLOCATED(njcoth)) ALLOCATE(njcoth(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmzjac)) ALLOCATE(nmzjac(jpcspf))
IF (.NOT. ALLOCATED(nzjac1)) ALLOCATE(nzjac1(jpnr,jpcspf))

IF (.NOT. ALLOCATED(njcoss)) ALLOCATE(njcoss(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmsjac)) ALLOCATE(nmsjac(jpcspf))
IF (.NOT. ALLOCATED(nsjac1)) ALLOCATE(nsjac1(jpnr,jpcspf))

! the following had save attribs.
IF (.NOT. ALLOCATED(ilcf)) ALLOCATE(ilcf(jpspec))
IF (.NOT. ALLOCATED(ilss)) ALLOCATE(ilss(jpspec))
IF (.NOT. ALLOCATED(ilct)) ALLOCATE(ilct(jpspec))
IF (.NOT. ALLOCATED(ilftr)) ALLOCATE(ilftr(jpspec))
IF (.NOT. ALLOCATED(ilft)) ALLOCATE(ilft(jpspec))
IF (.NOT. ALLOCATED(ilstmin)) ALLOCATE(ilstmin(jpspec))


! Set integration method (1 = IMPACT; 3 = N-R solver; 5 = Backward-Euler)
method = ukca_config%ukca_int_method

! Initialize variables that may be changed in cinit
nitnr   = 10            ! Iterations in ftoy
nitfg   = 10            ! Max number of iterations in ftoy

! Initialise arrays
njcoth(:,:) = 0

IF (method == int_method_NR) THEN
  IF (.NOT. ALLOCATED(nonzero_map_unordered))                                  &
    ALLOCATE(nonzero_map_unordered(jpcspf, jpcspf))
  IF (.NOT. ALLOCATED(nonzero_map))  ALLOCATE(nonzero_map(jpcspf, jpcspf))
  IF (.NOT. ALLOCATED(reorder))  ALLOCATE(reorder(jpcspf))

  ! allocate arrays required by solver. These should only be done once, and not
  ! deallocated (hence no matching deallocate statements).
  ! Set size of arrays for NR solver. Needs to be larger for CRI-Strat
  ! mechanism, should be smaller for all others to keep runtime down.
  ! Should be tuned for each mechanism to improve speed
  IF (ukca_config%l_ukca_cristrat) THEN
    spfjsize_max = 3000
    maxterms     = 345
    maxfterms    = 300
  ELSE
    spfjsize_max = 1000
    maxterms     = 160
    maxfterms    = 100
  END IF

  ! Allocate production and loss Jacobian arrays for NR solver
  IF (.NOT. ALLOCATED(nposterms))   ALLOCATE(nposterms(spfjsize_max))
  IF (.NOT. ALLOCATED(nnegterms))   ALLOCATE(nnegterms(spfjsize_max))
  IF (.NOT. ALLOCATED(nfracterms))  ALLOCATE(nfracterms(spfjsize_max))
  IF (.NOT. ALLOCATED(posterms))    ALLOCATE(posterms(spfjsize_max, maxterms))
  IF (.NOT. ALLOCATED(negterms))    ALLOCATE(negterms(spfjsize_max, maxterms))
  IF (.NOT. ALLOCATED(fracterms))   ALLOCATE(fracterms(spfjsize_max, maxfterms))
  IF (.NOT. ALLOCATED(base_tracer)) ALLOCATE(base_tracer(spfjsize_max))
  IF (.NOT. ALLOCATED(ffrac))       ALLOCATE(ffrac(spfjsize_max, maxfterms))
END IF

n_parallel = 1
!$ n_parallel = omp_get_max_teams()

IF (.NOT. ALLOCATED(asad_state)) ALLOCATE(asad_state(n_parallel))

! The arrays which have copies on all threads
DO i = 1, n_parallel
  
  ! pd is a TARGET
  IF (.NOT. ALLOCATED(asad_state(i)%pd))                                       &
    ALLOCATE(asad_state(i)%pd(n_points,2*jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%co3)) ALLOCATE(asad_state(i)%co3(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%deriv))                                    &
    ALLOCATE(asad_state(i)%deriv(n_points,4,4))
  IF (.NOT. ALLOCATED(asad_state(i)%dpd))                                      &
    ALLOCATE(asad_state(i)%dpd(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%dpw))                                      &
    ALLOCATE(asad_state(i)%dpw(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%ej))                                       &
    ALLOCATE(asad_state(i)%ej(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%emr))                                      &
    ALLOCATE(asad_state(i)%emr(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%f))                                        &
    ALLOCATE(asad_state(i)%f(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%fdot))                                     &
    ALLOCATE(asad_state(i)%fdot(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%fj))                                       &
    ALLOCATE(asad_state(i)%fj(n_points,jpcspf,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%fpsc1))                                    &
    ALLOCATE(asad_state(i)%fpsc1(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%fpsc2))                                    &
    ALLOCATE(asad_state(i)%fpsc2(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%ftilde))                                   &
    ALLOCATE(asad_state(i)%ftilde(n_points, jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%ipa))                                      &
    ALLOCATE(asad_state(i)%ipa(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%lati))                                     &
    ALLOCATE(asad_state(i)%lati(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%linfam))                                   &
    ALLOCATE(asad_state(i)%linfam(n_points,0:jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%p))                                        &
    ALLOCATE(asad_state(i)%p(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%pmintnd))                                  &
    ALLOCATE(asad_state(i)%pmintnd(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%prk))                                      &
    ALLOCATE(asad_state(i)%prk(n_points,jpnr))
  IF (.NOT. ALLOCATED(asad_state(i)%qa))                                       &
    ALLOCATE(asad_state(i)%qa(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%ratio))                                    &
    ALLOCATE(asad_state(i)%ratio(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%rk))                                       &
    ALLOCATE(asad_state(i)%rk(n_points,jpnr))
  IF (.NOT. ALLOCATED(asad_state(i)%sh2o))                                     &
    ALLOCATE(asad_state(i)%sh2o(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%shno3))                                    &
    ALLOCATE(asad_state(i)%shno3(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%sph2o))                                    &
    ALLOCATE(asad_state(i)%sph2o(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%sphno3))                                   &
    ALLOCATE(asad_state(i)%sphno3(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%t)) ALLOCATE(asad_state(i)%t(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%t300))                                     &
    ALLOCATE(asad_state(i)%t300(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%tnd)) ALLOCATE(asad_state(i)%tnd(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%wp)) ALLOCATE(asad_state(i)%wp(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%co2)) ALLOCATE(asad_state(i)%co2(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%y))                                        &
    ALLOCATE(asad_state(i)%y(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%ydot))                                     &
    ALLOCATE(asad_state(i)%ydot(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%za)) ALLOCATE(asad_state(i)%za(n_points))
  
  IF (method == int_method_NR) THEN
    IF (.NOT. ALLOCATED(asad_state(i)%modified_map))                           &
      ALLOCATE(asad_state(i)%modified_map(jpcspf, jpcspf))
    IF (.NOT. ALLOCATED(asad_state(i)%spfj))                                   &
      ALLOCATE(asad_state(i)%spfj(n_points,spfjsize_max))
  END IF
  
  asad_state(i)%prod => asad_state(i)%pd(:,1:jpspec)
  asad_state(i)%slos => asad_state(i)%pd(:,jpspec+1:2*jpspec)
  
  ! Initialise arrays
  asad_state(i)%deriv(:,:,:) = 1.0   ! Temp fix for deriv being uninitialised
                                     ! in first solver iteration

  ! Initialise chemistry timestep info
  asad_state(i)%cdt = init_cdt
  asad_state(i)%cdt_diag = init_cdt_diag
  asad_state(i)%interval = init_interval
  asad_state(i)%ncsteps = init_ncsteps
  asad_state(i)%ncsteps_factor = init_ncsteps_factor
  
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mod_init

! ######################################################################

SUBROUTINE asad_mod_init_spatial_vars(n_points)

! To allocate and initialise 2D/3D ASAD arrays and variables
! Called every chemistry timestep

!$ USE omp_lib

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_points

INTEGER :: n_parallel, i
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ASAD_MOD_INIT_SPATIAL_VARS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n_parallel = 1
!$ n_parallel = omp_get_max_threads()

! Set integration method (1 = IMPACT; 3 = N-R solver; 5 = Backward-Euler)
method = ukca_config%ukca_int_method

IF (.NOT. ALLOCATED(asad_state)) ALLOCATE(asad_state(n_parallel))

DO i = 1, n_parallel
  ! prod and slos are pointers
  IF (ASSOCIATED(asad_state(i)%prod))       NULLIFY(asad_state(i)%prod)
  IF (ASSOCIATED(asad_state(i)%slos))       NULLIFY(asad_state(i)%slos)
  
  ! pd is a TARGET
  IF (.NOT. ALLOCATED(asad_state(i)%pd))                                      &
    ALLOCATE(asad_state(i)%pd(n_points,2*jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%co3)) ALLOCATE(asad_state(i)%co3(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%deriv)) &
    ALLOCATE(asad_state(i)%deriv(n_points,4,4))
  IF (.NOT. ALLOCATED(asad_state(i)%dpd)) &
    ALLOCATE(asad_state(i)%dpd(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%dpw)) &
    ALLOCATE(asad_state(i)%dpw(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%ej)) &
    ALLOCATE(asad_state(i)%ej(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%emr)) &
    ALLOCATE(asad_state(i)%emr(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%f))                                        &
    ALLOCATE(asad_state(i)%f(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%fdot))                                     &
    ALLOCATE(asad_state(i)%fdot(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%fj))                                       &
    ALLOCATE(asad_state(i)%fj(n_points,jpcspf,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%fpsc1))                                    &
    ALLOCATE(asad_state(i)%fpsc1(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%fpsc2))                                    &
    ALLOCATE(asad_state(i)%fpsc2(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%ftilde))                                   &
    ALLOCATE(asad_state(i)%ftilde(n_points, jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%ipa))                                      &
    ALLOCATE(asad_state(i)%ipa(n_points,jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%lati))                                     &
    ALLOCATE(asad_state(i)%lati(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%linfam))                                   &
    ALLOCATE(asad_state(i)%linfam(n_points,0:jpcspf))
  IF (.NOT. ALLOCATED(asad_state(i)%p)) ALLOCATE(asad_state(i)%p(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%pmintnd))                                  &
    ALLOCATE(asad_state(i)%pmintnd(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%prk))                                      &
    ALLOCATE(asad_state(i)%prk(n_points,jpnr))
  IF (.NOT. ALLOCATED(asad_state(i)%qa))                                       &
    ALLOCATE(asad_state(i)%qa(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%ratio))                                    &
    ALLOCATE(asad_state(i)%ratio(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%rk))                                       &
    ALLOCATE(asad_state(i)%rk(n_points,jpnr))
  IF (.NOT. ALLOCATED(asad_state(i)%sh2o))                                     &
    ALLOCATE(asad_state(i)%sh2o(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%shno3))                                    &
    ALLOCATE(asad_state(i)%shno3(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%sph2o))                                    &
    ALLOCATE(asad_state(i)%sph2o(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%sphno3))                                   &
    ALLOCATE(asad_state(i)%sphno3(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%t)) ALLOCATE(asad_state(i)%t(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%t300)) ALLOCATE(asad_state(i)%t300(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%tnd)) ALLOCATE(asad_state(i)%tnd(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%wp)) ALLOCATE(asad_state(i)%wp(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%co2)) ALLOCATE(asad_state(i)%co2(n_points))
  IF (.NOT. ALLOCATED(asad_state(i)%y))                                        &
    ALLOCATE(asad_state(i)%y(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%ydot))                                     &
    ALLOCATE(asad_state(i)%ydot(n_points,jpspec))
  IF (.NOT. ALLOCATED(asad_state(i)%za))                                       &
    ALLOCATE(asad_state(i)%za(n_points))
  
  IF (method == int_method_NR) THEN
    IF (.NOT. ALLOCATED(asad_state(i)%spfj))                                   &
      ALLOCATE(asad_state(i)%spfj(n_points,spfjsize_max))
    asad_state(i)%spfj(:,:) = 0.0
  END IF
  
  asad_state(i)%prod => asad_state(i)%pd(:,1:jpspec)
  asad_state(i)%slos => asad_state(i)%pd(:,jpspec+1:2*jpspec)
  
  ! (re-)initialise DERIV array to 1.0 before each call to ASAD_CDRIVE
  ! to ensure bit-comparability when changing domain decomposition
  asad_state(i)%deriv(:,:,:) = 1.0
  
  !     Clear the species arrays
  asad_state(i)%linfam(:,:) = .FALSE.
  
  asad_state(i)%co2(:)      = 0.0
  asad_state(i)%co3(:)      = 0.0
  asad_state(i)%dpd(:,:)    = 0.0
  asad_state(i)%dpw(:,:)    = 0.0
  asad_state(i)%ej(:,:)     = 0.0
  asad_state(i)%emr(:,:)    = 0.0
  asad_state(i)%f(:,:)      = 0.0
  asad_state(i)%fdot(:,:)   = 0.0
  asad_state(i)%fj(:,:,:)   = 0.0
  asad_state(i)%fpsc1(:)    = 0.0
  asad_state(i)%fpsc2(:)    = 0.0
  asad_state(i)%ftilde(:,:) = 0.0
  asad_state(i)%ipa(:,:)    = 0
  asad_state(i)%lati(:)     = 0.0
  asad_state(i)%p(:)        = 0.0
  asad_state(i)%pd(:,:)     = 0.0
  asad_state(i)%pmintnd(:)  = 0.0
  asad_state(i)%prod(:,:)   = 0.0
  asad_state(i)%qa(:,:)     = 0.0
  asad_state(i)%ratio(:,:)  = 0.0
  asad_state(i)%sh2o(:)     = 0.0
  asad_state(i)%shno3(:)    = 0.0
  asad_state(i)%slos(:,:)   = 0.0
  asad_state(i)%sph2o(:)    = 0.0
  asad_state(i)%sphno3(:)   = 0.0
  asad_state(i)%t(:)        = 0.0
  asad_state(i)%t300(:)     = 0.0
  asad_state(i)%tnd(:)      = 0.0
  asad_state(i)%wp(:)       = 0.0
  asad_state(i)%y(:,:)      = 0.0
  asad_state(i)%ydot(:,:)   = 0.0
  asad_state(i)%za(:)       = 0.0
  
  !     Clear the rates and index arrays
  asad_state(i)%rk(:,:)     = 0.0
  asad_state(i)%prk(:,:)    = 0.0
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mod_init_spatial_vars

! ######################################################################

SUBROUTINE asad_mod_dealloc_spatial_vars()
! ----------------------------------------------------------------------
! Description:
! Deallocate the persistent spatial arrays in asad_mod
! ----------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'ASAD_MOD_DEALLOC_SPATIAL_VARS'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle
INTEGER :: i

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(asad_state)) RETURN

! Set integration method (1 = IMPACT; 3 = N-R solver; 5 = Backward-Euler)
method = ukca_config%ukca_int_method

! Deallocate asad mod variables (in reverse order
! to which they are initially allocated in ukca_mod.F90)...
DO i = 1, SIZE(asad_state)
  IF (method == int_method_NR) THEN ! sparse_vars
    IF (ALLOCATED(asad_state(i)%spfj)) DEALLOCATE(asad_state(i)%spfj)
  END IF
  
  IF (ALLOCATED(asad_state(i)%za)) DEALLOCATE(asad_state(i)%za)
  IF (ALLOCATED(asad_state(i)%ydot)) DEALLOCATE(asad_state(i)%ydot)
  IF (ALLOCATED(asad_state(i)%y)) DEALLOCATE(asad_state(i)%y)
  IF (ALLOCATED(asad_state(i)%co2)) DEALLOCATE(asad_state(i)%co2)
  IF (ALLOCATED(asad_state(i)%wp)) DEALLOCATE(asad_state(i)%wp)
  IF (ALLOCATED(asad_state(i)%tnd)) DEALLOCATE(asad_state(i)%tnd)
  IF (ALLOCATED(asad_state(i)%t300)) DEALLOCATE(asad_state(i)%t300)
  IF (ALLOCATED(asad_state(i)%t)) DEALLOCATE(asad_state(i)%t)
  IF (ALLOCATED(asad_state(i)%sphno3)) DEALLOCATE(asad_state(i)%sphno3)
  IF (ALLOCATED(asad_state(i)%sph2o)) DEALLOCATE(asad_state(i)%sph2o)
  IF (ALLOCATED(asad_state(i)%shno3)) DEALLOCATE(asad_state(i)%shno3)
  IF (ALLOCATED(asad_state(i)%sh2o)) DEALLOCATE(asad_state(i)%sh2o)
  IF (ALLOCATED(asad_state(i)%rk)) DEALLOCATE(asad_state(i)%rk)
  IF (ALLOCATED(asad_state(i)%ratio)) DEALLOCATE(asad_state(i)%ratio)
  IF (ALLOCATED(asad_state(i)%qa)) DEALLOCATE(asad_state(i)%qa)
  IF (ALLOCATED(asad_state(i)%prk)) DEALLOCATE(asad_state(i)%prk)
  IF (ALLOCATED(asad_state(i)%pmintnd)) DEALLOCATE(asad_state(i)%pmintnd)
  IF (ALLOCATED(asad_state(i)%p)) DEALLOCATE(asad_state(i)%p)
  IF (ALLOCATED(asad_state(i)%linfam)) DEALLOCATE(asad_state(i)%linfam)
  IF (ALLOCATED(asad_state(i)%lati)) DEALLOCATE(asad_state(i)%lati)
  IF (ALLOCATED(asad_state(i)%ipa)) DEALLOCATE(asad_state(i)%ipa)
  IF (ALLOCATED(asad_state(i)%ftilde)) DEALLOCATE(asad_state(i)%ftilde)
  IF (ALLOCATED(asad_state(i)%fpsc2)) DEALLOCATE(asad_state(i)%fpsc2)
  IF (ALLOCATED(asad_state(i)%fpsc1))  DEALLOCATE(asad_state(i)%fpsc1)
  IF (ALLOCATED(asad_state(i)%fj)) DEALLOCATE(asad_state(i)%fj)
  IF (ALLOCATED(asad_state(i)%fdot)) DEALLOCATE(asad_state(i)%fdot)
  IF (ALLOCATED(asad_state(i)%f)) DEALLOCATE(asad_state(i)%f)
  IF (ALLOCATED(asad_state(i)%emr)) DEALLOCATE(asad_state(i)%emr)
  IF (ALLOCATED(asad_state(i)%ej)) DEALLOCATE(asad_state(i)%ej)
  IF (ALLOCATED(asad_state(i)%dpw)) DEALLOCATE(asad_state(i)%dpw)
  IF (ALLOCATED(asad_state(i)%dpd)) DEALLOCATE(asad_state(i)%dpd)
  IF (ALLOCATED(asad_state(i)%deriv)) DEALLOCATE(asad_state(i)%deriv)
  IF (ALLOCATED(asad_state(i)%co3)) DEALLOCATE(asad_state(i)%co3)
  IF (ALLOCATED(asad_state(i)%pd)) DEALLOCATE(asad_state(i)%pd)
END DO
DEALLOCATE(asad_state)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE asad_mod_dealloc_spatial_vars

END MODULE asad_mod
