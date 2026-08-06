!#PYFT transfo: --buildModi
!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
      SUBROUTINE TURB(CST,CSTURB,BUCONF,TURBN,NEBN,D,TLES,            &
              & KRR,KRRL,KRRI,HLBCX,HLBCY,KGRADIENTSLEO,              &
              & KGRADIENTSGOG,KHALO,                                  &
              & KSPLIT, OCLOUDMODIFLM, KSV,KSV_LGBEG,KSV_LGEND,       &
              & KSV_LIMA_NR, KSV_LIMA_NS, KSV_LIMA_NG, KSV_LIMA_NH,   &
              & O2D,ONOMIXLG,OFLAT,OCOUPLES,OBLOWSNOW,OIBM,OFLYER,    &
              & OCOMPUTE_SRC, PRSNOW,                                 &
              & OOCEAN,ODEEPOC,ODIAG_IN_RUN,                          &
              & HTURBLEN_CL,HCLOUD,HELEC,                             &
              & PTSTEP,TPFILE,                                        &
              & PDXX,PDYY,PDZZ,PDZX,PDZY,PZZ,                         &
              & PDIRCOSXW,PDIRCOSYW,PDIRCOSZW,PCOSSLOPE,PSINSLOPE,    &
              & PRHODJ,PTHVREF,PHGRADLEO,PHGRADGOG,PZS,               &
              & PSFTH,PSFRV,PSFSV,PSFU,PSFV,                          &
              & PSEA_UCU,PSEA_VCU,                                    &
              & PPABST,PUT,PVT,PWT,PTKET,PSVT,PSRCT,                  &
              & PLENGTHM,PLENGTHH,MFMOIST,                            &
              & PBL_DEPTH,PSBL_DEPTH,                                 &
              & PCEI,PCEI_MIN,PCEI_MAX,PCOEF_AMPL_SAT,                &
              & PTHLT,PRT,                                            &
              & PRUS,PRVS,PRWS,PRTHLS,PRRS,PRSVS,PRTKES,              &
              & PSIGS,                                                &
              & PFLXZTHVMF, PFLXZUMF, PFLXZVMF,                       &
              & PWTH,PWRC,PWSV,PDP,PTP,PTDIFF,PTDISS,      &
              & TBUDGETS, KBUDGETS,                                   &
              & PEDR,PLEM,PRTKEMS,PDPMF,PTPMF,                        &
              & PDRUS_TURB,PDRVS_TURB,                                &
              & PDRTHLS_TURB,PDRRTS_TURB,PDRSVS_TURB,PTR,PDISS,       &
              & PIBM_LS, PIBM_XMUT,                                   &
              & PCURRENT_TKE_DISS, PSSTFL, PSSTFL_C, PSSRFL_C,        &
              & PSSUFL_C, PSSVFL_C,PSSUFL,PSSVFL                      )

!$ACDC singlecolumn --nocreate-interface

!     #################################################################
!
!
!!****  *TURB* - computes the turbulent source terms for the prognostic
!!               variables.
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute the source terms in
!!    the evolution equations due to the turbulent mixing.
!!      The source term is computed as the divergence of the turbulent fluxes.
!!    The cartesian fluxes are obtained by a one and a half order closure, based
!!    on a prognostic equation for the Turbulence Kinetic Energy( TKE ). The
!!    system is closed by prescribing a turbulent mixing length. Different
!!    choices are available for this length.
!
!*      0. DECLARATIONS
!          ------------
!
USE MODE_SHUMAN_PHY, ONLY: MZF_PHY,MXF_PHY,MYF_PHY
USE YOMHOOK ,   ONLY: LHOOK, DR_HOOK, JPHOOK
!
USE MODD_BUDGET,     ONLY:  NBUDGET_U,  NBUDGET_V,  NBUDGET_W,  NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, &
                            NBUDGET_RI, NBUDGET_SV1, NBUDGET_RG, NBUDGET_RH, NBUDGET_RR, NBUDGET_RS, &
                            TBUDGETDATA_PTR, TBUDGETCONF_T
USE MODD_CST,        ONLY: CST_t
USE MODD_CTURB,      ONLY: CSTURB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_FIELD,      ONLY: TFIELDMETADATA, TYPEREAL
USE MODD_IO,         ONLY: TFILEDATA
USE MODD_LES,        ONLY: TLES_t
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB, XUNDEF
USE MODD_TURB_n,     ONLY: TURB_t
USE MODD_NEB_n,      ONLY: NEB_t
!
USE MODE_BL89,                ONLY: BL89
USE MODE_CLOUD_MODIF_LM,      ONLY: CLOUD_MODIF_LM
USE MODE_COMPUTE_FUNCTION_THERMO_NEW_STAT, ONLY: COMPUTE_FUNCTION_THERMO_NEW_STAT
USE MODE_COMPUTE_FUNCTION_THERMO, ONLY: COMPUTE_FUNCTION_THERMO
USE MODE_DEAR,                ONLY: DEAR
USE MODE_DELT,                ONLY: DELT
USE MODE_GRADIENT_U_PHY,      ONLY: GZ_U_UW_PHY
USE MODE_GRADIENT_V_PHY,      ONLY: GZ_V_VW_PHY
USE MODE_GRADIENT_W_PHY,      ONLY: GZ_W_M_PHY
USE MODE_GRADIENT_M_PHY,      ONLY: GZ_M_W_PHY
USE MODE_IBM_MIXINGLENGTH,    ONLY: IBM_MIXINGLENGTH
USE MODE_IO_FIELD_WRITE_PHY,      ONLY: IO_FIELD_WRITE_PHY
USE MODE_RMC01,               ONLY: RMC01
USE MODE_ROTATE_WIND,         ONLY: ROTATE_WIND, UPDATE_ROTATE_WIND
USE MODE_SBL_PHY,             ONLY: LMO
USE MODE_SOURCES_NEG_CORRECT, ONLY: SOURCES_NEG_CORRECT_PHY
USE MODE_TM06,                ONLY: TM06
USE MODE_TKE_EPS_SOURCES,     ONLY: TKE_EPS_SOURCES
USE MODE_TURB_HOR_SPLT,       ONLY: TURB_HOR_SPLT
USE MODE_TURB_VER,            ONLY: TURB_VER
USE MODE_UPDATE_LM,           ONLY: UPDATE_LM
USE MODE_MSG,                 ONLY: PRINT_MSG, NVERB_FATAL
!
USE MODI_LES_MEAN_SUBGRID_PHY, ONLY: LES_MEAN_SUBGRID_PHY
USE MODI_SECOND_MNH,          ONLY: SECOND_MNH
!
! These macro are handled by pft_tool.py --craybyPassDOCONCURRENT applied on Cray Rules
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
!$mnh_undef(OPENACC)
#endif
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D             ! PHYEX variables dimensions structure
TYPE(CST_t),            INTENT(IN)   :: CST           ! modd_cst general constant structure
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB        ! modd_csturb turb constant structure
TYPE(TBUDGETCONF_t),    INTENT(IN)   :: BUCONF        ! budget structure
TYPE(TURB_t),           INTENT(IN)   :: TURBN         ! modn_turbn (turb namelist) structure
TYPE(NEB_t),            INTENT(IN)   :: NEBN          ! modd_nebn structure
TYPE(TLES_t),           INTENT(INOUT)   :: TLES          ! modd_les structure
INTEGER,                INTENT(IN)   :: KGRADIENTSLEO ! Number of stored horizontal gradients for Moeng scheme
INTEGER,                INTENT(IN)   :: KGRADIENTSGOG ! Number of stored horizontal gradients for Goger scheme
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL          ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice water var.
INTEGER,                INTENT(IN)   :: KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
INTEGER,                INTENT(IN)   :: KSV_LIMA_NR,KSV_LIMA_NS,KSV_LIMA_NG,KSV_LIMA_NH
CHARACTER(LEN=4),DIMENSION(2),INTENT(IN):: HLBCX, HLBCY  ! X- and Y-direc LBC
INTEGER,                INTENT(IN)   :: KSPLIT        ! number of time-splitting
LOGICAL,                INTENT(IN)   :: OCLOUDMODIFLM ! cloud mixing length modifications
INTEGER,                INTENT(IN)   ::  KHALO        ! Size of the halo for parallel distribution
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ODEEPOC      ! activates sfc forcing for ideal ocean deep conv
LOGICAL,                INTENT(IN)   ::  OFLYER       ! MesoNH flyer diagnostic
LOGICAL,                INTENT(IN)   ::  OFLAT        ! Logical for zero ororography
LOGICAL,                INTENT(IN)   ::  OCOUPLES     ! switch to activate atmos-ocean LES version 
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
LOGICAL,                INTENT(IN)   ::  ODIAG_IN_RUN ! switch to activate online diagnostics (mesonh)
LOGICAL,                INTENT(IN)   ::  OIBM         ! switch to modity mixing length near building with IBM
CHARACTER(LEN=4),       INTENT(IN)   ::  HTURBLEN_CL  ! kind of cloud mixing length
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
CHARACTER (LEN=4),      INTENT(IN)   ::  HELEC        ! Kind of cloud electricity scheme
REAL,                   INTENT(IN)   ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL,                   INTENT(IN)   ::  PTSTEP       ! timestep
TYPE(TFILEDATA),        INTENT(INOUT)   ::  TPFILE       ! Output file
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PDXX,PDYY,PDZZ,PDZX,PDZY
                                        ! metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PZZ       !  physical distance
! between 2 succesive grid points along the K direction
REAL, DIMENSION(D%NIJT),   INTENT(IN)      ::  PDIRCOSXW, PDIRCOSYW, PDIRCOSZW
! Director Cosinus along x, y and z directions at surface w-point
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PCOSSLOPE       ! cosinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PSINSLOPE       ! sinus of the angle
                                 ! between i and the slope vector
REAL, DIMENSION(D%NIJT),   INTENT(IN)   ::  PZS ! orography (for LEONARD terms)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PRHODJ    ! dry density * Grid size
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  MFMOIST ! moist mass flux dual scheme
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PTHVREF   ! Virtual Potential
                                        ! Temperature of the reference state
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%LLEONARD),MERGE(D%NKT,0,TURBN%LLEONARD), &
                MERGE(KGRADIENTSLEO,0,TURBN%LLEONARD)),   INTENT(IN) ::  PHGRADLEO  ! horizontal gradients in Moeng
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%LGOGER),MERGE(D%NKT,0,TURBN%LGOGER),     &
                MERGE(KGRADIENTSGOG,0,TURBN%LGOGER)),     INTENT(IN) ::  PHGRADGOG  ! horizontal gradients in Goger
!
REAL, DIMENSION(D%NIJT),   INTENT(IN)      ::  PSFTH,PSFRV,   &
! normal surface fluxes of theta and Rv
                                            PSFU,PSFV
! normal surface fluxes of (u,v) parallel to the orography
REAL, DIMENSION(D%NIJT,KSV), INTENT(IN)      ::  PSFSV
! normal surface fluxes of Scalar var.
!
REAL, DIMENSION(D%NIJT), INTENT(IN)      ::  PSEA_UCU,PSEA_VCU !Sea surface currents
!
!    prognostic variables at t- deltat
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PPABST      ! Pressure at time t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PUT,PVT,PWT ! wind components
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN) ::  PTKET       ! TKE
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVT        ! passive scal. var.
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)),   INTENT(IN) ::  PSRCT       ! Second-order flux
                      ! s'rc'/2Sigma_s2 at time t-1 multiplied by Lambda_3
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%CTOM=='TM06')),INTENT(INOUT) :: PBL_DEPTH  ! BL height for TOMS
REAL, DIMENSION(MERGE(D%NIJT,0,TURBN%LRMC01))      ,INTENT(INOUT) :: PSBL_DEPTH ! SBL depth for RMC01
!
!    variables for cloud mixing length
REAL, DIMENSION(MERGE(D%NIJT,0,OCLOUDMODIFLM),&
                MERGE(D%NKT,0,OCLOUDMODIFLM)),INTENT(IN)      ::  PCEI
                                                 ! Cloud Entrainment instability
                                                 ! index to emphasize localy
                                                 ! turbulent fluxes
REAL, INTENT(IN)      ::  PCEI_MIN ! minimum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCEI_MAX ! maximum threshold for the instability index CEI
REAL, INTENT(IN)      ::  PCOEF_AMPL_SAT ! saturation of the amplification coefficient
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PTHLT       ! conservative pot. temp.
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRT         ! water var.  where
                             ! PRT(:,:,:,1) is the conservative mixing ratio
!
! sources of momentum, conservative potential temperature, Turb. Kin. Energy,
! TKE dissipation
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) ::  PRUS,PRVS,PRWS,PRTHLS,PRTKES
! Source terms for all water kinds, PRRS(:,:,:,1) is used for the conservative
! mixing ratio
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN),OPTIONAL    ::  PRTKEMS
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(INOUT) ::  PRRS
! Source terms for all passive scalar variables
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(INOUT) ::  PRSVS
! Sigma_s at time t+1 : square root of the variance of the deviation to the
! saturation
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(OUT)     ::  PSIGS
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRUS_TURB   ! evolution of rhoJ*U   by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRVS_TURB   ! evolution of rhoJ*V   by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRTHLS_TURB ! evolution of rhoJ*thl by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL     ::  PDRRTS_TURB  ! evolution of rhoJ*rt  by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(OUT),OPTIONAL ::  PDRSVS_TURB  ! evolution of rhoJ*Sv  by turbulence only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)      ::  PFLXZTHVMF
!                                           MF contribution for vert. turb. transport
!                                           used in the buoy. prod. of TKE
REAL, DIMENSION(MERGE(D%NIJT,0,OFLYER),MERGE(D%NKT,0,OFLYER)), INTENT(OUT)  :: PWTH       ! heat flux
REAL, DIMENSION(MERGE(D%NIJT,0,OFLYER),MERGE(D%NKT,0,OFLYER)), INTENT(OUT)  :: PWRC       ! cloud water flux
REAL, DIMENSION(MERGE(D%NIJT,0,OFLYER),MERGE(D%NKT,0,OFLYER),MERGE(KSV,0,OFLYER)),INTENT(OUT) :: PWSV       ! scalar flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PFLXZUMF   ! MF contribution for vert. turb. transport (dyn prod)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PFLXZVMF   ! MF contribution for vert. turb. transport (dyn prod)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTP        ! Thermal TKE production MassFlux + turb
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL  :: PTPMF      ! Thermal TKE production MassFlux Only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PDP        ! Dynamic TKE production MassFlux + turb
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT),OPTIONAL  :: PDPMF      ! Dynamic TKE production MassFlux Only
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDIFF     ! Diffusion TKE term
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: PTDISS     ! Dissipation TKE term
!
TYPE(TBUDGETDATA_PTR), DIMENSION(KBUDGETS), INTENT(INOUT) :: TBUDGETS
INTEGER, INTENT(IN) :: KBUDGETS
!
LOGICAL, INTENT(IN) :: ONOMIXLG          ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL, INTENT(IN) :: O2D               ! Logical for 2D model version (modd_conf)
!
! length scale from vdfexcu
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PLENGTHM, PLENGTHH
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL  :: PEDR  ! EDR
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT), OPTIONAL  :: PLEM  ! Mixing length
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PTR          ! Transport prod. of TKE
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT), OPTIONAL  ::  PDISS        ! Dissipation of TKE
REAL, DIMENSION(MERGE(D%NIJT,0,ODIAG_IN_RUN),MERGE(D%NKT,0,ODIAG_IN_RUN)), INTENT(INOUT), OPTIONAL  ::  PCURRENT_TKE_DISS ! if ODIAG_IN_RUN in mesonh
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSTFL        ! Time evol Flux of T at sea surface (LOCEAN)
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSTFL_C  ! O-A interface flux for theta(LOCEAN and LCOUPLES)
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSRFL_C  ! O-A interface flux for vapor (LOCEAN and LCOUPLES) 
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSUFL_C        ! Time evol Flux of U at sea surface (LOCEAN)
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSVFL_C  !
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSUFL   
REAL, DIMENSION(MERGE(D%NIJT,0,OCOUPLES)), INTENT(IN),OPTIONAL   ::  PSSVFL  !
!
REAL, DIMENSION(MERGE(D%NIJT,0,OIBM),MERGE(D%NKT,0,OIBM)), INTENT(OUT), OPTIONAL :: PIBM_XMUT ! IBM turbulent viscosity
REAL, DIMENSION(MERGE(D%NIJT,0,OIBM),MERGE(D%NKT,0,OIBM)), INTENT(IN), OPTIONAL  :: PIBM_LS ! IBM Level-set function
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT) ::     &
          ZCP,                        &  ! Cp at t-1
          ZTEMP_BUD
REAL, DIMENSION(D%NIJT,D%NKT,KRR) :: ZRM ! initial mixing ratio
INTEGER :: ISV
!
ISV=1
END SUBROUTINE TURB
