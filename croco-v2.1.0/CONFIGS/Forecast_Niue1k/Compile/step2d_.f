      subroutine step2d (tile)
      implicit none
      integer*4 tile, trd
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=650,   MMm0=690,   N=50)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=NPP)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Ntides
      parameter (Ntides=8)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=LLm, Mm=MMm)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      integer*4   ntrc_mld
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_mld=0)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_mld)
      parameter (NTot=NT)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,9,0:NPP-1)
     &    ,A3dHz(N3dHz,4,0:NPP-1)
      integer*4 B2d(N2d,0:NPP-1)
      common/private_scratch/ A2d,A3d,A3dHz
      common/private_scratch_bis/ B2d
C$    integer*4 omp_get_thread_num
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lm+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lm)/2
      chunk_size_E=(Mm+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mm)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lm)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mm)
      trd=0
C$    trd=omp_get_thread_num()
      call step2D_FB_tile ( Istr,Iend,Jstr,Jend, A2d(1,1,trd)
     &                    , A2d(1, 2,trd), A2d(1, 3,trd), A2d(1, 4,trd)
     &                    , A2d(1, 5,trd), A2d(1, 6,trd), A2d(1, 7,trd)
     &                    , A2d(1, 8,trd), A2d(1, 9,trd)
     &                    , A2d(1,10,trd), A2d(1,11,trd)
     &                    , A2d(1,12,trd), A2d(1,13,trd)
     &                    )
      return
      end
      subroutine step2D_FB_tile (Istr,Iend,Jstr,Jend, zeta_new,
     &                           Dnew,rubar,rvbar,
     &                           Drhs, UFx,UFe,
     &                           VFx,VFe
     &                          ,urhs,vrhs
     &                          ,DUon,DVom
     &                          )
      USE debug
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=650,   MMm0=690,   N=50)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 N_sl
      parameter (N_sl=0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=NPP)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Ntides
      parameter (Ntides=8)
      integer*4 stdout, Np, NpHz, padd_X,padd_E
      parameter (stdout=6)
      parameter (Np=N+1)
      parameter (NpHz=(N+1+N_sl))
      parameter (Lm=LLm, Mm=MMm)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d,N3dHz, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      parameter (N3dHz=size_XI*size_ETA*NpHz)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, NTA, itemp, NTot
      integer*4   ntrc_temp, ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      integer*4   ntrc_subs, ntrc_substot
      integer*4   ntrc_mld
      parameter (itemp=1)
      parameter (ntrc_temp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_mld=0)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_subs=0, ntrc_substot=0)
      parameter (ntrc_sed=0)
      parameter (NTA=itemp+ntrc_salt)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed+ntrc_mld)
      parameter (NTot=NT)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_diapv
      integer*4   ntrc_diaeddy, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_diapv=0)
      parameter (ntrc_diaeddy=0)
      parameter (ntrc_surf=0)
      integer*4 Istr,Iend,Jstr,Jend, i,j, kbak, kold,
     &        imin,imax,jmin,jmax
      integer*4,dimension(2) :: ijmax
      real sum_c
      real    VMAX,VMAXL
      real zeta_new(Istr-2:Iend+2,Jstr-2:Jend+2),  cff,
     &         Dnew(Istr-2:Iend+2,Jstr-2:Jend+2),  cff0,
     &        rubar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff1,
     &        rvbar(Istr-2:Iend+2,Jstr-2:Jend+2),  cff2,
     &         Drhs(Istr-2:Iend+2,Jstr-2:Jend+2),  cff3,
     &          UFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &          UFe(Istr-2:Iend+2,Jstr-2:Jend+2),  DUnew,
     &          VFx(Istr-2:Iend+2,Jstr-2:Jend+2),  DVnew,
     &          VFe(Istr-2:Iend+2,Jstr-2:Jend+2)
      real     urhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         vrhs(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         DUon(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         DVom(Istr-2:Iend+2,Jstr-2:Jend+2)
      real h(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real hinv(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real f(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real fomn(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real angler(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_angler/angler
      real latr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real lonr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real latu(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real lonu(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real latv(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real lonv(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /grid_latr/latr /grid_lonr/lonr
      common /grid_latu/latu /grid_lonu/lonu
      common /grid_latv/latv /grid_lonv/lonv
      real pm(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real om_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real on_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pm_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pm_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pn_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real dmde(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real dndx(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx
      real pmon_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmon_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmon_u(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pnom_v(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real grdscl(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real umask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real vmask(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real pmask2(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2
      real zob(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /Z0B_VAR/zob
      real zeta(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      real ubar(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      real vbar(0:Lm+1+padd_X,0:Mm+1+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
      real u(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3)
      real v(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3)
      real t(0:Lm+1+padd_X,0:Mm+1+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real Hz_bak(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real z_r(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real z_w(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real Huon(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real Hvom(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      real rho1(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real rho(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      real rhoA(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rhoS(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rvfrc(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rufrc_bak(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real rvfrc_bak(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real DU_avg1(0:Lm+1+padd_X,0:Mm+1+padd_E,5)
      real DV_avg1(0:Lm+1+padd_X,0:Mm+1+padd_E,5)
      real DU_avg2(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real DV_avg2(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
      real sustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real svstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real sustrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real svstrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer*4 itsms, sms_ncycle, sms_rec, lsusgrd
      integer*4 lsvsgrd,sms_tid, susid, svsid
      real    sms_origin_date_in_sec
      common /smsdat1/ sustrp, svstrp, sms_time
      common /smsdat2/ sms_origin_date_in_sec
      common /smsdat3/ sms_cycle, sms_scale
      common /smsdat4/ itsms, sms_ncycle, sms_rec, lsusgrd
      common /smsdat5/ lsvsgrd,sms_tid, susid, svsid
      integer*4 lwgrd, wid
      common /smsdat5/ lwgrd, wid
      real bustr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real bvstr(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real bvstrg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen, bms_tstart, bms_tend, tsbms, sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /forces_stflx/stflx
      real shflx_rsw(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /frc_shflx_rsw/shflx_rsw
      real shflx_rlw(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /frc_shflx_rlw/shflx_rlw
      real shflx_lat(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /frc_shflx_lat/shflx_lat
      real shflx_sen(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /frc_shflx_sen/shflx_sen
      real stflxg(0:Lm+1+padd_X,0:Mm+1+padd_E,2,NT)
      common /stfdat_stflxg/stflxg
      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer*4 itstf(NT), stf_ncycle(NT), stf_rec(NT)
      integer*4 lstfgrd(NT), stf_tid(NT), stf_id(NT)
      REAL(kind=8) :: stf_origin_date_in_sec
      common /stfdat1/ stflxp,  stf_time, stf_cycle, stf_scale
      common /stfdat2/ stf_origin_date_in_sec
      common /stfdat3/ itstf, stf_ncycle, stf_rec, lstfgrd
      common /stfdat4/  stf_tid, stf_id
      real btflx(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /forces_btflx/btflx
      real tair(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real rhum(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real prate(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real radlw(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real radsw(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real wspd(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real uwnd(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real vwnd(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /bulk_tair/ tair
      common /bulk_rhum/ rhum
      common /bulk_prate/ prate
      common /bulk_radlw/ radlw
      common /bulk_radsw/ radsw
      common /bulk_wspd/ wspd
      common /bulk_uwnd/ uwnd
      common /bulk_vwnd/ vwnd
      real tairg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real rhumg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real prateg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real radlwg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real radswg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real uwndg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real vwndg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real wspdg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /bulkdat_tairg/tairg
      common /bulkdat_rhumg/rhumg
      common /bulkdat_prateg/prateg
      common /bulkdat_radlwg/radlwg
      common /bulkdat_radswg/radswg
      common /bulk_uwndg/uwndg
      common /bulk_vwndg/vwndg
      common /bulkdat_wspdg/wspdg
      real    tairp(2),rhump(2),pratep(2),radlwp(2),radswp(2)
      real    uwndp(2),vwndp(2)
      real    bulk_time(2), bulk_cycle
      integer*4 tair_id,rhum_id,prate_id,radlw_id,radsw_id
      integer*4 ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      REAL(kind=8) :: blk_origin_date_in_sec
      integer*4 uwnd_id,vwnd_id,luwndgrd,lvwndgrd
      integer*4 itbulk,bulk_ncycle,bulk_rec,bulk_tid
      integer*4 bulkunused
      common /bulkdat1_for/ tair_id,rhum_id,prate_id,radlw_id,radsw_id
      common /bulkdat1_grd
     &                 / ltairgrd,lrhumgrd,lprategrd,lradlwgrd,lradswgrd
      common /bulkdat1_tim/ itbulk, bulk_ncycle, bulk_rec, bulk_tid
      common /bulkdat1_uns/ bulkunused
      common /bulkdat1_wnd/ uwnd_id,vwnd_id,luwndgrd,lvwndgrd
      common /bulkdat2_for/ tairp,rhump,pratep,radlwp,radswp
      common /bulkdat2_tim/ bulk_time, bulk_cycle, 
     &                                            blk_origin_date_in_sec
      common /bulkdat2_wnd/ uwndp,vwndp
      real srflx(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /forces_srflx/srflx
      real srflxg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /srfdat_srflxg/srflxg
      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer*4 itsrf, srf_ncycle, srf_rec
      integer*4 lsrfgrd, srf_tid, srf_id
      REAL(kind=8) :: srf_origin_date_in_sec
      common /srfdat1/ srflxp, srf_time, srf_cycle, srf_scale
      common /srfdat2/ srf_origin_date_in_sec
      common /srfdat3/ itsrf,srf_ncycle,srf_rec,lsrfgrd,srf_tid,srf_id
      real visc2_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real visc2_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real visc2_sponge_r(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real visc2_sponge_p(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real diff2(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real Akv(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real Akt(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real bvf(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      common /mixing_bvf/ bvf
      real ustar(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /lmd_kpp_ustar/ustar
      integer*4 kbl(0:Lm+1+padd_X,0:Mm+1+padd_E)
      integer*4 kbbl(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real hbbl(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /lmd_kpp_kbl/ kbl
      common /lmd_kpp_hbbl/ hbbl
      common /lmd_kpp_kbbl/ kbbl
      real hbls(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /lmd_kpp_hbl/ hbls
      real ghats(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      common /lmd_kpp_ghats/ghats
      real ssh(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /climat_ssh/ssh
      real Znudgcof(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /climat_Znudgcof/Znudgcof
      real sshg(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /climat_sshg/sshg
      real    ssh_time(2)
      real    ssh_cycle
      integer*4 itssh, ssh_ncycle, ssh_rec, ssh_tid, ssh_id
      REAL(kind=8) :: ssh_origin_date_in_sec
      common /climat_zdat1/ ssh_time, ssh_origin_date_in_sec
      common /climat_zdat2/ ssh_cycle
      common /climat_zdat3/
     &        itssh, ssh_ncycle, ssh_rec, ssh_tid, ssh_id
      real tclm(0:Lm+1+padd_X,0:Mm+1+padd_E,N,NT)
      common /climat_tclm/tclm
      real Tnudgcof(0:Lm+1+padd_X,0:Mm+1+padd_E,N,NT)
      common /climat_Tnudgcof/Tnudgcof
      real tclima(0:Lm+1+padd_X,0:Mm+1+padd_E,N,2,NT)
      common /climat_tclima/tclima
      real tclm_time(2,NT)
      real tclm_cycle(NT)
      integer*4 ittclm(NT), tclm_ncycle(NT), tclm_rec(NT),
     &        tclm_tid(NT), tclm_id(NT)
      logical got_tclm(NT)
      REAL(kind=8) :: tclm_origin_date_in_sec
      common /climat_tdat/  tclm_time,       tclm_cycle,
     &        ittclm,       tclm_ncycle,     tclm_rec,
     &                      tclm_tid,        tclm_id,
     &                                       got_tclm,
     &                        tclm_origin_date_in_sec
      real ubclm(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real vbclm(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /climat_ubclm/ubclm /climat_vbclm/vbclm
      real uclm(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real vclm(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /climat_uclm/uclm /climat_vclm/vclm
      real M2nudgcof(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /climat_M2nudgcof/M2nudgcof
      real ubclima(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      real vbclima(0:Lm+1+padd_X,0:Mm+1+padd_E,2)
      common /climat_ubclima/ubclima /climat_vbclima/vbclima
      real M3nudgcof(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /climat_M3nudgcof/M3nudgcof
      real uclima(0:Lm+1+padd_X,0:Mm+1+padd_E,N,2)
      real vclima(0:Lm+1+padd_X,0:Mm+1+padd_E,N,2)
      common /climat_uclima/uclima /climat_vclima/vclima
      real     uclm_time(2)
      real     uclm_cycle
      integer*4 ituclm, uclm_ncycle, uclm_rec, uclm_tid,
     &        ubclm_id, vbclm_id, uclm_id, vclm_id
      REAL(kind=8) :: uclm_origin_date_in_sec
      common /climat_udat1/  uclm_time, uclm_origin_date_in_sec
      common /climat_udat2/  uclm_cycle
      common /climat_udat3/
     &             ituclm,   uclm_ncycle, uclm_rec,
     &             uclm_tid, ubclm_id,    vbclm_id,
     &             uclm_id,  vclm_id
      real dt, dtfast, time, time2, time_start, tdays, start_time
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &     ndtfast, iic, kstp, krhs, knew, next_kstp,
     &     start_time,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zobt
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      logical ldefhis
      logical got_tini(NT)
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zobt,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1
     &           ,       tnu2,    tnu4
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , got_tini
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 Cu_Adv3d,  Cu_W, Cu_Nbq_X, Cu_Nbq_Y, Cu_Nbq_Z
      integer*4 i_cx_max, j_cx_max, k_cx_max
      common /diag_vars/ Cu_Adv3d,  Cu_W,
     &        i_cx_max, j_cx_max, k_cx_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
     &        , bc_flux, ubar_xs
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
     &        , bc_flux,  ubar_xs
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, Erotation, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  Erotation=7.292115090D-5,
     &           day2sec=86400.D0, sec2day=1.D0/86400.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0,
     &           jul_off=2440000.D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      real Tperiod(Ntides)
      common /tides_Tperiod/ Tperiod
      real SSH_Tamp(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_SSH_Tamp/ SSH_Tamp
      real SSH_Tphase(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_SSH_Tphase/ SSH_Tphase
      real UV_Tangle(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_UV_Tangle/ UV_Tangle
      real UV_Tmajor(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_UV_Tmajor/ UV_Tmajor
      real UV_Tminor(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_UV_Tminor/ UV_Tminor
      real UV_Tphase(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_UV_Tphase/ UV_Tphase
      real POT_Tamp(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_POT_Tamp/ POT_Tamp
      real POT_Tphase(0:Lm+1+padd_X,0:Mm+1+padd_E,Ntides)
      common /tides_POT_Tphase/ POT_Tphase
      real PTide(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /tides_Ptide/ Ptide
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
      if (istr.eq.1) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (iend.eq.Lm) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (jstr.eq.1) then
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (jend.eq.Mm) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      if (iif.eq.1) then
        kbak=kstp
        kold=kstp
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iif.eq.1+1) then
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      else
        kbak=kstp-1
        if (kbak.lt.1) kbak=4
        kold=kbak-1
        if (kold.lt.1) kold=4
        cff1= 1.781105D0
        cff2=-1.06221D0
        cff3= 0.281105D0
      endif
      imin=IstrU-2
      imax=Iend+1
      jmin=JstrV-2
      jmax=Jend+1
!$acc kernels if(compute_on_device) default(present) async(1)
      do j=jmin,jmax
        do i=imin,imax
          Drhs(i,j)=cff1*zeta(i,j,kstp)+cff2*zeta(i,j,kbak)
     &                                 +cff3*zeta(i,j,kold)
     &                                             + h(i,j)
        enddo
      enddo
      do j=Jstr-1,Jend+1
        do i=imin+1,imax
          urhs(i,j)=cff1*ubar(i,j,kstp) +cff2*ubar(i,j,kbak)
     &                                  +cff3*ubar(i,j,kold)
          DUon(i,j)=0.5D0*(Drhs(i,j)+Drhs(i-1,j))*on_u(i,j)*( urhs(i,j)
     &                                                              )
        enddo
      enddo
      do j=jmin+1,jmax
        do i=Istr-1,Iend+1
          vrhs(i,j)=cff1*vbar(i,j,kstp) +cff2*vbar(i,j,kbak)
     &                                  +cff3*vbar(i,j,kold)
          DVom(i,j)=0.5D0*(Drhs(i,j)+Drhs(i,j-1))*om_v(i,j)*( vrhs(i,j)
     &                                                              )
        enddo
      enddo
      call set_DUV_bc_tile (Istr,Iend,Jstr,Jend, Drhs, DUon,DVom)
      if (iif.eq.1) then
        cff0=0.D0
        cff1=1.D0
        cff2=0.D0
        cff3=0.D0
      elseif (iif.eq.1+1) then
        cff0= 1.0833333333333D0
        cff1=-0.1666666666666D0
        cff2= 0.0833333333333D0
        cff3= 0.D0
      else
        cff0=0.614D0
        cff1=0.285D0
        cff2=0.088D0
        cff3=0.013D0
      endif
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta(i,j,kstp) + dtfast*pm(i,j)*pn(i,j)
     &                                   *(DUon(i,j)-DUon(i+1,j  )
     &                                    +DVom(i,j)-DVom(i  ,j+1))
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta_new(i,j) + dtfast*Znudgcof(i,j)
     &                                 *(ssh(i,j)-zeta_new(i,j))
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          zeta_new(i,j)=zeta_new(i,j)*rmask(i,j)
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          UFx(i,j)=cff0*zeta_new(i,j) +cff1*zeta(i,j,kstp)
     &             +cff2*zeta(i,j,kbak)+cff3*zeta(i,j,kold)
          UFe(i,j)=(1.D0+rhoS(i,j))*UFx(i,j)
          VFe(i,j)=UFe(i,j)*UFx(i,j)
          VFx(i,j)=UFx(i,j)*(rhoS(i,j)-rhoA(i,j))
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          Dnew(i,j)=zeta_new(i,j)+h(i,j)
          zeta(i,j,knew)=zeta_new(i,j)
        enddo
      enddo
!$acc end kernels
      call zetabc_tile (Istr,Iend,Jstr,Jend)
      cff1=weight(1,iif)
      cff2=weight(2,iif)
!$acc kernels if(compute_on_device) default(present) async(1)
      if (iif.eq.1) then
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=cff1*zeta(i,j,knew)
            DU_avg1(i,j,nnew)=0.D0
            DV_avg1(i,j,nnew)=0.D0
            DU_avg2(i,j)=cff2*DUon(i,j)
            DV_avg2(i,j)=cff2*DVom(i,j)
          enddo
        enddo
      else
        do j=JstrR,JendR
          do i=IstrR,IendR
            Zt_avg1(i,j)=Zt_avg1(i,j)+cff1*zeta(i,j,knew)
            DU_avg2(i,j)=DU_avg2(i,j)+cff2*DUon(i,j)
            DV_avg2(i,j)=DV_avg2(i,j)+cff2*DVom(i,j)
          enddo
        enddo
      endif
      cff=0.5D0*g
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=cff*on_u(i,j)*(
     &                         (h(i-1,j)+h(i,j))*(UFe(i-1,j)
     &                        -UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &              +(h(i-1,j)-h(i,j))*( VFx(i-1,j)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i-1,j)-rhoA(i,j))
     &                                      *(UFx(i-1,j)-UFx(i,j)))
     &                                                              )
          rvbar(i,j)=cff*om_v(i,j)*(
     &            (h(i,j-1)+h(i,j))*(UFe(i,j-1)
     &                        -UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &              +(h(i,j-1)-h(i,j))*( VFx(i,j-1)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i,j-1)-rhoA(i,j))
     &                                      *(UFx(i,j-1)-UFx(i,j)))
     &                                                              )
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr-1,Iend
          UFx(i,j)=0.25D0*(DUon(i,j)+DUon(i+1,j))
     &                 *(urhs(i,j)+urhs(i+1,j))
          VFx(i+1,j)=0.25D0*(DUon(i+1,j)+DUon(i+1,j-1))
     &                   *(vrhs(i+1,j)+vrhs(i,j))
     &                                 *pmask(i+1,j)
        enddo
      enddo
      do j=Jstr-1,Jend
        do i=Istr,Iend
          VFe(i,j)=0.25D0*(DVom(i,j)+DVom(i,j+1))
     &                 *(vrhs(i,j)+vrhs(i,j+1))
          UFe(i,j+1)=0.25D0*(DVom(i,j+1)+DVom(i-1,j+1))
     &                   *(urhs(i,j+1)+urhs(i,j))
     &                                 *pmask(i,j+1)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          rubar(i,j)=rubar(i,j)-UFx(i,j)+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i,j)
          rvbar(i,j)=rvbar(i,j)-VFx(i+1,j)+VFx(i,j)
     &                         -VFe(i,j)+VFe(i,j-1)
        enddo
      enddo
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          cff=Drhs(i,j)*(
     &                   fomn(i,j)
     &          +0.5D0*( dndx(i,j)*(vrhs(i,j)+vrhs(i,j+1))
     &                -dmde(i,j)*(urhs(i,j)+urhs(i+1,j)))
     &                   )
          UFx(i,j)=cff*(vrhs(i,j)+vrhs(i,j+1))
          VFe(i,j)=cff*(urhs(i,j)+urhs(i+1,j))
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          rubar(i,j)=rubar(i,j)+0.25D0*(UFx(i,j)+UFx(i-1,j))
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          rvbar(i,j)=rvbar(i,j)-0.25D0*(VFe(i,j)+VFe(i,j-1))
        enddo
      enddo
      if (rdrg2.gt.0.D0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25D0*( vbar(i  ,j,kstp)+vbar(i  ,j+1,kstp)
     &                +vbar(i-1,j,kstp)+vbar(i-1,j+1,kstp))
            rubar(i,j)=rubar(i,j) - ubar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(ubar(i,j,kstp)*ubar(i,j,kstp)+cff*cff)
     &                               )*om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25D0*( ubar(i,j  ,kstp)+ubar(i+1,j  ,kstp)
     &                +ubar(i,j-1,kstp)+ubar(i+1,j-1,kstp))
            rvbar(i,j)=rvbar(i,j) - vbar(i,j,kstp)*( rdrg+rdrg2
     &              *sqrt(cff*cff+vbar(i,j,kstp)*vbar(i,j,kstp))
     &                               )*om_v(i,j)*on_v(i,j)
          enddo
        enddo
      else if (rdrg.gt.0.0D0) then
        do j=Jstr,Jend
          do i=IstrU,Iend
            rubar(i,j)=rubar(i,j) - rdrg*ubar(i,j,kstp)
     &                             *om_u(i,j)*on_u(i,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rvbar(i,j)=rvbar(i,j) - rdrg*vbar(i,j,kstp)
     &                             *om_v(i,j)*on_v(i,j)
          enddo
        enddo
      endif
!$acc end kernels
      if (iif.eq.1) then
        if (iic.eq.ntstart) then
          cff3=0.D0
          cff2=0.D0
          cff1=1.D0
        elseif (iic.eq.ntstart+1) then
          cff3=0.D0
          cff2=-0.5D0
          cff1=1.5D0
        else
          cff3=0.281105D0
          cff2=-0.5D0-2.D0*cff3
          cff1=1.5D0+cff3
        endif
!$acc kernels if(compute_on_device) default(present) async(1)
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=rufrc(i,j)-rubar(i,j)
            rufrc(i,j)=cff1*cff + cff2*rufrc_bak(i,j,3-nstp)
     &                          + cff3*rufrc_bak(i,j,nstp)
            rufrc_bak(i,j,nstp)=cff
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=rvfrc(i,j)-rvbar(i,j)
            rvfrc(i,j)=cff1*cff + cff2*rvfrc_bak(i,j,3-nstp)
     &                          + cff3*rvfrc_bak(i,j,nstp)
            rvfrc_bak(i,j,nstp)=cff
          enddo
        enddo
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            UFx(i,j)=zeta_new(i,j)-zeta(i,j,kstp)
            UFe(i,j)=(1.D0+rhoS(i,j))*UFx(i,j)
            VFe(i,j)=UFe(i,j)*(zeta_new(i,j)+zeta(i,j,kstp))
            VFx(i,j)=UFx(i,j)*(rhoS(i,j)-rhoA(i,j))
          enddo
        enddo
        cff=0.5D0*g
        do j=Jstr,Jend
          do i=Istr,Iend
            rubar(i,j)=rubar(i,j) +cff*on_u(i,j)*( (h(i-1,j)+h(i,j))
     &          *(UFe(i-1,j)-UFe(i,j)) +VFe(i-1,j)-VFe(i,j)
     &              +(h(i-1,j)-h(i,j))*( VFx(i-1,j)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i-1,j)-rhoA(i,j))
     &                                     *(UFx(i-1,j)-UFx(i,j)) )
     &                                                              )
            rvbar(i,j)=rvbar(i,j) +cff*om_v(i,j)*( (h(i,j-1)+h(i,j))
     &          *(UFe(i,j-1)-UFe(i,j)) +VFe(i,j-1)-VFe(i,j)
     &              +(h(i,j-1)-h(i,j))*( VFx(i,j-1)+VFx(i,j)
     &                        +0.333333333333D0*(rhoA(i,j-1)-rhoA(i,j))
     &                                     *(UFx(i,j-1)-UFx(i,j)) )
     &                                                              )
          enddo
        enddo
!$acc end kernels
      endif
!$acc kernels if(compute_on_device) default(present) async(1)
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          DUon(i,j)=zeta(i,j,kstp)+h(i,j)
        enddo
      enddo
      cff=0.5D0*dtfast
      cff1=0.5D0*weight(1,iif)
      do j=Jstr,Jend
        do i=IstrU,Iend
          DUnew=( (DUon(i,j)+DUon(i-1,j))*ubar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
     &                             *(rubar(i,j)+rufrc(i,j))
     &                                                    )
     &                                         *umask(i,j)
          ubar(i,j,knew)=DUnew/(Dnew(i,j)+Dnew(i-1,j))
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*on_u(i,j)*( DUnew
     &                                                  )
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          DVnew=( (DUon(i,j)+DUon(i,j-1))*vbar(i,j,kstp)
     &        +cff*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
     &                             *(rvbar(i,j)+rvfrc(i,j))
     &                                                    )
     &                                         *vmask(i,j)
          vbar(i,j,knew)=DVnew/(Dnew(i,j)+Dnew(i,j-1))
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*om_v(i,j)*(DVnew
     &                                                   )
        enddo
      enddo
      do j=Jstr,Jend
        do i=IstrU,Iend
          DUnew = dtfast*M2nudgcof(i,j)*(ubclm(i,j)-ubar(i,j,knew))
     &                 * umask(i,j)
          ubar(i,j,knew)=ubar(i,j,knew) + DUnew
          DU_avg1(i,j,nnew)=DU_avg1(i,j,nnew) +cff1*DUnew*
     &                         (Dnew(i,j)+Dnew(i-1,j))*on_u(i,j)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          DVnew = dtfast*M2nudgcof(i,j)*(vbclm(i,j)-vbar(i,j,knew))
     &                 * vmask(i,j)
          vbar(i,j,knew)=vbar(i,j,knew) + DVnew
          DV_avg1(i,j,nnew)=DV_avg1(i,j,nnew) +cff1*DVnew*
     &                         (Dnew(i,j)+Dnew(i,j-1))*om_v(i,j)
        enddo
      enddo
!$acc end kernels
      call u2dbc_tile (Istr,Iend,Jstr,Jend, UFx)
      call v2dbc_tile (Istr,Iend,Jstr,Jend, UFx)
      call obc_flux_tile (Istr,Iend,Jstr,Jend)
!$acc kernels if(compute_on_device) default(present) async(1)
      if (istr.eq.1) then
        do j=Jstr-1,JendR
          Dnew(Istr-1,j)=h(Istr-1,j)+zeta(Istr-1,j,knew)
        enddo
      endif
      if (iend.eq.Lm) then
        do j=Jstr-1,JendR
          Dnew(Iend+1,j)=h(Iend+1,j)+zeta(Iend+1,j,knew)
        enddo
      endif
      if (jstr.eq.1) then
        do i=Istr-1,IendR
          Dnew(i,Jstr-1)=h(i,Jstr-1)+zeta(i,Jstr-1,knew)
        enddo
      endif
      if (jend.eq.Mm) then
        do i=Istr-1,IendR
          Dnew(i,Jend+1)=h(i,Jend+1)+zeta(i,Jend+1,knew)
        enddo
      endif
!$acc end kernels
      cff1=0.5D0*weight(1,iif)
!$acc kernels if(compute_on_device) default(present) async(1)
      if (istr.eq.1) then
        do j=JstrR,JendR
          DU_avg1(IstrU-1,j,nnew)=DU_avg1(IstrU-1,j,nnew)
     &         +cff1*(Dnew(IstrU-1,j)
     &         +Dnew(IstrU-2,j))*(ubar(IstrU-1,j,knew)
     &                                             )*on_u(IstrU-1,j)
        enddo
        do j=JstrV,Jend
          DV_avg1(Istr-1,j,nnew)=DV_avg1(Istr-1,j,nnew)
     &       +cff1*(Dnew(Istr-1,j)
     &       +Dnew(Istr-1,j-1) )*(vbar(Istr-1,j,knew)
     &                                              )*om_v(Istr-1,j)
        enddo
      endif
      if (iend.eq.Lm) then
        do j=JstrR,JendR
          DU_avg1(Iend+1,j,nnew)=DU_avg1(Iend+1,j,nnew)
     &            +cff1*( Dnew(Iend+1,j)
     &            +Dnew(Iend,j) )*(ubar(Iend+1,j,knew)
     &                                              )*on_u(Iend+1,j)
        enddo
        do j=JstrV,Jend
          DV_avg1(Iend+1,j,nnew)=DV_avg1(Iend+1,j,nnew)
     &        +cff1*( Dnew(Iend+1,j)
     &        +Dnew(Iend+1,j-1) )*(vbar(Iend+1,j,knew)
     &                                              )*om_v(Iend+1,j)
        enddo
      endif
      if (jstr.eq.1) then
        do i=IstrU,Iend
          DU_avg1(i,Jstr-1,nnew)=DU_avg1(i,Jstr-1,nnew)
     &        +cff1*( Dnew(i,Jstr-1)
     &        +Dnew(i-1,Jstr-1) )*(ubar(i,Jstr-1,knew)
     &                                              )*on_u(i,Jstr-1)
        enddo
        do i=IstrR,IendR
          DV_avg1(i,JstrV-1,nnew)=DV_avg1(i,JstrV-1,nnew)
     &         +cff1*(Dnew(i,JstrV-1)
     &         +Dnew(i,JstrV-2))*(vbar(i,JstrV-1,knew)
     &                                              )*om_v(i,JstrV-1)
        enddo
      endif
      if (jend.eq.Mm) then
        do i=IstrU,Iend
          DU_avg1(i,Jend+1,nnew)=DU_avg1(i,Jend+1,nnew)
     &        +cff1*( Dnew(i,Jend+1)
     &        +Dnew(i-1,Jend+1) )*(ubar(i,Jend+1,knew)
     &                                               )*on_u(i,Jend+1)
        enddo
        do i=IstrR,IendR
          DV_avg1(i,Jend+1,nnew)=DV_avg1(i,Jend+1,nnew)
     &            +cff1*( Dnew(i,Jend+1)
     &            +Dnew(i,Jend) )*(vbar(i,Jend+1,knew)
     &                                               )*om_v(i,Jend+1)
        enddo
      endif
!$acc end kernels
!$acc kernels if(compute_on_device) default(present) async(1)
      VMAXL=100.D0
      VMAX=0.D0
      do j=Jstr,Jend
        do i=Istr,Iend
          VMAX=max(VMAX,abs(ubar(i,j,knew)))
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          VMAX=max(VMAX,abs(vbar(i,j,knew)))
        enddo
      enddo
!$acc end kernels
      IF (VMAX.GT.VMAXL) THEN
!$acc update if(compute_on_device) host(ubar(Istr:Iend,Jstr:Jend,knew))
!$acc update if(compute_on_device) host(vbar(Istr:Iend,Jstr:Jend,knew))
        cff1=maxval(abs(ubar(Istr:Iend,Jstr:Jend,knew)))
        cff2=maxval(abs(vbar(Istr:Iend,Jstr:Jend,knew)))
        IF (cff1>cff2) THEN
          ijmax= maxloc(abs(ubar(Istr:Iend,Jstr:Jend,knew)))
        ELSE
          ijmax= maxloc(abs(vbar(Istr:Iend,Jstr:Jend,knew)))
        ENDIF
        imax=ijmax(1)
        jmax=ijmax(2)
        write(stdout,'(9(A/))')
     &     '                                         ',
     &     '                                         ',
     &     ' ======================================= ',
     &     ' =                                     = ',
     &     ' =   STEP2D:   ABNORMAL JOB END        = ',
     &     ' =                 BLOW UP             = ',
     &     ' =                                     = ',
     &     ' ======================================= ',
     &     '                                         '
        if (VMAX.eq.666.D0) then
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =   NaN'
        else
          write(stdout,'(A,F10.2)')
     &                                            '  VMAX (M/S) =',VMAX
        endif
        write(stdout,'(A,2I6)')
     &                                       '  IMAX JMAX  =',imax,jmax
        write(stdout,'(A,2I6/)')
     &                                         '  IINT IEXT  =',iic,iif
        may_day_flag=1
        stop
      ENDIF
      if (iif.eq.nfast) then
!$acc wait
      endif
      return
      end
