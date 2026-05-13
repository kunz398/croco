      subroutine lmd_bkpp_tile (Istr,Iend,Jstr,Jend, Kv,Kt,Ks,
     &                     Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS,
     &                        wm,ws, my_hbbl,my_kbbl,wrk, Rib)
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
      integer*4 Istr,Iend,Jstr,Jend, i,j,k, ka,ku,ksave,
     &        imin,imax,jmin,jmax
      real    Kv(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &        Kt(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &        Ks(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &           Gm1(Istr-2:Iend+2,Jstr-2:Jend+2),
     &        dGm1dS(Istr-2:Iend+2,Jstr-2:Jend+2),
     &           Gt1(Istr-2:Iend+2,Jstr-2:Jend+2),
     &        dGt1dS(Istr-2:Iend+2,Jstr-2:Jend+2),
     &           Gs1(Istr-2:Iend+2,Jstr-2:Jend+2),
     &        dGs1dS(Istr-2:Iend+2,Jstr-2:Jend+2),
     &            wm(Istr-2:Iend+2,Jstr-2:Jend+2),
     &            ws(Istr-2:Iend+2,Jstr-2:Jend+2),
     &       my_hbbl(Istr-2:Iend+2,Jstr-2:Jend+2),
     &       my_kbbl(Istr-2:Iend+2,Jstr-2:Jend+2),
     &           wrk(Istr-2:Iend+2,Jstr-2:Jend+2),
     &         Rib(Istr-2:Iend+2,Jstr-2:Jend+2,0:N)
      real Vtc,    hekman,  dVsq,   Vtsq,
     &     sig,    Kv_bl,   Kt_bl,  Ks_bl,
     &     cff,    lmd_a1,  dKv_bl, dKt_bl, dKs_bl,
     &     cff_up, lmd_a2,  Gm,     Gt,     Gs,
     &     cff_dn, Ritop,   lmd_a3,
     &     zbl,    zsbl,    eps
      parameter (eps=1.D-20)
      real lmd_cs, lmd_Cv, Ric, lmd_betaT, lmd_epsilon,
     &     lmd_cekman, lmd_nu0c,
     &     Kv0, Kt0, Ks0
      parameter (
     &   lmd_cs=98.96D0,
     &   lmd_Cv=1.8D0,
     &   Ric=0.3D0,
     &   lmd_betaT=-0.2D0,
     &   lmd_epsilon=0.1D0,
     &   lmd_cekman=0.7D0,
     &   lmd_nu0c=0.1D0
     &                                                           )
      if (istr.eq.1) then
        imin=istr
      else
        imin=istr-1
      endif
      if (iend.eq.Lm) then
        imax=iend
      else
        imax=iend+1
      endif
      if (jstr.eq.1) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (jend.eq.Mm) then
        jmax=jend
      else
        jmax=jend+1
      endif
      Vtc=lmd_Cv*sqrt(-lmd_betaT)/( sqrt(lmd_cs*lmd_epsilon)
     &                                        *Ric*vonKar*vonKar )
!$acc kernels if(compute_on_device) default(present)
      do j=jmin,jmax
        do i=imin,imax
          ustar(i,j)=sqrt(sqrt( (0.5D0*(bustr(i,j)+bustr(i+1,j)))**2
     &                         +(0.5D0*(bvstr(i,j)+bvstr(i,j+1)))**2))
        enddo
      enddo
      do j=jmin,jmax
        do i=imin,imax
          wm(i,j)=vonKar*ustar(i,j)
          ws(i,j)=wm(i,j)
        enddo
      enddo
      ka=1
      ku=2
      do j=jmin,jmax
        do i=imin,imax
          my_hbbl(i,j)=z_r(i,j,N)-z_w(i,j,0)
          my_kbbl(i,j)=N
          Rib(i,j,1)=0.D0
        enddo
      enddo
      do k=2,N
        cff=g/rho0
        do j=jmin,jmax
          do i=imin,imax
            Ritop=-cff*(rho1(i,j,k)-rho1(i,j,1))
     &                           *(z_r(i,j,k)-z_r(i,j,1))
            dVsq=0.25D0*( (u(i  ,j,k,nstp)-u(i  ,j,1,nstp)+
     &                   u(i+1,j,k,nstp)-u(i+1,j,1,nstp))**2
     &                 +(v(i,j  ,k,nstp)-v(i,j  ,1,nstp)+
     &                   v(i,j+1,k,nstp)-v(i,j+1,1,nstp))**2)
            Vtsq=Vtc*(z_r(i,j,k)-z_r(i,j,1))*ws(i,j)
     &          *sqrt(max(0.D0,0.5D0*(bvf(i,j,k)+bvf(i,j,k-1))))
            Rib(i,j,k)=Ritop/(dVsq+Vtsq+eps)
          enddo
        enddo
      enddo
      do k=2,N
        do j=jmin,jmax
          do i=imin,imax
            if (my_kbbl(i,j).eq.N .and. Rib(i,j,k).gt.Ric) then
              zbl=z_r(i,j,k)-(z_r(i,j,k)-z_r(i,j,k-1))*
     &           (Ric-Rib(i,j,k))/(Rib(i,j,k-1)-Rib(i,j,k))
              my_hbbl(i,j)=zbl-z_w(i,j,0)
              my_kbbl(i,j)=k
            endif
          enddo
        enddo
      enddo
      do j=jmin,jmax
        do i=imin,imax
          hekman=lmd_cekman*ustar(i,j)/max(abs(f(i,j)),1.D-6)
          my_hbbl(i,j)=min(my_hbbl(i,j),hekman)
        enddo
      enddo
      if (istr.eq.1) then
        do j=jmin,jmax
          my_hbbl(Istr-1,j)=my_hbbl(Istr,j)
        enddo
      endif
      if (iend.eq.Lm) then
        do j=jmin,jmax
          my_hbbl(Iend+1,j)=my_hbbl(Iend,j)
        enddo
      endif
      if (jstr.eq.1) then
        do i=imin,imax
          my_hbbl(i,Jstr-1)=my_hbbl(i,Jstr)
        enddo
      endif
      if (jend.eq.Mm) then
        do i=imin,imax
          my_hbbl(i,Jend+1)=my_hbbl(i,Jend)
        enddo
      endif
      if (istr.eq.1.and.jstr.eq.1) then
        my_hbbl(Istr-1,Jstr-1)=my_hbbl(Istr,Jstr)
      endif
      if (istr.eq.1.and.jend.eq.Mm) then
        my_hbbl(Istr-1,Jend+1)=my_hbbl(Istr,Jend)
      endif
      if (iend.eq.Lm.and.jstr.eq.1) then
        my_hbbl(Iend+1,Jstr-1)=my_hbbl(Iend,Jstr)
      endif
      if (iend.eq.Lm.and.jend.eq.Mm) then
        my_hbbl(Iend+1,Jend+1)=my_hbbl(Iend,Jend)
      endif
      do j=Jstr,Jend+1
        do i=Istr,Iend+1
          wrk(i,j)=0.25D0*(my_hbbl(i,j)  +my_hbbl(i-1,j)
     &                  +my_hbbl(i,j-1)+my_hbbl(i-1,j-1))
     &                   *pmask2(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          cff=0.25D0*(pmask2(i,j)   +pmask2(i+1,j)
     &             +pmask2(i,j+1) +pmask2(i+1,j+1))
          my_hbbl(i,j)=(1.D0-cff)*my_hbbl(i,j)+
     &              0.25D0*(wrk(i,j)  +wrk(i+1,j)
     &                   +wrk(i,j+1)+wrk(i+1,j+1))
          my_hbbl(i,j)=my_hbbl(i,j)*rmask(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          hbbl(i,j)=min(my_hbbl(i,j),z_w(i,j,N)-z_w(i,j,0))
     &                                          *rmask(i,j)
          kbbl(i,j)=N
        enddo
      enddo
      do k=N,1,-1
        do j=Jstr,Jend
          do i=Istr,Iend
            if (z_r(i,j,k)-z_w(i,j,0).gt.hbbl(i,j)) then
              kbbl(i,j)=k
            endif
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          zbl=z_w(i,j,0)+hbbl(i,j)
          k=kbbl(i,j)
          if (zbl.lt.z_w(i,j,k-1)) k=k-1
          cff=1.D0/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(zbl-z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)-zbl)
          Kv_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dKv_bl=-cff*(Kv(i,j,k)-Kv(i,j,k-1))
          Gm1(i,j)=Kv_bl/(hbbl(i,j)*wm(i,j)+eps)
          dGm1dS(i,j)=min(0.D0,dKv_bl/(wm(i,j)+eps))
          Kt_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dKt_bl=-cff*(Kt(i,j,k)-Kt(i,j,k-1))
          Gt1(i,j)=Kt_bl/(hbbl(i,j)*ws(i,j)+eps)
          dGt1dS(i,j)=min(0.D0,dKt_bl/(ws(i,j)+eps))
          Ks_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dKs_bl=-cff*(Ks(i,j,k)-Ks(i,j,k-1))
          Gs1(i,j)=Ks_bl/(hbbl(i,j)*ws(i,j)+eps)
          dGs1dS(i,j)=min(0.D0,dKs_bl/(ws(i,j)+eps))
        enddo
      enddo
      do k=1,N-1
        do j=Jstr,Jend
          do i=Istr,Iend
            if (k.lt.kbbl(i,j)) then
              sig=min((z_w(i,j,k)-z_w(i,j,0))/(hbbl(i,j)+eps),1.D0)
              sig=sig*rmask(i,j)
              lmd_a1=sig-2.D0
              lmd_a2=3.D0-2.D0*sig
              lmd_a3=sig-1.D0
              Gm=lmd_a1+lmd_a2*Gm1(i,j)+lmd_a3*dGm1dS(i,j)
              Gt=lmd_a1+lmd_a2*Gt1(i,j)+lmd_a3*dGt1dS(i,j)
              Gs=lmd_a1+lmd_a2*Gs1(i,j)+lmd_a3*dGs1dS(i,j)
              Kv0=hbbl(i,j)*wm(i,j)*sig*(1.D0+sig*Gm)
              Kt0=hbbl(i,j)*ws(i,j)*sig*(1.D0+sig*Gt)
              Ks0=hbbl(i,j)*ws(i,j)*sig*(1.D0+sig*Gs)
              zsbl=z_w(i,j,N)-hbls(i,j,3-nstp)
              if (z_w(i,j,k).gt.zsbl) then
                Kv0=max(Kv0,Kv(i,j,k))
                Kt0=max(Kt0,Kt(i,j,k))
                Ks0=max(Ks0,Ks(i,j,k))
              endif
              Kv(i,j,k)=Kv0
              Kt(i,j,k)=Kt0
              Ks(i,j,k)=Ks0
            else
              if (bvf(i,j,k).lt.0.D0) then
             zsbl=z_w(i,j,N)-hbls(i,j,3-nstp)
                if (z_w(i,j,k).lt.zsbl) then
                  Kv(i,j,k)=Kv(i,j,k)+lmd_nu0c
                  Kt(i,j,k)=Kt(i,j,k)+lmd_nu0c
                  Ks(i,j,k)=Ks(i,j,k)+lmd_nu0c
                endif
              endif
            endif
          enddo
        enddo
      enddo
      if (istr.eq.1) then
        do j=jstr,jend
          hbbl(istr-1,j)=hbbl(istr,j)
        enddo
      endif
      if (iend.eq.Lm) then
        do j=jstr,jend
          hbbl(iend+1,j)=hbbl(iend,j)
        enddo
      endif
      if (jstr.eq.1) then
        do i=istr,iend
          hbbl(i,jstr-1)=hbbl(i,jstr)
        enddo
      endif
      if (jend.eq.Mm) then
        do i=istr,iend
          hbbl(i,jend+1)=hbbl(i,jend)
        enddo
      endif
      if (istr.eq.1 .and. jstr.eq.1) then
        hbbl(istr-1,jstr-1)=hbbl(istr,jstr)
      endif
      if (istr.eq.1 .and. jend.eq.Mm) then
        hbbl(istr-1,jend+1)=hbbl(istr,jend)
      endif
      if (iend.eq.Lm .and. jstr.eq.1) then
        hbbl(iend+1,jstr-1)=hbbl(iend,jstr)
      endif
      if (iend.eq.Lm .and. jend.eq.Mm) then
        hbbl(iend+1,jend+1)=hbbl(iend,jend)
      endif
!$acc end kernels
      return
      end
