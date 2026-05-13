      subroutine init_arrays (tile)
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
      integer*4 tile, trd
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
      call init_arrays_tile (Istr,Iend,Jstr,Jend)
      return
      end
      subroutine init_arrays_tile (Istr,Iend,Jstr,Jend)
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
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diags_vrt, filetype_diags_vrt_avg
     &       ,filetype_diags_ek, filetype_diags_ek_avg
     &       ,filetype_diags_pv, filetype_diags_pv_avg
     &       ,filetype_diags_eddy_avg
     &       ,filetype_surf, filetype_surf_avg
     &       ,filetype_diabio, filetype_diabio_avg
     &       ,filetype_abl, filetype_abl_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diags_vrt=7, filetype_diags_vrt_avg=8,
     &           filetype_diags_ek=9, filetype_diags_ek_avg=10,
     &           filetype_diags_pv=11, filetype_diags_pv_avg=12,
     &           filetype_diags_eddy_avg=17,
     &           filetype_surf=13, filetype_surf_avg=14,
     &           filetype_diabio=15,filetype_diabio_avg=16,
     &           filetype_abl=18, filetype_abl_avg=19)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV
      parameter (indxU=6, indxV=7)
      integer*4 indxT
      parameter (indxT=indxV+1)
      integer*4 indxS
      parameter (indxS=indxV+ntrc_temp+1)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxV+ntrc_temp+ntrc_salt+ntrc_mld+ntrc_pas+
     &           ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxV+ntrc_temp+ntrc_salt
     &                      +ntrc_mld+ntrc_pas+ntrc_bio
     &                      +ntrc_sed+ntrc_substot
     &           +ntrc_diats+ntrc_diauv+ntrc_diavrt+ntrc_diaek
     &           +ntrc_diapv+ntrc_diaeddy+ntrc_surf+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxAks
      parameter (indxAks=indxAkv+ntrc_temp+4)
      integer*4 indxHbl
      parameter (indxHbl=indxAkv+ntrc_temp+5)
      integer*4 indxHbbl
      parameter (indxHbbl=indxAkv+ntrc_temp+6)
      integer*4 indxSSH
      parameter (indxSSH=indxAkv+ntrc_temp+12)
      integer*4 indxbvf
      parameter (indxbvf=indxSSH+1)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+2, indxSVSTR=indxSSH+3)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+4)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSUSTR+5)
      integer*4 indxSwflx
      parameter (indxSwflx=indxShflx+1, indxShflx_rsw=indxShflx+2)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWSPD,indxTAIR,indxRHUM,indxRADLW,indxRADSW,
     &        indxPRATE,indxUWND,indxVWND,indxPATM
      parameter (indxWSPD=indxSST+3,  indxTAIR=indxSST+4,
     &           indxRHUM=indxSST+5,  indxRADLW=indxSST+6,
     &           indxRADSW=indxSST+7, indxPRATE=indxSST+8,
     &           indxUWND=indxSST+9,  indxVWND=indxSST+10,
     &           indxPATM=indxSST+11)
      integer*4 indxShflx_rlw,indxShflx_lat,indxShflx_sen
      parameter (indxShflx_rlw=indxSST+12,
     &           indxShflx_lat=indxSST+13, indxShflx_sen=indxSST+14)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+23)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+24)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+25)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+26)
      integer*4 indxBustr, indxBvstr
      parameter (indxBustr=indxSUSTR+27,  indxBvstr=indxBustr+1)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB,indxWED,indxWER
      parameter (indxWWA=indxSUSTR+42, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=11, b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms
     &     , ncidqbar, ncidbtf
     &     , ntsrf,  ntssh,  ntsst, ntsss, ntuclm
     &     , ntbulk, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
     &       , ntbtf(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4 rstHbl
      integer*4 rstHbbl
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw, hisBhflx, hisBwflx
     &      , hisUnbq, hisVnbq, hisWnbq, hisRnbq, hisCnbq
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
     &      , hisbvf
     &      , hisShflx_rlw
     &      , hisShflx_lat,   hisShflx_sen
      integer*4 hisT(NT)
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw, avgBhflx, avgBwflx
     &      , avgU,   avgV,   avgR,    avgHbl, avgHbbl
     &      , avgO,   avgW,   avgVisc, avgDiff
     &      , avgAkv, avgAkt, avgAks
     &      , avgbvf
      integer*4 avgT(NT)
      integer*4 avgShflx_rlw
     &      , avgShflx_lat,   avgShflx_sen
      logical wrthis(1000+NT)
     &      , wrtavg(1000+NT)
      common/incscrum/
     &     ncidfrc, ncidbulk,ncidclm, ncidqbar, ncidbtf
     &     , ntsms, ntsrf, ntssh, ntsst
     &     , ntuclm, ntsss, ntbulk, ntqbar, ntww
     &     ,  nttclm, ntstf, nttsrc, ntbtf
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
     & ,   rstT
     &      , rstHbl
     &      , rstHbbl
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisBustr, hisBvstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisBhflx, hisBwflx
     &      , hisUnbq, hisVnbq, hisWnbq, hisRnbq, hisCnbq
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , hisbvf
     &      , hisShflx_rlw
     &      , hisShflx_lat, hisShflx_sen
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , avgBustr, avgBvstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgBhflx, avgBwflx
     &      , avgU,    avgV
     &      ,     avgT
     &      ,     avgR
     &      , avgO,    avgW,     avgVisc,  avgDiff
     &      , avgAkv,  avgAkt,   avgAks
     &      , avgHbl,  avgHbbl
     &      , avgbvf
     &      , avgShflx_rlw
     &      , avgShflx_lat, avgShflx_sen
     &      , wrthis
     &      , wrtavg
      character*80 date_str, title
      character*80 origin_date, start_date_run, xios_origin_date
      integer*4      start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
      REAL(kind=8) :: origin_date_in_sec, xios_origin_date_in_sec
      character*180 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname
     &                                ,  avgname
     &                                ,   clmname
     &                                ,   bry_file
      character*75  vname(20, 1000)
      common /cncscrum/   date_str,   title
     &         ,   origin_date, start_date_run
     &         ,   xios_origin_date
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &         ,   btfname, origin_date_in_sec
     &         ,   xios_origin_date_in_sec
     &         ,   start_day, start_month, start_year
     &         ,   start_hour, start_minute, start_second
     &         ,   origin_day, origin_month, origin_year
     &         ,   origin_hour, origin_minute, origin_second
     &                                ,  avgname
     &                                ,   clmname
     &                                ,   bry_file
     &                                ,   vname
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
      real zeta_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real ubar_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real vbar_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_zeta/zeta_avg
     &       /avg_ubar/ubar_avg
     &       /avg_vbar/vbar_avg
      real bostr_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_bostr/bostr_avg
      real bustr_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_bustr/bustr_avg
      real bvstr_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_bvstr/bvstr_avg
      real wstr_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_wstr/wstr_avg
      real sustr_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_sustr/sustr_avg
      real svstr_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_svstr/svstr_avg
      real srflx_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_srflx/srflx_avg
      real u_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real v_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real t_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,N,NT)
      real rho_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      real bvf_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real omega_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real w_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,N)
      common /avg_u/u_avg /avg_v/v_avg /avg_t/t_avg
     &       /avg_rho/rho_avg /avg_omega/omega_avg
     &       /avg_bvf/bvf_avg
     &       /avg_w/w_avg
      real stflx_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /avg_stflx/stflx_avg
      real btflx_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,NT)
      common /avg_btflx/btflx_avg
      real hbl_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_hbl/hbl_avg
      real hbbl_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_hbbl/hbbl_avg
      real shflx_rsw_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real shflx_rlw_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real shflx_lat_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      real shflx_sen_avg(0:Lm+1+padd_X,0:Mm+1+padd_E)
      common /avg_shflx_rsw/shflx_rsw_avg
      common /avg_shflx_rlw/shflx_rlw_avg
      common /avg_shflx_lat/shflx_lat_avg
      common /avg_shflx_sen/shflx_sen_avg
      real Akv_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N)
      real Akt_avg(0:Lm+1+padd_X,0:Mm+1+padd_E,0:N,2)
      common /avg_Akv/Akv_avg /avg_Akt/Akt_avg
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
      logical got_tbry(NT)
      common /bry_logical/ got_tbry
      real bry_time(2)
      REAL(kind=8) :: bry_origin_date_in_sec
      common /bry_indices_array/ bry_time,
     &        bry_origin_date_in_sec
      real bry_cycle
      common /bry_indices_real/ bry_cycle
      integer*4 bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, ntbry
      common /bry_indices_integer/ bry_id, bry_time_id, bry_ncycle,
     &                             bry_rec, itbry, ntbry
      integer*4 zetabry_west_id
      common /zeta_west_id/ zetabry_west_id
      integer*4 ubarbry_west_id, vbarbry_west_id
      common /ubar_west_id/ ubarbry_west_id, vbarbry_west_id
      integer*4 ubry_west_id, vbry_west_id
      common /u_west_id/ ubry_west_id, vbry_west_id
      integer*4 tbry_west_id(NT)
      common /t_west_id/ tbry_west_id
      integer*4 zetabry_east_id
      common /zeta_east_id/ zetabry_east_id
      integer*4 ubarbry_east_id, vbarbry_east_id
      common /ubar_east_id/ ubarbry_east_id, vbarbry_east_id
      integer*4 ubry_east_id, vbry_east_id
      common /u_east_id/ ubry_east_id, vbry_east_id
      integer*4 tbry_east_id(NT)
      common /t_east_id/ tbry_east_id
      integer*4 zetabry_south_id
      common /zeta_south_id/ zetabry_south_id
      integer*4 ubarbry_south_id, vbarbry_south_id
      common /ubar_south_id/ ubarbry_south_id, vbarbry_south_id
      integer*4 ubry_south_id, vbry_south_id
      common /u_south_id/ ubry_south_id, vbry_south_id
      integer*4 tbry_south_id(NT)
      common /t_south_id/ tbry_south_id
      integer*4 zetabry_north_id
      common /zeta_north_id/ zetabry_north_id
      integer*4 ubarbry_north_id, vbarbry_north_id
      common /ubar_north_id/ ubarbry_north_id, vbarbry_north_id
      integer*4 ubry_north_id, vbry_north_id
      common /u_north_id/ ubry_north_id, vbry_north_id
      integer*4 tbry_north_id(NT)
      common /t_north_id/ tbry_north_id
      real zetabry_west(0:Mm+1+padd_E),
     &    zetabry_west_dt(0:Mm+1+padd_E,2)
      common /bry_zeta_west/ zetabry_west, zetabry_west_dt
      real ubarbry_west(0:Mm+1+padd_E),
     &    ubarbry_west_dt(0:Mm+1+padd_E,2)
     &    ,vbarbry_west(0:Mm+1+padd_E),
     &    vbarbry_west_dt(0:Mm+1+padd_E,2)
      common /bry_ubar_west/ ubarbry_west, ubarbry_west_dt,
     &                       vbarbry_west, vbarbry_west_dt
      real ubry_west(0:Mm+1+padd_E,N),
     &    ubry_west_dt(0:Mm+1+padd_E,N,2)
     &    ,vbry_west(0:Mm+1+padd_E,N),
     &    vbry_west_dt(0:Mm+1+padd_E,N,2)
      common /bry_u_west/ ubry_west, ubry_west_dt,
     &                    vbry_west, vbry_west_dt
      real tbry_west(0:Mm+1+padd_E,N,NT),
     &     tbry_west_dt(0:Mm+1+padd_E,N,2,NT)
      common /bry_t_west/ tbry_west, tbry_west_dt
      real zetabry_east(0:Mm+1+padd_E),
     &    zetabry_east_dt(0:Mm+1+padd_E,2)
      common /bry_zeta_east/ zetabry_east, zetabry_east_dt
      real ubarbry_east(0:Mm+1+padd_E),
     &    ubarbry_east_dt(0:Mm+1+padd_E,2)
     &    ,vbarbry_east(0:Mm+1+padd_E),
     &    vbarbry_east_dt(0:Mm+1+padd_E,2)
      common /bry_ubar_east/ ubarbry_east, ubarbry_east_dt,
     &                       vbarbry_east, vbarbry_east_dt
      real ubry_east(0:Mm+1+padd_E,N),
     &    ubry_east_dt(0:Mm+1+padd_E,N,2)
     &    ,vbry_east(0:Mm+1+padd_E,N),
     &    vbry_east_dt(0:Mm+1+padd_E,N,2)
      common /bry_u_east/ ubry_east, ubry_east_dt,
     &                    vbry_east, vbry_east_dt
      real tbry_east(0:Mm+1+padd_E,N,NT),
     &    tbry_east_dt(0:Mm+1+padd_E,N,2,NT)
      common /bry_t_east/ tbry_east, tbry_east_dt
      real zetabry_south(0:Lm+1+padd_X),
     &    zetabry_south_dt(0:Lm+1+padd_X,2)
      common /bry_zeta_south/ zetabry_south, zetabry_south_dt
      real ubarbry_south(0:Lm+1+padd_X),
     &    ubarbry_south_dt(0:Lm+1+padd_X,2)
     &    ,vbarbry_south(0:Lm+1+padd_X),
     &    vbarbry_south_dt(0:Lm+1+padd_X,2)
      common /bry_ubar_south/ ubarbry_south, ubarbry_south_dt,
     &                        vbarbry_south, vbarbry_south_dt
      real ubry_south(0:Lm+1+padd_X,N),
     &    ubry_south_dt(0:Lm+1+padd_X,N,2)
     &    ,vbry_south(0:Lm+1+padd_X,N),
     &    vbry_south_dt(0:Lm+1+padd_X,N,2)
      common /bry_u_south/ ubry_south, ubry_south_dt,
     &                     vbry_south, vbry_south_dt
      real tbry_south(0:Lm+1+padd_X,N,NT),
     &    tbry_south_dt(0:Lm+1+padd_X,N,2,NT)
      common /bry_t_south/ tbry_south, tbry_south_dt
      real zetabry_north(0:Lm+1+padd_X),
     &    zetabry_north_dt(0:Lm+1+padd_X,2)
      common /bry_zeta_north/ zetabry_north, zetabry_north_dt
      real ubarbry_north(0:Lm+1+padd_X),
     &    ubarbry_north_dt(0:Lm+1+padd_X,2)
     &    ,vbarbry_north(0:Lm+1+padd_X),
     &    vbarbry_north_dt(0:Lm+1+padd_X,2)
      common /bry_ubar_north/ ubarbry_north, ubarbry_north_dt,
     &                        vbarbry_north, vbarbry_north_dt
      real ubry_north(0:Lm+1+padd_X,N),
     &    ubry_north_dt(0:Lm+1+padd_X,N,2)
     &    ,vbry_north(0:Lm+1+padd_X,N),
     &    vbry_north_dt(0:Lm+1+padd_X,N,2)
      common /bry_u_north/ ubry_north, ubry_north_dt,
     &                     vbry_north, vbry_north_dt
      real tbry_north(0:Lm+1+padd_X,N,NT),
     &    tbry_north_dt(0:Lm+1+padd_X,N,2,NT)
      common /bry_t_north/ tbry_north, tbry_north_dt
      integer*4 Istr,Iend,Jstr,Jend, i,j
     &       ,k,itrc
     &       ,itide
      real init
      parameter (init=0.D0)
      integer*4 IstrR,IendR,JstrR,JendR
      if (Istr.eq.1) then
        IstrR=Istr-1
      else
        IstrR=Istr
      endif
      if (Iend.eq.Lm) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (Jstr.eq.1) then
        JstrR=Jstr-1
      else
        JstrR=Jstr
      endif
      if (Jend.eq.Mm) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      do j=JstrR,JendR
        do i=IstrR,IendR
          zeta(i,j,1)=0.D0
          zeta(i,j,2)=init
          zeta(i,j,3)=init
          zeta(i,j,4)=init
          ubar(i,j,1)=init
          ubar(i,j,2)=init
          ubar(i,j,3)=init
          ubar(i,j,4)=init
          vbar(i,j,1)=init
          vbar(i,j,2)=init
          vbar(i,j,3)=init
          vbar(i,j,4)=init
          zeta_avg(i,j)=init
          ubar_avg(i,j)=init
          vbar_avg(i,j)=init
          rufrc(i,j)=init
          rufrc(i,j)=init
          rufrc_bak(i,j,1)=init
          rvfrc_bak(i,j,1)=init
          rufrc_bak(i,j,2)=init
          rvfrc_bak(i,j,2)=init
          Zt_avg1(i,j)=0.D0
          DU_avg1(i,j,1)=0.D0
          DV_avg1(i,j,1)=0.D0
          DU_avg1(i,j,2)=0.D0
          DV_avg1(i,j,2)=0.D0
          DU_avg2(i,j)=0.D0
          DV_avg2(i,j)=0.D0
        enddo
      enddo
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            u(i,j,k,1)=init
            u(i,j,k,2)=init
            v(i,j,k,1)=init
            v(i,j,k,2)=init
            rho(i,j,k) =init
            rho_avg(i,j,k)=init
            u_avg(i,j,k)=init
            v_avg(i,j,k)=init
            w_avg(i,j,k)=init
          enddo
        enddo
      enddo
      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
             We(i,j,k)=init
            omega_avg(i,j,k)=init
            bvf_avg(i,j,k)=init
          enddo
        enddo
      enddo
      do itrc=1,NT
        got_tini(itrc)=.false.
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              t(i,j,k,1,itrc)=init
              t(i,j,k,2,itrc)=init
              t_avg(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          sustr(i,j)=init
          svstr(i,j)=init
          sustrg(i,j,1)=init
          svstrg(i,j,1)=init
          sustrg(i,j,2)=init
          svstrg(i,j,2)=init
          bustr(i,j)=init
          bvstr(i,j)=init
          bostr_avg(i,j)=init
          wstr_avg(i,j)=init
          sustr_avg(i,j)=init
          svstr_avg(i,j)=init
          bustrg(i,j,1)=init
          bvstrg(i,j,1)=init
          bustrg(i,j,2)=init
          bvstrg(i,j,2)=init
        enddo
      enddo
      do itrc=1,NT
        do j=JstrR,JendR
          do i=IstrR,IendR
            stflx(i,j,itrc)=init
            stflx_avg(i,j,itrc)=init
            stflxg(i,j,1,itrc)=init
            stflxg(i,j,2,itrc)=init
            btflx(i,j,itrc)=init
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          tair (i,j)=init
          rhum (i,j)=init
          prate (i,j)=init
          radlw (i,j)=init
          radsw (i,j)=init
          wspd(i,j)=init
          uwnd(i,j)=init
          vwnd(i,j)=init
          tairg(i,j,1)=init
          rhumg (i,j,1)=init
          prateg (i,j,1)=init
          radlwg (i,j,1)=init
          radswg (i,j,1)=init
          tairg(i,j,2)=init
          rhumg (i,j,2)=init
          prateg (i,j,2)=init
          radlwg (i,j,2)=init
          radswg (i,j,2)=init
          uwndg(i,j,1)=init
          vwndg(i,j,1)=init
          uwndg(i,j,2)=init
          vwndg(i,j,2)=init
          shflx_rsw(i,j)=init
          shflx_rlw(i,j)=init
          shflx_lat(i,j)=init
          shflx_sen(i,j)=init
          shflx_rsw_avg(i,j)=init
          shflx_rlw_avg(i,j)=init
          shflx_lat_avg(i,j)=init
          shflx_sen_avg(i,j)=init
          srflx(i,j)=init
          srflx_avg(i,j)=init
          srflxg(i,j,1)=init
          srflxg(i,j,2)=init
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          Znudgcof(i,j)=init
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          ssh(i,j)=init
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          sshg(i,j,1)=init
          sshg(i,j,2)=init
        enddo
      enddo
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              Tnudgcof(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              tclm(i,j,k,itrc)=init
            enddo
          enddo
        enddo
      enddo
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              tclima(i,j,k,1,itrc)=init
              tclima(i,j,k,2,itrc)=init
            enddo
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          ubclm(i,j)=init
          vbclm(i,j)=init
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          ubclima(i,j,1)=init
          ubclima(i,j,2)=init
          vbclima(i,j,1)=init
          vbclima(i,j,2)=init
        enddo
      enddo
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            uclima(i,j,k,1)=init
            uclima(i,j,k,2)=init
            vclima(i,j,k,1)=init
            vclima(i,j,k,2)=init
          enddo
        enddo
      enddo
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            uclm(i,j,k)=init
            vclm(i,j,k)=init
          enddo
        enddo
      enddo
      if (istr.eq.1) then
        do j=JstrR,JendR
          zetabry_west(j)=init
        enddo
        do j=JstrR,JendR
          ubarbry_west(j)=init
          vbarbry_west(j)=init
        enddo
        do k=1,N
          do j=JstrR,JendR
            ubry_west(j,k)=init
            vbry_west(j,k)=init
            do itrc=1,NT
              tbry_west(j,k,itrc)=init
            enddo
          enddo
        enddo
      endif
      if (iend.eq.Lm) then
        do j=JstrR,JendR
          zetabry_east(j)=init
        enddo
        do j=JstrR,JendR
          ubarbry_east(j)=init
          vbarbry_east(j)=init
        enddo
        do k=1,N
          do j=JstrR,JendR
            ubry_east(j,k)=init
            vbry_east(j,k)=init
            do itrc=1,NT
              tbry_east(j,k,itrc)=init
            enddo
          enddo
        enddo
      endif
      if (jstr.eq.1) then
        do i=IstrR,IendR
          zetabry_south(i)=init
        enddo
        do i=IstrR,IendR
          ubarbry_south(i)=init
          vbarbry_south(i)=init
        enddo
        do k=1,N
          do i=IstrR,IendR
            ubry_south(i,k)=init
            vbry_south(i,k)=init
            do itrc=1,NT
              tbry_south(i,k,itrc)=init
            enddo
          enddo
        enddo
      endif
      if (jend.eq.Mm) then
        do i=IstrR,IendR
          zetabry_north(i)=init
        enddo
        do i=IstrR,IendR
          ubarbry_north(i)=init
          vbarbry_north(i)=init
        enddo
        do k=1,N
          do i=IstrR,IendR
            ubry_north(i,k)=init
            vbry_north(i,k)=init
            do itrc=1,NT
              tbry_north(i,k,itrc)=init
            enddo
          enddo
        enddo
      endif
      do itide=1,Ntides
        Tperiod(itide)=init
        do j=JstrR,JendR
          do i=IstrR,IendR
            SSH_Tamp(i,j,itide)=init
            SSH_Tphase(i,j,itide)=init
            UV_Tangle(i,j,itide)=init
            UV_Tphase(i,j,itide)=init
            UV_Tmajor(i,j,itide)=init
            UV_Tminor(i,j,itide)=init
          enddo
        enddo
      enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
            visc2_r(i,j)=visc2
            visc2_p(i,j)=visc2
            visc2_sponge_r(i,j)=init
            visc2_sponge_p(i,j)=init
          enddo
        enddo
        do itrc=1,NT
          do j=JstrR,JendR
            do i=IstrR,IendR
              diff2(i,j,itrc)=tnu2(itrc)
            enddo
          enddo
        enddo
          do j=JstrR,JendR
            do i=IstrR,IendR
              diff2_sponge(i,j)=init
            enddo
          enddo
      do k=0,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            Akv(i,j,k)=init
            bvf(i,j,k)=init
          enddo
        enddo
        do j=JstrR,JendR
          do i=IstrR,IendR
              Akt(i,j,k,itemp)=init
              Akt(i,j,k,isalt)=init
          enddo
        enddo
      enddo
      do k=1,N
        do j=JstrR,JendR
          do i=IstrR,IendR
            ghats(i,j,k)=init
          enddo
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          hbls(i,j,1)=init
          hbls(i,j,2)=init
          kbl(i,j)=init
          hbl_avg(i,j)=init
        enddo
      enddo
      do j=JstrR,JendR
        do i=IstrR,IendR
          hbbl(i,j)=init
          hbbl_avg(i,j)=init
        enddo
      enddo
      return
      end
