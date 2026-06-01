# 0 "tools_fort.F"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/home/claudiaa/miniconda3/envs/croco_pyenv/x86_64-conda-linux-gnu/sysroot/usr/include/stdc-predef.h" 1 3 4
# 0 "<command-line>" 2
# 1 "tools_fort.F"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROMS ROUTINES
!!
!! copied from actual ROMS scripts
!!
!! compile with:
!! "cpp R_tools_fort.F R_tools_fort.f"
!! "f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m R_tools_fort R_tools_fort.f" for python use
!!
!! print R_tools_fort.rho_eos.__doc__
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 1 "sigma_to_z_intr.F" 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     & z_lev, var_zlv, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev.
!


      implicit none

      integer Lm,Mm,N, nz, imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(imin:Lm+1,jmin:Mm+1,kmin:N),
     & z_r(0:Lm+1,0:Mm+1,N), rmask(0:Lm+1,0:Mm+1),
     & z_w(0:Lm+1,0:Mm+1,0:N), z_lev(imin:Lm+1,jmin:Mm+1,nz),
     & FillValue, var_zlv(imin:Lm+1,jmin:Mm+1,nz),
     & zz(0:Lm+1,0:N+1), dpth

     & , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max( trd *chunk_size -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv
# 54 "sigma_to_z_intr.F"
      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,j,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,j,0)
              zz(i,N+1)=z_w(i,j,N)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
              zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     & +z_r(i,j-1,k)+z_r(i-1,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     & +z_w(i,j-1,0)+z_w(i-1,j-1,0))

              zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     & +z_w(i,j-1,N)+z_w(i-1,j-1,N))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,j,k)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     & +z_w(i,j-1,k)+z_w(i-1,j-1,k))
              enddo
            enddo
          endif
        endif

        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1 !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff.gt.0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     & *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo

        do m=1,nz


          if (kmin.eq.0) then !
            do i=imin,imax !
              dpth=zz(i,N)-zz(i,0)
              if (rmask(i,j).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)
              if (rmask(i,j).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N+1)).gt.0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2 !<-- below bottom
              elseif (dpth*(zz(i,1)-z_lev(i,j,m)).gt.0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(i,k+1)-z_lev(i,j,m))*(z_lev(i,j,m)-zz(i,k))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(i,j,m)=0. !<-- masked out
            elseif (km(i).eq.-2) then
              var_zlv(i,j,m)=FillValue !<-- below bottom
            elseif (km(i).eq.N+2) then
              var_zlv(i,j,m)=-FillValue !<-- above surface
            elseif (km(i).eq.N) then
              var_zlv(i,j,m)=var(i,j,N) !-> R-point, above z_r(N)

     & +FC(i,N)*(z_lev(i,j,m)-zz(i,N))




            elseif (km(i).eq.kmin-1) then !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin) ! but above bottom

     & -FC(i,kmin)*(zz(i,kmin)-z_lev(i,j,m))




            else
              k=km(i)
              !write(*,*) k,km

              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=z_lev(i,j,m)-zz(i,k)
              q=zz(i,k+1)-z_lev(i,j,m)

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     & -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     & +p*FC(i,k+1) -q*FC(i,k) )
     & )







            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m
# 255 "sigma_to_z_intr.F"
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
# 15 "tools_fort.F" 2
# 1 "sigma_to_z_intr_bot.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bot (Lm,Mm,N,nz,z_r,z_w,rmask,var,
     & z_lev, var_zlv, below, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev.
!


      implicit none

      integer Lm,Mm,N, nz, imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(kmin:N,jmin:Mm+1,imin:Lm+1),
     & z_r(N,0:Mm+1,0:Lm+1), rmask(0:Mm+1,0:Lm+1),
     & z_w(0:N,0:Mm+1,0:Lm+1), z_lev(nz,jmin:Mm+1,imin:Lm+1),
     & FillValue, var_zlv(nz,jmin:Mm+1,imin:Lm+1),
     & zz(0:N+1,0:Lm+1), dpth, below

     & , dz(kmin-1:N,0:Lm+1), FC(kmin-1:N,0:Lm+1), p,q,cff

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max( trd *chunk_size -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, below, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv
# 53 "sigma_to_z_intr_bot.F"
      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=z_r(k,j,i)
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=z_w(0,j,i)
              zz(N+1,i)=z_w(N,j,i)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_r(k,j,i)+z_r(k,j,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i-1)+z_w(0,j,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i-1)+z_w(N,j,i))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5*(z_r(k,j,i)+z_r(k,j-1,i))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i)+z_w(0,j-1,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i)+z_w(N,j-1,i))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_r(k,j,i)+z_r(k,j,i-1)
     & +z_r(k,j-1,i)+z_r(k,j-1,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.25D0*( z_w(0,j,i)+z_w(0,j,i-1)
     & +z_w(0,j-1,i)+z_w(0,j-1,i-1))

              zz(N+1,i)=0.25D0*( z_w(N,j,i)+z_w(N,j,i-1)
     & +z_w(N,j-1,i)+z_w(N,j-1,i-1))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=z_w(k,j,i)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_w(k,j,i)+z_w(k,j,i-1))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5*(z_w(k,j,i)+z_w(k,j-1,i))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_w(k,j,i)+z_w(k,j,i-1)
     & +z_w(k,j-1,i)+z_w(k,j-1,i-1))
              enddo
            enddo
          endif
        endif

        do k=kmin,N-1
          do i=imin,imax
            dz(k,i)=zz(k+1,i)-zz(k,i)
            FC(k,i)=var(k+1,j,i)-var(k,j,i)
          enddo
        enddo
        do i=imin,imax
          dz(kmin-1,i)=dz(kmin,i)
          FC(kmin-1,i)=FC(kmin,i)

          dz(N,i)=dz(N-1,i)
          FC(N,i)=FC(N-1,i)
        enddo
        do k=N,kmin,-1 !--> irreversible
          do i=imin,imax
            cff=FC(k,i)*FC(k-1,i)
            if (cff.gt.0.D0) then
              FC(k,i)=cff*(dz(k,i)+dz(k-1,i))/( (FC(k,i)+FC(k-1,i))
     & *dz(k,i)*dz(k-1,i) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo

        do m=1,nz


          if (kmin.eq.0) then !
            do i=imin,imax !
              dpth=zz(N,i)-zz(0,i)
              if (rmask(j,i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(N+1,i)-zz(0,i)
              if (rmask(j,i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N+1,i)).gt.0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(0,i)-below-z_lev(m,j,i)).gt.0.) then
                km(i)=-3 !<-- below bottom
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2 !<-- below bottom
              elseif (dpth*(zz(1,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(k+1,i)-z_lev(m,j,i))*(z_lev(m,j,i)-zz(k,i))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(m,j,i)=FillValue !<-- masked out
            elseif (km(i).eq.-2) then
# 213 "sigma_to_z_intr_bot.F"
              var_zlv(m,j,i)=FillValue !<-- below bottom

            elseif (km(i).eq.N+2) then
# 225 "sigma_to_z_intr_bot.F"
              var_zlv(m,j,i)=-FillValue !<-- above surface

            elseif (km(i).eq.N) then
              var_zlv(m,j,i)=var(N,j,i) !-> R-point, above z_r(N)

     & +FC(N,i)*(z_lev(m,j,i)-zz(N,i))




            elseif (km(i).eq.kmin-1) then !-> R-point below z_r(1),
              var_zlv(m,j,i)=var(kmin,j,i) ! but above bottom

     & -FC(kmin,i)*(zz(kmin,i)-z_lev(m,j,i))




            else
              k=km(i)
              !write(*,*) k,km

              cff=1.D0/(zz(k+1,i)-zz(k,i))
              p=z_lev(m,j,i)-zz(k,i)
              q=zz(k+1,i)-z_lev(m,j,i)

              var_zlv(m,j,i)=cff*( q*var(k,j,i) + p*var(k+1,j,i)
     & -cff*p*q*( cff*(q-p)*(var(k+1,j,i)-var(k,j,i))
     & +p*FC(k+1,i) -q*FC(k,i) )
     & )







            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m
# 276 "sigma_to_z_intr_bot.F"
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
# 16 "tools_fort.F" 2
# 1 "sigma_to_z_intr_bot_2d.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bot_2d (Lm,N, nz, z_r, z_w, rmask, var,
     & z_lev, var_zlv, below,imin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev.
!


      implicit none

      integer Lm,Mm,N, nz, imin,imax, kmin, i,k,m

      integer km(0:Lm+1)

      real*8 var(kmin:N,imin:Lm+1),
     & z_r(N,0:Lm+1), rmask(0:Lm+1),
     & z_w(0:N,0:Lm+1), z_lev(nz,imin:Lm+1),
     & FillValue, var_zlv(nz,imin:Lm+1),
     & zz(0:N+1,0:Lm+1), dpth, below

     & , dz(kmin-1:N,0:Lm+1), FC(kmin-1:N,0:Lm+1), p,q,cff

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, below, imin,kmin, FillValue
Cf2py intent(out) var_zlv



        if (kmin.eq.1) then
          if (imin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=z_r(k,i)
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=z_w(0,i)
              zz(N+1,i)=z_w(N,i)
            enddo
          elseif (imin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_r(k,i)+z_r(k,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,i-1)+z_w(0,i))
              zz(N+1,i)=0.5D0*(z_w(N,i-1)+z_w(N,i))
            enddo
          endif
        else
          if (imin.eq.0 ) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=z_w(k,i)
              enddo
            enddo
          elseif (imin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_w(k,i)+z_w(k,i-1))
              enddo
            enddo
          endif
        endif

        do k=kmin,N-1
          do i=imin,imax
            dz(k,i)=zz(k+1,i)-zz(k,i)
            FC(k,i)=var(k+1,i)-var(k,i)
          enddo
        enddo
        do i=imin,imax
          dz(kmin-1,i)=dz(kmin,i)
          FC(kmin-1,i)=FC(kmin,i)

          dz(N,i)=dz(N-1,i)
          FC(N,i)=FC(N-1,i)
        enddo
        do k=N,kmin,-1 !--> irreversible
          do i=imin,imax
            cff=FC(k,i)*FC(k-1,i)
            if (cff.gt.0.D0) then
              FC(k,i)=cff*(dz(k,i)+dz(k-1,i))/( (FC(k,i)+FC(k-1,i))
     & *dz(k,i)*dz(k-1,i) )
            else
              FC(k,i)=0.D0
            endif
          enddo
        enddo

        do m=1,nz


          if (kmin.eq.0) then !
            do i=imin,imax !
              dpth=zz(N,i)-zz(0,i)
              if (rmask(i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,i)-zz(N,i)).gt.0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(0,i)-z_lev(m,i)).gt.0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(N+1,i)-zz(0,i)
              if (rmask(i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,i)-zz(N+1,i)).gt.0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(z_lev(m,i)-zz(N,i)).gt.0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(0,i)-below-z_lev(m,i)).gt.0.) then
                km(i)=-3 !<-- below bottom
              elseif (dpth*(zz(0,i)-z_lev(m,i)).gt.0.) then
                km(i)=-2 !<-- below bottom but close
              elseif (dpth*(zz(1,i)-z_lev(m,i)).gt.0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(k+1,i)-z_lev(m,i))*(z_lev(m,i)-zz(k,i))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(m,i)=FillValue !<-- masked out
            elseif (km(i).eq.-2) then
# 159 "sigma_to_z_intr_bot_2d.F"
              var_zlv(m,i)=FillValue !<-- below bottom

            elseif (km(i).eq.N+2) then
# 171 "sigma_to_z_intr_bot_2d.F"
              var_zlv(m,i)=-FillValue !<-- above surface

            elseif (km(i).eq.N) then
              var_zlv(m,i)=var(N,i) !-> R-point, above z_r(N)

     & +FC(N,i)*(z_lev(m,i)-zz(N,i))




            elseif (km(i).eq.kmin-1) then !-> R-point below z_r(1),
              var_zlv(m,i)=var(kmin,i) ! but above bottom

     & -FC(kmin,i)*(zz(kmin,i)-z_lev(m,i))




            else
              k=km(i)
              !write(*,*) k,km

              cff=1.D0/(zz(k+1,i)-zz(k,i))
              p=z_lev(m,i)-zz(k,i)
              q=zz(k+1,i)-z_lev(m,i)

              var_zlv(m,i)=cff*( q*var(k,i) + p*var(k+1,i)
     & -cff*p*q*( cff*(q-p)*(var(k+1,i)-var(k,i))
     & +p*FC(k+1,i) -q*FC(k,i) )
     & )







            !write(*,*) 'bof',i,k,zz(i,k), zz(i,k+1), z_lev(i,m), m



            endif
          enddo
        enddo ! <-- m


      return
      end
# 17 "tools_fort.F" 2
# 1 "sigma_to_z_intr_bounded.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bounded (Lm,Mm,N, nz, z_r, z_w,
     & rmask, var,z_lev, var_zlv,
     & imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev.
!


      implicit none

      integer Lm,Mm,N, nz, imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(kmin:N,jmin:Mm+1,imin:Lm+1),
     & z_r(N,0:Mm+1,0:Lm+1), rmask(0:Mm+1,0:Lm+1),
     & z_w(0:N,0:Mm+1,0:Lm+1), z_lev(nz,jmin:Mm+1,imin:Lm+1),
     & FillValue, var_zlv(nz,jmin:Mm+1,imin:Lm+1),
     & zz(0:N+1,0:Lm+1), dpth



      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max( trd *chunk_size -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv
# 54 "sigma_to_z_intr_bounded.F"
      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=z_r(k,j,i)
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=z_w(0,j,i)
              zz(N+1,i)=z_w(N,j,i)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_r(k,j,i)+z_r(k,j,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i-1)+z_w(0,j,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i-1)+z_w(N,j,i))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5*(z_r(k,j,i)+z_r(k,j-1,i))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i)+z_w(0,j-1,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i)+z_w(N,j-1,i))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_r(k,j,i)+z_r(k,j,i-1)
     & +z_r(k,j-1,i)+z_r(k,j-1,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.25D0*( z_w(0,j,i)+z_w(0,j,i-1)
     & +z_w(0,j-1,i)+z_w(0,j-1,i-1))

              zz(N+1,i)=0.25D0*( z_w(N,j,i)+z_w(N,j,i-1)
     & +z_w(N,j-1,i)+z_w(N,j-1,i-1))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=z_w(k,j,i)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_w(k,j,i)+z_w(k,j,i-1))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5*(z_w(k,j,i)+z_w(k,j-1,i))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_w(k,j,i)+z_w(k,j,i-1)
     & +z_w(k,j-1,i)+z_w(k,j-1,i-1))
              enddo
            enddo
          endif
        endif
# 155 "sigma_to_z_intr_bounded.F"
        do m=1,nz


          if (kmin.eq.0) then !
            do i=imin,imax !
              dpth=zz(N,i)-zz(0,i)
              if (rmask(j,i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(N+1,i)-zz(0,i)
              if (rmask(j,i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N+1,i)).gt.0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2 !<-- below bottom
              elseif (dpth*(zz(1,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(k+1,i)-z_lev(m,j,i))*(z_lev(m,j,i)-zz(k,i))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(m,j,i)=FillValue !<-- masked out
            elseif (km(i).eq.-2) then
              var_zlv(m,j,i)=var(1,j,i) !<-- below bottom
            elseif (km(i).eq.N+2) then
              var_zlv(m,j,i)=var(N,j,i) !<-- above surface
            elseif (km(i).eq.N) then
              var_zlv(m,j,i)=var(N,j,i) !-> R-point, above z_r(N)
            elseif (km(i).eq.kmin-1) then !-> R-point below z_r(1),
              var_zlv(m,j,i)=var(1,j,i) !<-- below bottom
            else
              k=km(i)
              !write(*,*) k,km
# 223 "sigma_to_z_intr_bounded.F"
              var_zlv(m,j,i)=( var(k,j,i)*(zz(k+1,i)-z_lev(m,j,i))
     & +var(k+1,j,i)*(z_lev(m,j,i)-zz(k,i))
     & )/(zz(k+1,i)-zz(k,i))



            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m
# 243 "sigma_to_z_intr_bounded.F"
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
# 18 "tools_fort.F" 2
# 1 "sigma_to_z_intr_sfc.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_sfc (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     & z_lev, var_zlv, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev.
!


      implicit none

      integer Lm,Mm,N, nz, imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(kmin:N,jmin:Mm+1,imin:Lm+1),
     & z_r(N,0:Mm+1,0:Lm+1), rmask(0:Mm+1,0:Lm+1),
     & z_w(0:N,0:Mm+1,0:Lm+1), z_lev(nz,jmin:Mm+1,imin:Lm+1),
     & FillValue, var_zlv(nz,jmin:Mm+1,imin:Lm+1),
     & zz(0:N+1,0:Lm+1), dpth

     & , dz(kmin-1:N,0:Lm+1), FC(kmin-1:N,0:Lm+1), p,q,cff

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max( trd *chunk_size -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv
# 53 "sigma_to_z_intr_sfc.F"
      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=z_r(k,j,i)
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=z_w(0,j,i)
              zz(N+1,i)=z_w(N,j,i)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_r(k,j,i)+z_r(k,j,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i-1)+z_w(0,j,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i-1)+z_w(N,j,i))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.5*(z_r(k,j,i)+z_r(k,j-1,i))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.5D0*(z_w(0,j,i)+z_w(0,j-1,i))
              zz(N+1,i)=0.5D0*(z_w(N,j,i)+z_w(N,j-1,i))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_r(k,j,i)+z_r(k,j,i-1)
     & +z_r(k,j-1,i)+z_r(k,j-1,i-1))
              enddo
            enddo
            do i=imin,imax
              zz(0,i)=0.25D0*( z_w(0,j,i)+z_w(0,j,i-1)
     & +z_w(0,j-1,i)+z_w(0,j-1,i-1))

              zz(N+1,i)=0.25D0*( z_w(N,j,i)+z_w(N,j,i-1)
     & +z_w(N,j-1,i)+z_w(N,j-1,i-1))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=z_w(k,j,i)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5D0*(z_w(k,j,i)+z_w(k,j,i-1))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.5*(z_w(k,j,i)+z_w(k,j-1,i))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(k,i)=0.25D0*( z_w(k,j,i)+z_w(k,j,i-1)
     & +z_w(k,j-1,i)+z_w(k,j-1,i-1))
              enddo
            enddo
          endif
        endif

        do k=kmin,N-1
          do i=imin,imax
            dz(k,i)=zz(k+1,i)-zz(k,i)
            FC(k,i)=var(k+1,j,i)-var(k,j,i)
          enddo
        enddo
        do i=imin,imax
          dz(kmin-1,i)=dz(kmin,i)
          FC(kmin-1,i)=FC(kmin,i)

          dz(N,i)=dz(N-1,i)
          FC(N,i)=FC(N-1,i)
        enddo
        do k=N,kmin,-1 !--> irreversible
          do i=imin,imax
            cff=FC(k,i)*FC(k-1,i)
            if (cff.gt.0.D0) then
              FC(k,i)=cff*(dz(k,i)+dz(k-1,i))/( (FC(k,i)+FC(k-1,i))
     & *dz(k,i)*dz(k-1,i) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo

        do m=1,nz


          if (kmin.eq.0) then !
            do i=imin,imax !
              dpth=zz(N,i)-zz(0,i)
              if (rmask(j,i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(N+1,i)-zz(0,i)
              if (rmask(j,i).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(m,j,i)-zz(N+1,i)).gt.0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(z_lev(m,j,i)-zz(N,i)).gt.0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(0,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=-2 !<-- below bottom
              elseif (dpth*(zz(1,i)-z_lev(m,j,i)).gt.0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(k+1,i)-z_lev(m,j,i))*(z_lev(m,j,i)-zz(k,i))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(m,j,i)=0. !<-- masked out
            elseif (km(i).eq.-2) then

              var_zlv(m,j,i)=FillValue !<-- below bottom

            elseif (km(i).eq.N+2) then
# 215 "sigma_to_z_intr_sfc.F"
              var_zlv(m,j,i)=-FillValue !<-- above surface

            elseif (km(i).eq.N) then
              var_zlv(m,j,i)=var(N,j,i) !-> R-point, above z_r(N)

     & +FC(N,i)*(z_lev(m,j,i)-zz(N,i))




            elseif (km(i).eq.kmin-1) then !-> R-point below z_r(1),
              var_zlv(m,j,i)=var(kmin,j,i) ! but above bottom

     & -FC(kmin,i)*(zz(kmin,i)-z_lev(m,j,i))




            else
              k=km(i)
              !write(*,*) k,km

              cff=1.D0/(zz(k+1,i)-zz(k,i))
              p=z_lev(m,j,i)-zz(k,i)
              q=zz(k+1,i)-z_lev(m,j,i)

              var_zlv(m,j,i)=cff*( q*var(k,j,i) + p*var(k+1,j,i)
     & -cff*p*q*( cff*(q-p)*(var(k+1,j,i)-var(k,j,i))
     & +p*FC(k+1,i) -q*FC(k,i) )
     & )







            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m
# 266 "sigma_to_z_intr_sfc.F"
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
# 19 "tools_fort.F" 2
# 1 "single_connect.F" 1
      subroutine single_connect(Lm,Mm,i0,j0,rmask,outmask)

! Purpose: to enforce single-domain connectedness of unmasked area.
!-------------------------------------------------------------------
! The problem to be addressed is as follows: the landmask is generated
! from topography by declaring that land is where depth is less than
! user-specified minimum depth, this may leave holes in the landmask
! -- lakes which are not connected to the main body of water.
! This program takes a user-specified initial point, indices (i0,j0),
! which is assumed to belong to the main body of water and, starting
! from this point starts filling the adjacent non-masked points as
! "water" (so for each non-masked point on the grid it looks at all the
! immediate neighboring points and if at least one of them is already
! labelled as "water", then the point itself becomes "water").
! Note that at any time during the process there are three types of
! points: 1-masked ("land" according to the initial landmask mask),
! which remain unchanged; 2-unmasked, but non-labelled as "water",
! (initially all unmasked points); and 3-unmasked labelled as "water"
! (initially there is only one such point, i0,j0). The labelling
! continues until the number of "water" points no longer grows, which
! is used as termination signal. At the end the new landmask is
! defined as the points which are not "water" (i.e., either initially
! masked points, or unmasked, but not labelled as "water").




      implicit none
      real(kind=8), dimension(Lm+2,Mm+2) :: rmask,outmask
      integer(kind=2), dimension(:,:), allocatable :: imask, mss,mss2
      integer nargs, iargc, ncgrd, i0,j0, Lm,Mm, i,j, ierr, lstr, lgrd
      logical show_changes

      include "netcdf.inc"

cf2py intent(in) Lm,Mm,i0,j0,rmask
cf2py intent(out) outmask
      show_changes=.false.
        if (rmask(i0,j0) > 0.5) then
          allocate(imask(0:Lm+1, 0:Mm+1))
          allocate(mss (-1:Lm+2,-1:Mm+2))
          allocate(mss2(-1:Lm+2,-1:Mm+2))

C$OMP PARALLEL SHARED(Lm,Mm, i0,j0,show_changes, rmask,imask,mss,mss2)
          call sin_con_thread(Lm,Mm, i0,j0, show_changes, rmask,imask,
     & mss,mss2)
          outmask=rmask
C$OMP END PARALLEL

        else
          write(*,'(/1x,2A,i4,2x,A,i4,1x,A/)') '### ERROR: selected ',
     & 'point i =', i0, 'j =', j0, 'is on land. Try another point.'
        endif

      end subroutine single_connect

      subroutine sin_con_thread(Lm,Mm, i0,j0,show_changes, rmask,
     & imask, mss,mss2)
      implicit none
      integer :: trd_count, wtr_pts, wtr_pts_bak
      integer Lm,Mm, i0,j0
      logical show_changes
      real(kind=8) rmask(0:Lm+1,0:Mm+1)
      integer(kind=2) imask(0:Lm+1,0:Mm+1), mss(-1:Lm+2,-1:Mm+2),
     & mss2(-1:Lm+2,-1:Mm+2)
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last,
     & range, tile, istr,iend,jstr,jend, my_sum, iter, i,j

     & , m

C$ integer omp_get_num_threads, omp_get_thread_num

      integer(kind=4) iclk_start, iclk_end, clk_rate, clk_max
      call system_clock(iclk_start, clk_rate, clk_max)

      trd_count=0
      wtr_pts=0
      wtr_pts_bak=-1
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()

      call set_tiles(Lm,Mm, nsub_x,nsub_y)

c** nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      if (trd==0) write(*,'(/1x,2A,2(I4,1x,A))') 'Enforce that all ',
     & 'water points can be reached from point  i =',i0, 'j =',j0,'...'

      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

        if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
        if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

        do j=jstr,jend
          do i=istr,iend
            if (rmask(i,j) > 0.5) then
              imask(i,j)=1
            else
              imask(i,j)=0
            endif
          enddo
        enddo

        if (istr==0) istr=istr-1 ; if (iend==Lm+1) iend=iend+1
        if (jstr==0) jstr=jstr-1 ; if (jend==Mm+1) jend=jend+1

        do j=jstr,jend
          do i=istr,iend
            mss(i,j)=0
            mss2(i,j)=0
          enddo
        enddo
        if (istr<=i0 .and. i0<=iend .and.
     & jstr<=j0 .and. j0<=jend) then
          mss(i0,j0)=1
        endif
      enddo !<-- tile
C$OMP BARRIER
                                          ! The following while() loop
      iter=0 ! body consists of two nearly
      do while (wtr_pts /= wtr_pts_bak) ! identical segments which
        iter=iter+1 ! differ only by switching
        my_sum=0 ! arrays "mss" and "mss2".
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
          if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

          do j=jstr,jend
            do i=istr,iend
              if ( imask(i,j)>0 .and. ( mss(i,j )>0 .or.
     & mss(i+1,j)>0 .or. mss(i,j+1)>0 .or.
     & mss(i-1,j)>0 .or. mss(i,j-1)>0 )) then
                my_sum=my_sum+1
                mss2(i,j)=1
              endif
            enddo
          enddo
        enddo !<-- tile
C$OMP CRITICAL(cr_region)
        if (trd_count==0) then
          wtr_pts_bak=wtr_pts
          wtr_pts=my_sum
        else
          wtr_pts=wtr_pts+my_sum
        endif
        trd_count=trd_count+1
        if (trd_count==numthreads) then
          trd_count=0
          if (mod(iter,50)==0 .or. wtr_pts==wtr_pts_bak) then
            write(*,*) 'iter =', iter, '  wtr_pts =', wtr_pts,
     & '  changes =', wtr_pts-wtr_pts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER

        iter=iter+1
        my_sum=0
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
          if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

          do j=jstr,jend
            do i=istr,iend
              if ( imask(i,j)>0 .and. ( mss2(i,j )>0 .or.
     & mss2(i+1,j)>0 .or. mss2(i,j+1)>0 .or.
     & mss2(i-1,j)>0 .or. mss2(i,j-1)>0 )) then
                my_sum=my_sum+1
                mss(i,j)=1
              endif
            enddo
          enddo
        enddo !<-- tile
C$OMP CRITICAL(cr_region)
        if (trd_count==0) then
          wtr_pts_bak=wtr_pts
          wtr_pts=my_sum
        else
          wtr_pts=wtr_pts+my_sum
        endif
        trd_count=trd_count+1
        if (trd_count==numthreads) then
          trd_count=0
          if (mod(iter,50)==0 .or. wtr_pts==wtr_pts_bak) then
            write(*,*) 'iter =', iter, '  wtr_pts =', wtr_pts,
     & '  changes =', wtr_pts-wtr_pts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo !<-- while()


      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

        if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
        if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

        if (show_changes) then
          do j=jstr,jend
            do i=istr,iend
              if (mss(i,j)==0 .and. imask(i,j)>0) imask(i,j)=-1
            enddo
          enddo
        else
          do j=jstr,jend
            do i=istr,iend
              if (mss(i,j)==0) imask(i,j)=0
            enddo
          enddo
        endif
      enddo
C$OMP BARRIER


! The following part is to identify and close unresolved features
! of land mask, such as single-point bays. For each water point the
! immediately adjacent surroundings is inspected and a weighting value
! is calculated by adding 8 for each water point adjacent on the east,
! west, north, and south sides, while adding 1 for each water point
! adjacent in diagonal direction. If the resultant value falls below
! a pre-set threshold, the point is set to land. The procedure is
! repeated iteratively until it is detected that further iterations
! do not result in any progress. In the three examples below "."
! means water, "x" land.
!
! . . x x x x . . x x x x x . . x x x x x x x .
! . . . . x x . . . x x x x . . . x x x x . . .
! . . . . x x . . . . o o x . . . . . . . . . .
! . . . x x x . . x x x x x . . x x x x . . . .
! . . x x x x . . x x x x x . . x x x x x x . .
!
! all water and initially water narrow passage: all
! land points will points "o" will points will be kept
! be kept as is be turned to land as they are



C$OMP CRITICAL(cr_region)
      if (trd_count==0) then
        wtr_pts_bak=wtr_pts+1
      endif

      trd_count=trd_count+1

      if (trd_count==numthreads) then
        trd_count=0
        write(*,'(/1x,A)') 'Closing unresolved narrow bays ...'
      endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER

      iter=0
      do while (wtr_pts /= wtr_pts_bak)
        iter=iter+1

        do tile=my_first,my_last
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          do j=jstr,jend
            do i=istr,iend
              m=0
              if (imask(i+1,j )>0) m=m+8
              if (imask(i+1,j+1)>0) m=m+1
              if (imask(i ,j+1)>0) m=m+8
              if (imask(i-1,j+1)>0) m=m+1
              if (imask(i-1,j )>0) m=m+8
              if (imask(i-1,j-1)>0) m=m+1
              if (imask(i ,j-1)>0) m=m+8
              if (imask(i+1,j-1)>0) m=m+1
              mss(i,j)=m
            enddo
          enddo

          if (istr==1) then
            do j=jstr,jend
              m=0
              if (imask(istr-1,j-1)>0) m=m+8
              if (imask(istr ,j-1)>0) m=m+1
              if (imask(istr ,j )>0) m=m+8
              if (imask(istr ,j+1)>0) m=m+1
              if (imask(istr-1,j+1)>0) m=m+8
              mss(istr-1,j)=m
            enddo
          endif

          if (iend==Lm) then
            do j=jstr,jend
              m=0
              if (imask(iend+1,j+1)>0) m=m+8
              if (imask(iend ,j+1)>0) m=m+1
              if (imask(iend ,j )>0) m=m+8
              if (imask(iend ,j-1)>0) m=m+1
              if (imask(iend+1,j-1)>0) m=m+8
              mss(iend+1,j)=m
            enddo
          endif

          if (jstr==1) then
            do i=istr,iend
              m=0
              if (imask(i+1,jstr-1)>0) m=m+8
              if (imask(i+1,jstr )>0) m=m+1
              if (imask(i ,jstr )>0) m=m+8
              if (imask(i-1,jstr )>0) m=m+1
              if (imask(i-1,jstr-1)>0) m=m+8
              mss(i,jstr-1)=m
            enddo
          endif

          if (jend==Mm) then
            do i=istr,iend
              m=0
              if (imask(i-1,jend+1)>0) m=m+8
              if (imask(i-1,jend )>0) m=m+1
              if (imask(i ,jend )>0) m=m+8
              if (imask(i+1,jend )>0) m=m+1
              if (imask(i+1,jend+1)>0) m=m+8
              mss(i,jend+1)=m
            enddo
          endif

          if (istr==1 .and. jstr==1) then
            m=0
            if (imask(istr ,jstr-1)>0) m=m+8
            if (imask(istr ,jstr )>0) m=m+1
            if (imask(istr-1,jstr )>0) m=m+8
            mss(istr-1,jstr-1)=m
          endif

          if (istr==1 .and. jend==Mm) then
            m=0
            if (imask(istr-1,jend )>0) m=m+8
            if (imask(istr ,jend )>0) m=m+1
            if (imask(istr ,jend+1)>0) m=m+8
            mss(istr-1,jend+1)=m
          endif

          if (iend==Lm .and. jstr==1) then
            m=0
            if (imask(iend+1,jstr )>0) m=m+8
            if (imask(iend ,jstr )>0) m=m+1
            if (imask(iend ,jstr-1)>0) m=m+8
            mss(iend+1,jstr-1)=m
          endif

          if (iend==Lm .and. jend==Mm) then
            m=0
            if (imask(iend ,jend+1)>0) m=m+8
            if (imask(iend ,jend )>0) m=m+1
            if (imask(iend+1,jend )>0) m=m+8
            mss(iend+1,jend+1)=m
          endif
        enddo !<-- tile
C$OMP BARRIER

        my_sum=0
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
          if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

          if (show_changes) then
            do j=jstr,jend
              do i=istr,iend
                if (mss(i,j)<12 .and. imask(i,j)>0) imask(i,j)=-2
                if (imask(i,j)>0) my_sum=my_sum+1
              enddo
            enddo
          else
            do j=jstr,jend
              do i=istr,iend
                if (mss(i,j)<12) imask(i,j)=0
                if (imask(i,j)>0) my_sum=my_sum+1
              enddo
            enddo
          endif
        enddo

C$OMP CRITICAL(cr_region)
        if (trd_count==0) then
          wtr_pts_bak=wtr_pts
          wtr_pts=my_sum
        else
          wtr_pts=wtr_pts+my_sum
        endif

        trd_count=trd_count+1

        if (trd_count==numthreads) then
          trd_count=0
          if (mod(iter,1)==0 .or. wtr_pts==wtr_pts_bak) then
            write(*,*) 'iter =', iter, '  wtr_pts =', wtr_pts,
     & '  changes =', wtr_pts-wtr_pts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo !<-- while()



      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

        if (istr==1) istr=istr-1 ; if (iend==Lm) iend=iend+1
        if (jstr==1) jstr=jstr-1 ; if (jend==Mm) jend=jend+1

        if (show_changes) then
          do j=jstr,jend
            do i=istr,iend
              if (imask(i,j) > 0) then
                rmask(i,j)=1.D0 !<-- water stays water
              elseif (imask(i,j)==-1) then
                rmask(i,j)=-1.D0 !<-- water turned into land
              elseif (imask(i,j)==-2) then
                rmask(i,j)=-2.D0 !<-- narrow bays closed
              else
                rmask(i,j)=0.D0
              endif
            enddo
          enddo
        else
          do j=jstr,jend
            do i=istr,iend
              if (imask(i,j) > 0) then
                rmask(i,j)=1.D0
              else
                rmask(i,j)=0.D0
              endif
            enddo
          enddo
        endif
      enddo
C$OMP BARRIER


      call system_clock(iclk_end, clk_rate, clk_max)
      if (clk_rate > 0) then
C$OMP CRITICAL(cr_region)
        trd_count=trd_count+1
        if (trd_count==numthreads) then
          trd_count=0
          write(*,'(/ /1x,A,F8.2,1x,A,I4,1x,A/10x,A/)')
     & 'Wall clock time spent in computational part',
     & (iclk_end-iclk_start)/dble(clk_rate),
     & 'sec  running', numthreads, 'threads.',
     & '[reading and writing files excluded]'
        endif
C$OMP END CRITICAL(cr_region)
      endif

      end subroutine sin_con_thread
# 20 "tools_fort.F" 2
!!! diverse !!!
# 1 "roms_read_write.F" 1
! This package contains a complete set of basic operators for reading
! and writing ROMS-style netCDF data files. The "standard" implies that
! horizontal dimensions are named as "xi_","eta_", vertical "s_" with
! corresponding suffix, "rho", "u", and "v" for horizontal dimensions;
! "rho" and "w" for vertical consistently with grid staggering rules
! within ROMS code. Time dimension (whether or not is "unlimited" from
! netCDF point of view) has its name ending with "time".

! Other than the spatial grid staggering rules, all other aspects
! related to netCDF file structure are expected to follow so called
! "CF conventions" as closely as possible. As the result, all get_*
! and put_* routines from this package are known to work for files
! other than ROMS-standard (leaving only 4 of them "read_roms_grid",
! "write_roms_grid", "roms_find_dims", and "roms_check_dims" be
! strictly ROMS specific).

! It should be noted that somewhat similar functionality for reading
! and writing netCDF files can be found in "nc_read_write.F", however
! the distinction between there and the routines in this package is
! that all the ones below having argument "ncid" are expected to have
! netCDF file in open state, while access to a specific variable is
! done by name, hence it is "file-by-ID -- var-by-name" semantics,
! while "nc_read_write.F" uses "file-by-name -- var-by-name". For this
! reason argument "ncid" is always placed before the filename, and the
! latter is used only to write error messages, but has no effect other
! than that.

! Another note is that FORTRAN 2003 standard mandates that the rank
! of argument (scalar vs. array) should be the same for both calling
! routine and the callee, even in the trivial case where the array
! consists of just a single element. Thus, it is formally illegal
! (thought works correctly in practice) to pass a scalar as an
! argument to a routine expecting an array of size 1. The fact
! that size is equal to 1 in known only during runtime, but not at
! compiling, so the compiler instrumented to verify F2003 compliance
! issues an error message and quits. It is for this and only this
! reason routines containing _sclr_ in their names in the list below
! were introduced, even thought their functionality may seem to be
! redundant (below "value" is scalar, while "var" is array).

! The content is:

! init_time(ncid, fname, tname, nrecs, init_year, ierr)

! read_roms_grid (fname, Lm,Mm)
! write_roms_grid (fname, Lm,Mm)
! roms_find_dims (ncid, fname, Lm,Mm,N)
! roms_check_dims (ncid, fname, Lm,Mm,N)

! put_sclr_by_name_real (ncid, vname, value)
! get_sclr_by_name_real (ncid, vname, value)
! put_sclr_by_name_double (ncid, vname, value)
! get_sclr_by_name_double (ncid, vname, value)

! put_var_by_name_real (ncid, vname, var)
! get_var_by_name_real (ncid, vname, var)
! put_var_by_name_double (ncid, vname, var)
! get_var_by_name_double (ncid, vname, var)

! put_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! put_sclr_rec_by_name_double (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_double (ncid, fname, vname, rec, value)

! put_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! put_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)

! put_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! put_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)

! All are subroutines designed to provide sufficient diagnostic
! messages and terminate the execution if something goes wrong
! rather than functions returning non-zero status.

! With the exception of "read_roms_grid" which opens the named file,
! creates its netCDF file ID as an internal variable, reads all the
! relevant data, and closes it after that, all the above procedures
! imply that the file is in opened state, hence input argument "ncid"
! has meaningful value at entry, while argument "fname" is used only
! for error messages id something goes wrong.






! The following module is designed to be completely initialized by
! "read_roms_grid", which includes both allocation arrays with proper
! dimensions matching the actual grid file and filling them with data.
! Note that "angle" is no longer part of the module because all what
! is needed in most cases when grid file is read is cos and sin of
! angle to rotate vector componets, but not angle itself.

      subroutine init_time(ncid, fname, tname, nrecs, units,
     & init_year, ierr)

! Takes netCDF ID "ncid" of a file in open state and name of timing
! variable "tname" which is expected to have a single dimension (not
! necessarily having the same named as the variable itself (as in CF-
! compliant case) and finds the number of records "nrecs" available
! in the file, time units (days, hours, seconds, etc), and the initial
! year from which the time is counted. Basically it expects the timing
! variable to have attribute "units" which looks like
!
! swf_time:units = "days since 2009-01-01 00:00:00" ;
!
! and just reads it.

      implicit none
      character(len=*) fname, tname, units
      integer ncid, nrecs, init_year, ierr, varid, vartype, vardims,
     & vardimids(8), varatts, i,is,ie, lfnm, ltnm,lunt,lstr
      character(len=64) str, varname, dimname

      include "netcdf.inc"

      nrecs=-1 ; call lenstr(fname,lfnm)
      init_year=-1 ; call lenstr(tname,ltnm)

      ierr=nf_inq_varid(ncid, tname(1:ltnm), varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, varname, vartype,
     & vardims, vardimids, varatts)
        if (ierr == nf_noerr) then
          if (vardims == 1) then
            ierr=nf_inq_dim(ncid, vardimids(vardims), dimname, nrecs)
            if (ierr == nf_noerr) then
              ierr=nf_get_att_text(ncid, varid, 'units', str)
              if (ierr == nf_noerr) then

! Decode 'units' attribute: the first word is expected to be units
! as such: days, hours, seconds, the initial year is expected to be
! the first set of four digits.

                call lenstr(str,lstr)
                i=1
                do while(str(i:i) /= ' ' .and. i < lstr)
                  i=i+1
                enddo
                if (str(i:i) == ' ') i=i-1
                units=str(1:i) !<-- the first word in the string
                lunt=i

                do while(i<lstr .and. (str(i:i)<'0'.or.'9'<str(i:i)))
                  i=i+1
                enddo
                if ('0' <= str(i:i) .and. str(i:i) <= '9') then
                  is=i ; ie=0
                  do while('0'<=str(i:i).and.str(i:i)<='9'.and.i<lstr)
                    i=i+1
                  enddo
                  ie=i-1
                  if (ie == is+3) then
                    init_year=0
                    do i=is,ie
                      init_year = 10*init_year +ichar(str(i:i))-48
                    enddo
                  endif
                endif

                write(*,'(1x,A,2(I6,1x,3A))') 'init_time :: found',
     & nrecs, 'records, units = ''', units(1:lunt),
     & ''' starting from year', init_year, 'in ''',
     & fname(1:lfnm), '''.'
                if (init_year < 0) then
                  write(*,'(/1x,6A/)') '### ERROR: init_time :: ',
     & 'Cannot find initial year segment within attribute ',
     & '''units'' for variable ''', tname(1:ltnm),
     & ''' in file ''', fname(1:lfnm), '''.'

                  write(*,*) 'units =''', str(1:lstr), ''''
                  ierr=ierr-999
                endif
              else
                write(*,'(/1x,6A/12x,A/)') '### ERROR: init_time ',
     & ':: Cannot find attribute ''units'' for variable ''',
     & tname(1:ltnm), ''' in file ''', fname(1:lfnm),
     & '''.', nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'init_time :: Cannot determine name and size of ',
     & 'dimension #', vardimids(vardims), ''' in file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)
            endif
          elseif (vardims > 1) then
            write(*,'(/1x,6A,I3/)') '### ERROR: init_time :: ',
     & 'Ambiguous: variable  ''', tname(1:ltnm), ''' in ''',
     & fname(1:lfnm),''' has more than one dimension:',vardims
            ierr=ierr+1
          else
            write(*,'(/1x,6A/)') '### ERROR: init_time :: No ',
     & 'dimension found for variable ''', tname(1:ltnm),
     & ''' in ''', fname(1:lfnm), '''.'
            ierr=ierr+1
          endif
        else
         write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: init_time ',
     & ':: Cannot make general inquiry for variable #', varid,
     & 'in ''', fname(1:lfnm), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A/)') '### ERROR: init_time :: Cannot ',
     & 'find variable ''', tname(1:ltnm), ''' in file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)
      endif
      end


      subroutine read_roms_grid(fname, Lm,Mm)

! Open ROMS grid netCDF file (argument "fname", input), read its
! dimensions (arguments Lm,Mm used in output mode), allocate all the
! arrays defined in the module above, read their values, and compute
! csA,snA, which are cos and sin of the angle between the geographical
! east direction and XI-direction of the curvilinear grid.
! For universality all the variables are allocated and read regardless
! whether they are actually needed or not.

      use roms_grid_vars
! use mod_io_size_acct
      implicit none
      character(len=*) fname
      integer Lm,Mm, ncgrd, i,j, ierr, lfnm
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc

      include "netcdf.inc"

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0
      if (.not.allocated(lon_r)) then
        call lenstr(fname,lfnm)
        ierr=nf_open(fname, nf_nowrite, ncgrd)
        if (ierr == nf_noerr) then
          write(*,'(/1x,4A)') 'Reading ROMS grid file ''',
     & fname(1:lfnm), '''...'
          call roms_find_dims(ncgrd, fname, Lm,Mm, i)

          allocate( lon_r(0:Lm+1,0:Mm+1), lon_p(1:Lm+1,1:Mm+1),
     & lat_r(0:Lm+1,0:Mm+1), lat_p(1:Lm+1,1:Mm+1),
     & pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     & csA(0:Lm+1,0:Mm+1), snA(0:Lm+1,0:Mm+1),
     & f(0:Lm+1,0:Mm+1), rmask(0:Lm+1,0:Mm+1),
     & stat=ierr )
          if (ierr /= 0) then
            write(*,'(/1x,2A/)') '### ERROR: read_roms_grid ',
     & ':: Cannot allocate memory.'
            stop
          endif

          call get_var_by_name_double(ncgrd, 'lon_rho', lon_r)
          call get_var_by_name_double(ncgrd, 'lat_rho', lat_r)
          call get_var_by_name_double(ncgrd, 'pm', pm)
          call get_var_by_name_double(ncgrd, 'pn', pn)
          call get_var_by_name_double(ncgrd, 'f', f)
          call get_var_by_name_double(ncgrd, 'mask_rho', rmask)
          sz_read_acc=sz_read_acc + 8*6*(Lm+2)*(Mm+2)

          ierr=nf_inq_varid(ncgrd, 'lon_psi', i)
          if (ierr == nf_noerr) then
            call get_var_by_name_double(ncgrd, 'lon_psi', lon_p)
            call get_var_by_name_double(ncgrd, 'lat_psi', lat_p)
            sz_read_acc=sz_read_acc + 8*2*(Lm+1)*(Mm+1)
          else
            do j=1,Mm+1
              do i=1,Lm+1
                lon_p(i,j)=0.25D0*( lon_r(i,j)+lon_r(i-1,j)
     & +lon_r(i,j-1) +lon_r(i-1,j-1))
                lat_p(i,j)=0.25D0*( lat_r(i,j)+lat_r(i-1,j)
     & +lat_r(i,j-1) +lat_r(i-1,j-1))
              enddo
            enddo
          endif

          call read_angle(ncgrd,fname, 1,1,Lm+2,Mm+2, csA,snA)
          curv_grid=.false.
          do j=0,Mm+1
            do i=0,Lm+1
              if (abs(snA(i,j)) > 1.D-12) curv_grid=.true.
            enddo
          enddo
          if (curv_grid) then
            write(*,'(1x,2A/)') 'Curvilinear grid discovered. ',
     & 'Vector components will be rotated.'
          endif
          ierr=nf_close (ncgrd)
          return !---> successful return
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: read_roms_grid ',
     & ':: Cannot open netCDF file ''', fname(1:lfnm),
     & ''' for reading,', nf_strerror(ierr)
          stop
        endif
      else
        write(*,*) 'WARNING: Grid is already initialized.'
      endif !<-- .not.allocated(lon_r)
      end


      subroutine read_angle(ncid,fname, iwest,jsouth, isize,jsize,
     & csA,snA)

! Read a rectangular portion (subdomain) of netCDF variable "angle"
! from ROMS grid file, which the angle between local XI-coordinate of
! ROMS grid and true EAST direction and compute its sine and cosine,
! csA=cos(angle) and snA=sin(angle). Arguments iwest,jsouth specify
! the western and southern edges, while isize,jsize -- the sizes of
! the subdomain within the i- and j-dimensions (with this respect it
! is conformal to "get_patch_by_name_double" routine defined below
! within in this file). The angle may be in either radians or degrees,
! which is determined by checking attribute "units". NetCDF file is
! expected to be in open state ("ncid" is a valid netCDF file ID;
! argument "fname" is merely to write error messages on the screen).

      implicit none
      character(len=*) fname
      integer iwest,jsouth, isize,jsize, ncid,varid, ierr,
     & start(2),count(2), i,j, lfnm,lstr
      real(kind=8) csA(isize,jsize),snA(isize,jsize), cff
      character(len=16) str

      include "phys_const.h"
      include "netcdf.inc"

      call lenstr(fname,lfnm)
      ierr=nf_inq_varid(ncid, 'angle', varid)
      if (ierr == nf_noerr) then
        start(1)=iwest ; count(1)=isize
        start(2)=jsouth ; count(2)=jsize
        ierr=nf_get_vara_double(ncid, varid, start,count, csA)
        if (ierr == nf_noerr) then
          ierr=nf_get_att_text(ncid, varid, 'units', str)
          if (ierr == nf_noerr) then
            call lenstr(str,lstr)
            write(*,'(2x,5A)') 'retrieved east angle from ''',
     & fname(1:lfnm), ''', units = ''',str(1:lstr),'''.'
            if (str(1:6) == 'degree') then
C$OMP PARALLEL SHARED(isize,jsize, csA,snA) PRIVATE(i,j, cff)
C$OMP DO
              do j=1,jsize
                do i=1,isize
                  cff=deg2rad*csA(i,j)
                  csA(i,j)=cos(cff)
                  snA(i,j)=sin(cff)
                enddo
              enddo
C$OMP END DO
C$OMP END PARALLEL
            elseif (str(1:6) == 'radian') then
C$OMP PARALLEL SHARED(isize,jsize, csA,snA) PRIVATE(i,j, cff)
C$OMP DO
              do j=1,jsize
                do i=1,isize
                  cff=csA(i,j)
                  csA(i,j)=cos(cff)
                  snA(i,j)=sin(cff)
                enddo
              enddo
C$OMP END DO
C$OMP END PARALLEL
            else
              write(*,'(/1x,4A/)') '### ERROR: Unknown units for ',
     & 'variable ''angle'' in ''', fname(1:lfnm), '''.'
              stop
            endif
          else
            write(*,'(/1x,4A/)') '### ERROR: Cannot read attribute ',
     & '''units'' for variable ''angle'' in ''',fname(1:lfnm),'''.'
          endif
        else
          write(*,'(/1x,4A/)') '### ERROR: Cannot read variable ',
     & '''angle'' from ''', fname(1:lfnm), '''.'
        endif
      else
        write(*,'(/1x,4A/)') '### ERROR: Cannot find variable ',
     & '''angle'' in ''', fname(1:lfnm), '''.'
      endif
      if (ierr /= nf_noerr) stop
      end subroutine read_angle


      subroutine write_roms_grid(fname, Lm,Mm)

! Create ROMS grid file and write all the variables associated with
! horizontal curvilinear coordinates (these are stored inside module
! "roms_grid" where all the arrays are expected to be allocated with
! their dimensions specified by Lm,Mm and assigned meaningful values).
! This routine also creates netcdf variables for model topography and
! mask, "hraw", "h", and "mask_rho", however they are left
! uninitialized to be filled in later.

      use roms_grid_vars
      use roms_grid_params
      implicit none
      character(len=*) fname
      integer Lm,Mm, xi_rho,eta_rho, xi_u,eta_v, ncgrd, varid, ierr,
     & r2dgrd(2), p2dgrd(2), lfnm, lstt, lstr
      character(len=32) str
      character(len=256) settings

      include "netcdf.inc"

      xi_rho=Lm+2 ; eta_rho=Mm+2
      xi_u=xi_rho-1 ; eta_v=eta_rho-1

! Create netcdf file.
!------- ------ -----

      call lenstr(fname,lfnm)
      ierr=nf_create(fname(1:lfnm), nf_netcdf4, ncgrd)
      if (ierr == nf_noerr) then
        ierr=nf_def_dim(ncgrd, 'xi_rho', xi_rho, r2dgrd(1))
        ierr=nf_def_dim(ncgrd, 'xi_u', xi_u, p2dgrd(1))
        ierr=nf_def_dim(ncgrd, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=nf_def_dim(ncgrd, 'eta_v', eta_v, p2dgrd(2))

! Grid type switch: Spherical or Cartesian.

        ierr=nf_def_var(ncgrd, 'spherical', nf_char, 0, 0, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',24,
     & 'grid type logical switch')
        ierr=nf_put_att_text(ncgrd,varid, 'option_T', 9, 'spherical')
        ierr=nf_put_att_text(ncgrd,varid, 'option_F', 9, 'cartesian')

! Longitude/latitude at RHO-points.

        ierr=nf_def_var(ncgrd, 'lon_rho', nf_double, 2, r2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 23,
     & 'longitude of RHO-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 11, 'degree_east')


        ierr=nf_def_var(ncgrd, 'lat_rho', nf_double, 2, r2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid,'long_name',22,
     & 'latitude of RHO-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 12,'degree_north')

! Longitude/latitude at PSI-points.

        ierr=nf_def_var(ncgrd, 'lon_psi', nf_double, 2, p2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 23,
     & 'longitude of PSI-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 11, 'degree_east')


        ierr=nf_def_var(ncgrd, 'lat_psi', nf_double, 2, p2dgrd,varid)
        ierr=nf_put_att_text(ncgrd,varid,'long_name',22,
     & 'latitude of PSI-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 12,'degree_north')

! Curvilinear coordinate metric coefficients pm,pn.

        ierr=nf_def_var(ncgrd, 'pm', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 35,
     & 'curvilinear coordinate metric in XI')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 7, 'meter-1')

        ierr=nf_def_var(ncgrd, 'pn', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd,varid, 'long_name', 36,
     & 'curvilinear coordinate metric in ETA')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 7, 'meter-1')

! Angle between direction to the EAST and XI-axis, at RHO-points

        ierr=nf_def_var(ncgrd, 'angle', nf_double, 2, r2dgrd,varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',30,
     & 'angle between EAST and XI-axis')
        ierr=nf_put_att_text(ncgrd, varid, 'units', 7, 'degrees')

! Coriolis Parameter

        ierr=nf_def_var(ncgrd, 'f', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name', 32,
     & 'Coriolis parameter at RHO-points')
        ierr=nf_put_att_text(ncgrd,varid, 'units', 8, 'second-1')


! Land-Sea mask at RHO-points.

        ierr=nf_def_var(ncgrd,'mask_rho',nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name',18,
     & 'mask on RHO-points')
        ierr=nf_put_att_text(ncgrd, varid, 'option_0', 4, 'land' )
        ierr=nf_put_att_text(ncgrd, varid, 'option_1', 5, 'water')

! Raw and smoothed bathymetry.

        ierr=nf_def_var(ncgrd, 'hraw', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name', 28,
     & 'raw bathymetry at RHO-points')
        ierr=nf_put_att_text(ncgrd, varid, 'units', 5, 'meter')

        ierr=nf_def_var(ncgrd, 'h', nf_double, 2, r2dgrd, varid)
        ierr=nf_put_att_text(ncgrd, varid, 'long_name', 24,
     & 'bathymetry at RHO-points')
        ierr=nf_put_att_text(ncgrd, varid, 'units', 5, 'meter')

        if (allocated(orterr)) then
          ierr=nf_def_var(ncgrd,'ort_error', nf_double, 2,r2dgrd,varid)
          ierr=nf_put_att_text(ncgrd, varid, 'long_name', 19,
     & 'orthogonality error')
          ierr=nf_put_att_text(ncgrd, varid, 'units', 7, 'degrees')
        endif

! Create signature containing parameters used for generating grid and
! save it a global attribute so the grid can be reproduced if needed.

        write(str,*) nx ; call lenstr(str,lstr)
        settings='nx='/ /str(1:lstr); call lenstr(settings,lstt)
        write(str,*) ny ; call lenstr(str,lstr)
        settings=settings(1:lstt)/ /' ny='/ /str(1:lstr)
        call lenstr(settings,lstt)

        if (lat_max > lat_min .or. lon_max > lon_min) then

          write(str,*) lat_min ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' lat_min='/ /str(1:lstr)
          call lenstr(settings,lstt)

          write(str,*) lat_max ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' lat_max='/ /str(1:lstr)
          call lenstr(settings,lstt)

          write(str,*) lon_min ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' lon_min='/ /str(1:lstr)
          call lenstr(settings,lstt)

          write(str,*) lon_max ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' lon_max='/ /str(1:lstr)
          call lenstr(settings,lstt)

        elseif (size_x>0.D0 .or. size_y>0.D0) then !--> convert to km

          write(str,*) size_x * 1.0D-3 ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' size_x='/ /str(1:lstr)
          call lenstr(settings,lstt)

          write(str,*) size_y * 1.0D-3 ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' size_y='/ /str(1:lstr)
          call lenstr(settings,lstt)
        endif

        if (cent_lat /= 0.D0) then
          write(str,*) cent_lat ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' cent_lat='/ /str(1:lstr)
          call lenstr(settings,lstt)
        endif

        if (psi0 /= 0.D0) then
          write(str,*) psi0 ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' Lon='/ /str(1:lstr)
          call lenstr(settings,lstt)
        endif

        if (theta0 /= 0.D0) then
          write(str,*) theta0 ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' Lat='/ /str(1:lstr)
          call lenstr(settings,lstt)
        endif

        if (alpha /= 0.D0) then
          write(str,*) alpha ; call lenstr(str,lstr)
          settings=settings(1:lstt)/ /' rotate='/ /str(1:lstr)
          call lenstr(settings,lstt)
        endif

        write(str,*) flip_xy; call lenstr(str,lstr)
        settings=settings(1:lstt)/ /' flip_xy='/ /str(1:lstr)
        call lenstr(settings,lstt)

        ierr=nf_put_att_text(ncgrd, nf_global, 'Settings',
     & lstt, settings)

! Leave definition mode.
! ----- ---------- -----

        ierr=nf_enddef(ncgrd)

! Grid type switch: ALWAYS SPHERICAL

        ierr=nf_inq_varid (ncgrd, 'spherical', varid)
        ierr=nf_put_var1_text (ncgrd, varid, 1, 'T')

! Longitude/latitude at RHO- and PSI-points.

        call put_var_by_name_double(ncgrd, 'lon_rho', lon_r)
        call put_var_by_name_double(ncgrd, 'lat_rho', lat_r)
        call put_var_by_name_double(ncgrd, 'lon_psi', lon_p)
        call put_var_by_name_double(ncgrd, 'lat_psi', lat_p)

! Curvilinear coordinate metric coefficients pm,pn.

        call put_var_by_name_double(ncgrd, 'pm', pm)
        call put_var_by_name_double(ncgrd, 'pn', pn)

! Angle between XI-axis and EAST at RHO-points

c* call put_var_by_name_double(ncgrd, 'angle', angle)

! Coriolis Parameter.

        call put_var_by_name_double(ncgrd, 'f', f)

        if (allocated(orterr)) then
          call put_var_by_name_double(ncgrd, 'ort_error', orterr)
        endif

! Close netCDF file
! ----- ------ ----

        ierr=nf_close(ncgrd)
      else
        write(*,'(/1x,4A/12x,A)') '### ERROR: Cannot create netCDF ',
     & 'file ''', fname(1:lfnm), '''.', nf_strerror(ierr)
      endif
      end subroutine write_roms_grid


! The following two routines are either to find dimensions from already
! opened netCDF file, or to check whether dimensions in the file match
! the supplied values. The second variant is needed merely as a check
! point to detect the mismatch and to stop the execution if it happens.
! In both cases the presence of both horizontal dimensions is mandatory
! while vertical is optional (i.e., "find_dims" does not touch its
! argument N, if no dimension, so N retains its previous value or
! remains uninitialized, if not assigned by the caller; "check_dims"
! ignores vertical mismatch, if input value of "N_ck" is zero, or if
! vertical dimension cannot be determined from the file. Both routines
! are designed to terminate execution if something goes wrong.


      subroutine roms_find_dims(ncgrd, fname, Lm, Mm, N)





      implicit none
      character(len=*) fname
      integer ncgrd, Lm,Mm,N, xi_rho,xi_u, eta_rho,eta_v, s_rho,s_w,
     & ndims, size, id, i,is, ierr, lvar, lfnm



      character(len=16) dname
      character(len=128) string

      include "netcdf.inc"

      xi_rho=0 ; xi_u=0 ; s_rho=0
      eta_rho=0 ; eta_v=0 ; s_w=0

      call lenstr(fname,lfnm)
      ierr=nf_inq_ndims(ncgrd, ndims)
      if (ierr == nf_noerr) then
        do id=1,ndims
          dname='                '
          ierr=nf_inq_dim (ncgrd, id, dname, size)
          if (ierr == nf_noerr) then
            call lenstr(dname,lvar)
            if (lvar == 6 .and. dname(1:lvar) == 'xi_rho') then
              xi_rho=size
            elseif (lvar == 4 .and. dname(1:lvar) == 'xi_u') then
              xi_u=size

            elseif (lvar == 7 .and. dname(1:lvar) == 'eta_rho') then
              eta_rho=size
            elseif (lvar == 5 .and. dname(1:lvar) == 'eta_v') then
              eta_v=size

            elseif (lvar == 5 .and. dname(1:lvar) == 's_rho') then
              s_rho=size
            elseif (lvar == 3 .and. dname(1:lvar) == 's_w') then
              s_w=size

            elseif (lvar == 5 .and. dname(1:lvar) == 'depth') then
              s_rho=size
            elseif (lvar == 7 .and. dname(1:lvar) == 'rho_ntr') then
              s_rho=size
            endif
          else
            write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     & 'determine name and size of dimension #', id,
     & 'in ''', fname(1:lfnm), '''.', nf_strerror(ierr)
          endif
        enddo

        write(string,'(A,6(1x,A,I4))')

     & ' roms_find_dims ::',



     & 'xi_rho=',xi_rho, 'xi_u=',xi_u, 'eta_rho=',eta_rho,
     & 'eta_v=',eta_v, 's_rho=',s_rho, 's_w=',s_w
        call lenstr(string,lvar)
        i=0 ! Write dimensions into
        do while(i < lvar) ! character string first,
          i=i+1 ! and then suppress blank
          if (string(i:i) == '=') then ! characters after =sign.
            i=i+1 ! This is merely to make
            if (string(i:i) == ' ') then
              is=1
              do while(string(i+is:i+is) == ' ' .and. i+is < lvar)
                is=is+1
              enddo
              string(i:lvar-is)=string(i+is:lvar) ; lvar=lvar-is
            endif
          endif
        enddo ! a narrower printout
        write(*,'(2x,A)') string(1:lvar) ! on the screen.

        ierr=0
        if (xi_rho > 0) then
          Lm=xi_rho-2
        elseif (xi_u > 0) then
          Lm=xi_u-1
        else
          write(*,'(/1x,4A/)') '### ERROR: Cannot determine size ',
     & 'of horizontal XI-dimension in netCDF file ''',
     & fname(1:lfnm), '''.'
          ierr=ierr+1
        endif
        if (eta_rho > 0) then
          Mm=eta_rho-2
        elseif (eta_v > 0) then
          Mm=eta_rho-1
        else
          write(*,'(/1x,4A/)') '### ERROR: Cannot determine size ',
     & 'of horizontal ETA-dimension in netCDF file ''',
     & fname(1:lfnm), '''.'
          ierr=ierr+1
        endif ! The policy here is that vertical



        if (s_rho > 0) then ! it is filled up only if found, and
          N=s_rho ! "not touched" otherwise (hence it
        elseif (s_w > 0) then ! is possible to call this function
          N=s_w-1 ! while passing a constant, if no
        endif ! vertical dimension is expected to
                                   ! exist in the file.
# 774 "roms_read_write.F"
        if (ierr /= 0) stop !--> ERROR
      else
        write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot determine ',
     & 'number of dimensions in netCDF file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      end


! The rest are the reading-writting subroutins generated by CPP
! from the same source code package (in fact, quadrupled: X2 due
! to read/write functionality and another X2 due to single/double
! precision version).
# 844 "roms_read_write.F"
! The following eight routines are just instrumented wrappers around
! the standard sequence of netCDF calls which (1) inquire variable ID
! and (2) put/get the ENTIRE variable into/from netCDF file. These
! wrapper is needed solely to write error messages if something goes
! wrong. These are
!
! put/get_sclr/var_by_name_TYPE (ncid, vname, value/var)
!
! where get/put and TYPE=real/double occur in all permutations (hence
! it adds up to a total of eight). Because of semantically identical
! code real/double is implemented by CPP-redefinition of basic netCDF
! functions using the same source code.




      subroutine get_sclr_by_name_real(ncid, vname, value)


! These two routines are for reading or writing just a single number
! which may exist in netCDF file either as a variable or a global
! attribute containing just a single number of the proper type.
! Selection between variable or attribute is by the file, while no
! attempt to change is format is made here. For this reason is either
! variable of attribute with given name must pre-exist in order for
! these operations to succeed.

      implicit none
      integer ncid, varid, type, size, ierr, lvar
      character(len=*) vname
      real(kind=4) value
      include "netcdf.inc"

      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname(1:lvar), varid)
      if (ierr == nf_noerr) then
# 889 "roms_read_write.F"
        ierr=nf_get_var_real(ncid, varid, value)
        if (ierr == nf_noerr) then
          write(*,'(9x,3A)') 'read ''', vname(1:lvar), ''''
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read netCDF ',
     & 'variable ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif

      else
        ierr=nf_inq_att(ncid, nf_global, vname(1:lvar), type, size)
        if (ierr == nf_noerr) then




          if (size == 1 .and. type == nf_real) then
# 931 "roms_read_write.F"
            ierr=nf_get_att_real(ncid, nf_global, vname(1:lvar),
     & value)
            if (ierr == nf_noerr) then
              write(*,'(9x,4A)') 'read ''', vname(1:lvar),
     & ''' as global attribute'
            else
              write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read ',
     & 'global attribute ''', vname(1:lvar),
     & ''' from netCDF file.', nf_strerror(ierr)
            endif

          else
            if (size /= 1) then
              write(*,'(/1x,5A,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong size', size, 'instead of 1.'
              ierr=ierr-1
            endif



            if (type /= nf_real) then

              write(*,'(/1x,5A,I4,1x,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong type', type, 'instead of',



     & nf_real, 'which is real.'

              ierr=ierr-1
            endif
          endif
        else
          write(*,'(/1x,4A/)') '### ERROR: Neither variable, nor ',
     & 'global attribute named ''', vname,
     & ''' is present in netCDF file.'
        endif
      endif
      if (ierr /= nf_noerr) stop
      end
# 981 "roms_read_write.F"
      subroutine get_var_by_name_real(ncid, vname, var)

      implicit none
      integer ncid, ierr, varid, lvar
      character(len=*) vname
      real(kind=4) var(*)
      include "netcdf.inc"
      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then





        ierr=nf_get_var_real(ncid, varid, var)
        if (ierr == nf_noerr) then
          write(*,'(9x,3A)') 'read ''', vname(1:lvar), ''''

        else



          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read ',

     & 'netCDF variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot find ',
     & 'netCDF variable ID for ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
      endif
      if (ierr /= nf_noerr) stop
      end


! The following routines are to put/get just a single number into
! a specified record of one-dimensional netCDF variable. The netCDF
! file "fname" is expected to be open (hence input argument "ncid" is
! a valid file ID, while the name of the file is needed only to write
! error messages if something goes wrong, but otherwise is not used);
! "vname" is name of the variable (to be translated inside into netCDF
! ID with error message if not found), and "rec" is record number;
! "var" is input for put_ (output for get_) is just a scalar (single
! number). The corresponding netCDF variable is expected to be either
! a one-dimensional array (having more dimensions results in error
! message) or a scalar. If array it puts/gets the value at location
! "rec" with performing necessary error checking; if scalar it takes
! the only value, while argument "rec" is not used. Again, get/put
! and TYPE=real/double occur in all four permutations.






      subroutine get_sclr_rec_by_name_real(ncid, fname, vname,
     & rec, value)

      implicit none
      integer ncid, rec
      character(len=*) fname, vname
      real(kind=4) value

      character(len=16) name
      integer varid, vtype, ndims, natts, dimid(8), size,
     & start(4), count(4), ierr, lfnm, lvar
      include "netcdf.inc"

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then
          if (ndims == 1) then
            ierr=nf_inq_dimlen(ncid, dimid(1), size)
            if (ierr == nf_noerr) then
              start(1)=rec ; count(1)=1
# 1074 "roms_read_write.F"
              if (0 < rec .and. rec <= size) then
                ierr=nf_get_vara_real(ncid,varid, start,count, value)
                if (ierr == nf_noerr) then
                  write(*,'(6x,A,I5,1x,5A)') 'read rec', rec,
     & 'of scalar ''', vname(1:lvar), ''' from ''',
     & fname(1:lfnm), '''.'
                  return !---> successful return
                else
                  write(*,'(/1x,7A/12x,A)') '### ERROR: Cannot read ',
     & 'variable ''', vname(1:lvar), ''' from netCDF ',
     & 'file ''',fname(1:lfnm), ''':', nf_strerror(ierr)
                endif
              else
                write(*,'(/1x,2A,I4,1x,6A,I4/)') '### ERROR: ',
     & 'Requested record number ', rec, 'for scalar ',
     & 'variable ''', vname(1:lvar), ''' in file ''',
     & fname(1:lfnm), ''' exceeds dimension bound', size
              endif

            else
              write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     & 'Cannot determine size of dimension #', dimid(1),
     & 'in file ''', fname(1:lfnm), '''.', nf_strerror(ierr)
            endif
          elseif (ndims == 0) then
# 1112 "roms_read_write.F"
            ierr=nf_get_var_real(ncid, varid, value)
            if (ierr == nf_noerr) then
              write(*,'(7x,5A)') 'read scalar ''', vname(1:lvar),
     & ''' from ''', fname(1:lfnm), '''.'
              return !---> successful return
            else
              write(*,'(/1x,5A/12x,A)') '### ERROR: Cannot read ',
     & 'scalar variable ''', vname(1:lvar), ''' from ''',
     & fname(1:lfnm),''':', nf_strerror(ierr)
            endif

          else
            write(*,'(/1x,5A,I4/)') '### ERROR: Variable ''',
     & vname(1:lvar), ''' from file ''', fname(1:lfnm),
     & ''' has more then one dimension, ndims =', ndims
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! A more sophisticated functions, which read/write a portion of data
! into/from an existing netCDF variable accessing it by name. In doing
! so, it checks whether the first one, two, or three dimensions for
! that variable as defined in the netCDF file are consistent with the
! shape of array "var" specified as n1,n2,n3. If not, it complains
! about the error and quits. The variable "var" may be one, two, or
! three-dimensional, and, in addition to that may have record dimension
! (not necessarily unlimited). If the variable is two dimensional,
! then the calling routine must be specify n3=0 to avoid checking of
! the non-existing dimension. Similar policy applies for n1 and n2
! (e.g., n1,n2,n3=0,0,0 means that the variable is a scalar). However,
! if the variable in netCDF file has record dimension (identified as
! either unlimited netCDF dimension, or as the last dimension of the
! variable with its name ending with "...time" AND the actual number
! of spatial dimensions of netCDF variable in the file is LESS that
! the number of non-zero arguments THEN the excess spatial dimensions
! (n3, or both n3 and n2) will be ignored, while "rec" will be
! interpreted normally.

! Note: in the code below arguments n1,n2,n3 are used EXCLUSIVELY for
! checking, while the starts and counts are computed from the the
! dimensions of the variable in netCDF file.





      subroutine get_rec_by_name_real(ncid, fname, vname,
     & n1,n2,n3, rec, var)

! use mod_io_size_acct
      implicit none
      integer ncid, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=4), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0 ! Determine the number of spatial
        count(id)=0 ! dimensions "nspc" for the variable.
      enddo
      nspc=0 ! Note that "nspc" found here should be
      if (n1 > 0) then ! either
        nspc=nspc+1 ! less by 1 than the actual number of
        count(nspc)=n1 ! dimensions "ndims" of the variable
      endif ! stored in the file -- in this case
      if (n2 > 0) then ! the extra dimension will be treated
        nspc=nspc+1 ! as record dimension;
        count(nspc)=n2 ! or
      endif ! be the same as "ndims" -- in this
      if (n3 > 0) then ! case there is no record dimension,
        nspc=nspc+1 ! and argument "rec" is ignored.
        count(nspc)=n3
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then






          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1244 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  start(id)=1 ; rec_size=rec_size*size

                  if (count(id) /= size) then
                    call lenstr(name,ldim)
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',
# 1271 "roms_read_write.F"
     & 'get_rec_by_name_real',


     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1331 "roms_read_write.F"
              endif







              ierr=nf_get_vara_real(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                write(*,'(6x,A,I5,1x,5A)') 'read rec', rec, 'of ''',
     & vname(1:lvar), ''' from ''', fname(1:lfnm), '''.'
                sz_read_acc = sz_read_acc + rec_size * 4


                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max



                read_clk = read_clk + inc_clk


                return !---> successful return
              else





                write(*,'(/1x,7A/12x,A)') '### ERROR: Cannot read ',
     & 'variable ''',vname(1:lvar),''' from netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)

              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! The following pair is essentially an instrumented versions of
! put_vara and get_vara to read a single record of a 2D-subdomain
! [iwest : iwest+n1-1] x [jsouth : jsouth+n2-1]
! within the (i,j)-index-space of netCDF array. The third dimension
! (if exists) is treated as the whole: n3 is checked against the
! actual dimension of netCDF variable and the mismatch is treated
! as an error; n1,n2 are checked only for upper-bound overrun, e.g.,
! iwest+n1-1 exceeds the actual size of the first dimension
! resulting in error.





      subroutine get_patch_by_name_real(ncid, fname, vname,
     & iwest,jsouth, n1,n2,n3, rec, var)

c--#define VERBOSE
! use mod_io_size_acct
      implicit none
      integer ncid, iwest,jsouth, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=4), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0
      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0
        count(id)=0
      enddo ! Determine the number of spatial
      nspc=0 ! dimensions "nspc" for the variable.
      if (n1 > 0) then
        nspc=nspc+1 ! Note that "nspc" found here should be
        start(nspc)=iwest !
        count(nspc)=n1 ! either
      endif ! less by 1 than the actual number of
      if (n2 > 0) then ! dimensions "ndims" of the variable
        nspc=nspc+1 ! stored in the file -- in this case
        start(nspc)=jsouth ! the extra dimension will be treated
        count(nspc)=n2 ! as record dimension,
      endif ! or
      if (n3 > 0) then ! be the same as "ndims" -- in this
        nspc=nspc+1 ! case there is no record dimension,
        count(nspc)=n3 ! and argument "rec" is ignored.
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then





          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1478 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  rec_size=rec_size*count(id)
                  if (start(id) == 0) then
                    start(id)=1
                    if (count(id) /= size) then
                      write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',
# 1505 "roms_read_write.F"
     & 'get_patch_by_name_real',


     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                      stop
                    endif
                  elseif (start(id)+count(id)-1 > size) then
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,3(I5,1x,A)/)')
     & '### ERROR: ',
# 1527 "roms_read_write.F"
     & 'get_patch_by_name_real',


     & ' :: Overrun dimension bound #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', start(id), '+',
     & count(id), '-1 >', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1587 "roms_read_write.F"
              endif
# 1600 "roms_read_write.F"
              ierr=nf_get_vara_real(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                if (start(ndims) == rec .and. count(ndims) == 1) then
                  write(*,'(6x,A,I5,1x,5A)') 'read rec', rec,'of ''',
     & vname(1:lvar), ''' from ''', fname(1:lfnm), ''''
                else
                  write(*,'(6x,5A)') 'read ''', vname(1:lvar),
     & ''' from ''', fname(1:lfnm), ''''
                endif
                sz_read_acc = sz_read_acc + rec_size * 4


                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max



                read_clk = read_clk + inc_clk


                return !---> successful return
              else





                write(*,'(/1x,7A/12x,A)') '### ERROR: Cannot read ',
     & 'variable ''',vname(1:lvar),''' from netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)

              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end




# 1 "roms_read_write.F" 1
! This package contains a complete set of basic operators for reading
! and writing ROMS-style netCDF data files. The "standard" implies that
! horizontal dimensions are named as "xi_","eta_", vertical "s_" with
! corresponding suffix, "rho", "u", and "v" for horizontal dimensions;
! "rho" and "w" for vertical consistently with grid staggering rules
! within ROMS code. Time dimension (whether or not is "unlimited" from
! netCDF point of view) has its name ending with "time".

! Other than the spatial grid staggering rules, all other aspects
! related to netCDF file structure are expected to follow so called
! "CF conventions" as closely as possible. As the result, all get_*
! and put_* routines from this package are known to work for files
! other than ROMS-standard (leaving only 4 of them "read_roms_grid",
! "write_roms_grid", "roms_find_dims", and "roms_check_dims" be
! strictly ROMS specific).

! It should be noted that somewhat similar functionality for reading
! and writing netCDF files can be found in "nc_read_write.F", however
! the distinction between there and the routines in this package is
! that all the ones below having argument "ncid" are expected to have
! netCDF file in open state, while access to a specific variable is
! done by name, hence it is "file-by-ID -- var-by-name" semantics,
! while "nc_read_write.F" uses "file-by-name -- var-by-name". For this
! reason argument "ncid" is always placed before the filename, and the
! latter is used only to write error messages, but has no effect other
! than that.

! Another note is that FORTRAN 2003 standard mandates that the rank
! of argument (scalar vs. array) should be the same for both calling
! routine and the callee, even in the trivial case where the array
! consists of just a single element. Thus, it is formally illegal
! (thought works correctly in practice) to pass a scalar as an
! argument to a routine expecting an array of size 1. The fact
! that size is equal to 1 in known only during runtime, but not at
! compiling, so the compiler instrumented to verify F2003 compliance
! issues an error message and quits. It is for this and only this
! reason routines containing _sclr_ in their names in the list below
! were introduced, even thought their functionality may seem to be
! redundant (below "value" is scalar, while "var" is array).

! The content is:

! init_time(ncid, fname, tname, nrecs, init_year, ierr)

! read_roms_grid (fname, Lm,Mm)
! write_roms_grid (fname, Lm,Mm)
! roms_find_dims (ncid, fname, Lm,Mm,N)
! roms_check_dims (ncid, fname, Lm,Mm,N)

! put_sclr_by_name_real (ncid, vname, value)
! get_sclr_by_name_real (ncid, vname, value)
! put_sclr_by_name_double (ncid, vname, value)
! get_sclr_by_name_double (ncid, vname, value)

! put_var_by_name_real (ncid, vname, var)
! get_var_by_name_real (ncid, vname, var)
! put_var_by_name_double (ncid, vname, var)
! get_var_by_name_double (ncid, vname, var)

! put_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! put_sclr_rec_by_name_double (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_double (ncid, fname, vname, rec, value)

! put_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! put_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)

! put_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! put_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)

! All are subroutines designed to provide sufficient diagnostic
! messages and terminate the execution if something goes wrong
! rather than functions returning non-zero status.

! With the exception of "read_roms_grid" which opens the named file,
! creates its netCDF file ID as an internal variable, reads all the
! relevant data, and closes it after that, all the above procedures
! imply that the file is in opened state, hence input argument "ncid"
! has meaningful value at entry, while argument "fname" is used only
! for error messages id something goes wrong.
# 783 "roms_read_write.F"
! The rest are the reading-writting subroutins generated by CPP
! from the same source code package (in fact, quadrupled: X2 due
! to read/write functionality and another X2 due to single/double
! precision version).
# 844 "roms_read_write.F"
! The following eight routines are just instrumented wrappers around
! the standard sequence of netCDF calls which (1) inquire variable ID
! and (2) put/get the ENTIRE variable into/from netCDF file. These
! wrapper is needed solely to write error messages if something goes
! wrong. These are
!
! put/get_sclr/var_by_name_TYPE (ncid, vname, value/var)
!
! where get/put and TYPE=real/double occur in all permutations (hence
! it adds up to a total of eight). Because of semantically identical
! code real/double is implemented by CPP-redefinition of basic netCDF
! functions using the same source code.


      subroutine put_sclr_by_name_real(ncid, vname, value)




! These two routines are for reading or writing just a single number
! which may exist in netCDF file either as a variable or a global
! attribute containing just a single number of the proper type.
! Selection between variable or attribute is by the file, while no
! attempt to change is format is made here. For this reason is either
! variable of attribute with given name must pre-exist in order for
! these operations to succeed.

      implicit none
      integer ncid, varid, type, size, ierr, lvar
      character(len=*) vname
      real(kind=4) value
      include "netcdf.inc"

      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname(1:lvar), varid)
      if (ierr == nf_noerr) then

        ierr=nf_put_var_real(ncid, varid, value)
        if (ierr == nf_noerr) then
          write(*,'(8x,3A)') 'wrote ''', vname(1:lvar), ''''
        else
          write(*,'(/1x,4A,1x,A/)') '### ERROR: Cannot write ',
     & 'netCDF variable ''', vname, ''':', nf_strerror(ierr)
        endif
# 897 "roms_read_write.F"
      else
        ierr=nf_inq_att(ncid, nf_global, vname(1:lvar), type, size)
        if (ierr == nf_noerr) then




          if (size == 1 .and. type == nf_real) then


            ierr=nf_redef(ncid)
            if (ierr == nf_noerr) then
              ierr=nf_put_att_real(ncid, nf_global, vname(1:lvar),
     & value)
              if (ierr == nf_noerr) then
                ierr=nf_enddef(ncid)
                if (ierr == nf_noerr) then
                  write(*,'(8x,3A)') 'wrote ''', vname(1:lvar),
     & ''' as global attribute'
                else
                  write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot ',
     & 'close redefinition mode for netCDF file.',
     & nf_strerror(ierr)
                endif
              else
                write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot put ',
     & 'global attribute ''', vname(1:lvar),
     & ''' into netCDF file.', nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot switch ',
     & 'netCDF file into redefinition mode.', nf_strerror(ierr)
            endif
# 942 "roms_read_write.F"
          else
            if (size /= 1) then
              write(*,'(/1x,5A,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong size', size, 'instead of 1.'
              ierr=ierr-1
            endif



            if (type /= nf_real) then

              write(*,'(/1x,5A,I4,1x,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong type', type, 'instead of',



     & nf_real, 'which is real.'

              ierr=ierr-1
            endif
          endif
        else
          write(*,'(/1x,4A/)') '### ERROR: Neither variable, nor ',
     & 'global attribute named ''', vname,
     & ''' is present in netCDF file.'
        endif
      endif
      if (ierr /= nf_noerr) stop
      end






      subroutine put_var_by_name_real(ncid, vname, var)



      implicit none
      integer ncid, ierr, varid, lvar
      character(len=*) vname
      real(kind=4) var(*)
      include "netcdf.inc"
      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then

        ierr=nf_put_var_real(ncid, varid, var)
        if (ierr == nf_noerr) then
          write(*,'(8x,3A)') 'wrote ''', vname(1:lvar), ''''





        else

          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot write ',



     & 'netCDF variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot find ',
     & 'netCDF variable ID for ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
      endif
      if (ierr /= nf_noerr) stop
      end


! The following routines are to put/get just a single number into
! a specified record of one-dimensional netCDF variable. The netCDF
! file "fname" is expected to be open (hence input argument "ncid" is
! a valid file ID, while the name of the file is needed only to write
! error messages if something goes wrong, but otherwise is not used);
! "vname" is name of the variable (to be translated inside into netCDF
! ID with error message if not found), and "rec" is record number;
! "var" is input for put_ (output for get_) is just a scalar (single
! number). The corresponding netCDF variable is expected to be either
! a one-dimensional array (having more dimensions results in error
! message) or a scalar. If array it puts/gets the value at location
! "rec" with performing necessary error checking; if scalar it takes
! the only value, while argument "rec" is not used. Again, get/put
! and TYPE=real/double occur in all four permutations.



      subroutine put_sclr_rec_by_name_real(ncid, fname, vname,
     & rec, value)




      implicit none
      integer ncid, rec
      character(len=*) fname, vname
      real(kind=4) value

      character(len=16) name
      integer varid, vtype, ndims, natts, dimid(8), size,
     & start(4), count(4), ierr, lfnm, lvar
      include "netcdf.inc"

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then
          if (ndims == 1) then
            ierr=nf_inq_dimlen(ncid, dimid(1), size)
            if (ierr == nf_noerr) then
              start(1)=rec ; count(1)=1

              ierr=nf_put_vara_real(ncid,varid, start,count, value)
              if (ierr == nf_noerr) then
                write(*,'(5x,A,I5,1x,5A)') 'wrote rec', rec,
     & 'of scalar ''', vname(1:lvar), ''' into ''',
     & fname(1:lfnm), '''.'
                return !---> successful return
              else
                write(*,'(/1x,5A/12x,A)') '### ERROR: Cannot write ',
     & 'variable ''', vname(1:lvar), ''' into ''',
     & fname(1:lfnm),''':', nf_strerror(ierr)
              endif
# 1093 "roms_read_write.F"
            else
              write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     & 'Cannot determine size of dimension #', dimid(1),
     & 'in file ''', fname(1:lfnm), '''.', nf_strerror(ierr)
            endif
          elseif (ndims == 0) then


            ierr=nf_put_var_real(ncid, varid, value)
            if (ierr == nf_noerr) then
              write(*,'(6x,5A)') 'wrote scalar ''', vname(1:lvar),
     & ''' into ''', fname(1:lfnm), '''.'
              return !---> successful return
            else
              write(*,'(/1x,5A/12x,A)') '### ERROR: Cannot write ',
     & 'scalar variable ''', vname(1:lvar), ''' into ''',
     & fname(1:lfnm),''':', nf_strerror(ierr)
            endif
# 1123 "roms_read_write.F"
          else
            write(*,'(/1x,5A,I4/)') '### ERROR: Variable ''',
     & vname(1:lvar), ''' from file ''', fname(1:lfnm),
     & ''' has more then one dimension, ndims =', ndims
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! A more sophisticated functions, which read/write a portion of data
! into/from an existing netCDF variable accessing it by name. In doing
! so, it checks whether the first one, two, or three dimensions for
! that variable as defined in the netCDF file are consistent with the
! shape of array "var" specified as n1,n2,n3. If not, it complains
! about the error and quits. The variable "var" may be one, two, or
! three-dimensional, and, in addition to that may have record dimension
! (not necessarily unlimited). If the variable is two dimensional,
! then the calling routine must be specify n3=0 to avoid checking of
! the non-existing dimension. Similar policy applies for n1 and n2
! (e.g., n1,n2,n3=0,0,0 means that the variable is a scalar). However,
! if the variable in netCDF file has record dimension (identified as
! either unlimited netCDF dimension, or as the last dimension of the
! variable with its name ending with "...time" AND the actual number
! of spatial dimensions of netCDF variable in the file is LESS that
! the number of non-zero arguments THEN the excess spatial dimensions
! (n3, or both n3 and n2) will be ignored, while "rec" will be
! interpreted normally.

! Note: in the code below arguments n1,n2,n3 are used EXCLUSIVELY for
! checking, while the starts and counts are computed from the the
! dimensions of the variable in netCDF file.


      subroutine put_rec_by_name_real(ncid, fname, vname,
     & n1,n2,n3, rec, var)




! use mod_io_size_acct
      implicit none
      integer ncid, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=4), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0 ! Determine the number of spatial
        count(id)=0 ! dimensions "nspc" for the variable.
      enddo
      nspc=0 ! Note that "nspc" found here should be
      if (n1 > 0) then ! either
        nspc=nspc+1 ! less by 1 than the actual number of
        count(nspc)=n1 ! dimensions "ndims" of the variable
      endif ! stored in the file -- in this case
      if (n2 > 0) then ! the extra dimension will be treated
        nspc=nspc+1 ! as record dimension;
        count(nspc)=n2 ! or
      endif ! be the same as "ndims" -- in this
      if (n3 > 0) then ! case there is no record dimension,
        nspc=nspc+1 ! and argument "rec" is ignored.
        count(nspc)=n3
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then






          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1244 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  start(id)=1 ; rec_size=rec_size*size

                  if (count(id) /= size) then
                    call lenstr(name,ldim)
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',




     & 'put_rec_by_name_real',
# 1274 "roms_read_write.F"
     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1331 "roms_read_write.F"
              endif

              ierr=nf_put_vara_real(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                write(*,'(5x,A,I5,1x,5A)') 'wrote rec', rec, 'of ''',
     & vname(1:lvar), ''' into ''', fname(1:lfnm), '''.'
                sz_write_acc = sz_write_acc + rec_size * 4
# 1346 "roms_read_write.F"
                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max

                write_clk = write_clk + inc_clk




                return !---> successful return
              else

                write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot write ',
     & 'variable ''',vname(1:lvar),''' into netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)





              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! The following pair is essentially an instrumented versions of
! put_vara and get_vara to read a single record of a 2D-subdomain
! [iwest : iwest+n1-1] x [jsouth : jsouth+n2-1]
! within the (i,j)-index-space of netCDF array. The third dimension
! (if exists) is treated as the whole: n3 is checked against the
! actual dimension of netCDF variable and the mismatch is treated
! as an error; n1,n2 are checked only for upper-bound overrun, e.g.,
! iwest+n1-1 exceeds the actual size of the first dimension
! resulting in error.


      subroutine put_patch_by_name_real(ncid, fname, vname,
     & iwest,jsouth, n1,n2,n3, rec, var)




c--#define VERBOSE
! use mod_io_size_acct
      implicit none
      integer ncid, iwest,jsouth, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=4), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0
      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0
        count(id)=0
      enddo ! Determine the number of spatial
      nspc=0 ! dimensions "nspc" for the variable.
      if (n1 > 0) then
        nspc=nspc+1 ! Note that "nspc" found here should be
        start(nspc)=iwest !
        count(nspc)=n1 ! either
      endif ! less by 1 than the actual number of
      if (n2 > 0) then ! dimensions "ndims" of the variable
        nspc=nspc+1 ! stored in the file -- in this case
        start(nspc)=jsouth ! the extra dimension will be treated
        count(nspc)=n2 ! as record dimension,
      endif ! or
      if (n3 > 0) then ! be the same as "ndims" -- in this
        nspc=nspc+1 ! case there is no record dimension,
        count(nspc)=n3 ! and argument "rec" is ignored.
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then





          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1478 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  rec_size=rec_size*count(id)
                  if (start(id) == 0) then
                    start(id)=1
                    if (count(id) /= size) then
                      write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',




     & 'put_patch_by_name_real',
# 1508 "roms_read_write.F"
     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                      stop
                    endif
                  elseif (start(id)+count(id)-1 > size) then
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,3(I5,1x,A)/)')
     & '### ERROR: ',




     & 'put_patch_by_name_real',
# 1530 "roms_read_write.F"
     & ' :: Overrun dimension bound #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', start(id), '+',
     & count(id), '-1 >', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1587 "roms_read_write.F"
              endif

              ierr=nf_put_vara_real(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                if (start(ndims) == rec .and. count(ndims) == 1) then
                  write(*,'(5x,A,I5,1x,5A)') 'wrote rec',rec,'of ''',
     & vname(1:lvar), ''' into ''', fname(1:lfnm), ''''
                else
                  write(*,'(5x,5A)') 'wrote''', vname(1:lvar),
     & ''' into ''', fname(1:lfnm), ''''
                endif
                sz_write_acc = sz_write_acc + rec_size * 4
# 1612 "roms_read_write.F"
                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max

                write_clk = write_clk + inc_clk




                return !---> successful return
              else

                write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot write ',
     & 'variable ''',vname(1:lvar),''' into netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)





              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end
# 1657 "roms_read_write.F" 2


# 1 "roms_read_write.F" 1
! This package contains a complete set of basic operators for reading
! and writing ROMS-style netCDF data files. The "standard" implies that
! horizontal dimensions are named as "xi_","eta_", vertical "s_" with
! corresponding suffix, "rho", "u", and "v" for horizontal dimensions;
! "rho" and "w" for vertical consistently with grid staggering rules
! within ROMS code. Time dimension (whether or not is "unlimited" from
! netCDF point of view) has its name ending with "time".

! Other than the spatial grid staggering rules, all other aspects
! related to netCDF file structure are expected to follow so called
! "CF conventions" as closely as possible. As the result, all get_*
! and put_* routines from this package are known to work for files
! other than ROMS-standard (leaving only 4 of them "read_roms_grid",
! "write_roms_grid", "roms_find_dims", and "roms_check_dims" be
! strictly ROMS specific).

! It should be noted that somewhat similar functionality for reading
! and writing netCDF files can be found in "nc_read_write.F", however
! the distinction between there and the routines in this package is
! that all the ones below having argument "ncid" are expected to have
! netCDF file in open state, while access to a specific variable is
! done by name, hence it is "file-by-ID -- var-by-name" semantics,
! while "nc_read_write.F" uses "file-by-name -- var-by-name". For this
! reason argument "ncid" is always placed before the filename, and the
! latter is used only to write error messages, but has no effect other
! than that.

! Another note is that FORTRAN 2003 standard mandates that the rank
! of argument (scalar vs. array) should be the same for both calling
! routine and the callee, even in the trivial case where the array
! consists of just a single element. Thus, it is formally illegal
! (thought works correctly in practice) to pass a scalar as an
! argument to a routine expecting an array of size 1. The fact
! that size is equal to 1 in known only during runtime, but not at
! compiling, so the compiler instrumented to verify F2003 compliance
! issues an error message and quits. It is for this and only this
! reason routines containing _sclr_ in their names in the list below
! were introduced, even thought their functionality may seem to be
! redundant (below "value" is scalar, while "var" is array).

! The content is:

! init_time(ncid, fname, tname, nrecs, init_year, ierr)

! read_roms_grid (fname, Lm,Mm)
! write_roms_grid (fname, Lm,Mm)
! roms_find_dims (ncid, fname, Lm,Mm,N)
! roms_check_dims (ncid, fname, Lm,Mm,N)

! put_sclr_by_name_real (ncid, vname, value)
! get_sclr_by_name_real (ncid, vname, value)
! put_sclr_by_name_double (ncid, vname, value)
! get_sclr_by_name_double (ncid, vname, value)

! put_var_by_name_real (ncid, vname, var)
! get_var_by_name_real (ncid, vname, var)
! put_var_by_name_double (ncid, vname, var)
! get_var_by_name_double (ncid, vname, var)

! put_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! put_sclr_rec_by_name_double (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_double (ncid, fname, vname, rec, value)

! put_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! put_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)

! put_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! put_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)

! All are subroutines designed to provide sufficient diagnostic
! messages and terminate the execution if something goes wrong
! rather than functions returning non-zero status.

! With the exception of "read_roms_grid" which opens the named file,
! creates its netCDF file ID as an internal variable, reads all the
! relevant data, and closes it after that, all the above procedures
! imply that the file is in opened state, hence input argument "ncid"
! has meaningful value at entry, while argument "fname" is used only
! for error messages id something goes wrong.
# 645 "roms_read_write.F"
      subroutine roms_check_dims(ncgrd, fname, Lm_ck, Mm_ck, N_ck)



      implicit none
      character(len=*) fname
      integer ncgrd, Lm,Mm,N, xi_rho,xi_u, eta_rho,eta_v, s_rho,s_w,
     & ndims, size, id, i,is, ierr, lvar, lfnm

     & , Lm_ck, Mm_ck, N_ck

      character(len=16) dname
      character(len=128) string

      include "netcdf.inc"

      xi_rho=0 ; xi_u=0 ; s_rho=0
      eta_rho=0 ; eta_v=0 ; s_w=0

      call lenstr(fname,lfnm)
      ierr=nf_inq_ndims(ncgrd, ndims)
      if (ierr == nf_noerr) then
        do id=1,ndims
          dname='                '
          ierr=nf_inq_dim (ncgrd, id, dname, size)
          if (ierr == nf_noerr) then
            call lenstr(dname,lvar)
            if (lvar == 6 .and. dname(1:lvar) == 'xi_rho') then
              xi_rho=size
            elseif (lvar == 4 .and. dname(1:lvar) == 'xi_u') then
              xi_u=size

            elseif (lvar == 7 .and. dname(1:lvar) == 'eta_rho') then
              eta_rho=size
            elseif (lvar == 5 .and. dname(1:lvar) == 'eta_v') then
              eta_v=size

            elseif (lvar == 5 .and. dname(1:lvar) == 's_rho') then
              s_rho=size
            elseif (lvar == 3 .and. dname(1:lvar) == 's_w') then
              s_w=size

            elseif (lvar == 5 .and. dname(1:lvar) == 'depth') then
              s_rho=size
            elseif (lvar == 7 .and. dname(1:lvar) == 'rho_ntr') then
              s_rho=size
            endif
          else
            write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     & 'determine name and size of dimension #', id,
     & 'in ''', fname(1:lfnm), '''.', nf_strerror(ierr)
          endif
        enddo

        write(string,'(A,6(1x,A,I4))')



     & 'roms_check_dims ::',

     & 'xi_rho=',xi_rho, 'xi_u=',xi_u, 'eta_rho=',eta_rho,
     & 'eta_v=',eta_v, 's_rho=',s_rho, 's_w=',s_w
        call lenstr(string,lvar)
        i=0 ! Write dimensions into
        do while(i < lvar) ! character string first,
          i=i+1 ! and then suppress blank
          if (string(i:i) == '=') then ! characters after =sign.
            i=i+1 ! This is merely to make
            if (string(i:i) == ' ') then
              is=1
              do while(string(i+is:i+is) == ' ' .and. i+is < lvar)
                is=is+1
              enddo
              string(i:lvar-is)=string(i+is:lvar) ; lvar=lvar-is
            endif
          endif
        enddo ! a narrower printout
        write(*,'(2x,A)') string(1:lvar) ! on the screen.

        ierr=0
        if (xi_rho > 0) then
          Lm=xi_rho-2
        elseif (xi_u > 0) then
          Lm=xi_u-1
        else
          write(*,'(/1x,4A/)') '### ERROR: Cannot determine size ',
     & 'of horizontal XI-dimension in netCDF file ''',
     & fname(1:lfnm), '''.'
          ierr=ierr+1
        endif
        if (eta_rho > 0) then
          Mm=eta_rho-2
        elseif (eta_v > 0) then
          Mm=eta_rho-1
        else
          write(*,'(/1x,4A/)') '### ERROR: Cannot determine size ',
     & 'of horizontal ETA-dimension in netCDF file ''',
     & fname(1:lfnm), '''.'
          ierr=ierr+1
        endif ! The policy here is that vertical

        N=0 !<-- initialize ! dimension is optional, therefore

        if (s_rho > 0) then ! it is filled up only if found, and
          N=s_rho ! "not touched" otherwise (hence it
        elseif (s_w > 0) then ! is possible to call this function
          N=s_w-1 ! while passing a constant, if no
        endif ! vertical dimension is expected to
                                   ! exist in the file.

        if (Lm /= Lm_ck) then
          write(*,'(/1x,2A,I4,1x,3A,I4/)') '### ERROR: Size of XI-',
     & 'dimension Lm =', Lm, 'from file ''', fname(1:lfnm),
     & ''' does not match the previous size ', Lm_ck
          ierr=ierr+1
        endif
        if (Mm /= Mm_ck) then
          write(*,'(/1x,2A,I4,1x,3A,I4/)') '### ERROR: Size of ETA-',
     & 'dimension Mm =', Mm, 'from file ''', fname(1:lfnm),
     & ''' does not match the previous size ', Mm_ck
          ierr=ierr+1
        endif
        if (N > 0 .and. N_ck > 0 .and. N /= N_ck) then
          write(*,'(/1x,2A,I4,1x,3A,I4/)') '### ERROR: Size of ',
     & 'vertical dimension N =', N, 'from file ''', fname(1:lfnm),
     & ''' does not match the previous size ', N_ck
          ierr=ierr+1
        endif

        if (ierr /= 0) stop !--> ERROR
      else
        write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot determine ',
     & 'number of dimensions in netCDF file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      end


! The rest are the reading-writting subroutins generated by CPP
! from the same source code package (in fact, quadrupled: X2 due
! to read/write functionality and another X2 due to single/double
! precision version).
# 844 "roms_read_write.F"
! The following eight routines are just instrumented wrappers around
! the standard sequence of netCDF calls which (1) inquire variable ID
! and (2) put/get the ENTIRE variable into/from netCDF file. These
! wrapper is needed solely to write error messages if something goes
! wrong. These are
!
! put/get_sclr/var_by_name_TYPE (ncid, vname, value/var)
!
! where get/put and TYPE=real/double occur in all permutations (hence
! it adds up to a total of eight). Because of semantically identical
! code real/double is implemented by CPP-redefinition of basic netCDF
! functions using the same source code.




      subroutine get_sclr_by_name_double(ncid, vname, value)


! These two routines are for reading or writing just a single number
! which may exist in netCDF file either as a variable or a global
! attribute containing just a single number of the proper type.
! Selection between variable or attribute is by the file, while no
! attempt to change is format is made here. For this reason is either
! variable of attribute with given name must pre-exist in order for
! these operations to succeed.

      implicit none
      integer ncid, varid, type, size, ierr, lvar
      character(len=*) vname
      real(kind=8) value
      include "netcdf.inc"

      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname(1:lvar), varid)
      if (ierr == nf_noerr) then
# 889 "roms_read_write.F"
        ierr=nf_get_var_double(ncid, varid, value)
        if (ierr == nf_noerr) then
          write(*,'(9x,3A)') 'read ''', vname(1:lvar), ''''
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read netCDF ',
     & 'variable ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif

      else
        ierr=nf_inq_att(ncid, nf_global, vname(1:lvar), type, size)
        if (ierr == nf_noerr) then


          if (size == 1 .and. type == nf_double) then
# 931 "roms_read_write.F"
            ierr=nf_get_att_double(ncid, nf_global, vname(1:lvar),
     & value)
            if (ierr == nf_noerr) then
              write(*,'(9x,4A)') 'read ''', vname(1:lvar),
     & ''' as global attribute'
            else
              write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read ',
     & 'global attribute ''', vname(1:lvar),
     & ''' from netCDF file.', nf_strerror(ierr)
            endif

          else
            if (size /= 1) then
              write(*,'(/1x,5A,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong size', size, 'instead of 1.'
              ierr=ierr-1
            endif

            if (type /= nf_double) then



              write(*,'(/1x,5A,I4,1x,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong type', type, 'instead of',

     & nf_double, 'which is double precision.'



              ierr=ierr-1
            endif
          endif
        else
          write(*,'(/1x,4A/)') '### ERROR: Neither variable, nor ',
     & 'global attribute named ''', vname,
     & ''' is present in netCDF file.'
        endif
      endif
      if (ierr /= nf_noerr) stop
      end
# 981 "roms_read_write.F"
      subroutine get_var_by_name_double(ncid, vname, var)

      implicit none
      integer ncid, ierr, varid, lvar
      character(len=*) vname
      real(kind=8) var(*)
      include "netcdf.inc"
      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then





        ierr=nf_get_var_double(ncid, varid, var)
        if (ierr == nf_noerr) then
          write(*,'(9x,3A)') 'read ''', vname(1:lvar), ''''

        else



          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read ',

     & 'netCDF variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot find ',
     & 'netCDF variable ID for ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
      endif
      if (ierr /= nf_noerr) stop
      end


! The following routines are to put/get just a single number into
! a specified record of one-dimensional netCDF variable. The netCDF
! file "fname" is expected to be open (hence input argument "ncid" is
! a valid file ID, while the name of the file is needed only to write
! error messages if something goes wrong, but otherwise is not used);
! "vname" is name of the variable (to be translated inside into netCDF
! ID with error message if not found), and "rec" is record number;
! "var" is input for put_ (output for get_) is just a scalar (single
! number). The corresponding netCDF variable is expected to be either
! a one-dimensional array (having more dimensions results in error
! message) or a scalar. If array it puts/gets the value at location
! "rec" with performing necessary error checking; if scalar it takes
! the only value, while argument "rec" is not used. Again, get/put
! and TYPE=real/double occur in all four permutations.






      subroutine get_sclr_rec_by_name_double(ncid, fname, vname,
     & rec, value)

      implicit none
      integer ncid, rec
      character(len=*) fname, vname
      real(kind=8) value

      character(len=16) name
      integer varid, vtype, ndims, natts, dimid(8), size,
     & start(4), count(4), ierr, lfnm, lvar
      include "netcdf.inc"

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then
          if (ndims == 1) then
            ierr=nf_inq_dimlen(ncid, dimid(1), size)
            if (ierr == nf_noerr) then
              start(1)=rec ; count(1)=1
# 1074 "roms_read_write.F"
              if (0 < rec .and. rec <= size) then
                ierr=nf_get_vara_double(ncid,varid, start,count, value)
                if (ierr == nf_noerr) then
                  write(*,'(6x,A,I5,1x,5A)') 'read rec', rec,
     & 'of scalar ''', vname(1:lvar), ''' from ''',
     & fname(1:lfnm), '''.'
                  return !---> successful return
                else
                  write(*,'(/1x,7A/12x,A)') '### ERROR: Cannot read ',
     & 'variable ''', vname(1:lvar), ''' from netCDF ',
     & 'file ''',fname(1:lfnm), ''':', nf_strerror(ierr)
                endif
              else
                write(*,'(/1x,2A,I4,1x,6A,I4/)') '### ERROR: ',
     & 'Requested record number ', rec, 'for scalar ',
     & 'variable ''', vname(1:lvar), ''' in file ''',
     & fname(1:lfnm), ''' exceeds dimension bound', size
              endif

            else
              write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     & 'Cannot determine size of dimension #', dimid(1),
     & 'in file ''', fname(1:lfnm), '''.', nf_strerror(ierr)
            endif
          elseif (ndims == 0) then
# 1112 "roms_read_write.F"
            ierr=nf_get_var_double(ncid, varid, value)
            if (ierr == nf_noerr) then
              write(*,'(7x,5A)') 'read scalar ''', vname(1:lvar),
     & ''' from ''', fname(1:lfnm), '''.'
              return !---> successful return
            else
              write(*,'(/1x,5A/12x,A)') '### ERROR: Cannot read ',
     & 'scalar variable ''', vname(1:lvar), ''' from ''',
     & fname(1:lfnm),''':', nf_strerror(ierr)
            endif

          else
            write(*,'(/1x,5A,I4/)') '### ERROR: Variable ''',
     & vname(1:lvar), ''' from file ''', fname(1:lfnm),
     & ''' has more then one dimension, ndims =', ndims
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! A more sophisticated functions, which read/write a portion of data
! into/from an existing netCDF variable accessing it by name. In doing
! so, it checks whether the first one, two, or three dimensions for
! that variable as defined in the netCDF file are consistent with the
! shape of array "var" specified as n1,n2,n3. If not, it complains
! about the error and quits. The variable "var" may be one, two, or
! three-dimensional, and, in addition to that may have record dimension
! (not necessarily unlimited). If the variable is two dimensional,
! then the calling routine must be specify n3=0 to avoid checking of
! the non-existing dimension. Similar policy applies for n1 and n2
! (e.g., n1,n2,n3=0,0,0 means that the variable is a scalar). However,
! if the variable in netCDF file has record dimension (identified as
! either unlimited netCDF dimension, or as the last dimension of the
! variable with its name ending with "...time" AND the actual number
! of spatial dimensions of netCDF variable in the file is LESS that
! the number of non-zero arguments THEN the excess spatial dimensions
! (n3, or both n3 and n2) will be ignored, while "rec" will be
! interpreted normally.

! Note: in the code below arguments n1,n2,n3 are used EXCLUSIVELY for
! checking, while the starts and counts are computed from the the
! dimensions of the variable in netCDF file.





      subroutine get_rec_by_name_double(ncid, fname, vname,
     & n1,n2,n3, rec, var)

! use mod_io_size_acct
      implicit none
      integer ncid, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=8), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0 ! Determine the number of spatial
        count(id)=0 ! dimensions "nspc" for the variable.
      enddo
      nspc=0 ! Note that "nspc" found here should be
      if (n1 > 0) then ! either
        nspc=nspc+1 ! less by 1 than the actual number of
        count(nspc)=n1 ! dimensions "ndims" of the variable
      endif ! stored in the file -- in this case
      if (n2 > 0) then ! the extra dimension will be treated
        nspc=nspc+1 ! as record dimension;
        count(nspc)=n2 ! or
      endif ! be the same as "ndims" -- in this
      if (n3 > 0) then ! case there is no record dimension,
        nspc=nspc+1 ! and argument "rec" is ignored.
        count(nspc)=n3
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then






          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1244 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  start(id)=1 ; rec_size=rec_size*size

                  if (count(id) /= size) then
                    call lenstr(name,ldim)
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',
# 1269 "roms_read_write.F"
     & 'get_rec_by_name_double',




     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1331 "roms_read_write.F"
              endif







              ierr=nf_get_vara_double(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                write(*,'(6x,A,I5,1x,5A)') 'read rec', rec, 'of ''',
     & vname(1:lvar), ''' from ''', fname(1:lfnm), '''.'
                sz_read_acc = sz_read_acc + rec_size * 8


                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max



                read_clk = read_clk + inc_clk


                return !---> successful return
              else





                write(*,'(/1x,7A/12x,A)') '### ERROR: Cannot read ',
     & 'variable ''',vname(1:lvar),''' from netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)

              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! The following pair is essentially an instrumented versions of
! put_vara and get_vara to read a single record of a 2D-subdomain
! [iwest : iwest+n1-1] x [jsouth : jsouth+n2-1]
! within the (i,j)-index-space of netCDF array. The third dimension
! (if exists) is treated as the whole: n3 is checked against the
! actual dimension of netCDF variable and the mismatch is treated
! as an error; n1,n2 are checked only for upper-bound overrun, e.g.,
! iwest+n1-1 exceeds the actual size of the first dimension
! resulting in error.





      subroutine get_patch_by_name_double(ncid, fname, vname,
     & iwest,jsouth, n1,n2,n3, rec, var)

c--#define VERBOSE
! use mod_io_size_acct
      implicit none
      integer ncid, iwest,jsouth, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=8), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0
      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0
        count(id)=0
      enddo ! Determine the number of spatial
      nspc=0 ! dimensions "nspc" for the variable.
      if (n1 > 0) then
        nspc=nspc+1 ! Note that "nspc" found here should be
        start(nspc)=iwest !
        count(nspc)=n1 ! either
      endif ! less by 1 than the actual number of
      if (n2 > 0) then ! dimensions "ndims" of the variable
        nspc=nspc+1 ! stored in the file -- in this case
        start(nspc)=jsouth ! the extra dimension will be treated
        count(nspc)=n2 ! as record dimension,
      endif ! or
      if (n3 > 0) then ! be the same as "ndims" -- in this
        nspc=nspc+1 ! case there is no record dimension,
        count(nspc)=n3 ! and argument "rec" is ignored.
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then





          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1478 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  rec_size=rec_size*count(id)
                  if (start(id) == 0) then
                    start(id)=1
                    if (count(id) /= size) then
                      write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',
# 1503 "roms_read_write.F"
     & 'get_patch_by_name_double',




     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                      stop
                    endif
                  elseif (start(id)+count(id)-1 > size) then
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,3(I5,1x,A)/)')
     & '### ERROR: ',
# 1525 "roms_read_write.F"
     & 'get_patch_by_name_double',




     & ' :: Overrun dimension bound #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', start(id), '+',
     & count(id), '-1 >', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1587 "roms_read_write.F"
              endif
# 1600 "roms_read_write.F"
              ierr=nf_get_vara_double(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                if (start(ndims) == rec .and. count(ndims) == 1) then
                  write(*,'(6x,A,I5,1x,5A)') 'read rec', rec,'of ''',
     & vname(1:lvar), ''' from ''', fname(1:lfnm), ''''
                else
                  write(*,'(6x,5A)') 'read ''', vname(1:lvar),
     & ''' from ''', fname(1:lfnm), ''''
                endif
                sz_read_acc = sz_read_acc + rec_size * 8


                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max



                read_clk = read_clk + inc_clk


                return !---> successful return
              else





                write(*,'(/1x,7A/12x,A)') '### ERROR: Cannot read ',
     & 'variable ''',vname(1:lvar),''' from netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)

              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end
# 1660 "roms_read_write.F" 2

# 1 "roms_read_write.F" 1
! This package contains a complete set of basic operators for reading
! and writing ROMS-style netCDF data files. The "standard" implies that
! horizontal dimensions are named as "xi_","eta_", vertical "s_" with
! corresponding suffix, "rho", "u", and "v" for horizontal dimensions;
! "rho" and "w" for vertical consistently with grid staggering rules
! within ROMS code. Time dimension (whether or not is "unlimited" from
! netCDF point of view) has its name ending with "time".

! Other than the spatial grid staggering rules, all other aspects
! related to netCDF file structure are expected to follow so called
! "CF conventions" as closely as possible. As the result, all get_*
! and put_* routines from this package are known to work for files
! other than ROMS-standard (leaving only 4 of them "read_roms_grid",
! "write_roms_grid", "roms_find_dims", and "roms_check_dims" be
! strictly ROMS specific).

! It should be noted that somewhat similar functionality for reading
! and writing netCDF files can be found in "nc_read_write.F", however
! the distinction between there and the routines in this package is
! that all the ones below having argument "ncid" are expected to have
! netCDF file in open state, while access to a specific variable is
! done by name, hence it is "file-by-ID -- var-by-name" semantics,
! while "nc_read_write.F" uses "file-by-name -- var-by-name". For this
! reason argument "ncid" is always placed before the filename, and the
! latter is used only to write error messages, but has no effect other
! than that.

! Another note is that FORTRAN 2003 standard mandates that the rank
! of argument (scalar vs. array) should be the same for both calling
! routine and the callee, even in the trivial case where the array
! consists of just a single element. Thus, it is formally illegal
! (thought works correctly in practice) to pass a scalar as an
! argument to a routine expecting an array of size 1. The fact
! that size is equal to 1 in known only during runtime, but not at
! compiling, so the compiler instrumented to verify F2003 compliance
! issues an error message and quits. It is for this and only this
! reason routines containing _sclr_ in their names in the list below
! were introduced, even thought their functionality may seem to be
! redundant (below "value" is scalar, while "var" is array).

! The content is:

! init_time(ncid, fname, tname, nrecs, init_year, ierr)

! read_roms_grid (fname, Lm,Mm)
! write_roms_grid (fname, Lm,Mm)
! roms_find_dims (ncid, fname, Lm,Mm,N)
! roms_check_dims (ncid, fname, Lm,Mm,N)

! put_sclr_by_name_real (ncid, vname, value)
! get_sclr_by_name_real (ncid, vname, value)
! put_sclr_by_name_double (ncid, vname, value)
! get_sclr_by_name_double (ncid, vname, value)

! put_var_by_name_real (ncid, vname, var)
! get_var_by_name_real (ncid, vname, var)
! put_var_by_name_double (ncid, vname, var)
! get_var_by_name_double (ncid, vname, var)

! put_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_real (ncid, fname, vname, rec, value)
! put_sclr_rec_by_name_double (ncid, fname, vname, rec, value)
! get_sclr_rec_by_name_double (ncid, fname, vname, rec, value)

! put_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_real (ncid, fname, vname, n1,n2,n3,rec, var)
! put_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)
! get_rec_by_name_double (ncid, fname, vname, n1,n2,n3,rec, var)

! put_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_real (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! put_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)
! get_patch_by_name_double (ncid, fname, vname, iwest,jsouth,
! n1,n2,n3,rec, var)

! All are subroutines designed to provide sufficient diagnostic
! messages and terminate the execution if something goes wrong
! rather than functions returning non-zero status.

! With the exception of "read_roms_grid" which opens the named file,
! creates its netCDF file ID as an internal variable, reads all the
! relevant data, and closes it after that, all the above procedures
! imply that the file is in opened state, hence input argument "ncid"
! has meaningful value at entry, while argument "fname" is used only
! for error messages id something goes wrong.
# 783 "roms_read_write.F"
! The rest are the reading-writting subroutins generated by CPP
! from the same source code package (in fact, quadrupled: X2 due
! to read/write functionality and another X2 due to single/double
! precision version).
# 844 "roms_read_write.F"
! The following eight routines are just instrumented wrappers around
! the standard sequence of netCDF calls which (1) inquire variable ID
! and (2) put/get the ENTIRE variable into/from netCDF file. These
! wrapper is needed solely to write error messages if something goes
! wrong. These are
!
! put/get_sclr/var_by_name_TYPE (ncid, vname, value/var)
!
! where get/put and TYPE=real/double occur in all permutations (hence
! it adds up to a total of eight). Because of semantically identical
! code real/double is implemented by CPP-redefinition of basic netCDF
! functions using the same source code.


      subroutine put_sclr_by_name_double(ncid, vname, value)




! These two routines are for reading or writing just a single number
! which may exist in netCDF file either as a variable or a global
! attribute containing just a single number of the proper type.
! Selection between variable or attribute is by the file, while no
! attempt to change is format is made here. For this reason is either
! variable of attribute with given name must pre-exist in order for
! these operations to succeed.

      implicit none
      integer ncid, varid, type, size, ierr, lvar
      character(len=*) vname
      real(kind=8) value
      include "netcdf.inc"

      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname(1:lvar), varid)
      if (ierr == nf_noerr) then

        ierr=nf_put_var_double(ncid, varid, value)
        if (ierr == nf_noerr) then
          write(*,'(8x,3A)') 'wrote ''', vname(1:lvar), ''''
        else
          write(*,'(/1x,4A,1x,A/)') '### ERROR: Cannot write ',
     & 'netCDF variable ''', vname, ''':', nf_strerror(ierr)
        endif
# 897 "roms_read_write.F"
      else
        ierr=nf_inq_att(ncid, nf_global, vname(1:lvar), type, size)
        if (ierr == nf_noerr) then


          if (size == 1 .and. type == nf_double) then




            ierr=nf_redef(ncid)
            if (ierr == nf_noerr) then
              ierr=nf_put_att_double(ncid, nf_global, vname(1:lvar),
     & value)
              if (ierr == nf_noerr) then
                ierr=nf_enddef(ncid)
                if (ierr == nf_noerr) then
                  write(*,'(8x,3A)') 'wrote ''', vname(1:lvar),
     & ''' as global attribute'
                else
                  write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot ',
     & 'close redefinition mode for netCDF file.',
     & nf_strerror(ierr)
                endif
              else
                write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot put ',
     & 'global attribute ''', vname(1:lvar),
     & ''' into netCDF file.', nf_strerror(ierr)
              endif
            else
              write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot switch ',
     & 'netCDF file into redefinition mode.', nf_strerror(ierr)
            endif
# 942 "roms_read_write.F"
          else
            if (size /= 1) then
              write(*,'(/1x,5A,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong size', size, 'instead of 1.'
              ierr=ierr-1
            endif

            if (type /= nf_double) then



              write(*,'(/1x,5A,I4,1x,I4,1x,A/)') '### ERROR: Global ',
     & 'attribute ''', vname(1:lvar), ''' is present, ',
     & 'but has wrong type', type, 'instead of',

     & nf_double, 'which is double precision.'



              ierr=ierr-1
            endif
          endif
        else
          write(*,'(/1x,4A/)') '### ERROR: Neither variable, nor ',
     & 'global attribute named ''', vname,
     & ''' is present in netCDF file.'
        endif
      endif
      if (ierr /= nf_noerr) stop
      end






      subroutine put_var_by_name_double(ncid, vname, var)



      implicit none
      integer ncid, ierr, varid, lvar
      character(len=*) vname
      real(kind=8) var(*)
      include "netcdf.inc"
      call lenstr(vname,lvar)
      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then

        ierr=nf_put_var_double(ncid, varid, var)
        if (ierr == nf_noerr) then
          write(*,'(8x,3A)') 'wrote ''', vname(1:lvar), ''''





        else

          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot write ',



     & 'netCDF variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot find ',
     & 'netCDF variable ID for ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
      endif
      if (ierr /= nf_noerr) stop
      end


! The following routines are to put/get just a single number into
! a specified record of one-dimensional netCDF variable. The netCDF
! file "fname" is expected to be open (hence input argument "ncid" is
! a valid file ID, while the name of the file is needed only to write
! error messages if something goes wrong, but otherwise is not used);
! "vname" is name of the variable (to be translated inside into netCDF
! ID with error message if not found), and "rec" is record number;
! "var" is input for put_ (output for get_) is just a scalar (single
! number). The corresponding netCDF variable is expected to be either
! a one-dimensional array (having more dimensions results in error
! message) or a scalar. If array it puts/gets the value at location
! "rec" with performing necessary error checking; if scalar it takes
! the only value, while argument "rec" is not used. Again, get/put
! and TYPE=real/double occur in all four permutations.



      subroutine put_sclr_rec_by_name_double(ncid, fname, vname,
     & rec, value)




      implicit none
      integer ncid, rec
      character(len=*) fname, vname
      real(kind=8) value

      character(len=16) name
      integer varid, vtype, ndims, natts, dimid(8), size,
     & start(4), count(4), ierr, lfnm, lvar
      include "netcdf.inc"

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then
          if (ndims == 1) then
            ierr=nf_inq_dimlen(ncid, dimid(1), size)
            if (ierr == nf_noerr) then
              start(1)=rec ; count(1)=1

              ierr=nf_put_vara_double(ncid,varid, start,count, value)
              if (ierr == nf_noerr) then
                write(*,'(5x,A,I5,1x,5A)') 'wrote rec', rec,
     & 'of scalar ''', vname(1:lvar), ''' into ''',
     & fname(1:lfnm), '''.'
                return !---> successful return
              else
                write(*,'(/1x,5A/12x,A)') '### ERROR: Cannot write ',
     & 'variable ''', vname(1:lvar), ''' into ''',
     & fname(1:lfnm),''':', nf_strerror(ierr)
              endif
# 1093 "roms_read_write.F"
            else
              write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: Cannot ',
     & 'Cannot determine size of dimension #', dimid(1),
     & 'in file ''', fname(1:lfnm), '''.', nf_strerror(ierr)
            endif
          elseif (ndims == 0) then


            ierr=nf_put_var_double(ncid, varid, value)
            if (ierr == nf_noerr) then
              write(*,'(6x,5A)') 'wrote scalar ''', vname(1:lvar),
     & ''' into ''', fname(1:lfnm), '''.'
              return !---> successful return
            else
              write(*,'(/1x,5A/12x,A)') '### ERROR: Cannot write ',
     & 'scalar variable ''', vname(1:lvar), ''' into ''',
     & fname(1:lfnm),''':', nf_strerror(ierr)
            endif
# 1123 "roms_read_write.F"
          else
            write(*,'(/1x,5A,I4/)') '### ERROR: Variable ''',
     & vname(1:lvar), ''' from file ''', fname(1:lfnm),
     & ''' has more then one dimension, ndims =', ndims
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! A more sophisticated functions, which read/write a portion of data
! into/from an existing netCDF variable accessing it by name. In doing
! so, it checks whether the first one, two, or three dimensions for
! that variable as defined in the netCDF file are consistent with the
! shape of array "var" specified as n1,n2,n3. If not, it complains
! about the error and quits. The variable "var" may be one, two, or
! three-dimensional, and, in addition to that may have record dimension
! (not necessarily unlimited). If the variable is two dimensional,
! then the calling routine must be specify n3=0 to avoid checking of
! the non-existing dimension. Similar policy applies for n1 and n2
! (e.g., n1,n2,n3=0,0,0 means that the variable is a scalar). However,
! if the variable in netCDF file has record dimension (identified as
! either unlimited netCDF dimension, or as the last dimension of the
! variable with its name ending with "...time" AND the actual number
! of spatial dimensions of netCDF variable in the file is LESS that
! the number of non-zero arguments THEN the excess spatial dimensions
! (n3, or both n3 and n2) will be ignored, while "rec" will be
! interpreted normally.

! Note: in the code below arguments n1,n2,n3 are used EXCLUSIVELY for
! checking, while the starts and counts are computed from the the
! dimensions of the variable in netCDF file.


      subroutine put_rec_by_name_double(ncid, fname, vname,
     & n1,n2,n3, rec, var)




! use mod_io_size_acct
      implicit none
      integer ncid, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=8), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0

      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0 ! Determine the number of spatial
        count(id)=0 ! dimensions "nspc" for the variable.
      enddo
      nspc=0 ! Note that "nspc" found here should be
      if (n1 > 0) then ! either
        nspc=nspc+1 ! less by 1 than the actual number of
        count(nspc)=n1 ! dimensions "ndims" of the variable
      endif ! stored in the file -- in this case
      if (n2 > 0) then ! the extra dimension will be treated
        nspc=nspc+1 ! as record dimension;
        count(nspc)=n2 ! or
      endif ! be the same as "ndims" -- in this
      if (n3 > 0) then ! case there is no record dimension,
        nspc=nspc+1 ! and argument "rec" is ignored.
        count(nspc)=n3
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then






          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1244 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  start(id)=1 ; rec_size=rec_size*size

                  if (count(id) /= size) then
                    call lenstr(name,ldim)
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',


     & 'put_rec_by_name_double',
# 1274 "roms_read_write.F"
     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1331 "roms_read_write.F"
              endif

              ierr=nf_put_vara_double(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                write(*,'(5x,A,I5,1x,5A)') 'wrote rec', rec, 'of ''',
     & vname(1:lvar), ''' into ''', fname(1:lfnm), '''.'
                sz_write_acc = sz_write_acc + rec_size * 8
# 1346 "roms_read_write.F"
                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max

                write_clk = write_clk + inc_clk




                return !---> successful return
              else

                write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot write ',
     & 'variable ''',vname(1:lvar),''' into netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)





              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end


! The following pair is essentially an instrumented versions of
! put_vara and get_vara to read a single record of a 2D-subdomain
! [iwest : iwest+n1-1] x [jsouth : jsouth+n2-1]
! within the (i,j)-index-space of netCDF array. The third dimension
! (if exists) is treated as the whole: n3 is checked against the
! actual dimension of netCDF variable and the mismatch is treated
! as an error; n1,n2 are checked only for upper-bound overrun, e.g.,
! iwest+n1-1 exceeds the actual size of the first dimension
! resulting in error.


      subroutine put_patch_by_name_double(ncid, fname, vname,
     & iwest,jsouth, n1,n2,n3, rec, var)




c--#define VERBOSE
! use mod_io_size_acct
      implicit none
      integer ncid, iwest,jsouth, n1,n2,n3, rec
      character(len=*) fname, vname
      real(kind=8), dimension(*) :: var

      character(len=16) name
      integer varid, vtype, nspc, ndims, natts, dimid(8),
     & rec_dimid, rec_size, size, start(4),count(4),
     & id, ierr, lfnm, lvar, ldim
      logical matched_dims
      integer(kind=8), save :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      include "netcdf.inc"

      integer iclk1, iclk2, clk_rate, clk_max, inc_clk
      call system_clock(iclk1, clk_rate, clk_max)

      read_clk=0
      sz_read_acc=0
      write_clk=0
      sz_write_acc=0
      call lenstr(fname,lfnm) ; call lenstr(vname,lvar)
      do id=1,4
        start(id)=0
        count(id)=0
      enddo ! Determine the number of spatial
      nspc=0 ! dimensions "nspc" for the variable.
      if (n1 > 0) then
        nspc=nspc+1 ! Note that "nspc" found here should be
        start(nspc)=iwest !
        count(nspc)=n1 ! either
      endif ! less by 1 than the actual number of
      if (n2 > 0) then ! dimensions "ndims" of the variable
        nspc=nspc+1 ! stored in the file -- in this case
        start(nspc)=jsouth ! the extra dimension will be treated
        count(nspc)=n2 ! as record dimension,
      endif ! or
      if (n3 > 0) then ! be the same as "ndims" -- in this
        nspc=nspc+1 ! case there is no record dimension,
        count(nspc)=n3 ! and argument "rec" is ignored.
      endif

      ierr=nf_inq_varid(ncid, vname, varid)
      if (ierr == nf_noerr) then
        ierr=nf_inq_var(ncid, varid, name, vtype, ndims, dimid, natts)
        if (ierr == nf_noerr) then





          if (nspc == ndims .or. nspc == ndims-1) then
            matched_dims=.true.

! Note: in the code below there are three ways how record dimension is
! identified: it is either
! (1) unlimited dimension, or
! (2) dimension with name ending as "...time", or
! if (3) ndims=nspc+1.

            rec_dimid=-1
            ierr=nf_inq_unlimdim(ncid, rec_dimid)
# 1478 "roms_read_write.F"
            rec_size=1
            do id=1,ndims
              ierr=nf_inq_dim(ncid, dimid(id), name, size)
              call lenstr(name,ldim)
              if (ldim > 3 .and. id == ndims) then
                if (name(ldim-3:ldim) == 'time') rec_dimid=dimid(id)
              endif
              if (ierr == nf_noerr) then
                if (dimid(id) == rec_dimid) then
                  start(id)=rec ; count(id)=1
                else
                  rec_size=rec_size*count(id)
                  if (start(id) == 0) then
                    start(id)=1
                    if (count(id) /= size) then
                      write(*,'(/1x,3A,I2,1x,3A/12x,3A,2(I5,1x,A)/)')
     & '### ERROR: ',


     & 'put_patch_by_name_double',
# 1508 "roms_read_write.F"
     & ' :: Mismatch of dimension #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', count(id),
     & 'instead of', size, 'in the netCDF file.'
                      stop
                    endif
                  elseif (start(id)+count(id)-1 > size) then
                    write(*,'(/1x,3A,I2,1x,3A/12x,3A,3(I5,1x,A)/)')
     & '### ERROR: ',


     & 'put_patch_by_name_double',
# 1530 "roms_read_write.F"
     & ' :: Overrun dimension bound #', id, 'named ''',
     & name(1:ldim), '''', 'for variable ''',
     & vname(1:lvar), ''': attempted', start(id), '+',
     & count(id), '-1 >', size, 'in the netCDF file.'
                    stop
                  endif
                endif
              else
                write(*,'(/1x,2A,I3,1x,3A/12x,A/)') '### ERROR: ',
     & 'Cannot get name and size of dimension #', id,
     & 'for variable ''', vname(1:lvar), '''.',
     & nf_strerror(ierr)
                stop
              endif
            enddo !<-- ndims

            if (matched_dims) then
              if (ndims == nspc+1) then
                start(ndims)=rec ; count(ndims)=1
# 1587 "roms_read_write.F"
              endif

              ierr=nf_put_vara_double(ncid,varid, start,count, var)
              if (ierr == nf_noerr) then
                if (start(ndims) == rec .and. count(ndims) == 1) then
                  write(*,'(5x,A,I5,1x,5A)') 'wrote rec',rec,'of ''',
     & vname(1:lvar), ''' into ''', fname(1:lfnm), ''''
                else
                  write(*,'(5x,5A)') 'wrote''', vname(1:lvar),
     & ''' into ''', fname(1:lfnm), ''''
                endif
                sz_write_acc = sz_write_acc + rec_size * 8
# 1612 "roms_read_write.F"
                call system_clock(iclk2, clk_rate, clk_max)
                inc_clk=iclk2-iclk1
                if (inc_clk < 0) inc_clk=inc_clk+clk_max

                write_clk = write_clk + inc_clk




                return !---> successful return
              else

                write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot write ',
     & 'variable ''',vname(1:lvar),''' into netCDF file ''',
     & fname(1:lfnm), '''.', nf_strerror(ierr)





              endif
            endif
          else
            write(*,'(/1x,4A,I2,A,I2,A/)') '### ERROR: Wrong number ',
     & 'of dimensions for variable ''', vname(1:lvar),
     & ''': requested ', nspc, '[+1], but found in file is',
     & ndims, '.'
          endif
        else
          write(*,'(/1x,2A,I3,1x,2A/12x,A)') '### ERROR: Cannot ',
     & 'make general inquiry for variable ID =', varid,
     & 'named ''', vname(1:lvar), '''.', nf_strerror(ierr)
        endif
      else
        write(*,'(/1x,6A/12x,A)') '### ERROR: Cannot get netCDF ID ',
     & 'for variable ''', vname(1:lvar), ''' from file ''',
     & fname(1:lfnm), ''':', nf_strerror(ierr)
      endif
      stop
      end
# 1662 "roms_read_write.F" 2
# 22 "tools_fort.F" 2
# 1 "srtopo.F" 1
      subroutine srtopo(Lm,Mm,srcdir,lon_r,lat_r,pm,pn,radius,hraw)

! Purpose: Read "raw" topography from SRTM30 dataset -- a directory
! containing 33 tiles in netCDF format -- and interpolate or
! coarsen it onto ROMS model grid. The outcome is written as
! netCDF variable "hraw".

! NOTES: (1) the content of this file is just a driver: the actual
! interpolation/coarsening routine is in "compute_hraw.F".

! (2) nothing is done here about setting mask, which needs to
! be processed separately via "rmask.F". This is to avoid
! accidental overwriting of land mask field "rmask" which
! may be subject to elaborate hand editing.

! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com

! use roms_grid_vars
! use comm_vars_hraw

      implicit none
      real(kind=8) radius
      integer nx_alc,ny_alc, nx,ny, nx_lon,ny_lat, nxy_pts,
     & nargs, Lm,Mm, ncid, ndims, varid, tile, id, size,
     & i,j, ierr, lsrc, lgrd, ldim,llat,llon
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: lon_r,lat_r,pm,pn
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: hraw
      real(kind=8), allocatable :: x(:),y(:), xlon(:),ylat(:)
      real(kind=8), allocatable :: htopo(:,:)
      integer(kind=2), allocatable :: bffr(:)
!! roms_grid_vars
      logical curv_grid
!! comm_vars_hraw
      real(kind=8) :: lon_min,lon_max, lat_min,lat_max
      integer :: comm_j1_indx, comm_j2_indx

      integer, parameter :: ntiles=33
      character(len=7), dimension(ntiles) :: srtile = (/
     & 'w180s60','w120s60','w060s60','w000s60','e060s60','e120s60',
     & 'w180s10','w140s10','w100s10','w060s10','w020s10',
     & 'e020s10','e060s10','e100s10','e140s10',
     & 'w180n40','w140n40','w100n40','w060n40','w020n40',
     & 'e020n40','e060n40','e100n40','e140n40',
     & 'w180n90','w140n90','w100n90','w060n90','w020n90',
     & 'e020n90', 'e060n90','e100n90','e140n90'/)

      character(len=7) curr_srtile
      real(kind=8), dimension(ntiles) :: xW,xE, yS,yN
      real(kind=8) xmin,xmax,ymin,ymax
      integer, dimension(ntiles) :: tndx, istr,iend,jstr,jend
      integer, dimension(ntiles) :: ewshft, iwest,ieast,jsouth
      integer, dimension(ntiles) :: jnorth
      integer :: ntls, ishft,jshft,imin,imax,jmin,jmax

      character(len=160) :: srcdir
      character(len=16) :: str, dname, lon_name, lat_name
      real(kind=8), parameter :: spv=-99999.D0

      include "netcdf.inc"

cf2py intent(in) Lm,Mm,srcdir,lon_r,lat_r,pm,pn,radius
cf2py intent(out) hraw


      comm_j1_indx=0
      comm_j2_indx=0

      call lenstr(srcdir,lsrc)
      if (srcdir(lsrc:lsrc) /= '/') then
        lsrc=lsrc+1
        srcdir(lsrc:lsrc)='/'
      endif


! Reset everything: ! "nx,ny_alc" are the actual allocated sizes
! ----------------- ! of tile-sized buffer arrays [they will grow
      nx_alc=0 ; ny_alc=0 ! if the tiles are of different sizes];
      ntls=0 ! number of active tiles;
      do tile=1,ntiles
        tndx(tile)=0
        iwest(tile)=0 ; ieast(tile)=0
        jsouth(tile)=0 ; jnorth(tile)=0
      enddo

! Open and read ROMS grid file. Find geographical limits of the grid.

      call roms_grid_geo_bounds(lon_r,lat_r, Lm,Mm, radius,
     & lon_min,lon_max,lat_min,lat_max)

      write(*,'(1x,2A,2F16.8/22x,A,2F16.8 )') 'roms grid extremes: ',
     & 'longitude:', lon_min, lon_max, 'latitude:', lat_min,lat_max

! Open and scan topography data files.

      do tile=1,ntiles
        ierr=nf_open(srcdir(1:lsrc)/ /srtile(tile)/ /'.nc',
     & nf_nowrite, ncid)
        if (ierr == nf_noerr) then
          ierr=nf_inq_ndims (ncid, ndims)
          if (ierr == nf_noerr) then
            nx=0 ; ny=0
            do id=1,ndims
              dname='       '
              ierr=nf_inq_dim (ncid, id, dname, size)
              if (ierr == nf_noerr) then
                call lenstr(dname,ldim)
                if( (ldim==1 .and. dname(1:ldim)=='x') .or.
     & (ldim==3 .and. dname(1:ldim)=='nx') .or.
     & (ldim==9 .and. dname(1:ldim)=='longitude') ) then
                  lon_name=dname(1:ldim)
                  nx=size
                elseif( (ldim==1 .and.dname(1:ldim)=='y') .or.
     & (ldim==3 .and.dname(1:ldim)=='ny') .or.
     & (ldim==8 .and.dname(1:ldim)=='latitude') ) then
                  lat_name=dname(1:ldim)
                  ny=size
                endif
              else
                write(*,*) '### ERROR: dimension id =', id, '?'
              endif
            enddo
          else
            write(*,*) '### ERROR: nf_inq_ndims?'
          endif
!#ifdef VERBOSE
! write(*,'(1x,A,3(1x,A),I6,2x,2(1x,A),I6)') srtile(tile),
! & ': dimensions:', lon_name(1:lenstr(lon_name)), '=', nx,
! & lat_name(1:lenstr(lat_name)), '=', ny
!#endif

          if (nx > 0 .and. ny > 0) then
            if (nx > nx_alc) then ! NOTE: Since the data
              if (allocated(x)) deallocate(x) ! files are expected to
              allocate(x(nx)) ! be what is called
              nx_alc=nx ! "CF-compliant" its
            endif ! coordinate variables
                                               ! should have the same
            if (ny > ny_alc) then ! names as their
              if (allocated(y)) deallocate(y) ! corresponding
              allocate(y(ny)) ! dimensions.
              ny_alc=ny
            endif

            ierr=nf_inq_varid (ncid, lon_name, varid)
            if (ierr == nf_noerr) then
              ierr=nf_get_var_double (ncid, varid, x)
              if (ierr == nf_noerr) then
                ierr=nf_inq_varid (ncid, lat_name, varid)
                if (ierr == nf_noerr) then
                  ierr=nf_get_var_double (ncid, varid, y)
                  if (ierr == nf_noerr) then

! Note that fiction "indx_bound" returns 0 on "ny" if the test value
! "lon_mix/max" is outside the range of "y" for this tile, so the first
! if-condition indicates that that there is an overlap in y-coordinate
! (latitude of data) between ROMS grid and the tile.

                    jmin=indx_bound3(y,ny, lat_min)
                    jmax=indx_bound3(y,ny, lat_max)
                    if (jmin<ny .and. jmax>0) then
                      if (jmin == 0) jmin=1
                      if (jmin > 1 ) jmin=jmin-1
                      if (jmax < ny) jmax=jmax+1

! Longitude coordinate is defined with 360 degree periodicity, so if
! the tile dies not fit right a way, try to shift it east of west by
! 360 degrees and then check again whether it has overlap.

                      ishft=0
                      imin=indx_bound3(x,nx, lon_min)
                      imax=indx_bound3(x,nx, lon_max)
                      if (imin==nx) then
                        ishft=+1
                        imin=indx_bound3(x,nx, lon_min-360.D0)
                        imax=indx_bound3(x,nx, lon_max-360.D0)
                      elseif (imax==0) then
                        ishft=-1
                        imin=indx_bound3(x,nx, lon_min+360.D0)
                        imax=indx_bound3(x,nx, lon_max+360.D0)
                      endif
                      if (imin<nx .and. imax>0) then
                        if (imin==0) imin=1
                        if (imin >1) imin=imin-1
                        if (imax<nx) imax=imax+1

! Once it passed all the logical check above, add the tile to the list
! of tiles to be read, and record all its attributes.

                        ntls=ntls+1
                        tndx(ntls)=tile ; ewshft(ntls)=ishft
                        iwest(ntls)=imin ; jsouth(ntls)=jmin
                        ieast(ntls)=imax ; jnorth(ntls)=jmax

                        if (ishft > 0) then
                          xW(ntls)=x(imin) +360.D0
                          xE(ntls)=x(imax) +360.D0
                        elseif (ishft < 0) then
                          xW(ntls)=x(imin) -360.D0
                          xE(ntls)=x(imax) -360.D0
                        else
                          xW(ntls)=x(imin)
                          xE(ntls)=x(imax)
                        endif
                        yS(ntls)=y(jsouth(ntls))
                        yN(ntls)=y(jnorth(ntls))

                      endif !<-- imin<nx .and. imax>0
                    endif !<-- jmin<ny .and. jmax>0

                  else
                    call lenstr(lat_name,llat)
                    write(*,'(/1x,6A/12x,A/)') '### ERROR: Cannot ',
     & 'read variable ''', lat_name(1:llat),
     & ''' from ''', srtile(tile), '.nc''.',
     & nf_strerror(ierr)
                  endif
                else
                  call lenstr(lat_name,llat)
                  write(*,'(/1x,6A/12x,A/)') '### ERROR: Cannot get ',
     & 'netCDF variable ID for ''',
     & lat_name(1:llat), ''' from  ''',
     & srtile(tile), '.nc''.', nf_strerror(ierr)
                endif
              else
                call lenstr(lon_name,llon)
                write(*,'(/1x,6A/12x,A/)') '### ERROR: Cannot read ',
     & 'variable ''', lon_name(1:llon),
     & ''' from ''', srtile(tile), '.nc''.', nf_strerror(ierr)
              endif
            else
              call lenstr(lon_name,llon)
              write(*,'(/1x,6A/12x,A/)') '### ERROR: Cannot get ',
     & 'netCDF variable ID for ''', lon_name(1:llon),
     & ''' from ''', srtile(tile), '.nc''.', nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,4A/)') '### ERROR: Cannot determine ',
     & 'dimension sizes for topography data file ''',
     & srtile(tile), '.nc''.'
          endif
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot open netCDF ',
     & 'file ''', srtile(tile), '.nc'' in read-only mode.',
     & nf_strerror(ierr)
        endif
        ierr=nf_close(ncid)
      enddo
      if (allocated(x)) deallocate(x)
      if (allocated(y)) deallocate(y)



! Summarize the results of tile scan.


      write(*,'(/ /1x,2A,I3,3x,A/1x,2A/)') 'SUMMARY: number of ',
     & 'selected tiles, ntls =', ntls, 'Their file names, bounding',
     & 'indices of actually used portions, and geographical limits ',
     & 'of used portions:'

      write(*,'(2x,A,1x,A,3x,A,1x,A,1x,A,1x,A,2x,A,11x,A,10x,A,9x,A)')
     & '#', 'file', 'iwst','iest','jsth','jnrth', 'xW','xE','yS','yN'
      write(*,'(1x,2A)') '---------------------------------------',
     & '---------------------------------------'
      do tile=1,ntls
        write(*,'(I3,1x,A,4I5,2F13.7,2F11.7)') tile, srtile(tndx(tile)),
     & iwest(tile),ieast(tile),jsouth(tile),jnorth(tile),
     & xW(tile),xE(tile),yS(tile),yN(tile)
      enddo

      xmin=xW(1) ; xmax=xE(1)
      ymin=YS(1) ; ymax=yN(1)
      do tile=2,ntls
        if (xW(tile)<xmin) xmin=xW(tile)
        if (xE(tile)>xmax) xmax=xE(tile)
        if (yS(tile)<ymin) ymin=yS(tile)
        if (yN(tile)>ymax) ymax=yN(tile)
      enddo

      write(*,'(/1x,2A,/1x,A,2F13.8,2x,A,2F12.8/)') 'Geographical ',
     & 'extremes of the used portions among all tiles:',
     & 'xmin,max =',xmin,xmax, 'ymin,max =',ymin,ymax


! The following code segment determines bounding indices of each tile
! (actually used portion thereof) within the index coordinates of would
! be global data array. This is needed for correct placement of data
! from each individual tile.
! The mapping is as follows:
! iwest:ieast --> istr:iend
! jsouth:jnorth --> jstr:jend

! The algorithm is essentially to convert geographical limits of each
! used tile portion into integer numbers implicitly assuming that the
! grid spacing is globally uniform among all the tiles. [Note that 0.5
! is added inside each int() to counter roundoff error - mathematically
! the result of division (ieast-iwest)/(xE-xW) yields the inverse of
! grid spacing while (xW)-xmin) should be an integer number of grid
! spaces, so the whole expression inside int() without 0.5 should be
! an integer number. However minute roundoff error may cause it to
! be 1 less than it should.]

! The dimensions of global data arrays are defined as the maximum
! iend and jend among all the tiles.

      nxy_pts=0 ; nx_lon=0 ; ny_lat=0

      do tile=1,ntls
        istr(tile)=1+int( 0.5D0 + dble( ieast(tile)-iwest(tile) )
     & *(xW(tile)-xmin)/(xE(tile)-xW(tile)) )
        iend(tile)=istr(tile) + ieast(tile)-iwest(tile)
        if (iend(tile)>nx_lon) nx_lon=iend(tile)

        jstr(tile)=1+int( 0.5D0 + dble( jnorth(tile)-jsouth(tile) )
     & *(yS(tile)-ymin)/(yN(tile)-yS(tile)) )
        jend(tile)=jstr(tile) +jnorth(tile)-jsouth(tile)
        if (jend(tile)>ny_lat) ny_lat=jend(tile)

        nxy_pts=nxy_pts + (iend(tile)-istr(tile)+1)
     & *(jend(tile)-jstr(tile)+1)
      enddo

      write(*,'(/2(1x,A,I7)/15x,A,I12/42x,A,I12/66x,A/47x,A,I12/)')
     & 'Required dimensions of array to hold the entire used data:',
     & nx_lon, 'x', ny_lat,
     & 'Aggregate number of points in all selected tiles:', nxy_pts,
     & 'product of dimensions:', nx_lon*ny_lat,
     & '-----------', 'their difference:', nxy_pts-nx_lon*ny_lat


      write(*,'(/1x,2A/1x,A/)') 'Starting and ending indices for ',
     & 'the used portions of each tile as defined',
     & 'within the logical coordinates of global data array:'

      write(*,'(10x,A,1x,A,2x,A/4x,2A)') '#', 'file', 'ew_shft',
     & '------------------------------------',
     & '------------------------------------'
      do tile=1,ntls
        write(*,'(8x,I3,1x,A,1x,I3,1x,A,2I7,2x,A,2I7)') tile,
     & srtile(tndx(tile)), ewshft(tile),
     & 'istr,iend =', istr(tile),iend(tile),
     & 'jstr,jend =', jstr(tile),jend(tile)
      enddo

      if (nxy_pts < nx_lon*ny_lat) then
        write(*,'(/1x,2A/)') '### ERROR: Available data tiles do ',
     & 'not cover the entire ROMS grid.'
        stop
      elseif (nxy_pts > nx_lon*ny_lat) then
        write(*,'(/1x,2A/)') '### ERROR: Possible overlapping tiles ',
     & 'or tile selection algorithm failure.'
        stop
      endif


      write(*,'(/1x,A/)') 'Reading topographic data..'

! Note that the tile-sized coordinate arrays x,y were deallocated above
! and are allocated again just below, but after the arrays covering the
! entire grid area. This is to avoid memory fragmentation by keeping
! the tile-sized x,y, and bffr at the end of allocated memory, are they
! will be deallocated as soon as reading of tiled data is complete, and
! a new array to hold interpolated topography will be allocated.

      allocate(xlon(nx_lon)) ; allocate(ylat(ny_lat))
      allocate(htopo(nx_lon,ny_lat))

      allocate(x(nx_alc)) ; allocate(y(ny_alc))
      allocate(bffr(nx_alc*ny_alc))

      do i=1,nx_lon ! Initialize coordinate arrays with a special
        xlon(i)=spv ! value. This is needed to check consistency of
      enddo ! coordinates stored in different tiles. Each
      do j=1,ny_lat ! coordinate value will be recorded only once,
        ylat(j)=spv ! and thereafter checked that the value from a
      enddo ! different tile matches the already recorded.

      do tile=1,ntls
        curr_srtile=srtile(tndx(tile))
        ierr=nf_open(srcdir(1:lsrc)/ /curr_srtile/ /'.nc',
     & nf_nowrite, ncid)
        if (ierr == nf_noerr) then
          ierr=nf_inq_dimid (ncid, lon_name, id)
          if (ierr == nf_noerr) then
            ierr=nf_inq_dimlen (ncid, id, nx)
            if (ierr == nf_noerr) then
              ierr=nf_inq_dimid (ncid, lat_name, id)
              if (ierr == nf_noerr) then
                ierr=nf_inq_dimlen (ncid, id, ny)
                if (ierr == nf_noerr) then
                  call lenstr(lon_name,llon)
                  call lenstr(lat_name,llat)
                  write(*,'(I3,4(1x,A),I6,2x,2(1x,A),I6)',advance='no')
     & tile, curr_srtile, 'dimensions:',
     & lon_name(1:llon), '=', nx,
     & lat_name(1:llat), '=', ny
                else
                  call lenstr(lat_name,llat)
                  write(*,'(1x,2A,I3,1x,5A/12x,A)') '### ERROR: ',
     & 'Cannot determine length of dimension ', id,
     & 'named ''', lat_name(1:llat),
     & ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
                endif
              else
                call lenstr(lat_name,llat)
                write(*,'(1x,6A/12x,A)') '### ERROR: Cannot get ',
     & 'dimension ID for ''', lat_name(1:llat),
     & ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
              endif
            else
              call lenstr(lon_name,llon)
              write(*,'(1x,2A,I3,1x,5A/12x,A)') '### ERROR: Cannot ',
     & 'determine length of dimension ', id, 'named ''',
     & lon_name(1:llon), ''' from ''',
     & curr_srtile, '''.', nf_strerror(ierr)
            endif
          else
            call lenstr(lon_name,llon)
            write(*,'(1x,6A/12x,A)') '### ERROR: Cannot get ',
     & 'dimension ID for ''', lon_name(1:llon),
     & ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
          endif

          if (ierr == nf_noerr) then
            ierr=nf_inq_varid (ncid, lon_name, varid)
            if (ierr == nf_noerr) then
              ierr=nf_get_var_double (ncid, varid, x)
              if (ierr == nf_noerr) then
                ierr=nf_inq_varid (ncid, lat_name, varid)
                if (ierr == nf_noerr) then
                  ierr=nf_get_var_double (ncid, varid, y)
                  if (ierr /= nf_noerr) then
                    call lenstr(lat_name,llat)
                    write(*,'(1x,2A,I3,1x,5A/12x,A)') '### ERROR: ',
     & 'Cannot read coordinate variable ', varid,
     & 'named ''', lat_name(1:llat),
     & ''' from ''', curr_srtile, '''.', nf_strerror(ierr)






                  endif
                else
                  call lenstr(lat_name,llat)
                  write(*,'(1x,6A/12x,A)') '### ERROR: Cannot get ',
     & 'variable ID for ''', lat_name(1:llat),
     & ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
                endif
              else
                call lenstr(lon_name,llon)
                write(*,'(1x,2A,I3,1x,5A/12x,A)') '### ERROR: Cannot ',
     & 'read coordinate variable ', varid, 'named ''',
     & lon_name(1:llon), ''' from ''',
     & curr_srtile, '''.', nf_strerror(ierr)
              endif
            else
              call lenstr(lon_name,llon)
              write(*,'(1x,6A/12x,A)') '### ERROR: Cannot determine ',
     & 'variable ID for ''', lon_name(1:llon),
     & ''' from ''', curr_srtile, '''.', nf_strerror(ierr)
            endif
          endif

          if (ierr == nf_noerr) then
            ierr=nf_inq_varid (ncid, 'topo', varid)
            if (ierr == nf_noerr) then
              write(*,'(12x,A)') 'found variable ''topo''.'
            else
              ierr=nf_inq_varid (ncid, 'elevation', varid)
              if (ierr == nf_noerr) then
                write(*,'(12x,A)') 'found variable ''elevation''.'
              else
                ierr=nf_inq_varid (ncid, 'z', varid)




              endif
            endif
          endif

          if (ierr == nf_noerr) then
            ierr=nf_get_var_int2 (ncid, varid, bffr)
            if (ierr == nf_noerr) then
              write(*,'(1x,A)') 'retrieved topographic data'
            else
              write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot read ',
     & 'topographic data from ''', curr_srtile, '''.',
     & nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot determine ',
     & 'netCDF variable ID for topography field in ''',
     & curr_srtile, '''.', nf_strerror(ierr)
          endif

          if (ierr == nf_noerr) then
            ishft=istr(tile)-iwest(tile)
            if (ewshft(tile)>0) then
              do i=iwest(tile),ieast(tile)
                x(i)=x(i) +360.D0
              enddo
            elseif (ewshft(tile)<0) then
              do i=iwest(tile),ieast(tile)
                x(i)=x(i) -360.D0
              enddo
            endif
            ierr=0
            do i=iwest(tile),ieast(tile)
              if (xlon(i+ishft) > spv) then
                if (xlon(i+ishft) /= x(i)) ierr=ierr+1
              else
                xlon(i+ishft)=x(i)
              endif
            enddo
            if (ierr>0) write(*,'(/1x,2A/)') '### ERROR: Conflicting ',
     & 'longitude coordinate data between tiles.'

            jshft=jstr(tile)-jsouth(tile)
            ierr=0
            do j=jsouth(tile),jnorth(tile)
              if (ylat(j+jshft) > spv) then
                if (ylat(j+jshft) /= y(j)) ierr=ierr+1
              else
                ylat(j+jshft)=y(j)
              endif
            enddo
            if (ierr>0) write(*,'(/1x,2A/)') '### ERROR: Conflicting ',
     & 'latitude coordinate data between tiles.'


            do j=jsouth(tile),jnorth(tile)
              do i=iwest(tile),ieast(tile)
                htopo(i+ishft,j+jshft)=bffr(i+nx*(j-1))
              enddo
            enddo
# 551 "srtopo.F"
          endif
        else
          write(*,'(/1x,4A/12x,A/)') '### ERROR: Cannot open netCDF ',
     & 'file ''', srtile(tile), '.nc'' in read-only mode.',
     & nf_strerror(ierr)
        endif
        ierr=nf_close(ncid)
      enddo

      deallocate(bffr) ; deallocate(y) ; deallocate(x)

! 3. Once both the input topography data and the target roms grid
! file are successfully opened and the relevant data is read, do the
! actual interpolation:

      if (ierr == nf_noerr) then
! allocate (hraw(0:Lm+1,0:Mm+1)) ; ierr=0

C$OMP PARALLEL SHARED(nx_lon,ny_lat, xlon,ylat, htopo, Lm,Mm,
C$OMP& lon_r,lat_r, pm,pn, hraw, radius, ierr)
        call compute_hraw(nx_lon,ny_lat, xlon,ylat, htopo,
     & Lm,Mm, lon_r,lat_r, pm,pn, hraw, radius)
C$OMP END PARALLEL
! if (ierr==0) call write_hraw(grid, srcdir,radius, Lm,Mm,hraw)
      endif

      contains
        integer function indx_bound3(x,n, x0)
        integer n, i ! bounded by
        real(kind=8) x(n), x0 ! x(i) <= x0 < x(i+1)

        if (x0 < x(1)) then
          i=0 ! if x0 is outside the full range
        elseif (x0 > x(n)) then ! of x(1) ... x(n), then return
          i=n ! i=0 or i=n.
        else
          i=int( ( x(n)-x0 +n*(x0-x(1)) )/(x(n)-x(1)) )
          if (x(i+1)<x0) then
            do while (x(i+1) < x0) ! This algorithm computes "i" as
              i=i+1 ! linear interpolation between x(1)
            enddo ! and x(n) which should yield the
          elseif (x(i) > x0) then ! correct value for "i" right a way
            do while (x(i) > x0) ! because array elements x(i) are
              i=i-1 ! equidistantly spaced. The while
            enddo ! loops are here merely to address
          endif ! possible roundoff errors.

          if (x(i+1)-x0 < 0 .or. x0-x(i) < 0) then
            write(*,'(1x,A,5F12.6)') '### ERROR: indx_bound :: ',
     & x(i), x0, x(i+1), x0-x(i), x(i+1)-x0
            stop
          endif
        endif
        indx_bound3=i
        end function indx_bound3


      end subroutine srtopo
# 23 "tools_fort.F" 2
!!! compute hraw !!!
# 1 "compute_hraw.F" 1
      subroutine compute_hraw (lon,lat, xlon,ylat, htopo, Lm,Mm,
     & lon_r,lat_r, pm,pn, h, radius)

! A universal interpolation/coarsening procedure to compute topography
! at model grid. It takes topography data "htopo" defined at source
! grid (xlon,ylat) and transfers it onto the target grid (lon_r,lat_r)
! using weighted averaging,
!
! weight(r) = [1-(r/width)^2]^2, r < width
!
! or, ideally, but computationally costly, using HAN WINDOW,
!
! weight(r) = 1/2 + 1/2*cos(pi*r/width), r < width
!
! or its polynomial fit
!
! weight(r) = 1/2 + 1/2*cos(pi*r/width), r < width
!
! or its polynomial fit
!
! weight(r)=[ 1- 1.228..*(r/width)^2 +0.228..*(r/width)^4 ]^2
!
! with coefficient chosen to make weight(width/2)= 1/2 (see below).
!
! where "r" is the distance between target point and data point, and
! "width" is specified averaging width. The goal here is to achieve
! isotropic response of the algorithm to unresolved spike in the data
! (delta function) in sense that "width", expressed as real distance
! (in meters) should be the same in all directions. Furthermore,
! "width" controlled by non-dimensional input parameter "radius" as
! multiple of local grid spacing
!
! width[meters] = radius/sqrt(0.5*(pm^2+pn^2))
!
! Translated into geographical lon-lat coordinates, this means that
! "wdthX" (width along lon, in degrees) and "wdthY" (along lat,
! degrees) are related as wdthY = wdthX*cos(lat), and,
!
! wdthY[degrees] = width[meters]*180/(pi*Eradius)
!
! where Eradius is radius of the earth. The weighted coarsening is
! activated by CPP-switch WEIGHTED_AVERAGING; otherwise a simple
! bi-linear interpolation is performed.
!
! PARALLELIZATION METHOD: This routine is expected to be called from
! inside an Open MP parallel region and it employs dynamic scheduling
! by processing one j-row of array at a time,
!
! do while(.true.)
! C$OMP CRITICAL
! j=comm_j_index <-- check out one j-row to process
! comm_j_index=j+1 <-- increment shared counter to signal
! C$OMP END CRITICAL to the other threads that this row
! if (j > Mm+1) --> DONE is already taken
! do i=....
! .... do useful work
! enddo
! enddo
!
! where "comm_j_index" is a shared counter initialized to zero used as
! a signalling variable. This method is preferred to more conventional
! tiling here because the input topographic data array "htopo" may be
! very large, so the algorithm here is designed to keep the threads
! working close to each other in memory space (perhaps synergetically
! benefit from shared shared cache of a multicore CPU by mutually
! reusing data brought into cache by different threads).

! Created and maintained by Alexander Shchepetkin old_galaxy@yahoo.com



c-#define VERBOSE

      implicit none
      integer:: lon,lat,Lm,Mm!, ierr
      real(kind=8), dimension(lon,lat):: htopo
      real(kind=8) :: radius, xlon(lon), ylat(lat)
      real(kind=8), dimension(0:Lm+1,0:Mm+1):: lon_r,lat_r,pm,pn
      real(kind=8), dimension(0:Lm+1,0:Mm+1):: h
      real(kind=8) :: lon_min,lon_max, lat_min,lat_max
      integer :: comm_j1_indx, comm_j2_indx

!--> internal variables

      integer, dimension(0:Lm+1) :: ic,jc
      real(kind=8), dimension(0:Lm+1) :: xi,eta
      real(kind=8):: xwest,xeast,ysouth,ynorth, dx_inv, dy_inv, xr,yr
      integer:: numthreads, trd, iwest,ieast,jsouth,jnorth,
     & i,j, ir,jr, iter
C$ integer omp_get_num_threads, omp_get_thread_num


      logical :: bounded
      integer :: istr,iend,jstr,jend, imin,imax,jmin,jmax, mwx,mwy,
     & mwx_alloc,mwy_alloc, irc,jrc, itm,jtm, iups,jups
      real(kind=8), dimension(:,:), allocatable :: wgt,wgt_new
      real(kind=8):: SmScle,wdthX,wdthY, sum,biasX,biasY, xc,yc,
     & cff,ctr, cx,cy, dcx,dcy, acx,acy


      integer(kind=4):: iclk_start, iclk_end, clk_rate, clk_max


cf2py intent(in) Lm,Mm,lon,lat,xlon,ylat,htopo,lon_r,lat_r,pm,pn,radius
cf2py intent(out) h

      include "phys_const.h"

      mwx_alloc=0 ; mwy_alloc=0
      dx_inv=dble(lon-1)/(xlon(lon)-xlon(1))
      dy_inv=dble(lat-1)/(ylat(lat)-ylat(1))

      comm_j1_indx=0
      comm_j2_indx=0
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads()
C$ trd=omp_get_thread_num()


      if(trd == 0) call system_clock(iclk_start, clk_rate, clk_max)

                     ! Check whether the data coverage is sufficient
! Preliminary step: ! for the entire ROMS grid. Complain about error
!------------------ ! and signal to quit, if anything goes wrong.

! Find extremes of lon-lat oordinates for topography data set

      xwest=xlon(1) ; xeast=xlon(lon)
      ysouth=ylat(1) ; ynorth=ylat(lat)

! Find extremes of lon-lat coordinates for croco data set
      call roms_grid_geo_bounds(lon_r,lat_r,Lm,Mm,radius,
     & lon_min,lon_max,lat_min,lat_max)

C$OMP MASTER
      write(*,'(/1x,2A,2(/4x,A,2F12.6))') 'Geographical extents of ',
     & 'the dataset:', 'Longitude min/max =', xwest, xeast,
     & 'Latitude  min/max =', ysouth, ynorth
      write(*,'(/1x,2A,2(/4x,A,2F12.6))') 'Geographical extents of ',
     & 'requested ROMS grid:', 'Longitude min/max =', lon_min,lon_max,
     & 'Latitude  min/max =', lat_min,lat_max
C$OMP END MASTER

      iter=0
      if (lon_min < xwest) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)') '### WARNING: ',
     & 'Western extent of model grid is beyond western bound of ',
     & 'dataset.', 'lon_min_data =', xwest, 'lon_min =', lon_min
C$OMP END MASTER
      endif
      if (lon_max > xeast) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)') '### WARNING: ',
     & 'Eastern extent of model grid is beyond eastern bound of ',
     & 'dataset.', 'lon_max_data =', xeast, 'lon_max =', lon_max
C$OMP END MASTER
      endif
      if (lat_min < ysouth) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)') '### WARNING: ',
     & 'Southern extent of model grid is beyond southern bound of ',
     & 'dataset.', 'lat_min_data =', ysouth, 'lat_min =', lat_min
C$OMP END MASTER
      endif
      if (lat_max > ynorth) then
        iter=iter+1
C$OMP MASTER
        write(*,'(/1x,3A/12x,A,F12.6,1x,A,F12.6/)') '### WARNING: ',
     & 'Northern extent of model grid is beyond northern bound of ',
     & 'dataset.', 'lat_max_data =', ynorth, 'lat_max =', lat_max
C$OMP END MASTER
      endif
      if (iter > 0) then ! Note: Because the error condition checked
C$OMP MASTER ! here occurs simultaneously on all threads,
! ierr=iter ! all exit simultaneously. The only reason
C$OMP END MASTER ! for setting "ierr" is to prevent writing
c return ! "hraw" into netCDF file in the case of
      endif ! error detected here.

      iwest=indx_bound2(xlon,lon, lon_min)
      ieast=indx_bound2(xlon,lon, lon_max)
      jsouth=indx_bound2(ylat,lat, lat_min)
      jnorth=indx_bound2(ylat,lat, lat_max)

C$OMP MASTER
      write(*,'(/1x,2A/2(1x,A,I5,I6,1x,A,I6))') 'Bounding indices ',
     & 'for the portion of dataset covering the entire model grid:',
     & 'iwest,ieast =', iwest,ieast, 'out of', lon,
     & 'jsouth,jnorth =',jsouth,jnorth, 'out of', lat

      if (iwest == 0) write(*,*) 'WARNING: compute_hraw :: ',
     & 'restricting iwest from 0 to 1'
      if (ieast == lon) write(*,*) 'WARNING: compute_hraw :: ',
     & 'restricting ieast from', lon, ' to', lon-1
      if (jsouth == 0) write(*,*) 'WARNING: compute_hraw :: ',
     & 'restricting jsouth from 0 to 1'
      if (jnorth == lat) write(*,*) 'WARNING: compute_hraw :: ',
     & 'restricting jnorth from', lat, ' to', lat-1




      write(*,'(/1x,A,F8.5,1x,A)') 'Averaging window width =',
     & radius, 'in grid spaces.'



C$OMP END MASTER
C$OMP BARRIER







! Start interpolation: For each point of model grid, "lon_r,lat_r"
!------ -------------- find indices "ic,jc" of the data grid such
! that
! xlon(ic) <= lon_r < xlon(ic+1)
! ylat(jc) <= lat_r < ylat(jc+1)
!
! after which compute fractional distances "xi,eta"
!
! 0 <= xi,eta < 1
!
! which identify position of point "lon_r,lat_r" relatively to the
! rectangle [xlon(ic):xlon(jc+1)]X[ylat(jc):ylat(jc+1)].
! This part is the same for weighted and bi-linear interpolation.
!
! The rest of this code is organized as as single giant j-loop lasting
! all the way to the end.
      do while(.true.)
C$OMP CRITICAL(progress_sync)
        j=comm_j2_indx ; comm_j2_indx=j+1




        if (mod(j,5) == 0) then
          if (mod(j,375) == 0) write(6,*)
          write(6,'(A)',advance='no') '.'
          flush(unit=6)
        endif

C$OMP END CRITICAL(progress_sync)
        if (j > Mm+1) exit !goto 98 !--> DONE

! Find bounding indices (ic,jc) for each point (i,j) of ROMS grid,
! such that ROMS grid points with coordinates lon_r(i,j),lat_r(i,j)
! falls between ic and ic+1, as well as jc and jc+1 on the data grid
! with coordinates xlon,ylat. Mathematically speaking the "while"
! searches in the longitudinal direction in the code below are not
! necessary because the data grid is always uniform in this direction.
! However, the dataset may store coordinates in single precision which
! may cause round-off-level violations of "xi,eta" bounds triggering
! error messages below. Latitudinal direction may theoretically use
! non-uniform grid spacing, although in most datasets it does not.

! The code below (starting from here and all the way to the end) is
! adapted to handle the situation that the dataset may not cover the
! entire ROMS grid, and even not to cover the entire rectangular
! portion [xwest,xeast] x [ysouth,ynorth] due to tiled approach for
! storing data in multiple netCDF files in such a way that individual
! tiles use lon,lat coordinates within, but mutial arrangement of the
! tiles is rather arbitrary. For this reason (ic,jc) are set first to
! unrealistic special values, which are then owerwritten with
! meaningful indices only if interpolation is possible.






        do i=0,Lm+1
          ic(i)=-1 ; jc(i)=-1 !<-- special values
          if (xwest < lon_r(i,j) .and. lon_r(i,j) < xeast) then
            ir=1+int(dx_inv*(lon_r(i,j)-xwest))
            ir=max(min(ir, lon-2), 2)
            do while(xlon(ir) > lon_r(i,j) .and. ir > 2)
              ir=ir-1
            enddo
            do while(xlon(ir+1) <= lon_r(i,j) .and. ir < lon-2)
              ir=ir+1
            enddo
            if ( xlon(ir) <= lon_r(i,j) .and. lon_r(i,j) < xlon(ir+1)
     & ) ic(i)=ir
          endif

          if (ysouth < lat_r(i,j) .and. lat_r(i,j) < ynorth) then
            jr=1+int(dy_inv*(lat_r(i,j)-ysouth))
            jr=max(min(jr, lat-2), 2)
            do while(ylat(jr)> lat_r(i,j) .and. jr > 2)
              jr=jr-1
            enddo
            do while(ylat(jr+1) <= lat_r(i,j) .and. jr < lat-2)
              jr=jr+1
            enddo
            if ( ylat(jr) <= lat_r(i,j) .and. lat_r(i,j) < ylat(jr+1)
     & ) jc(i)=jr
          endif


! Compute compute fractional distances "xi,eta" which locate the
! position of "lon_r,lar_r" within the 2D-cell [ir:ir+1]X[jr:jr+1] on
! the data grid. If everything goes correctly both "xi" and "eta" must
! be bounded as 0 <= xi,eta < 1.

          if (ic(i) > 0 .and. jc(i) > 0) then
            ir=ic(i) ; jr=jc(i)
            xi(i)= (lon_r(i,j)-xlon(ir))/(xlon(ir+1)-xlon(ir))
            eta(i)=(lat_r(i,j)-ylat(jr))/(ylat(jr+1)-ylat(jr))
          else
            xi(i)=0.D0 ; eta(i)=0.D0
          endif
        enddo
# 348 "compute_hraw.F"
! Find bounding indices "istr:iend" and "jstr:jend" enclosing all data
! grid points which lie within the distances "wdthX" (in longitudinal)
! and "wdthY" (in meridional direction) from the point "lon_r,lat_r"
! of the model grid.

        SmScle=radius *180.D0/(pi*Eradius*sqrt(0.5D0))

        do i=0,Lm+1
          bounded=.false.
          if (ic(i) > 0 .and. jc(i) > 0) then
            irc=ic(i) ; xc=lon_r(i,j) ; istr=irc ; iend=irc+1
            jrc=jc(i) ; yc=lat_r(i,j) ; jstr=jrc ; jend=jrc+1

            wdthY=SmScle/sqrt(pm(i,j)**2+pn(i,j)**2)
            wdthX=wdthY/cos(deg2rad*yc)

c* wdthX=0.1*wdthX !<-- checking
c* wdthY=0.1*wdthY !<-- isotropy


! Expand eastern and western bounding indices until the interval
! xlon(istr):xlon(iend) is wide enough to cover "xc-wdthX:xc+wdthX".
! Similarly for southern and northern bounding indices to make
! ylat(jstr):ylat(jend) cover "yc-wdthY:yc+wdthY".


            do while(xc-wdthX < xlon(istr) .and. istr > 3)
              istr=istr-1
            enddo
            do while(xlon(iend) < xc+wdthX .and. iend < lon-2)
              iend=iend+1
            enddo

            do while(yc-wdthY < ylat(jstr) .and. jstr > 3)
              jstr=jstr-1
            enddo
            do while(ylat(jend) < yc+wdthY .and. jend < lat-2)
              jend=jend+1
            enddo

            if (xlon(istr) < xc-wdthX .and. xc+wdthX < xlon(iend) .and.
     & ylat(jstr) < yc-wdthY .and. yc+wdthY < ylat(jend)) then
              bounded=.true.
            endif
          endif !<-- ic(i) > 0 .and. jc(i) > 0



          if (bounded) then
            mwx=max(iend-irc+1,irc-istr+1)
            mwy=max(jend-jrc+1,jrc-jstr+1)
            if (mwx > mwx_alloc .or. mwy > mwy_alloc) then
              mwx_alloc=mwx+1 ! add 1 merely to reduce
              mwy_alloc=mwy+1 ! probability of the need
              if (allocated(wgt)) then ! to allocate arrays of a
                deallocate(wgt) ! larger size later
                deallocate(wgt_new)
              endif
              allocate(wgt(-mwx_alloc:mwx_alloc,-mwy_alloc:mwy_alloc),
     & wgt_new(-mwx_alloc:mwx_alloc,-mwy_alloc:mwy_alloc))






            endif

! Construct interpolation/weighting coefficients "wgt":
!---------- ------------- --------- ------------ -------
! The intended shape of weight is a cosine function,
!
! f(x) = 1/2 + (1/2)*cos(pi*x) , |x|<1
!
! Below it is approximated by f(x)=(1- 1.228..*x^2 +0.228..*x^4)^2
! where 1.22876383367175=5-8*sqrt(2)/3 is chosen to make f(1/2)=1/2
! like it is for cosine function above. The difference between the
! two functions within the entire interval |x|<1 is between -1.1e-3
! to +5e-4.

! Note that in the case when width of averaging window is narrower
! than the distance between adjacent data points, there is a chance
! that all weights computed by the algorithm below are all zeros. In
! anticipation of this possibility weights are set initially to delta
! function at the closest discrete point, and are then only if the
! window is wide enough for averaging to occur (hence the assignment
! of "weight" occurs conditionally in the 2D ir-jr loop below).
! Subsequently, in the case when "weight" remain delta-function, the
! exact coefficients for bi-linear interpolation are recovered by the
! iterative iterative centering procedure.

            do jr=jstr-1,jend+1
              do ir=istr-1,iend+1
                wgt(ir-irc,jr-jrc)=0.D0 ; wgt_new(ir-irc,jr-jrc)=0.D0
              enddo
            enddo
            if (xi(i) < 0.5D0) then
              ir=0
            else
              ir=1
            endif
            if (eta(i) < 0.5D0) then
              jr=0
            else
              jr=1
            endif
            wgt(ir,jr)=1.D0 !<-- delta function

            dcx=1.D0/wdthX ; dcy=1.D0/wdthY
            sum=0.D0 ; biasX=0.D0 ; biasY=0.D0


            do jr=jstr,jend
              do ir=istr,iend
                xr=xlon(ir) ; yr=ylat(jr)
                cff=(dcx*(xr-xc))**2 + (dcy*(yr-yc))**2
                if (cff < 1.D0) then
                  wgt(ir-irc,jr-jrc)=( 1.D0-cff*( 1.22876383367175D0
     & -0.22876383367175D0*cff ))**2
                endif
                sum =sum + wgt(ir-irc,jr-jrc)
                biasX=biasX + wgt(ir-irc,jr-jrc)*xr
                biasY=biasY + wgt(ir-irc,jr-jrc)*yr
              enddo
            enddo

            if (sum > 0.D0) then
              cx=dx_inv*(xc-biasX/sum) ! initial estimate of
              cy=dy_inv*(yc-biasY/sum) ! off-centering errors

              if (abs(cx) > 1.D0 .or. abs(cy) > 1.D0) then
                 write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     & 'failure, cx =', cx, ' cy =',cy, ' istr,iend =',
     & istr,iend, ' jstr,jend =', jstr,jend
                 stop
              endif
            else
              write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     & 'failure, line 483.'
              stop
            endif

! Weight-centering iterations: Although the initial setting of "wgt"
!----------------- ----------- is centered already close to (xi,eta)
! defined as the fractional distance between points (ir:ir+1,jr:jr+1),
! this centering may not be accurate enough if too few points are
! participating in averaging. The iterative procedure below compares
! the actual location of the center of gravity with the desired one,
! (xi,eta), and modifies the weights using a first-order upstream
! advection scheme. In essence, it computes the fractional shift
! (cx,cy) which moves the center of gravity exactly to (xi,eta).

            do iter=1,6
              if (cx > 0.D0) then
                iups=+1 ; acx=cx ; imin=istr ; imax=iend+1
              else
                iups=-1 ; acx=-cx ; imin=istr-1 ; imax=iend
              endif
              if (cy > 0.D0) then
                jups=+1 ; acy=cy ; jmin=jstr ; jmax=jend+1
              else
                jups=-1 ; acy=-cy ; jmin=jstr-1 ; jmax=jend
              endif

              ctr=(1.D0-acx)*(1.D0-acy) ; dcx=acx*(1.D0-acy)
              dcy=(1.D0-acx)* acy ; cff=acx * acy

              sum=0.D0 ; biasX=0.D0 ; biasY=0.D0

              do jr=jmin,jmax
                do ir=imin,imax
                  itm=ir-irc ; jtm=jr-jrc

                  wgt_new(itm,jtm)=ctr*wgt(itm ,jtm )
     & +dcx*wgt(itm-iups,jtm )
     & +dcy*wgt(itm ,jtm-jups)
     & +cff*wgt(itm-iups,jtm-jups)

                  sum = sum + wgt_new(itm,jtm)
                  biasX=biasX + wgt_new(itm,jtm)*xlon(ir)
                  biasY=biasY + wgt_new(itm,jtm)*ylat(jr)
                enddo
              enddo

              if (sum > 0.D0) then
                dcx=dx_inv*(xc-biasX/sum) ! normalized residual
                dcy=dy_inv*(yc-biasY/sum) ! off-centering errors

c** write(*,'(I2,4F19.16)') iter,cx,dcx, cy,dcy

                cx=cx+dcx ; cy=cy+dcy !<-- adjust cx,cy

                if (abs(cx) > 1.D0 .or. abs(cy) > 1.D0) then
                   write(*,*) '### ERROR: compute_hraw :: cx =',cx,
     & ' cy =',cy
                   stop
                endif
              else
                write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     & 'failure at line 542.'
                stop
              endif ! After this adjustment
            enddo !<-- iter ! "wgt_new" is centered
                                               ! exactly at (xi,eta);
            sum=0.D0 ; biasX=0.D0 ; biasY=1.D0 ! use it to interpolate
            do jr=jmin,jmax ! htopo--> hraw.
              do ir=imin,imax
                if (htopo(ir,jr) < 1.E+32) then
                  sum = sum +wgt_new(ir-irc,jr-jrc)
                  biasX=biasX+wgt_new(ir-irc,jr-jrc)*dble(htopo(ir,jr))
                else
                  biasY=-1.D0
                endif
              enddo
            enddo
            if (biasY < 0.D0) then
              h(i,j)=1.D+33
            elseif (sum > 0.D0) then
              h(i,j)=biasX/sum
            else
              write(*,*) '### ERROR: compute_hraw :: Algorithm ',
     & 'failure at line 562.'
              stop
            endif

          else ! <-- .not.bounded
            h(i,j)=0.D0
          endif
        enddo ! <-- i
# 601 "compute_hraw.F"
        do i=0,Lm+1
            h(i,j)=-h(i,j) ! Convert into CROCO standart

                            ! but do not set h<0 to 0
                            ! as it can serve when doing wet_and_dry
                                         ! Convert topographic data
! if (h(i,j) < 0.D0) then ! into ROMS standard: ETOPO
! h(i,j)=-h(i,j) ! and SRTM conventions imply
! else ! that "htopo" is elevation,
! h(i,j)=0.D0 ! hence it is positive above
! endif ! and negative below sea
        enddo ! level. ROMS convention
      enddo !<-- j, while(.true.) ! anticipates positive depth.
! 98 continue
C$OMP BARRIER

      if (trd == 0) then
        call system_clock(iclk_end, clk_rate, clk_max)
        if (clk_rate>0) then
          write(*,'(/ /1x,2A,F8.2,1x,A,I4,1x,A)') 'Wall Clock time ',
     & 'spent to compute hraw', (iclk_end-iclk_start)/dble(clk_rate),
     & 'sec running', numthreads, 'threads.'
        endif
      endif


      contains
        integer function indx_bound2(x,n, x0)
        integer n, i ! bounded by
        real(kind=8) x(n), x0 ! x(i) <= x0 < x(i+1)

        if (x0 < x(1)) then
          i=0 ! if x0 is outside the full range
        elseif (x0 > x(n)) then ! of x(1) ... x(n), then return
          i=n ! i=0 or i=n.
        else
          i=int( ( x(n)-x0 +n*(x0-x(1)) )/(x(n)-x(1)) )
          if (x(i+1)<x0) then
            do while (x(i+1) < x0) ! This algorithm computes "i" as
              i=i+1 ! linear interpolation between x(1)
            enddo ! and x(n) which should yield the
          elseif (x(i) > x0) then ! correct value for "i" right a way
            do while (x(i) > x0) ! because array elements x(i) are
              i=i-1 ! equidistantly spaced. The while
            enddo ! loops are here merely to address
          endif ! possible roundoff errors.

          if (x(i+1)-x0 < 0 .or. x0-x(i) < 0) then
            write(*,'(1x,A,5F12.6)') '### ERROR: indx_bound :: ',
     & x(i), x0, x(i+1), x0-x(i), x(i+1)-x0
            stop
          endif
        endif
        indx_bound2=i
        end function indx_bound2

      end subroutine compute_hraw


      subroutine roms_grid_geo_bounds(lon_r,lat_r, Lm,Mm,radius,
     & lon_min,lon_max,lat_min,lat_max)

! Determine geographical limits of data needed to generate ROMS grid
! topography. Algorithmically this is done by going along the perimeter
! of the grid and slightly extrapolating nx-ny coordinates outside to
! accommodate for the data points needed for averaging to compute the
! topography for the points on the perimeter of ROMS grid. The distance
! for extrapolation depends on ROMS grid spacing and the smoothing
! "radius" (with an extra 50% margin for safety of the logic).
!
! input: lon_r,lat_r, Lm,Mm, radius,
! output: lon_min,lon_max,lat_min,lat_max

! use comm_vars_hraw
      implicit none
      integer Lm,Mm, i,j
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: lon_r,lat_r
      real(kind=8) radius, cff, lon_ext, lat_ext

      real(kind=8) :: lon_min,lon_max, lat_min,lat_max
      integer :: comm_j1_indx, comm_j2_indx

cf2py intent(in) lon_r,lat_r, Lm,Mm,radius
cf2py intent(out) lon_min,lon_max,lat_min,lat_max
      comm_j1_indx=0
      comm_j2_indx=0
      cff=1.5D0*radius !<-- safety margin

      lon_min=lon_r(0,0)+cff*(lon_r(0,0)-lon_r(1,1))
      lon_max=lon_min
      lat_min=lat_r(0,0)+cff*(lat_r(0,0)-lat_r(1,1))
      lat_max=lat_min

      do i=1,Lm
        lon_ext=lon_r(i,0)+cff*(lon_r(i,0)-lon_r(i,1))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(i,0)+cff*(lat_r(i,0)-lat_r(i,1))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      lon_ext=lon_r(Lm+1,0)+cff*(lon_r(Lm+1,0)-lon_r(Lm,1))
      if (lon_ext<lon_min) then
        lon_min=lon_ext
      elseif (lon_ext>lon_max) then
        lon_max=lon_ext
      endif

      lat_ext=lat_r(Lm+1,0)+cff*(lat_r(Lm+1,0)-lat_r(Lm,1))
      if (lat_ext<lat_min) then
        lat_min=lat_ext
      elseif (lat_ext>lat_max) then
        lat_max=lat_ext
      endif

      do j=1,Mm
        lon_ext=lon_r(Lm+1,j)+cff*(lon_r(Lm+1,j)-lon_r(Lm,j))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(Lm+1,j)+cff*(lat_r(Lm+1,j)-lat_r(Lm,j))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      lon_ext=lon_r(Lm+1,Mm+1)+cff*(lon_r(Lm+1,Mm+1)-lon_r(Lm,Mm))
      if (lon_ext<lon_min) then
        lon_min=lon_ext
      elseif (lon_ext>lon_max) then
        lon_max=lon_ext
      endif

      lat_ext=lat_r(Lm+1,Mm+1)+cff*(lat_r(Lm+1,Mm+1)-lat_r(Lm,Mm))
      if (lat_ext<lat_min) then
        lat_min=lat_ext
      elseif (lat_ext>lat_max) then
        lat_max=lat_ext
      endif

      do i=Lm,1,-1
        lon_ext=lon_r(i,Mm+1)+cff*(lon_r(i,Mm+1)-lon_r(i,Mm))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(i,Mm+1)+cff*(lat_r(i,Mm+1)-lat_r(i,Mm))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      lon_ext=lon_r(0,Mm+1)+cff*(lon_r(0,Mm+1)-lon_r(1,Mm))
      if (lon_ext<lon_min) then
        lon_min=lon_ext
      elseif (lon_ext>lon_max) then
        lon_max=lon_ext
      endif

      lat_ext=lat_r(0,Mm+1)+cff*(lat_r(0,Mm+1)-lat_r(1,Mm))
      if (lat_ext<lat_min) then
        lat_min=lat_ext
      elseif (lat_ext>lat_max) then
        lat_max=lat_ext
      endif

      do j=Mm,1,-1
        lon_ext=lon_r(0,j)+cff*(lon_r(0,j)-lon_r(1,j))
        if (lon_ext<lon_min) then
          lon_min=lon_ext
        elseif (lon_ext>lon_max) then
          lon_max=lon_ext
        endif

        lat_ext=lat_r(0,j)+cff*(lat_r(0,j)-lat_r(1,j))
        if (lat_ext<lat_min) then
          lat_min=lat_ext
        elseif (lat_ext>lat_max) then
          lat_max=lat_ext
        endif
      enddo

      end subroutine roms_grid_geo_bounds
# 25 "tools_fort.F" 2

!!! smooth !!!
# 1 "topo_smooth_subs.F" 1



      subroutine smooth_thread(Lm,Mm, msk, h, Lgh,hmin,hmax,r_max,
     & method,Lgh1)

      implicit none
      integer Lm,Mm
      character(len=64) method
      character(len=64), dimension(4) :: log_list

      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: h,Lgh,Lgh1
      real(kind=8), allocatable, dimension(:) :: FX,FE, FX1,FE1
      real(kind=8) :: hmin,lcl_hmin, hmax, r_max

      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      integer istr,iend,jstr,jend, i,j, iter,ifrst_call
      integer numthreads,trd, nsub_x,nsub_y, isize,jsize,
     & tile, my_first,my_last, range, size
      integer iters_cond, iters_lin
      integer(kind=2) :: transform_to_log

C$ integer omp_get_num_threads, omp_get_thread_num

      integer(kind=4) iclk_start, iclk_end, clk_rate, clk_max

      log_list = [character(len=64) :: 'LOG_SMOOTHING',
     & 'LEGACY_LOG_SMOOTH','LOG_SMOOTH_2','LOG_SMOOTH_1']
      transform_to_log=0
      iters_cond=0
      iters_lin=0
      numthreads=1 ; trd=0
      ifrst_call=1
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()

      if (trd == 0) then

        call system_clock(iclk_start, clk_rate, clk_max)

        iters_cond=500 !<-- set number of iterations
        write(*,'(1x,2A,F10.7,3(2x,A,I5))') 'enter smooth_thread ',
     & 'r_max =', r_max, 'iters_cond =', iters_cond,
     & 'iters_lin =', iters_lin
C$ & , 'numthreads =', numthreads
      endif

      call set_tiles(Lm,Mm, nsub_x,nsub_y)

c** nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      isize=(Lm+nsub_x-1)/nsub_x ; jsize=(Mm+nsub_y-1)/nsub_y
      size=(isize+4)*(jsize+4)
      allocate(FX(size),FE(size), FX1(size),FE1(size))

      write(*,'(5(2x,A,I4))') 'nsub_x =',nsub_x, 'nsub_y =',nsub_y,
     & 'isize =',isize, 'jsize =',jsize, 'trd =',trd

      do tile=my_first,my_last
        call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

        if (istr == 1) istr=istr-1 ; if (iend == Lm) iend=iend+1
        if (jstr == 1) jstr=jstr-1 ; if (jend == Mm) jend=jend+1

        do j=jstr,jend
          do i=istr,iend
            h(i,j)=min( hmax, max( hmin, h(i,j) ))
            if ( hmin<=0 ) then
                  h(i,j)=h(i,j)-hmin+0.1 ! raise topo when hmin<=0 to avoid log error
                  lcl_hmin=0.1
              else
                  lcl_hmin=hmin
              endif

            if ( any( log_list == method ) ) then
              transform_to_log=1
              if (h(i,j) > lcl_hmin) then
                Lgh(i,j)=log(h(i,j)/lcl_hmin)
              else
                Lgh(i,j)=0.D0
              endif
            endif
          enddo
        enddo
      enddo
C$OMP BARRIER

! Conditional smoothing or log-smoothing to reduce r-factors

      do iter=1,iters_cond
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          call rx_diag_tile(istr,iend,jstr,jend, Lm,Mm, h,msk,
     & iter-1, nsub_x*nsub_y,ifrst_call)
          ifrst_call=0
          if ( method == 'LOG_SMOOTHING') then
              call lsmooth_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh1, FX,FE,FX1,FE1)

          elseif (method == 'LOG_SMOOTH_2') then
            call lsmooth_2_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh1, FX,FE,FE1)

          elseif (method == 'LOG_SMOOTH_1') then
            call lsmooth_1_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh1, FX,FE)

          elseif (method == 'LEGACY_LOG_SMOOTH') then
            call lsmth_legacy_tile(istr,iend,jstr,jend, Lm,Mm,
     & r_max, Lgh,Lgh1, FX,FE,FX1,FE1)

          else
            call cond_smth_tile(istr,iend,jstr,jend, Lm,Mm,
     & r_max, msk, h, Lgh1, FX,FE,FX1,FE1)
          endif
        enddo
C$OMP BARRIER

        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          if (method == 'LOG_SMOOTHING') then
            call land_lsmth_tile (istr,iend,jstr,jend, Lm,Mm,
     & r_max, msk, Lgh1,Lgh, FX,FE,FX1,FE1)

          elseif (method == 'LOG_SMOOTH_2') then
            call land_lsmth_2_tile(istr,iend,jstr,jend, Lm,Mm,
     & r_max, msk, Lgh,Lgh1, FX,FE,FE1)

          elseif (method == 'LOG_SMOOTH_1') then
            call land_lsmth_1_tile(istr,iend,jstr,jend, Lm,Mm,
     & r_max, msk, Lgh,Lgh1, FX,FE)

          elseif (method == 'LEGACY_LOG_SMOOTH') then
            call lsmth_legacy_tile(istr,iend,jstr,jend, Lm,Mm,
     & r_max, Lgh1,Lgh, FX,FE,FX1,FE1)

          else
            call land_cnd_smth_tile(istr,iend,jstr,jend, Lm,Mm,
     & r_max, msk, Lgh1,h, FX,FE, FX1,FE1)
          endif

          if (transform_to_log == 1) then
            if (istr == 1) istr=istr-1 ; if (iend == Lm) iend=iend+1
            if (jstr == 1) jstr=jstr-1 ; if (jend == Mm) jend=jend+1

            do j=jstr,jend ! this backward conversion is
              do i=istr,iend ! needed because "rx_diag_tile"
                h(i,j)=lcl_hmin*exp(Lgh(i,j)) ! expects "h" to report "r_max"
              enddo ! values achieved at every
            enddo ! iteration.
          endif

        enddo !<-- tile
C$OMP BARRIER
      enddo !<-- iter

      if ( hmin<=0 ) then
        h=h+hmin-0.1 ! Deepens topo if hmin<0 as we previously raise it
      endif


! Final linear smoothing

      if (iters_lin > 0 .and. trd == 0) then
        write(*,*) 'Applying linear smoothing.'
      endif

      do iter=1,iters_lin
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          call smooth_tile(istr,iend,jstr,jend, Lm,Mm, msk,
     & h,Lgh1, FX,FE,FE1)
        enddo
C$OMP BARRIER

        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          call land_smth_tile(istr,iend,jstr,jend, Lm,Mm,
     & msk, Lgh1,h, FX,FE,FE1)
        enddo
C$OMP BARRIER

        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, Lm,Mm, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          call rx_diag_tile(istr,iend,jstr,jend, Lm,Mm,
     & h,msk, iter, nsub_x*nsub_y)

        enddo
C$OMP BARRIER
      enddo !<-- iter

      if (trd == 0) then
        call system_clock(iclk_end, clk_rate, clk_max)
        if (clk_rate > 0) then
          write(*,'(/1x,2A,F8.2,1x,A,I4,1x,A)') 'Wall Clock time ',
     & 'spent in smoothing', (iclk_end-iclk_start)/dble(clk_rate),
     & 'sec running', numthreads, 'threads.'
        endif
      endif

      end
# 28 "tools_fort.F" 2
# 1 "tools_topo.F" 1
! Land-mask-dependent smoothing operators imply Neumann (no-flux) b.c.
! on the coastline from the water side: that is any point adjacent to
! the coast line is allowed to be changed without regard of the value
! of its nearest neighbor inside the land. This is accomplished by
! applying the usual U- and V-type masking rules to the straight
! fluxes, FX and FE. Masking for diagonal fluxes FX1 and FE1 imposes
! an additional condition: its end points must be both water points,
! AND MUST NOT BE SEPARATED BY LAND, e.g.,
!
! FX1(i,j)=h(i,j)-h(i-1,j-1)
!
! is unmasked only if (1) both (i,j) and (i-1,j-1) are water points;
! AND (2) at least one of (i-1,j) or (i,j-1) is water
! as well -- thus there is exist a two-step
! 90-degree passage from (i,j) to (i-1,j-1)
! by water points.
!
! This explains the somewhat awkward-looking expressions like
!
! FX1(i,j)=(h(i,j)-h(i-1,j-1)) *dble( msk(i,j)*msk(i-1,j-1)
! & *max(msk(i-1,j),msk(i,j-1)) )

! Because ROMS relies to some degree on continuity of the topography
! across coastline (spatial interpolations may reach grid-box heights
! at least one point inside land), each routine below is accompanied
! by a land_smooth_something routine generated by CPP from nearly the
! same source code), which changes values inside land only.
! These routines use essentially Dirichlet-type boundary conditions:
! their fluxes are unmasked, and from the land-side point of view the
! points outside the coastline (water points nearest to the coastline)
! are fixed and used as source for Dirichlet side boundaries.
! So, as a rule of thumb, when processing water points h-values inside
! land have no influence; in contrast, points inside land are
! completely slaved to the adjacent water points.

! Without masking the stencil of the smoothing operators have the
! following weights, depending on parameter settings in the code below
!
! 1/8 1/32 1/8 1/32 1/16 1/8 1/16
!
! 1/8 1/2 1/8 1/8 3/8 1/8 1/8 1/4 1/8
!
! 1/8 1/32 1/8 1/32 1/16 1/8 1/16
!
! 5-point isotropic two-dimensional
! Laplacian Laplacian 1-2-1
! smoother smoother Hann window
!
! Note that all 5 or 9 coefficients add up to 1. All three operators
! suppress the checkerboard mode in just a single iteration, however,
! only the last one eliminates flat-front 2dx-mode in one iteration;
! the first and the second attenuate 1D 2dx-mode by factors of 1/2 and
! 1/4 respectively.

!!! LOG_SMOOTHING !!!

      subroutine lsmooth_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh_new, FX,FE,FX1,FE1)




      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: Lgh,Lgh_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) ::
     & FX,FE,FX1,FE1
      real(kind=8) r_max, lgr_max, lgr1_max, grad
      real(kind=8), parameter :: OneEights=0.125D0

      if (r_max > 0.D0) then ! Set threshold
        lgr_max=log((1.D0+r_max)/(1.D0-r_max)) ! for logarithmic
        lgr1_max=lgr_max*sqrt(2.D0) ! slopes
      else
        lgr_max=0.D0 ; lgr1_max=0.D0
      endif

      do j=jstr,jend
        do i=istr,iend+1

          grad=(Lgh(i,j)-Lgh(i-1,j))*dble(msk(i,j)*msk(i-1,j))



          if (grad > lgr_max) then
            FX(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FX(i,j)=grad +lgr_max
          else
            FX(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend

          grad=(Lgh(i,j)-Lgh(i,j-1))*dble(msk(i,j)*msk(i,j-1))



          if (grad > lgr_max) then
            FE(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FE(i,j)=grad +lgr_max
          else
            FE(i,j)=0.D0
          endif
        enddo
        do i=istr,iend+1

          grad=(Lgh(i,j)-Lgh(i-1,j-1)) *dble( msk(i,j)*msk(i-1,j-1)
     & *max(msk(i-1,j),msk(i,j-1)) )



          if (grad > lgr1_max) then
            FX1(i,j)=grad-lgr1_max
          elseif (grad < -lgr1_max) then
            FX1(i,j)=grad+lgr1_max
          else
            FX1(i,j)=0.D0
          endif
        enddo
        do i=istr-1,iend

          grad=(Lgh(i,j)-Lgh(i+1,j-1)) *dble( msk(i,j)*msk(i+1,j-1)
     & *max(msk(i+1,j),msk(i,j-1)) )



          if (grad > lgr1_max) then
            FE1(i,j)=grad-lgr1_max
          elseif (grad < -lgr1_max) then
            FE1(i,j)=grad+lgr1_max
          else
            FE1(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend

          if (msk(i,j) > 0) then !--> water points



            Lgh_new(i,j)=Lgh(i,j) + OneEights*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j)
     & +0.25D0*( FX1(i+1,j+1)-FX1(i,j)
     & +FE1(i-1,j+1)-FE1(i,j) ))
          else
            Lgh_new(i,j)=Lgh(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,Lgh_new)
      end

! Same as above, but instead of explicitly computing diagonal fluxes
! FX1,FE1 apply cross-averaging to 90-degree fluxes. Linear version of
! the resultant operator is equivalent to isotropic Laplacian, however
! r_max threshold is applied differently, and only to construct the
! initial 90-degree slops.



!!! LOG_SMOOTH_2 !!!

      subroutine lsmooth_2_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh_new, FX,FE,FE1)




      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: Lgh,Lgh_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: FX,FE,FE1
      real(kind=8) r_max, lgr_max, grad
      real(kind=8), parameter :: ThreeSixteenth=3.D0/16.D0,
     & OneTwelfth=1.D0/12.D0

      if (r_max > 0.D0) then ! Set threshold
        lgr_max=log((1.D0+r_max)/(1.D0-r_max)) ! for logarithmic
      else ! slopes
        lgr_max=0.D0
      endif

      do j=jstr-1,jend+1
        do i=istr,iend+1

          grad=(Lgh(i,j)-Lgh(i-1,j))*dble(msk(i,j)*msk(i-1,j))



          if (grad > lgr_max) then
            FX(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FX(i,j)=grad +lgr_max
          else
            FX(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1

          grad=(Lgh(i,j)-Lgh(i,j-1))*dble(msk(i,j)*msk(i,j-1))



          if (grad > lgr_max) then
            FE1(i,j)=grad-lgr_max
          elseif (grad < -lgr_max) then
            FE1(i,j)=grad+lgr_max
          else
            FE1(i,j)=0.D0
          endif
        enddo
        do i=istr,iend
          FE(i,j)=FE1(i,j) + OneTwelfth*( FX(i+1,j)+FX(i ,j-1)
     & -FX(i ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=FX(i,j) + OneTwelfth*( FE1(i,j+1)+FE1(i-1,j )
     & -FE1(i,j )-FE1(i-1,j+1))
        enddo
        do i=istr,iend

          if (msk(i,j) > 0) then !--> water points



            Lgh_new(i,j)=Lgh(i,j) +ThreeSixteenth*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j))
          else
            Lgh_new(i,j)=Lgh(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,Lgh_new)
      end


! The same as above, but without diagonal fluxes.

!!! LOG_SMOOTH_1 !!!

      subroutine lsmooth_1_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh_new, FX,FE)




      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: Lgh, Lgh_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: FX,FE
      real(kind=8) r_max, lgr_max, grad
      real(kind=8), parameter :: OneEighth=0.125D0

      if (r_max > 0.D0) then ! Set threshold
        lgr_max=log((1.D0+r_max)/(1.D0-r_max)) ! for logarithmic
      else ! slopes
        lgr_max=0.D0
      endif

      do j=jstr,jend+1
        do i=istr,iend

          grad=(Lgh(i,j)-Lgh(i,j-1))*dble(msk(i,j)*msk(i,j-1))



          if (grad > lgr_max) then
            FE(i,j)=grad-lgr_max
          elseif (grad < -lgr_max) then
            FE(i,j)=grad+lgr_max
          else
            FE(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1

          grad=(Lgh(i,j)-Lgh(i-1,j))*dble(msk(i,j)*msk(i-1,j))



          if (grad > lgr_max) then
            FX(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FX(i,j)=grad +lgr_max
          else
            FX(i,j)=0.D0
          endif
        enddo
        do i=istr,iend

          if (msk(i,j) > 0) then !--> water points



            Lgh_new(i,j)=Lgh(i,j) + OneEighth*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j))
          else
            Lgh_new(i,j)=Lgh(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm, Lgh_new)
      end




!!! LEGACY_LOG_SMOOTH !!!


! Mask-independent logarithmic r-factor limiting algorithm. Somewhat
! similar to the one from RomsTools by Pierrick Penven, the code below
! uses both straight and diagonal fluxes arranged into isotropic
! Laplacian operator to counter the diamond-shaped spreading patterns
! of single point spikes.

      subroutine lsmth_legacy_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & Lgh,Lgh_new, FX,FE, FX1,FE1)
      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: Lgh,Lgh_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) ::
     & FX,FE,FX1,FE1
      real(kind=8) r_max, lgr_max, lgr1_max, grad
      real(kind=8), parameter :: OneEights=0.125D0

      if (r_max > 0.D0) then ! Set threshold
        lgr_max=log((1.D0+r_max)/(1.D0-r_max)) ! for logarithmic
        lgr1_max=lgr_max*sqrt(2.D0) ! slopes
      else
        lgr_max=0.D0 ; lgr1_max=0.D0
      endif

      do j=jstr,jend
        do i=istr,iend+1
          grad=Lgh(i,j)-Lgh(i-1,j)
          if (grad > lgr_max) then
            FX(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FX(i,j)=grad +lgr_max
          else
            FX(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend
          grad=Lgh(i,j)-Lgh(i,j-1)
          if (grad > lgr_max) then
            FE(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FE(i,j)=grad +lgr_max
          else
            FE(i,j)=0.D0
          endif
        enddo
        do i=istr,iend+1
          grad=Lgh(i,j)-Lgh(i-1,j-1)
          if (grad > lgr1_max) then
            FX1(i,j)=grad-lgr1_max
          elseif (grad < -lgr1_max) then
            FX1(i,j)=grad+lgr1_max
          else
            FX1(i,j)=0.D0
          endif
        enddo
        do i=istr-1,iend
          grad=Lgh(i,j)-Lgh(i+1,j-1)
          if (grad > lgr1_max) then
            FE1(i,j)=grad-lgr1_max
          elseif (grad < -lgr1_max) then
            FE1(i,j)=grad+lgr1_max
          else
            FE1(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          Lgh_new(i,j)=Lgh(i,j) + OneEights*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j)
     & +0.25D0*( FX1(i+1,j+1)-FX1(i,j)
     & +FE1(i-1,j+1)-FE1(i,j) ))
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,Lgh_new)
      end


!!! ELSE !!!

! Physical-space (non-logarithmic) mask-dependent r-factor restriction
! algorithm. It is somewhat suggested in Mellor, Ezer, and Oey (1994)
! "conundrum" article (page 1130, Fig. 6, text on the right, however no
! detail about the specific procedure are provided). Algorithm below
! computes both 90-degree (FX,FE) and diagonal (FX1,FE1) "fluxes" as
! the excesses of gradient over the maximum allowed threshold, e.g.,
!
! h(i,j)-h(i-1,j)
! ----------------- > r_max [assuming that h(i,j)-h(i-1,j) > 0]
! h(i,j)+h(i-1,j)
!
! results in FX(i,j) = h(i,j)-h(i-1,j) -r_max*[h(i,j)+h(i-1,j)]
!
! while FX(i,j) = 0 if |h(i,j)-h(i-1,j)|/[h(i,j)+h(i-1,j)] < r_max
!
! so this procedure will not alter topography where the slopes are
! already within the limits.



      subroutine cond_smth_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, h,h_new, FX,FE, FX1,FE1)




      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: h,h_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2)
     & :: FX,FE, FX1,FE1
      real(kind=8) r_max, r1_max, grad, cr
      real(kind=8), parameter :: OneEights=0.125D0

      r1_max=r_max*sqrt(2.D0)

      do j=jstr,jend
        do i=istr,iend+1

          grad=(h(i,j)-h(i-1,j))*dble(msk(i,j)*msk(i-1,j))



          cr=r_max*(h(i,j)+h(i-1,j))
          if (grad > cr) then
            FX(i,j)=grad -cr
          elseif (grad < -cr) then
            FX(i,j)=grad +cr
          else
            FX(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend

          grad=(h(i,j)-h(i,j-1))*dble(msk(i,j)*msk(i,j-1))



          cr=r_max*(h(i,j)+h(i,j-1))
          if (grad > cr) then
            FE(i,j)=grad -cr
          elseif (grad < -cr) then
            FE(i,j)=grad +cr
          else
            FE(i,j)=0.D0
          endif
        enddo
        do i=istr,iend+1

          grad=(h(i,j)-h(i-1,j-1)) *dble( msk(i,j)*msk(i-1,j-1)
     & *max(msk(i-1,j),msk(i,j-1)) )



          cr=r1_max*(h(i,j)+h(i-1,j-1))
          if (grad > cr) then
            FX1(i,j)=grad -cr
          elseif (grad < -cr) then
            FX1(i,j)=grad +cr
          else
            FX1(i,j)=0.D0
          endif
        enddo
        do i=istr-1,iend

          grad=(h(i,j)-h(i+1,j-1)) *dble( msk(i,j)*msk(i+1,j-1)
     & *max(msk(i+1,j),msk(i,j-1)) )



          cr=r1_max*(h(i,j)+h(i+1,j-1))
          if (grad > cr) then
            FE1(i,j)=grad -cr
          elseif (grad < -cr) then
            FE1(i,j)=grad +cr
          else
            FE1(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend

          if (msk(i,j) > 0) then !--> water points



            h_new(i,j)=h(i,j) + OneEights*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j)
     & +0.25D0*( FX1(i+1,j+1)-FX1(i,j)
     & +FE1(i-1,j+1)-FE1(i,j) ))
          else
            h_new(i,j)=h(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,h_new)
      end
!#endif



      subroutine smooth_tile(istr,iend,jstr,jend, Lm,Mm, msk,
     & h,h_new, FX,FE,FE1)




      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: h,h_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2)
     & :: FX,FE,FE1
      real(kind=8), parameter :: ThreeSixteenth=3.D0/16.D0,
     & OneTwelfth=1.D0/12.D0

      do j=jstr-1,jend+1
        do i=istr,iend+1

          FX(i,j)=(h(i,j)-h(i-1,j))*dble(msk(i,j)*msk(i-1,j))



        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1

          FE1(i,j)=(h(i,j)-h(i,j-1))*dble(msk(i,j)*msk(i,j-1))



        enddo
        do i=istr,iend
          FE(i,j)=FE1(i,j) + OneTwelfth*( FX(i+1,j)+FX(i ,j-1)
     & -FX(i ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=FX(i,j) + OneTwelfth*( FE1(i,j+1)+FE1(i-1,j )
     & -FE1(i,j )-FE1(i-1,j+1))
        enddo
        do i=istr,iend

          if (msk(i,j) > 0) then !--> water points



            h_new(i,j)=h(i,j) + ThreeSixteenth*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j))
          else
            h_new(i,j)=h(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,h_new)
      end



# 1 "tools_topo.F" 1
! Land-mask-dependent smoothing operators imply Neumann (no-flux) b.c.
! on the coastline from the water side: that is any point adjacent to
! the coast line is allowed to be changed without regard of the value
! of its nearest neighbor inside the land. This is accomplished by
! applying the usual U- and V-type masking rules to the straight
! fluxes, FX and FE. Masking for diagonal fluxes FX1 and FE1 imposes
! an additional condition: its end points must be both water points,
! AND MUST NOT BE SEPARATED BY LAND, e.g.,
!
! FX1(i,j)=h(i,j)-h(i-1,j-1)
!
! is unmasked only if (1) both (i,j) and (i-1,j-1) are water points;
! AND (2) at least one of (i-1,j) or (i,j-1) is water
! as well -- thus there is exist a two-step
! 90-degree passage from (i,j) to (i-1,j-1)
! by water points.
!
! This explains the somewhat awkward-looking expressions like
!
! FX1(i,j)=(h(i,j)-h(i-1,j-1)) *dble( msk(i,j)*msk(i-1,j-1)
! & *max(msk(i-1,j),msk(i,j-1)) )

! Because ROMS relies to some degree on continuity of the topography
! across coastline (spatial interpolations may reach grid-box heights
! at least one point inside land), each routine below is accompanied
! by a land_smooth_something routine generated by CPP from nearly the
! same source code), which changes values inside land only.
! These routines use essentially Dirichlet-type boundary conditions:
! their fluxes are unmasked, and from the land-side point of view the
! points outside the coastline (water points nearest to the coastline)
! are fixed and used as source for Dirichlet side boundaries.
! So, as a rule of thumb, when processing water points h-values inside
! land have no influence; in contrast, points inside land are
! completely slaved to the adjacent water points.

! Without masking the stencil of the smoothing operators have the
! following weights, depending on parameter settings in the code below
!
! 1/8 1/32 1/8 1/32 1/16 1/8 1/16
!
! 1/8 1/2 1/8 1/8 3/8 1/8 1/8 1/4 1/8
!
! 1/8 1/32 1/8 1/32 1/16 1/8 1/16
!
! 5-point isotropic two-dimensional
! Laplacian Laplacian 1-2-1
! smoother smoother Hann window
!
! Note that all 5 or 9 coefficients add up to 1. All three operators
! suppress the checkerboard mode in just a single iteration, however,
! only the last one eliminates flat-front 2dx-mode in one iteration;
! the first and the second attenuate 1D 2dx-mode by factors of 1/2 and
! 1/4 respectively.

!!! LOG_SMOOTHING !!!




      subroutine land_lsmth_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh_new, FX,FE,FX1,FE1)

      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: Lgh,Lgh_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) ::
     & FX,FE,FX1,FE1
      real(kind=8) r_max, lgr_max, lgr1_max, grad
      real(kind=8), parameter :: OneEights=0.125D0

      if (r_max > 0.D0) then ! Set threshold
        lgr_max=log((1.D0+r_max)/(1.D0-r_max)) ! for logarithmic
        lgr1_max=lgr_max*sqrt(2.D0) ! slopes
      else
        lgr_max=0.D0 ; lgr1_max=0.D0
      endif

      do j=jstr,jend
        do i=istr,iend+1



          grad=Lgh(i,j)-Lgh(i-1,j)

          if (grad > lgr_max) then
            FX(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FX(i,j)=grad +lgr_max
          else
            FX(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend



          grad= Lgh(i,j)-Lgh(i,j-1)

          if (grad > lgr_max) then
            FE(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FE(i,j)=grad +lgr_max
          else
            FE(i,j)=0.D0
          endif
        enddo
        do i=istr,iend+1




          grad= Lgh(i,j)-Lgh(i-1,j-1)

          if (grad > lgr1_max) then
            FX1(i,j)=grad-lgr1_max
          elseif (grad < -lgr1_max) then
            FX1(i,j)=grad+lgr1_max
          else
            FX1(i,j)=0.D0
          endif
        enddo
        do i=istr-1,iend




          grad= Lgh(i,j)-Lgh(i+1,j-1)

          if (grad > lgr1_max) then
            FE1(i,j)=grad-lgr1_max
          elseif (grad < -lgr1_max) then
            FE1(i,j)=grad+lgr1_max
          else
            FE1(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend



          if (msk(i,j) < 1) then !--> on land only

            Lgh_new(i,j)=Lgh(i,j) + OneEights*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j)
     & +0.25D0*( FX1(i+1,j+1)-FX1(i,j)
     & +FE1(i-1,j+1)-FE1(i,j) ))
          else
            Lgh_new(i,j)=Lgh(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,Lgh_new)
      end

! Same as above, but instead of explicitly computing diagonal fluxes
! FX1,FE1 apply cross-averaging to 90-degree fluxes. Linear version of
! the resultant operator is equivalent to isotropic Laplacian, however
! r_max threshold is applied differently, and only to construct the
! initial 90-degree slops.



!!! LOG_SMOOTH_2 !!!




      subroutine land_lsmth_2_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh_new, FX,FE,FE1)

      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: Lgh,Lgh_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: FX,FE,FE1
      real(kind=8) r_max, lgr_max, grad
      real(kind=8), parameter :: ThreeSixteenth=3.D0/16.D0,
     & OneTwelfth=1.D0/12.D0

      if (r_max > 0.D0) then ! Set threshold
        lgr_max=log((1.D0+r_max)/(1.D0-r_max)) ! for logarithmic
      else ! slopes
        lgr_max=0.D0
      endif

      do j=jstr-1,jend+1
        do i=istr,iend+1



          grad= Lgh(i,j)-Lgh(i-1,j)

          if (grad > lgr_max) then
            FX(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FX(i,j)=grad +lgr_max
          else
            FX(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1



          grad= Lgh(i,j)-Lgh(i,j-1)

          if (grad > lgr_max) then
            FE1(i,j)=grad-lgr_max
          elseif (grad < -lgr_max) then
            FE1(i,j)=grad+lgr_max
          else
            FE1(i,j)=0.D0
          endif
        enddo
        do i=istr,iend
          FE(i,j)=FE1(i,j) + OneTwelfth*( FX(i+1,j)+FX(i ,j-1)
     & -FX(i ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=FX(i,j) + OneTwelfth*( FE1(i,j+1)+FE1(i-1,j )
     & -FE1(i,j )-FE1(i-1,j+1))
        enddo
        do i=istr,iend



          if (msk(i,j) < 1) then !--> on land only

            Lgh_new(i,j)=Lgh(i,j) +ThreeSixteenth*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j))
          else
            Lgh_new(i,j)=Lgh(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,Lgh_new)
      end


! The same as above, but without diagonal fluxes.

!!! LOG_SMOOTH_1 !!!




      subroutine land_lsmth_1_tile(istr,iend,jstr,jend, Lm,Mm, r_max,
     & msk, Lgh,Lgh_new, FX,FE)

      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: Lgh, Lgh_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2) :: FX,FE
      real(kind=8) r_max, lgr_max, grad
      real(kind=8), parameter :: OneEighth=0.125D0

      if (r_max > 0.D0) then ! Set threshold
        lgr_max=log((1.D0+r_max)/(1.D0-r_max)) ! for logarithmic
      else ! slopes
        lgr_max=0.D0
      endif

      do j=jstr,jend+1
        do i=istr,iend



          grad= Lgh(i,j)-Lgh(i,j-1)

          if (grad > lgr_max) then
            FE(i,j)=grad-lgr_max
          elseif (grad < -lgr_max) then
            FE(i,j)=grad+lgr_max
          else
            FE(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1



          grad= Lgh(i,j)-Lgh(i-1,j)

          if (grad > lgr_max) then
            FX(i,j)=grad -lgr_max
          elseif (grad < -lgr_max) then
            FX(i,j)=grad +lgr_max
          else
            FX(i,j)=0.D0
          endif
        enddo
        do i=istr,iend



          if (msk(i,j) < 1) then !--> on land only

            Lgh_new(i,j)=Lgh(i,j) + OneEighth*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j))
          else
            Lgh_new(i,j)=Lgh(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm, Lgh_new)
      end




!!! LEGACY_LOG_SMOOTH !!!
# 405 "tools_topo.F"
!!! ELSE !!!

! Physical-space (non-logarithmic) mask-dependent r-factor restriction
! algorithm. It is somewhat suggested in Mellor, Ezer, and Oey (1994)
! "conundrum" article (page 1130, Fig. 6, text on the right, however no
! detail about the specific procedure are provided). Algorithm below
! computes both 90-degree (FX,FE) and diagonal (FX1,FE1) "fluxes" as
! the excesses of gradient over the maximum allowed threshold, e.g.,
!
! h(i,j)-h(i-1,j)
! ----------------- > r_max [assuming that h(i,j)-h(i-1,j) > 0]
! h(i,j)+h(i-1,j)
!
! results in FX(i,j) = h(i,j)-h(i-1,j) -r_max*[h(i,j)+h(i-1,j)]
!
! while FX(i,j) = 0 if |h(i,j)-h(i-1,j)|/[h(i,j)+h(i-1,j)] < r_max
!
! so this procedure will not alter topography where the slopes are
! already within the limits.






      subroutine land_cnd_smth_tile(istr,iend,jstr,jend, Lm,Mm,
     & r_max, msk, h,h_new, FX,FE, FX1,FE1)

      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: h,h_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2)
     & :: FX,FE, FX1,FE1
      real(kind=8) r_max, r1_max, grad, cr
      real(kind=8), parameter :: OneEights=0.125D0

      r1_max=r_max*sqrt(2.D0)

      do j=jstr,jend
        do i=istr,iend+1



          grad= h(i,j)-h(i-1,j)

          cr=r_max*(h(i,j)+h(i-1,j))
          if (grad > cr) then
            FX(i,j)=grad -cr
          elseif (grad < -cr) then
            FX(i,j)=grad +cr
          else
            FX(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend+1
        do i=istr,iend



          grad= h(i,j)-h(i,j-1)

          cr=r_max*(h(i,j)+h(i,j-1))
          if (grad > cr) then
            FE(i,j)=grad -cr
          elseif (grad < -cr) then
            FE(i,j)=grad +cr
          else
            FE(i,j)=0.D0
          endif
        enddo
        do i=istr,iend+1




          grad= h(i,j)-h(i-1,j-1)

          cr=r1_max*(h(i,j)+h(i-1,j-1))
          if (grad > cr) then
            FX1(i,j)=grad -cr
          elseif (grad < -cr) then
            FX1(i,j)=grad +cr
          else
            FX1(i,j)=0.D0
          endif
        enddo
        do i=istr-1,iend




          grad=h(i,j)-h(i+1,j-1)

          cr=r1_max*(h(i,j)+h(i+1,j-1))
          if (grad > cr) then
            FE1(i,j)=grad -cr
          elseif (grad < -cr) then
            FE1(i,j)=grad +cr
          else
            FE1(i,j)=0.D0
          endif
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend



          if (msk(i,j) < 1) then !--> on land only

            h_new(i,j)=h(i,j) + OneEights*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j)
     & +0.25D0*( FX1(i+1,j+1)-FX1(i,j)
     & +FE1(i-1,j+1)-FE1(i,j) ))
          else
            h_new(i,j)=h(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,h_new)
      end
!#endif






      subroutine land_smth_tile(istr,iend,jstr,jend, Lm,Mm, msk,
     & h,h_new, FX,FE,FE1)

      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: h,h_new
      real(kind=8), dimension(istr-2:iend+2,jstr-2:jend+2)
     & :: FX,FE,FE1
      real(kind=8), parameter :: ThreeSixteenth=3.D0/16.D0,
     & OneTwelfth=1.D0/12.D0

      do j=jstr-1,jend+1
        do i=istr,iend+1



          FX(i,j)= h(i,j)-h(i-1,j)

        enddo
      enddo
      do j=jstr,jend+1
        do i=istr-1,iend+1



          FE1(i,j)=h(i,j)-h(i,j-1)

        enddo
        do i=istr,iend
          FE(i,j)=FE1(i,j) + OneTwelfth*( FX(i+1,j)+FX(i ,j-1)
     & -FX(i ,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend+1
          FX(i,j)=FX(i,j) + OneTwelfth*( FE1(i,j+1)+FE1(i-1,j )
     & -FE1(i,j )-FE1(i-1,j+1))
        enddo
        do i=istr,iend



          if (msk(i,j) < 1) then !--> on land only

            h_new(i,j)=h(i,j) + ThreeSixteenth*( FX(i+1,j)-FX(i,j)
     & +FE(i,j+1)-FE(i,j))
          else
            h_new(i,j)=h(i,j)
          endif
        enddo
      enddo
      call neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,h_new)
      end
# 593 "tools_topo.F" 2
# 29 "tools_fort.F" 2
# 1 "smooth.F" 1
! Purpose: to smooth or "log-smooth" or rx-condition model topography
!--------- in order to achieve acceptable r-factors,
!
! | h(i+1) - h(i) |
! rx = ------------------- < rx_max
! h(i+1) + h(i)
!
! It reads raw topography "hraw" from netCDF grid file calls the actual
! smoothing routine. Upon completion it puts the resultant smoothed
! topography into netCDF variable "hsmth" (if exists) or into "h" of
! the same file leaving the original "hraw" unchanged.

! There are multiple versions of such smoothing/rx-limiting procedures
! which can be subdivided into two major groups: (i) the ones which are
! applied directly to "h" (so called rx-"smooth"ers); and (ii) the ones
! with first transform "h" into its logarithm,
!
! Lgh(i,j) = log[h(i,j)/hmin]
!
! and smooth it trying to achieve condition
!
! | Lgh(i+1,j) - Lgh(i+1,j) | < log[(1 + r_max)/(1 - r_max)]
!
! which is equivalent to rx < rx_max condition above, then transform
! it back, h(i,j)=hmin*exp[Lgh(i,j)]. These are called log-smoothers,
! or "lmsooth". Besides this distinctions, there are others:
! handling boundary conditions at coast line; choice of discrete
! operator, etc. - resulting in a fair number of options.

! Usage: "smooth" or "lsmooth" take four arguments:
!
! lsmooth hmin hmax r_max file.nc
! where
! #1 hmin desired minimum depth limit, meters, positive;
! #2 hmax --/ /-- maximum --/ /--
! #3 r_max desired maximum r-factor, nondimensional;
! #4- name of ROMS grid netCDF file;
!
! This file is merely a driver which decodes command-line arguments,
! reads model grid topography, calls relevant subroutine selected by
! CPP-switches (defined from compiler line to support multiple
! versions), and creates appropriate signature in the resultant grid
! file so all the parameters of execution these operators are
! documented to allow exact reproduction of the result relying
! solely on what is stored in netCDF file.
# 58 "smooth.F"
      subroutine smooth (Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)


      integer Lm,Mm
      character(len=80) :: grid
      character(len=64) method
      character(len=16) str
      real(kind=8) :: hmin, hmax, r_max
      real(kind=8) :: hraw(1:Lm+2,1:Mm+2)
      real(kind=8) :: h(1:Lm+2,1:Mm+2)
      real(kind=8) :: wrk(1:Lm+2,1:Mm+2)
      real(kind=8), dimension(:,:), allocatable :: wrk1

      integer iters_cond, iters_lin
      integer(kind=2), allocatable, dimension(:,:) :: msk
      integer nargs, ncx,ncy, i,j, ierr, ncid, lgrd


      include "netcdf.inc"

cf2py intent(in) hraw,hmin,hmax,r_max,wrk
cf2py intent(out) h

      iters_cond=0
      iters_lin=0

! call lenstr(grid,lgrd)

! ierr=nf_open(grid(1:lgrd), nf_nowrite, ncid)

      ncx=Lm+2 ; ncy=Mm+2
      allocate( wrk1(ncx,ncy) ,msk(ncx,ncy) )

! call get_var_by_name_double(ncid, 'mask_rho', wrk)

      do j=1,ncy
        do i=1,ncx ! temporary use "wrk" just
          if (wrk(i,j) > 0.5D0) then ! to read mask from the file;
            msk(i,j)=1 ! from now on this code will
          else ! use integer(kind=2) mask
            msk(i,j)=0
          endif
          h(i,j)=max(hmin,min(hmax,hraw(i,j))) !<-- restrict
          wrk(i,j)=h(i,j) !<-- make another copy
        enddo
      enddo
# 116 "smooth.F"
      method='SMOOTH'



C$OMP PARALLEL SHARED(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,wrk1)
      call smooth_thread(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,method,wrk1)
C$OMP END PARALLEL






! call write_topo(ncid, grid, ncx,ncy, h)
      end
# 30 "tools_fort.F" 2

!!! lsmooth !!!

# 1 "smooth.F" 1
! Purpose: to smooth or "log-smooth" or rx-condition model topography
!--------- in order to achieve acceptable r-factors,
!
! | h(i+1) - h(i) |
! rx = ------------------- < rx_max
! h(i+1) + h(i)
!
! It reads raw topography "hraw" from netCDF grid file calls the actual
! smoothing routine. Upon completion it puts the resultant smoothed
! topography into netCDF variable "hsmth" (if exists) or into "h" of
! the same file leaving the original "hraw" unchanged.

! There are multiple versions of such smoothing/rx-limiting procedures
! which can be subdivided into two major groups: (i) the ones which are
! applied directly to "h" (so called rx-"smooth"ers); and (ii) the ones
! with first transform "h" into its logarithm,
!
! Lgh(i,j) = log[h(i,j)/hmin]
!
! and smooth it trying to achieve condition
!
! | Lgh(i+1,j) - Lgh(i+1,j) | < log[(1 + r_max)/(1 - r_max)]
!
! which is equivalent to rx < rx_max condition above, then transform
! it back, h(i,j)=hmin*exp[Lgh(i,j)]. These are called log-smoothers,
! or "lmsooth". Besides this distinctions, there are others:
! handling boundary conditions at coast line; choice of discrete
! operator, etc. - resulting in a fair number of options.

! Usage: "smooth" or "lsmooth" take four arguments:
!
! lsmooth hmin hmax r_max file.nc
! where
! #1 hmin desired minimum depth limit, meters, positive;
! #2 hmax --/ /-- maximum --/ /--
! #3 r_max desired maximum r-factor, nondimensional;
! #4- name of ROMS grid netCDF file;
!
! This file is merely a driver which decodes command-line arguments,
! reads model grid topography, calls relevant subroutine selected by
! CPP-switches (defined from compiler line to support multiple
! versions), and creates appropriate signature in the resultant grid
! file so all the parameters of execution these operators are
! documented to allow exact reproduction of the result relying
! solely on what is stored in netCDF file.


      subroutine lsmooth (Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
# 61 "smooth.F"
      integer Lm,Mm
      character(len=80) :: grid
      character(len=64) method
      character(len=16) str
      real(kind=8) :: hmin, hmax, r_max
      real(kind=8) :: hraw(1:Lm+2,1:Mm+2)
      real(kind=8) :: h(1:Lm+2,1:Mm+2)
      real(kind=8) :: wrk(1:Lm+2,1:Mm+2)
      real(kind=8), dimension(:,:), allocatable :: wrk1

      integer iters_cond, iters_lin
      integer(kind=2), allocatable, dimension(:,:) :: msk
      integer nargs, ncx,ncy, i,j, ierr, ncid, lgrd


      include "netcdf.inc"

cf2py intent(in) hraw,hmin,hmax,r_max,wrk
cf2py intent(out) h

      iters_cond=0
      iters_lin=0

! call lenstr(grid,lgrd)

! ierr=nf_open(grid(1:lgrd), nf_nowrite, ncid)

      ncx=Lm+2 ; ncy=Mm+2
      allocate( wrk1(ncx,ncy) ,msk(ncx,ncy) )

! call get_var_by_name_double(ncid, 'mask_rho', wrk)

      do j=1,ncy
        do i=1,ncx ! temporary use "wrk" just
          if (wrk(i,j) > 0.5D0) then ! to read mask from the file;
            msk(i,j)=1 ! from now on this code will
          else ! use integer(kind=2) mask
            msk(i,j)=0
          endif
          h(i,j)=max(hmin,min(hmax,hraw(i,j))) !<-- restrict
          wrk(i,j)=h(i,j) !<-- make another copy
        enddo
      enddo


      method='LOG_SMOOTHING'
# 120 "smooth.F"
C$OMP PARALLEL SHARED(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,wrk1)
      call smooth_thread(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,method,wrk1)
C$OMP END PARALLEL






! call write_topo(ncid, grid, ncx,ncy, h)
      end
# 34 "tools_fort.F" 2


!!! legacy_lsmooth !!!

# 1 "smooth.F" 1
! Purpose: to smooth or "log-smooth" or rx-condition model topography
!--------- in order to achieve acceptable r-factors,
!
! | h(i+1) - h(i) |
! rx = ------------------- < rx_max
! h(i+1) + h(i)
!
! It reads raw topography "hraw" from netCDF grid file calls the actual
! smoothing routine. Upon completion it puts the resultant smoothed
! topography into netCDF variable "hsmth" (if exists) or into "h" of
! the same file leaving the original "hraw" unchanged.

! There are multiple versions of such smoothing/rx-limiting procedures
! which can be subdivided into two major groups: (i) the ones which are
! applied directly to "h" (so called rx-"smooth"ers); and (ii) the ones
! with first transform "h" into its logarithm,
!
! Lgh(i,j) = log[h(i,j)/hmin]
!
! and smooth it trying to achieve condition
!
! | Lgh(i+1,j) - Lgh(i+1,j) | < log[(1 + r_max)/(1 - r_max)]
!
! which is equivalent to rx < rx_max condition above, then transform
! it back, h(i,j)=hmin*exp[Lgh(i,j)]. These are called log-smoothers,
! or "lmsooth". Besides this distinctions, there are others:
! handling boundary conditions at coast line; choice of discrete
! operator, etc. - resulting in a fair number of options.

! Usage: "smooth" or "lsmooth" take four arguments:
!
! lsmooth hmin hmax r_max file.nc
! where
! #1 hmin desired minimum depth limit, meters, positive;
! #2 hmax --/ /-- maximum --/ /--
! #3 r_max desired maximum r-factor, nondimensional;
! #4- name of ROMS grid netCDF file;
!
! This file is merely a driver which decodes command-line arguments,
! reads model grid topography, calls relevant subroutine selected by
! CPP-switches (defined from compiler line to support multiple
! versions), and creates appropriate signature in the resultant grid
! file so all the parameters of execution these operators are
! documented to allow exact reproduction of the result relying
! solely on what is stored in netCDF file.




      subroutine lsmooth_legacy(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
# 61 "smooth.F"
      integer Lm,Mm
      character(len=80) :: grid
      character(len=64) method
      character(len=16) str
      real(kind=8) :: hmin, hmax, r_max
      real(kind=8) :: hraw(1:Lm+2,1:Mm+2)
      real(kind=8) :: h(1:Lm+2,1:Mm+2)
      real(kind=8) :: wrk(1:Lm+2,1:Mm+2)
      real(kind=8), dimension(:,:), allocatable :: wrk1

      integer iters_cond, iters_lin
      integer(kind=2), allocatable, dimension(:,:) :: msk
      integer nargs, ncx,ncy, i,j, ierr, ncid, lgrd


      include "netcdf.inc"

cf2py intent(in) hraw,hmin,hmax,r_max,wrk
cf2py intent(out) h

      iters_cond=0
      iters_lin=0

! call lenstr(grid,lgrd)

! ierr=nf_open(grid(1:lgrd), nf_nowrite, ncid)

      ncx=Lm+2 ; ncy=Mm+2
      allocate( wrk1(ncx,ncy) ,msk(ncx,ncy) )

! call get_var_by_name_double(ncid, 'mask_rho', wrk)

      do j=1,ncy
        do i=1,ncx ! temporary use "wrk" just
          if (wrk(i,j) > 0.5D0) then ! to read mask from the file;
            msk(i,j)=1 ! from now on this code will
          else ! use integer(kind=2) mask
            msk(i,j)=0
          endif
          h(i,j)=max(hmin,min(hmax,hraw(i,j))) !<-- restrict
          wrk(i,j)=h(i,j) !<-- make another copy
        enddo
      enddo




      method='LEGACY_LOG_SMOOTH'
# 120 "smooth.F"
C$OMP PARALLEL SHARED(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,wrk1)
      call smooth_thread(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,method,wrk1)
C$OMP END PARALLEL






! call write_topo(ncid, grid, ncx,ncy, h)
      end
# 39 "tools_fort.F" 2


!!! lsmooth2 !!!

# 1 "smooth.F" 1
! Purpose: to smooth or "log-smooth" or rx-condition model topography
!--------- in order to achieve acceptable r-factors,
!
! | h(i+1) - h(i) |
! rx = ------------------- < rx_max
! h(i+1) + h(i)
!
! It reads raw topography "hraw" from netCDF grid file calls the actual
! smoothing routine. Upon completion it puts the resultant smoothed
! topography into netCDF variable "hsmth" (if exists) or into "h" of
! the same file leaving the original "hraw" unchanged.

! There are multiple versions of such smoothing/rx-limiting procedures
! which can be subdivided into two major groups: (i) the ones which are
! applied directly to "h" (so called rx-"smooth"ers); and (ii) the ones
! with first transform "h" into its logarithm,
!
! Lgh(i,j) = log[h(i,j)/hmin]
!
! and smooth it trying to achieve condition
!
! | Lgh(i+1,j) - Lgh(i+1,j) | < log[(1 + r_max)/(1 - r_max)]
!
! which is equivalent to rx < rx_max condition above, then transform
! it back, h(i,j)=hmin*exp[Lgh(i,j)]. These are called log-smoothers,
! or "lmsooth". Besides this distinctions, there are others:
! handling boundary conditions at coast line; choice of discrete
! operator, etc. - resulting in a fair number of options.

! Usage: "smooth" or "lsmooth" take four arguments:
!
! lsmooth hmin hmax r_max file.nc
! where
! #1 hmin desired minimum depth limit, meters, positive;
! #2 hmax --/ /-- maximum --/ /--
! #3 r_max desired maximum r-factor, nondimensional;
! #4- name of ROMS grid netCDF file;
!
! This file is merely a driver which decodes command-line arguments,
! reads model grid topography, calls relevant subroutine selected by
! CPP-switches (defined from compiler line to support multiple
! versions), and creates appropriate signature in the resultant grid
! file so all the parameters of execution these operators are
! documented to allow exact reproduction of the result relying
! solely on what is stored in netCDF file.






      subroutine lsmooth2(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)
# 61 "smooth.F"
      integer Lm,Mm
      character(len=80) :: grid
      character(len=64) method
      character(len=16) str
      real(kind=8) :: hmin, hmax, r_max
      real(kind=8) :: hraw(1:Lm+2,1:Mm+2)
      real(kind=8) :: h(1:Lm+2,1:Mm+2)
      real(kind=8) :: wrk(1:Lm+2,1:Mm+2)
      real(kind=8), dimension(:,:), allocatable :: wrk1

      integer iters_cond, iters_lin
      integer(kind=2), allocatable, dimension(:,:) :: msk
      integer nargs, ncx,ncy, i,j, ierr, ncid, lgrd


      include "netcdf.inc"

cf2py intent(in) hraw,hmin,hmax,r_max,wrk
cf2py intent(out) h

      iters_cond=0
      iters_lin=0

! call lenstr(grid,lgrd)

! ierr=nf_open(grid(1:lgrd), nf_nowrite, ncid)

      ncx=Lm+2 ; ncy=Mm+2
      allocate( wrk1(ncx,ncy) ,msk(ncx,ncy) )

! call get_var_by_name_double(ncid, 'mask_rho', wrk)

      do j=1,ncy
        do i=1,ncx ! temporary use "wrk" just
          if (wrk(i,j) > 0.5D0) then ! to read mask from the file;
            msk(i,j)=1 ! from now on this code will
          else ! use integer(kind=2) mask
            msk(i,j)=0
          endif
          h(i,j)=max(hmin,min(hmax,hraw(i,j))) !<-- restrict
          wrk(i,j)=h(i,j) !<-- make another copy
        enddo
      enddo






      method='LOG_SMOOTH_2'
# 120 "smooth.F"
C$OMP PARALLEL SHARED(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,wrk1)
      call smooth_thread(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,method,wrk1)
C$OMP END PARALLEL






! call write_topo(ncid, grid, ncx,ncy, h)
      end
# 44 "tools_fort.F" 2


!!! lsmooth1 !!!

# 1 "smooth.F" 1
! Purpose: to smooth or "log-smooth" or rx-condition model topography
!--------- in order to achieve acceptable r-factors,
!
! | h(i+1) - h(i) |
! rx = ------------------- < rx_max
! h(i+1) + h(i)
!
! It reads raw topography "hraw" from netCDF grid file calls the actual
! smoothing routine. Upon completion it puts the resultant smoothed
! topography into netCDF variable "hsmth" (if exists) or into "h" of
! the same file leaving the original "hraw" unchanged.

! There are multiple versions of such smoothing/rx-limiting procedures
! which can be subdivided into two major groups: (i) the ones which are
! applied directly to "h" (so called rx-"smooth"ers); and (ii) the ones
! with first transform "h" into its logarithm,
!
! Lgh(i,j) = log[h(i,j)/hmin]
!
! and smooth it trying to achieve condition
!
! | Lgh(i+1,j) - Lgh(i+1,j) | < log[(1 + r_max)/(1 - r_max)]
!
! which is equivalent to rx < rx_max condition above, then transform
! it back, h(i,j)=hmin*exp[Lgh(i,j)]. These are called log-smoothers,
! or "lmsooth". Besides this distinctions, there are others:
! handling boundary conditions at coast line; choice of discrete
! operator, etc. - resulting in a fair number of options.

! Usage: "smooth" or "lsmooth" take four arguments:
!
! lsmooth hmin hmax r_max file.nc
! where
! #1 hmin desired minimum depth limit, meters, positive;
! #2 hmax --/ /-- maximum --/ /--
! #3 r_max desired maximum r-factor, nondimensional;
! #4- name of ROMS grid netCDF file;
!
! This file is merely a driver which decodes command-line arguments,
! reads model grid topography, calls relevant subroutine selected by
! CPP-switches (defined from compiler line to support multiple
! versions), and creates appropriate signature in the resultant grid
! file so all the parameters of execution these operators are
! documented to allow exact reproduction of the result relying
! solely on what is stored in netCDF file.
# 54 "smooth.F"
      subroutine lsmooth1(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)






      integer Lm,Mm
      character(len=80) :: grid
      character(len=64) method
      character(len=16) str
      real(kind=8) :: hmin, hmax, r_max
      real(kind=8) :: hraw(1:Lm+2,1:Mm+2)
      real(kind=8) :: h(1:Lm+2,1:Mm+2)
      real(kind=8) :: wrk(1:Lm+2,1:Mm+2)
      real(kind=8), dimension(:,:), allocatable :: wrk1

      integer iters_cond, iters_lin
      integer(kind=2), allocatable, dimension(:,:) :: msk
      integer nargs, ncx,ncy, i,j, ierr, ncid, lgrd


      include "netcdf.inc"

cf2py intent(in) hraw,hmin,hmax,r_max,wrk
cf2py intent(out) h

      iters_cond=0
      iters_lin=0

! call lenstr(grid,lgrd)

! ierr=nf_open(grid(1:lgrd), nf_nowrite, ncid)

      ncx=Lm+2 ; ncy=Mm+2
      allocate( wrk1(ncx,ncy) ,msk(ncx,ncy) )

! call get_var_by_name_double(ncid, 'mask_rho', wrk)

      do j=1,ncy
        do i=1,ncx ! temporary use "wrk" just
          if (wrk(i,j) > 0.5D0) then ! to read mask from the file;
            msk(i,j)=1 ! from now on this code will
          else ! use integer(kind=2) mask
            msk(i,j)=0
          endif
          h(i,j)=max(hmin,min(hmax,hraw(i,j))) !<-- restrict
          wrk(i,j)=h(i,j) !<-- make another copy
        enddo
      enddo
# 112 "smooth.F"
      method='LOG_SMOOTH_1'







C$OMP PARALLEL SHARED(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,wrk1)
      call smooth_thread(Lm,Mm, msk, h, wrk,hmin,hmax,r_max,method,wrk1)
C$OMP END PARALLEL






! call write_topo(ncid, grid, ncx,ncy, h)
      end
# 49 "tools_fort.F" 2


!!! non_decreasing_rx_cond !!!

# 1 "cond_rx0_topo.F" 1
      subroutine cond_rx0_thread(ncx,ncy, msk, h,hmin, hmax, r_max,wrk)

! Conditions topography to restrict the maximum r-factor via
! deepening-only procedure: scan the grid to find points where
! pair-wise ratio
! | h(i+1) - h(i) |
! rx= -------------------
! h(i+1) + h(i)
!
! exceeds threshold "r_max" and adjust THE SMALLER value of
! depth toward THE LARGER; iterate until the condition is no
! longer exceeded anywhere.

      implicit none
      integer ncx,ncy
      integer(kind=2) msk(ncx,ncy)
      real(kind=8) h(ncx,ncy), wrk(ncx,ncy)
      integer :: trd_count=0, npts_all=1
      integer :: iters_cond, iters_lin
      real(kind=8) :: hmin, hmax, r_max
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     & range, istr,iend,jstr,jend, iter, my_npts

      iters_cond=0
      iters_lin=0
C$ integer omp_get_num_threads, omp_get_thread_num
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()

      call set_tiles(ncx,ncy, nsub_x,nsub_y)

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      iter=0
      do while(npts_all > 0)
        iter=iter+1 ; my_npts=0

        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          call cond_rx0_tile(istr,iend,jstr,jend, ncx,ncy,
     & msk, h,wrk, r_max, my_npts)
        enddo
C$OMP BARRIER

        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          call cond_rx0_tile(istr,iend,jstr,jend, ncx,ncy,
     & msk, wrk,h, r_max, my_npts)
        enddo
C$OMP CRITICAL(cond_rx0_smth)
        if (trd_count == 0) npts_all=0
        trd_count=trd_count+1 ; npts_all=npts_all+my_npts
        if (trd_count == numthreads) then
          trd_count=0
          write(*,'(8x,A,I7,2x,A,I10)') 'iter =', iter,
     & 'changes =', npts_all
        endif
C$OMP END CRITICAL(cond_rx0_smth)
C$OMP BARRIER
      enddo !<-- while
      end


      subroutine cond_rx0_tile(istr,iend,jstr,jend, ncx,ncy,
     & msk, src,targ, r_max, my_npts)

      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, my_npts
      integer(kind=2) msk(ncx,ncy)
      real(kind=8) src(ncx,ncy), targ(ncx,ncy), r_max
      integer i,j, iw,js,ie,jn
      real(kind=8) ratio, max_surr
      real(kind=8), parameter :: epsil=1.D-12
      ratio=(1.-r_max)/(1.+r_max)

      do j=jstr,jend
        js=max(j-1,1) ; jn=min(j+1,ncy)
        do i=istr,iend
          iw=max(i-1,1) ; ie=min(i+1,ncx)
          targ(i,j)=src(i,j) !<-- first copy then adjust

          if (msk(i,j) > 0) then
            max_surr=0. !<-- initialize
            if (msk(iw,j) > 0) max_surr=max(max_surr, src(iw,j))
            if (msk(i,js) > 0) max_surr=max(max_surr, src(i,js))
            if (msk(ie,j) > 0) max_surr=max(max_surr, src(ie,j))
            if (msk(i,jn) > 0) max_surr=max(max_surr, src(i,jn))

            max_surr=max_surr*ratio
            if (targ(i,j) < max_surr-epsil) then
              targ(i,j)=max_surr
              my_npts=my_npts+1
            endif
          endif !<-- msk(i,j)
        enddo !<-- i
      enddo !<-- j
      end
# 54 "tools_fort.F" 2
# 1 "smooth.F" 1
! Purpose: to smooth or "log-smooth" or rx-condition model topography
!--------- in order to achieve acceptable r-factors,
!
! | h(i+1) - h(i) |
! rx = ------------------- < rx_max
! h(i+1) + h(i)
!
! It reads raw topography "hraw" from netCDF grid file calls the actual
! smoothing routine. Upon completion it puts the resultant smoothed
! topography into netCDF variable "hsmth" (if exists) or into "h" of
! the same file leaving the original "hraw" unchanged.

! There are multiple versions of such smoothing/rx-limiting procedures
! which can be subdivided into two major groups: (i) the ones which are
! applied directly to "h" (so called rx-"smooth"ers); and (ii) the ones
! with first transform "h" into its logarithm,
!
! Lgh(i,j) = log[h(i,j)/hmin]
!
! and smooth it trying to achieve condition
!
! | Lgh(i+1,j) - Lgh(i+1,j) | < log[(1 + r_max)/(1 - r_max)]
!
! which is equivalent to rx < rx_max condition above, then transform
! it back, h(i,j)=hmin*exp[Lgh(i,j)]. These are called log-smoothers,
! or "lmsooth". Besides this distinctions, there are others:
! handling boundary conditions at coast line; choice of discrete
! operator, etc. - resulting in a fair number of options.

! Usage: "smooth" or "lsmooth" take four arguments:
!
! lsmooth hmin hmax r_max file.nc
! where
! #1 hmin desired minimum depth limit, meters, positive;
! #2 hmax --/ /-- maximum --/ /--
! #3 r_max desired maximum r-factor, nondimensional;
! #4- name of ROMS grid netCDF file;
!
! This file is merely a driver which decodes command-line arguments,
! reads model grid topography, calls relevant subroutine selected by
! CPP-switches (defined from compiler line to support multiple
! versions), and creates appropriate signature in the resultant grid
! file so all the parameters of execution these operators are
! documented to allow exact reproduction of the result relying
! solely on what is stored in netCDF file.
# 56 "smooth.F"
      subroutine cond_rx0_topo(Lm,Mm,hraw,hmin,hmax,r_max,wrk,h)




      integer Lm,Mm
      character(len=80) :: grid
      character(len=64) method
      character(len=16) str
      real(kind=8) :: hmin, hmax, r_max
      real(kind=8) :: hraw(1:Lm+2,1:Mm+2)
      real(kind=8) :: h(1:Lm+2,1:Mm+2)
      real(kind=8) :: wrk(1:Lm+2,1:Mm+2)
      real(kind=8), dimension(:,:), allocatable :: wrk1

      integer iters_cond, iters_lin
      integer(kind=2), allocatable, dimension(:,:) :: msk
      integer nargs, ncx,ncy, i,j, ierr, ncid, lgrd


      include "netcdf.inc"

cf2py intent(in) hraw,hmin,hmax,r_max,wrk
cf2py intent(out) h

      iters_cond=0
      iters_lin=0

! call lenstr(grid,lgrd)

! ierr=nf_open(grid(1:lgrd), nf_nowrite, ncid)

      ncx=Lm+2 ; ncy=Mm+2
      allocate( wrk1(ncx,ncy) ,msk(ncx,ncy) )

! call get_var_by_name_double(ncid, 'mask_rho', wrk)

      do j=1,ncy
        do i=1,ncx ! temporary use "wrk" just
          if (wrk(i,j) > 0.5D0) then ! to read mask from the file;
            msk(i,j)=1 ! from now on this code will
          else ! use integer(kind=2) mask
            msk(i,j)=0
          endif
          h(i,j)=max(hmin,min(hmax,hraw(i,j))) !<-- restrict
          wrk(i,j)=h(i,j) !<-- make another copy
        enddo
      enddo
# 114 "smooth.F"
      method='RX_COND'
# 124 "smooth.F"
C$OMP PARALLEL SHARED(ncx,ncy, msk, h,wrk)
      call cond_rx0_thread(ncx,ncy, msk, h,hmin,hmax,r_max,wrk)
C$OMP END PARALLEL


! call write_topo(ncid, grid, ncx,ncy, h)
      end
# 55 "tools_fort.F" 2


# 1 "r2r_match_topo.F" 1
      subroutine r2r_match_topo (OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH,
     & ic,jc,iip,jjp,wdth,xc,yc,hchd,xi,xp,yp,hp,hwg)


!xp lon_rho parent
!yp lat_rho parent
!hp topo parent
!
!xc lon_rho child
!yc lat_rho child
!hwg child topo
!xi mask_rho child
!wght weight child
!hc new htopo


! Modifies child-grid bottom topography in the vicinity of forced
! open boundaries to exactly match topography interpolated from the
! parent grid right at the boundary, while at the same time, time
! making smooth transition into the interior of child-grid domain
! where the topography is expected to be less smooth because of
! finer resolution.

! The method is essentially generating a 2D weighting function,
! wgt=wgt(i,j), 0 <= wgt <= 1, hence
!
! h_matched = (1-wgt) * h + wgt * h_parent
!
! where wgt=1 at open boundaries, and wgt=0 in the interior beyond
! merging zone, with smooth transition within - sounds quite simple,
! however is should be noted that merging should be avoided in the
! areas closed to the edge of the child domain, but blocked by the
! land mask. This in its turn leads to a rather complicated logical
! "etching" to determine the exact coverage of merging zone.
! Interpolation of parent-grid topography is by bi-cubic spline
! using Hermite polynomials.

! Created/maintained by Alexander Shchepetkin, old_galaxy@yahoo.com

c--#define LEGACY_METHOD

      implicit none
      character(len=4) merging_flags
      character(len=64) prntgrd, chldgrd, str
      integer :: ic,jc,iip,jjp
      real(kind=8), dimension(ic+2,jc+2) :: xc,yc,xi,eta
     & ,wgt,hwg,hchd
      real(kind=8), dimension(iip+2,jjp+2) :: hp,xp,yp

      real(kind=8), allocatable, dimension(:,:) :: hc,srX,srY,sXY,sYX
      real(kind=4), allocatable, dimension(:,:) :: mwgt
      integer(kind=4), allocatable, dimension(:,:) :: ip,jp
      integer(kind=2), allocatable, dimension(:,:) ::mask,mgz,ms1,ms2
      integer nargs, wdth, net_alloc_mem, nx,ny, ncx,ncy, ncsrc,nctarg,
     & ncid, ierr, ipvar,jpvar, xivar,etavar, hpvar,hvar, dhvar,
     & wgtvar,mwgvar, r2dgrd(2), i,j,k, iter, lprnt,lchld
      real(kind=8) xmin,xmax, xcmin,xcmax, cff
      logical OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH, rename_hvar

       integer rad,rad2, ict,jct

      include "netcdf.inc"

Cf2py intent(in) OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH,ic,jc,iip,jjp,wdth,xc,yc,hchd,xi,xp,yp,hp
Cf2py intent(out) hwg

      net_alloc_mem=0 !<-- allocated memory expressed in 4-byte numbers
      merging_flags='    '; rename_hvar=.false.

! Check whether all the arguments are specified correctly and both
! netCDF files can be opened (at least in read-only mode); diagnose
! errors and quite if an errors occurs; write help page and quit if
! the program is called with insufficient number of arguments.

        k=0
        if (OBC_WEST) then ! prepare signature
          k=k+1 ; merging_flags(k:k)='W' ! string to be saved
        endif ! as an attribute in
        if (OBC_EAST) then ! netCDF file
         k=k+1 ; merging_flags(k:k)='E'
        endif
        if (OBC_SOUTH) then
          k=k+1 ; merging_flags(k:k)='S'
        endif
        if (OBC_NORTH) then
          k=k+1 ; merging_flags(k:k)='N'
        endif

! Allocate arrays and read coordinates first for the target grid,
! then for the source. This leads to a more optimal memory use as
! xp,yp arrays for the source grid can be deallocated after computing
! ip,jp-indices and fractional offsets xi,eta, the xp,yp are allocated
! at the end to be deallocated to free memory for arrays associated
! with spline interpolation.

      ncx=ic+2 ; ncy=jc+2

      if (ierr == 0) allocate( mwgt(ncx,ncy), stat=ierr )
      if (ierr == 0) allocate( ip(ncx,ncy),jp(ncx,ncy), stat=ierr )
      if (ierr == 0) allocate( mask(ncx,ncy), mgz(ncx,ncy),
     & ms1(ncx,ncy),ms2(ncx,ncy), stat=ierr )
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem +ncx*ncy*(6*2 + 1*1 + 2*1)
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 1.'
        stop
      endif

      do j=1,ncy
        do i=1,ncx
          if (xi(i,j) > 0.5D0) then ! temporarily use array "xi"
            mask(i,j)=1 ! to read mask from the file;
          else ! thereafter this program uses
            mask(i,j)=0 ! only integer(kind=2) version
          endif ! of mask
        enddo
      enddo

      xcmin=xc(1,1) ; xcmax=xc(1,1)
      do j=1,ncy
        do i=1,ncx
          if (xc(i,j) < xcmin) then
            xcmin=xc(i,j)
          elseif (xc(i,j) > xcmax) then
            xcmax=xc(i,j)
          endif
        enddo
      enddo

! The same for parent grid

      nx=iip+2 ; ny=jjp+2

      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem+(3*2)*nx*ny
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 2.'
      stop
      endif

! Interpolate parent-grid topography onto child grid...

      xmin=xp(1,1) ; xmax=xp(1,1)
      do j=1,ny
        do i=1,nx
          if (xp(i,j) < xmin) then
            xmin=xp(i,j)
          elseif (xp(i,j) > xmax) then
            xmax=xp(i,j)
          endif
        enddo
      enddo
      write(*,*) 'Parent grid xmin,xmax =', xmin,xmax
      if (xmin > xcmax) then
        do j=1,ny
          do i=1,nx
            xp(i,j)=xp(i,j)-360.D0
          enddo
        enddo
        write(*,*) 'Adjusted to',xmin-360.D0,xmax-360.D0
      elseif (xmax < xcmin) then
        do j=1,ny
          do i=1,nx
            xp(i,j)=xp(i,j)+360.D0
          enddo
        enddo
        write(*,*) 'Adjusted to',xmin+360.D0, xmax+360.D0
      endif

C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp, xi,eta)
        call r2r_interp_init_thread(nx,ny, xp,yp, ncx,ncy, xc,yc,
     & ip,jp, xi,eta)
C$OMP END PARALLEL

      call check_search_indices(nx,ny,xp,yp, ncx,ncy, xc,yc, ip,jp)
      call check_offsets(nx,ny,xp,yp, ncx,ncy,xc,yc, ip,jp, xi,eta)

      net_alloc_mem=net_alloc_mem-4*nx*ny

      allocate(srX(nx,ny),srY(nx,ny), sXY(nx,ny),sYX(nx,ny))
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + (4*2)*nx*ny
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 3.'
        stop
      endif

      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)

      deallocate(sYX) ; net_alloc_mem=net_alloc_mem-2*nx*ny
      allocate(hc(ncx,ncy), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + (1*2)*ncx*ncy
        write(*,*) 'allocated', net_alloc_mem/262144, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 4.'
        stop
      endif

C$OMP PARALLEL SHARED(nx,ny, hp, srX,srY,sXY, ncx,ncy,
C$OMP& ip,jp, xi,eta, hc)
      call spln2d_interp_double(nx,ny, hp, srX,srY,sXY,
     & ncx,ncy, ip,jp, xi,eta, hc)
C$OMP END PARALLEL


! Start forming merging zone and weights: the two methods below differ
!------ --------------- -------- -------- only by handling land mask
! and produce identical results in the case of no land. Either way,
! the outcome is "mgz" assuming values 0 or 1 to define merging zone
! (unmasked water points only) and "wgt" initialized as constant-slope
! function wgt=wdth at each open boundary water points from which it
! descends linearly into the interior of the domain by decreasing by
! 1 for each row of grid points until it vanishes to 0. In the case
! of non-LEGACY_METHOD "wgt" also decreases by 1 every step when it
! goes around the corner of the coastline -- the minimal path from
! the costline to the given point is counted as the number of steps
! connecting unmasked water points by C-grid masking rules rather than
! along a straight line to the nearest boundary point, so the resultant
! weights and the extent of merging zone expected to be somewhat
! smaller in areas behind land.
# 257 "r2r_match_topo.F"
      do j=1,ncy ! Initialize "mgz" as a single
        do i=1,ncx ! row of perimeter points along
          mgz(i,j)=0 ! the unmasked parts of open
        enddo ! boundaries, then "etch" inward,
      enddo ! while obeying the connectivity
      if (OBC_WEST) then ! rules of land-water masking:
        do j=1,ncy ! the by the construction "mgz"
          mgz(1,j)=mask(1,j) ! has the property that each of
        enddo ! its points can be reached from
      endif ! the open boundary by water.
      if (OBC_SOUTH) then
        do i=1,ncx
          mgz(i,1)=mask(i,1)
        enddo
      endif
      if (OBC_EAST) then
        do j=1,ncy
          mgz(ncx,j)=mask(ncx,j)
        enddo
      endif
      if (OBC_NORTH) then
        do i=1,ncx
          mgz(i,ncy)=mask(i,ncy)
        enddo
      endif

C$OMP PARALLEL SHARED(ncx,ncy, mask,ms1,ms2, wdth)
      call etch_mgz_weights_thread(ncx,ncy, mask,mgz,ms2, wdth)
C$OMP END PARALLEL

      do j=1,ncy
        do i=1,ncx
          if (mgz(i,j) > 0) then
            wgt(i,j)=dble(mgz(i,j)) ; mgz(i,j)=1
          else
            wgt(i,j)=0.D0
          endif
        enddo
      enddo


! The following part makes round junctions between two adjacent open
! boundaries by extending merging zone inward. This is useful to avoid
! steep gradient near inner corner of merging zone when weighting
! finction "wgt" is subjected to Laplacian smoothing. Note that all
! loops below are non-reversible as they process points in row-by-row
! manner (both in i and j directions) with checking that the previous
! row is set to 1.


      rad=(wdth+1)/3 ; rad2=rad**2
      if (OBC_NORTH .and. OBC_WEST) then
        ict=wdth+rad+1 ; jct=ncy-wdth-rad !<-- center of the circle
        do j=jct+rad,jct,-1
          do i=ict-rad,ict,+1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     & .and. mgz(i-1,j) > 0 .and. mgz(i,j+1) > 0
     & .and. mgz(i-1,j+1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_WEST .and. OBC_SOUTH) then
        ict=wdth+rad+1 ; jct=wdth+rad+1
        do j=jct-rad,jct,+1
          do i=ict-rad,ict,+1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     & .and. mgz(i-1,j) > 0 .and. mgz(i,j-1) > 0
     & .and. mgz(i-1,j-1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_SOUTH .and. OBC_EAST) then
        ict=ncx-wdth-rad ; jct=wdth+rad+1 !<-- center of the circle
        do j=jct-rad,jct,+1
          do i=ict+rad,ict,-1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     & .and. mgz(i+1,j) > 0 .and. mgz(i,j-1) > 0
     & .and. mgz(i+1,j-1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif
      if (OBC_EAST .and. OBC_NORTH) then
        ict=ncx-wdth-rad ; jct=ncy-wdth-rad
        do j=jct+rad,jct,-1
          do i=ict+rad,ict,-1
            if ( (i-ict)**2+(j-jct)**2 > rad2 .and. mask(i,j) > 0
     & .and. mgz(i+1,j) > 0 .and. mgz(i,j+1) > 0
     & .and. mgz(i+1,j+1) > 0 ) mgz(i,j)=1
          enddo
        enddo
      endif


! Etch the merging zone into the land mask, in such a way that
! any possibility of reaching water behind mask is excluded.

C$OMP PARALLEL SHARED(ncx,ncy, mask, mgz,ms1,ms2)
      call etch_mgz_into_land_thread(ncx,ncy, mask, mgz,ms1,ms2)
C$OMP END PARALLEL

! At this stage "ms1 > 0" defines the area to which the merging zone
! should be allowed to grow. This area consists of all water points
! of the merging zone which can be reached by water starting from
! the boundary row of non-masked points, and the adjacent land area
! to which to merging zone can be expanded to ensure smoothness of
! the merged topography, but without interfering with water interior.
! The next stage is to etch the actual merging function into the
! land-masked area.

      do j=1,ncy ! Set "mgz" to be merging
        do i=1,ncx ! zone mask (positive "ms1"
          if (ms1(i,j) > 0) then ! area), which will remain
            if (mask(i,j) > 0) then ! unchanged from now on,
              ms1(i,j)=int(wgt(i,j)+0.5D0) ! then initialize both
            else ! "ms1" and "ms2" to the
              ms1(i,j)=0 ! constant-slope merging
            endif ! function (which later
            mgz(i,j)=1 ! be used as the argument
          else ! for the actual merging
            ms1(i,j)=0 ! function), using integer
            mgz(i,j)=0 ! numbers from the range
          endif ! from 0 to "wdth".
          ms2(i,j)=ms1(i,j)
        enddo
      enddo

C$OMP PARALLEL SHARED(wdth, ncx,ncy, mgz,ms1,ms2, OBC_WEST,OBC_EAST,
C$OMP& OBC_SOUTH,OBC_NORTH)
      call etch_weights_into_land_thread(wdth, ncx,ncy, mgz,ms1,ms2,
     & OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
C$OMP END PARALLEL

! At this moment if there would be no land "ms1" is constant-slope
! function: the same maximum value at the outermost row of points and
! descending linearly to zero within the band "wdth" into the domain -
! so smoothing iterations below (subject to Dirichlet B.C) would not
! actually cause any change to it with the exception of (i) in the
! vicinity of corners between adjacent open boundaries (resulting in
! sharp bend along the 45-degree row points where two side slopes
! come together), and (ii) irregular shapes caused by "etching" of
! "ms1" into the land.

      do j=1,ncy
        do i=1,ncx
          wgt(i,j)=dble(ms1(i,j)) ; hwg(i,j)=wgt(i,j)
        enddo
      enddo
      iter=4*wdth
C$OMP PARALLEL SHARED(iter, ncx,ncy, mgz,wgt,hwg)
      call smooth_wgt_thread(iter, ncx,ncy, mgz,wgt,hwg)
C$OMP END PARALLEL

! Convert the originally constant-slope --> etched --> smoothed "wgt"
! into merging function, then perform the actual merging of topography.
! Note that the incoming "wgt" is within the range of 0 <= wgt <= wdth
! nominally changing by 1 from one point to the next (nominally because
! masking-etching-smoothing may modify this). A slightly different
! rescaling,
! -1/(2*(wdth-1)) <= wgt <= (wdth-1/2)(wdth-1)
!
! subject to limiting 0 <= wgt <= 1 is used below with the rationale
! to take out discretely jagged inner border of merging zone (obtained
! by etching), and to make cos() function a bit closer to 1 at the
! second row of points near the open boundary.
      hwg=hchd

      write(*,*) 'Merging topography...'
      cff=1.D0/dble(wdth-1)
      do j=1,ncy
        do i=1,ncx
          if (mgz(i,j) > 0) then
            wgt(i,j)=min(1.D0, cff*(wgt(i,j)-0.5D0))
            if (wgt(i,j) < 0.D0) then
              wgt(i,j)=0.D0 ; mgz(i,j)=0
            endif
            wgt(i,j)=0.5D0-0.5D0*cos(3.14159265358979323D0*wgt(i,j))
            hwg(i,j)=hwg(i,j)*(1.D0-wgt(i,j))+hc(i,j)*wgt(i,j)
          endif

          if (mgz(i,j) > 0) then
            if (mask(i,j) > 0) then
              mwgt(i,j)=1.
            else ! Array "mwgt" is merely
              mwgt(i,j)=0.8 ! for 4-color illustration
            endif ! of the layout of merging
          else ! zone relatively to land
            wgt(i,j)=0.D0 ! mask:
            if (mask(i,j) > 0) then !
              mwgt(i,j)=0. ! mwgt=1 or 0.8 merge
            else ! mwgt=0 or 0.2 intact
              mwgt(i,j)=0.2 !
            endif ! mwgt=0 or 1 water
          endif ! mwgt= 0.8 or 0.2 land
        enddo
      enddo
      end subroutine
# 58 "tools_fort.F" 2


# 1 "r2r_init.F" 1
      subroutine r2r_init(ntrc,chldgrd,theta_s,theta_b,hc,N,prntgrd,
     & prnt_data,rec,tracer)

! Creates initial condition file for child-grid model by interpolating
! user-specified record of parent-grid solution. Should be used as
!
! r2r_init chld_grd.nc theta_s theta_b hc N prnt_grd.nc prnt_data.nc rec
!
! where theta_s theta_b hc N are S-coordinate parameters and number of
! vertical levels for the intended child-grid model; chld_grd.nc and
! prnt_grd.nc are netCDF grid files for child and parent respectively;
! prnt_data.nc is file name of parent-grid history file containing data
! for zeta,u,v,T,S (mandatory variables for the initial conditions),
! "rec" is record number within this file (only one record will be
! read); upon completion the interpolated fields are written into
! netCDF file "roms_init.nc" (always named this way, specified by
! parameter roms_init below) containing just a single record and an
! auxiliary file "r2r_init_diag.nc" for diagnostic purposes for this
! program itself -- there is no use for it other than to ncview it
! and verify sanity of the algorithms below.

! Method: 2D horizontal bi-cubic spline interpolation via Hermite basis
! functions followed by vertical splines by spline inversion (z-levels
! of child grid are translated into continuous k-index space of parent
! grid in such a way spline interpolation of the resultant kp=kp(z_r)
! values back into z-space of child grid yields the original child-grid
! z_r=z_r(k) exactly. in doing so vertical spline derivatives of both
! z_r and the field to be interpolated are constructed in k-index
! space). Prior to the horizontal interpolation data is extended to
! land-masked areas by etching algorithm.

! CPP-switch PARENT_GRID_SUBREGION make this program read only relevant
! subdomain within the parent grid as opposite to reading the whole
! data. The pros and cons for this are as follows: (1) all parent-grid
! arrays are allocated with smaller size; (2) computational savings
! mainly due to less data etching to land-masked areas (other savings
! as well, but not so dramatic); (3) as for the reading netCDF files
! it may or may not improve depending on the specific situation:
! reading sub-array from netCDF-3 file is often slower than reading
! the whole thing because it effectively reading many small records
! instead of reading in one single touch; reading compressed netCDF-4
! files may be faster because fewer blocks needs to be read (highly
! dependent on block structure of the source file).

! CPP-switch CONTOUR_CHILD_MASK reduces masking in the resultant
! file "roms_init.nc" to just one row of points along the coastline,
! while the internal pints inside land remain unmasked. This is useful
! to see what the result of etching procedure, but overall selecting
! this switch does not affect the usability of the file as the initial
! condition because ROMS code will mask all the land points on its own
! any way.

! CPP-switch WITH_REC_DIM makes all the variables the resultant file
! have time dimension which is unlimited dimension, but still having
! only one record written. Both versions can be used as the initial
! to start ROMS model.

! Parameter setting "btm_slp" and "btm_trc" control type of boundary
! condition at bottom for the purpose of computing spline derivatives:
! btm_slp=-1 is no-slip, meaning that u,v=0 exactly at the bottom, i.e.
! half-grid-interval below the k=1 value; This applies to velocities
! only; for tracers the meaningful b.c. are natural on Neumann (0 and
! +1 respectively). Top b.c. is always natural.

! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com




c-#define WITH_REC_DIM

c--#define VERBOSE

      implicit none
      character(len=16) :: VertCoordType, str
      character(len=160) :: chldgrd, prntgrd, prnt_data
      character(len=16) :: roms_init
      character(len=10) :: time_var_name
      integer :: ntrc,ipt_char_len
      character(len=20), dimension(ntrc) :: tracer
      real(kind=8) :: theta_s,theta_b,hc,hcp, xcmin,xcmax, cff,cff1,cff2
      real(kind=8), allocatable, dimension(:) :: Cs_w,Cs_r,Csp_w,Csp_r,
     & srX,srY,sXY,sYX, kprnt
      real(kind=8), allocatable, dimension(:,:) :: h,xc,yc, hp,xp,yp,
     & xi,eta, xpu,ypu,xiu,etau, xpv,ypv,xiv,etav, hprnt, csA,snA
      integer(kind=4), allocatable, dimension(:,:) :: ip,jp, ipu,jpu,
     & ipv,jpv
      integer(kind=2), allocatable, dimension(:,:) :: mskp, umsp, vmsp,
     & mask, umask,vmask
      real(kind=4), allocatable, dimension(:) :: wrk1,wrk2,wrk3,wrmx

      integer :: net_alloc_mem, nargs, rec, ncpgrd,nccgrd,ncsrc,nctarg,
     & nx,ny,Np, ncx,ncy,N, i,j,k, itrc,isrc,itrg, size, varid,
     & tvar_in,tvar_out, natts, ierr, lstr, lpgrd,lprnt,lcgrd
      integer :: btm_slp ! no-slip=-1; +1 free; 0 natural
     & , btm_trc ! +1 Neumann; 0 natural;

      integer(kind=8) :: read_clk, sz_read_acc,
     & write_clk, sz_write_acc
      real(kind=8) :: hmin

      integer :: iwestpg, jsouthpg
      integer :: imin,imax,jmin,jmax
      integer :: margn
# 113 "r2r_init.F"
      include "phys_const.h"
      include "netcdf.inc"
Cf2py intent(in) ntrc,chldgrd,theta_s,theta_b,hc,N,prntgrd,prnt_data,rec,tracer
# 126 "r2r_init.F"
      roms_init='croco_chd_ini.nc';time_var_name='scrum_time'
      margn=8
      ipt_char_len=20 ! nb of tracers char
      btm_slp=-1 ; btm_trc=+1

      call lenstr(chldgrd,lcgrd)
      ierr=nf_open(chldgrd(1:lcgrd), nf_nowrite, nccgrd)
      if (ierr /= nf_noerr) then
          write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #1 :: Cannot ',
     & 'open ''', chldgrd(1:lcgrd), '''.', nf_strerror(ierr)
          stop
      endif

      call lenstr(prntgrd,lpgrd)
      ierr=nf_open(prntgrd(1:lpgrd), nf_nowrite, ncpgrd)
      if (ierr /= nf_noerr) then
          write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #6 ::',
     & ' Cannot open ''', prntgrd(1:lpgrd), '''.',
     & nf_strerror(ierr)
          stop
      endif

      call lenstr(prnt_data,lprnt)
      ierr=nf_open(prnt_data(1:lprnt), nf_nowrite, ncsrc)
      if (ierr /= nf_noerr) then
          write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #7 ',
     & ':: Cannot open ''', prnt_data(1:lprnt),
     & '''.', nf_strerror(ierr)
          stop
      endif

      call roms_find_dims(nccgrd, chldgrd, i,j,k)
      ncx=i+2 ; ncy=j+2
      allocate(Cs_r(N), Cs_w(0:N))
      call set_scoord(theta_s,theta_b, N, Cs_r,Cs_w)
      call roms_find_dims(ncsrc, prnt_data, i,j,k)
      call roms_check_dims(ncpgrd, prntgrd, i,j,0)
      nx=i+2 ; ny=j+2 ; Np=k
      allocate(Csp_w(0:Np),Csp_r(Np))
      call read_scoord(ncsrc, Np, Csp_r,Csp_w,hcp, VertCoordType)
!_________________________

! Allocate arrays and read coordinates first for the child grid, then
! for the parent. This leads to a more optimal memory use because
! coordinate arrays xc,yc, xp,yp, xpu,ypu, xpv,ypv for both child and
! parent grid can be deallocated after computing ip,jp-indices and
! fractional offsets xi,eta. Keep track of the total amount of
! allocated memory "net_alloc_mem" expressed in 4-byte numbers.

      allocate( csA(ncx,ncy),snA(ncx,ncy), h(ncx,ncy), hprnt(ncx,ncy),
     & mask(ncx,ncy), umask(ncx-1,ncy), vmask(ncx,ncy-1),
     & stat=ierr )
      if (ierr == 0) then
        net_alloc_mem=2*(2*N+1) +2*(2*Np+1) +4*2*ncx*ncy +3*ncx*ncy/2
        allocate( xi(ncx,ncy),eta(ncx,ncy), xiu(ncx,ncy),etau(ncx,ncy),
     & xiv(ncx,ncy),etav(ncx,ncy),
     & stat=ierr )
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem +6*2*ncx*ncy
          allocate( ip(ncx,ncy),jp(ncx,ncy), ipu(ncx,ncy),jpu(ncx,ncy),
     & ipv(ncx,ncy),jpv(ncx,ncy),
     & stat=ierr )
          if (ierr == 0) then
            net_alloc_mem = net_alloc_mem + 6*ncx*ncy
            allocate( xc(ncx,ncy),yc(ncx,ncy), stat=ierr)
            if (ierr == 0) then
              net_alloc_mem = net_alloc_mem + 2*2*ncx*ncy
            endif
          endif
        endif
      endif
      if (ierr == 0) then
        write(*,'(10x,A,F9.2,1x,A)') 'allocated', dble(net_alloc_mem)
     & /262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 1.'
        stop
      endif

! Read child grid land mask and topography. Because land mask may or
! may not present check for it first using netCDF native functions; if
! there is no variable "mask_rho" in netCDF file, set it to water
! everywhere. Temporarily use array "xc" to read mask from the file,
! thereafter this program uses only integer(kind=2) version of mask.

      ierr=nf_inq_varid(nccgrd, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        ierr=nf_get_var_double(nccgrd, varid, xc) !<-- termorarily
        if (ierr == nf_noerr) then



          call set_mask(ncx,ncy, xc, mask,umask,vmask)

        else
          write(*,'(/1x,5A/)') '### ERROR: Cannot read ''mask_rho'' ',
     & 'from ''', chldgrd(1:lcgrd), ''', ', nf_strerror(ierr)
          stop
        endif
      else
        mask=1 ; umask=1 ; vmask=1
        write(*,'(9x,4A)') 'No land mask ''mask_rho'' is present ',
     & 'in ''', chldgrd(1:lcgrd), ''', assuming mask=1 everywhere.'
      endif
      call get_var_by_name_double(nccgrd, 'h', h)

! Read horizontal coordinates for the child grid...

      call get_var_by_name_double(nccgrd, 'lon_rho', xc)
      call get_var_by_name_double(nccgrd, 'lat_rho', yc)
      write(*,'(1x,A)',advance='no') 'child grid longitude '
      call compute_min_max(ncx,ncy, xc, xcmin,xcmax)



! Preliminary step: read horizontal coordinates for the parent
!------------ ----- grid, then pretend initializing parent --> child
! interpolation, but actually all what we need is (ip,jp)-indices:
! use them to find the smallest logically rectangular subdomain within
! the parent grid which encloses the child (the unmasked portion of it,
! to be more specific). This region is characterized by bounds
! imin,imax,jmin,jmax defined within the parent grid index space.
! Expand the region in all four directions by making "margn"-wide halo
! to allow more points to reduce the influence of artificial boundary
! conditions associated with spline interpolations, then save
! coordinates of south-west corner of the expanded region and redefine
! nx,ny consistently with its size. Thereafter deallocate coordinate
! arrays for both parent and child -- basically all what this step is
! needed for is just 4 integer numbers: illegal,illegal,nx,ny.

      allocate(xp(nx,ny), yp(nx,ny), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 2*2*nx*ny
      endif

      call get_var_by_name_double(ncpgrd, 'lon_rho', xp)
      call get_var_by_name_double(ncpgrd, 'lat_rho', yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)

C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp)
      call r2r_interp_search_thread( nx,ny, xp,yp, ncx,ncy, xc,yc,
     & ip,jp)
C$OMP END PARALLEL
      call check_search_indices(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp)

      call compute_index_bounds( ncx,ncy, ip,jp, mask, imin,imax,
     & jmin,jmax)

      write(*,'(1x,2A/2(4x,A,2I5,1x,A,I5))') 'minimal parent-grid ',
     & 'index bounds to accommodate child grid:',
     & 'imin,imax =', imin,imax, 'of nx =', nx,
     & 'jmin,jmax =', jmin,jmax, 'of ny =', ny
      imin=max(imin-margn,1) ; jmin=max(jmin-margn,1)
      imax=min(imax+margn,nx) ; jmax=min(jmax+margn,ny)
      write(*,'(1x,A/4x,A,2I5,17x,A,2I5)') 'adjusted to',
     & 'imin,imax =',imin,imax, 'jmin,jmax =',jmin,jmax

      iwestpg=imin ; nx=imax-imin+1 ; jsouthpg=jmin ; ny=jmax-jmin+1
      write(*,'(/2(2x,A,I5)/)') 'setting subdomain sizes  nx =', nx,
     & 'ny =', ny
      deallocate(xp,yp)
      net_alloc_mem=net_alloc_mem -2*2*nx*ny
      write(*,'(8x,A,F9.2,1x,A)') 'deallocated',
     & dble(2*2*nx*ny)/262144.D0,'MB'


! Note that the only outcome of the code segment above is 4 integer
! numbers, illegal,illegal, nx,ny, which are the indices of south-west
! corner of the subdomain within the parent grid, and the sizes of the
! subdomain. Everything else is discarded. Allocate subdomain-sized
! (within the parent grid) 2D arrays: note that xp(nx,ny),yp(nx,ny)
! were deallocated above and now are allocated again with a different
! (expected to be smaller) size).

      allocate( srX(nx*ny), srY(nx*ny), sXY(nx*ny), sYX(nx*ny),
     & stat=ierr )
      if (ierr == 0) then
        allocate(hp(nx,ny), mskp(nx,ny), umsp(nx-1,ny), vmsp(nx,ny-1),
     & xp(nx,ny), yp(nx,ny), xpu(nx-1,ny), xpv(nx,ny-1),
     & ypu(nx-1,ny), ypv(nx,ny-1),
     & stat=ierr )
      endif
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem +(4+7)*2*nx*ny +3*(nx*ny)/2
        write(*,'(10x,A,F9.2,1x,A)') 'allocated',dble(net_alloc_mem)
     & /262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 2.'
        stop
      endif

! Re-read (read first time if CPP is undefined)
! horizontal coordinates for parent grid. however this time only within
! the subdomain of parent grid. Add/subtract 360 degrees to/from
! longitude if necessary to be consistent with child-grid xcmin,xcmax
! range determined above.

      call get_patch_by_name_double(ncpgrd, prntgrd, 'lon_rho',
     & iwestpg,jsouthpg, nx,ny,0,0, xp)
      call get_patch_by_name_double(ncpgrd, prntgrd, 'lat_rho',
     & iwestpg,jsouthpg, nx,ny,0,0, yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)

! Initialize horizontal interpolation:
# 341 "r2r_init.F"
      write(*,'(/1x,A)')

     & 'initializing horizontal interpolation...'
C$OMP PARALLEL SHARED(nx,ny, xp,yp, ncx,ncy, xc,yc, ip,jp, xi,eta)
      call r2r_interp_init_thread(nx,ny, xp,yp, ncx,ncy, xc,yc,
     & ip,jp, xi,eta)
C$OMP END PARALLEL

      call check_search_indices(nx,ny,xp,yp, ncx,ncy, xc,yc, ip,jp)
      call check_offsets(nx,ny,xp,yp, ncx,ncy,xc,yc, ip,jp, xi,eta)

      do j=1,ny
        do i=1,nx-1
          xpu(i,j)=0.5D0*(xp(i,j)+xp(i+1,j))
          ypu(i,j)=0.5D0*(yp(i,j)+yp(i+1,j))
        enddo
      enddo
C$OMP PARALLEL SHARED(nx,ny, xpu,ypu, ncx,ncy,xc,yc, ipu,jpu,xiu,etau)
      call r2r_interp_init_thread(nx-1,ny, xpu,ypu, ncx,ncy, xc,yc,
     & ipu,jpu, xiu,etau)
C$OMP END PARALLEL
      call check_search_indices(nx-1,ny,xpu,ypu, ncx,ncy,xc,yc,ipu,jpu)
      call check_offsets(nx-1,ny,xpu,ypu, ncx,ncy,xc,yc,
     & ipu,jpu,xiu,etau)

      do j=1,ny-1
        do i=1,nx
          xpv(i,j)=0.5D0*(xp(i,j)+xp(i,j+1))
          ypv(i,j)=0.5D0*(yp(i,j)+yp(i,j+1))
        enddo
      enddo
C$OMP PARALLEL SHARED(nx,ny, xpv,ypv, ncx,ncy,xc,yc, ipv,jpv,xiv,etav)
      call r2r_interp_init_thread(nx,ny-1, xpv,ypv, ncx,ncy, xc,yc,
     & ipv,jpv, xiv,etav)
C$OMP END PARALLEL
      call check_search_indices(nx,ny-1,xpv,ypv, ncx,ncy,xc,yc,ipv,jpv)
      call check_offsets(nx,ny-1,xpv,ypv, ncx,ncy,xc,yc,
     & ipv,jpv,xiv,etav)
# 412 "r2r_init.F"
! Read parent-grid land mask and topography. Similarly to above,
! temporarily use use array "hp" to read "mask_rho" from the file,
! thereafter this program uses only integer(kind=2) version of mask.

      ierr=nf_inq_varid(ncpgrd, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        call get_patch_by_name_double(ncpgrd, prntgrd, 'mask_rho',
     & iwestpg,jsouthpg, nx,ny,0,0, hp)
        call set_mask(nx,ny, hp, mskp,umsp,vmsp)
      else
        mskp=1 ; umsp=1 ; vmsp=1
        write(*,'(9x,4A)') 'No land mask ''mask_rho'' is present ',
     & 'in ''', prntgrd(1:lpgrd), ''', assuming mask=1 everywhere.'
      endif
      call get_patch_by_name_double(ncpgrd, prntgrd, 'h',
     & iwestpg,jsouthpg, nx,ny,0,0, hp)

      if ( minval(h) <=0 ) then ! To prevent errors from
          hmin=minval(h)
          hp=hp-minval(h)+0.2 ! wet and drying, move
          h=h-minval(h)+0.2 ! topography so it is positive
      endif ! everywhere

! Interpolate parent grid topography onto child grid...

      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, hp, srX,srY,sXY, ncx,ncy,
C$OMP& ip,jp, xi,eta, hprnt)
      call spln2d_interp_double(nx,ny, hp, srX,srY,sXY,
     & ncx,ncy, ip,jp, xi,eta, hprnt)
C$OMP END PARALLEL
                                                       ! to prevent
      do j=1,ncy ! zeros left by
        do i=1,ncx ! horizontal
          if (hprnt(i,j) < 0.0001D0) hprnt(i,j)=h(i,j) ! interpolation
        enddo ! where parent
      enddo ! and child do
                                                       ! not overlap
! Read angle between true East and local direction of XI-coordinate
! of ROMS grid, then compute csA=cos(alpha) and snA=sin(alpha), first
! for child grid then for the parent. For the latter temporarily place
! the outcome into arrays xp,yp -- after all horizontal interpolations
! have been initialized the content arrays is no longer needed. Then
! interpolate cos(A) and sin(A) of the parent into child grid, then
! compute cos and sin of the child-parent difference of angles; these
! will be used to rotate velocity vector components.

      call read_angle(nccgrd, chldgrd, 1,1, ncx,ncy, csA,snA)
      call read_angle(ncpgrd, prntgrd, iwestpg,jsouthpg, nx,ny, xp,yp)

      call spln2d_double(nx,ny, xp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, xp, srX,srY,sXY, ncx,ncy,
C$OMP& ip,jp, xi,eta, xc)
      call spln2d_interp_double(nx,ny, xp, srX,srY,sXY,
     & ncx,ncy, ip,jp, xi,eta, xc)
C$OMP END PARALLEL

      call spln2d_double(nx,ny, yp, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, yp, srX,srY,sXY, ncx,ncy,
C$OMP& ip,jp, xi,eta, yc)
      call spln2d_interp_double(nx,ny, yp, srX,srY,sXY,
     & ncx,ncy, ip,jp, xi,eta, yc)
C$OMP END PARALLEL

! Convert child-grid angles cosA,sinA into cosA=cos(child-parent) and
! snA=sin(child-parent) where the parent-grid angle data is coming from
! above as arrays xc=cos(parent) and yc=sin(parent) interpolated into
! child grid. Because the intepolation does not preserve the property
! of having cos^2+sin^2=1 exactly, re-normalize them in the process.

      do j=1,ncy
        do i=1,ncx
          cff=1.D0/sqrt(xc(i,j)*xc(i,j)+yc(i,j)*yc(i,j))
          cff1=csA(i,j)*xc(i,j) +snA(i,j)*yc(i,j)
          cff2=snA(i,j)*xc(i,j) -csA(i,j)*yc(i,j)
          csA(i,j)=cff*cff1 ; snA(i,j)=cff*cff2
        enddo
      enddo

      deallocate(hp, xp,yp, xpu,ypu, xpv,ypv, xc,yc)
      net_alloc_mem=net_alloc_mem-(7*2*nx*ny+2*2*ncy*ncy)
      write(*,'(8x,A,F10.2,1x,A)') 'deallocated',
     & dble(7*2*nx*ny+2*2*ncy*ncy)/262144.D0,'MB'
      ierr=nf_close(nccgrd) ; ierr=nf_close(ncpgrd)


! Initialize vertical interpolation

      size=ncx*ncy*N ; allocate(kprnt(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 2*size
        write(*,'(/8x,A,F9.2,1x,A/)') 'allocated kprnt, reaching',
     & dble(net_alloc_mem)/262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 3.'
        stop
      endif
# 518 "r2r_init.F"
      write(*,'(/1x,A)')

     & 'initializing vertical interpolation...'


C$OMP PARALLEL SHARED(ncx,ncy, hprnt,Np,hcp,Csp_r, h,N,hc,Cs_r, kprnt)
      call r2r_init_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     & h, N, hc, Cs_r, kprnt)

      call r2r_check_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     & h, N, hc, Cs_r, kprnt)
C$OMP END PARALLEL
# 542 "r2r_init.F"
      call r2r_init_diag_file(ncx,ncy, N, ip,jp,xi,eta,
     & ipu,jpu,xiu,etau, ipv,jpv,xiv,etav,
     & csA,snA, h, hprnt, kprnt)

      write(*,'(1x,A/)') 'initialization complete'

c** stop !<-- to test initialization


! Allocate large 3D arrays: wrk1, wrk2, wrk3 must be of sufficient
! size to hold one time record of the largest-possible 3D field which
! may be either parent- or child-grid-size variable. Furthermore,
! horizontal interpolation creates an intermediate field with mixed
! dimensions: horizontally on child grid (ncx,ncy), but still of parent
! size vertically and therefore its total size may be bigger than
! either parent or child [say, if parent has fewer horizontal points
! but more vertical levels than the child]. So the forth array, wrmx,
! is allocated to a size sufficient for all three cases, parent, child,
! and intermediate.

      size=max(nx*ny*(Np+1), ncx*ncy*(N+1))
      allocate(wrk1(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem+size
        write(*,'(8x,A,F10.2,1x,A)') 'allocated wrk1, reaching',
     & dble(net_alloc_mem)/262144.D0, 'MB total'
        allocate(wrk2(size), stat=ierr)
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem+size
          write(*,'(8x,A,F10.2,1x,A)') 'allocated wrk2, reaching',
     & dble(net_alloc_mem)/262144.D0, 'MB total'
          allocate(wrk3(size), stat=ierr)
          if (ierr == 0) then
            net_alloc_mem=net_alloc_mem+size
            write(*,'(8x,A,F10.2,1x,A)') 'allocated wrk3, reaching',
     & dble(net_alloc_mem)/262144.D0, 'MB total'

            size=max(size, ncx*ncy*(Np+1)) !<-- mixed size array
            allocate(wrmx(size), stat=ierr)
            if (ierr == 0) then
              net_alloc_mem=net_alloc_mem+size
              write(*,'(8x,A,F10.2,1x,A/)')'allocated wrmx, reaching',
     & dble(net_alloc_mem)/262144.D0, 'MB total'
            endif
          endif
        endif
      endif
      if (ierr /= 0) then
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 4.'
        stop
      endif


! ***** ********* ****** ******* *********
! *** *** * *** * ** *** *** *** * *** *
! *** *** ** *** *** *** ***
! ***** *** *** *** *** ** ***
! *** *** ********* ****** ***
! *** *** *** *** *** *** ** ***
! ***** *** *** *** *** *** ***

! Create target file, save names of the parent grid and data source
! files, record number as global attributes so this target file can
! be reproduced later using exactly the same conditions...

      ierr=nf_create(roms_init, nf_netcdf4, nctarg)
      if (ierr == nf_noerr) then
        call def_roms_file(ntrc,nctarg, roms_init, ncx,ncy,N,
     & theta_s,theta_b, hc, Cs_w,Cs_r,tracer,ncsrc)
        ierr=nf_put_att_text(nctarg, nf_global, 'memo', 19,
     & 'created by r2r_init')
        if (ierr == nf_noerr) then
          ierr=nf_put_att_text(nctarg, nf_global, 'parent_grid',
     & lpgrd, prntgrd(1:lpgrd))
          if (ierr == nf_noerr) then
            ierr=nf_put_att_text(nctarg, nf_global,
     & 'parent_data_file', lprnt, prnt_data(1:lprnt))
            if (ierr == nf_noerr) then
              ierr=nf_put_att_int(nctarg, nf_global,
     & 'parent_data_file_record', nf_int, 1, rec)
              if (ierr == nf_noerr) then
                write(*,*) 'added r2r_init attributes'
              else
                write(*,*) '### ERROR 5: put_att_int'
              endif
            else
              write(*,*) '### ERROR 4: put_att_text'
            endif
          else
            write(*,*) '### ERROR 3: put_att_text'
          endif
        else
          write(*,*) '### ERROR 2: put_att_text'
        endif
      else
        write(*,*) '### ERROR 1: nf_create'
      endif
      if (ierr /= nf_noerr) stop

! Copy all attributes for time variable while the target file is still
! in redefinition mode, then switch it into input and copy time itself.

      ierr=nf_inq_varid(nctarg, time_var_name, tvar_out)
      if (ierr == nf_noerr) then
        ierr=nf_inq_varid(ncsrc, time_var_name, tvar_in)
        if (ierr == nf_noerr) then
          ierr=nf_inq_varnatts(ncsrc, tvar_in, natts)
          if (ierr == nf_noerr) then
            do i=1,natts
              ierr=nf_inq_attname(ncsrc, tvar_in, i, str)
              if (ierr == nf_noerr) then
                call lenstr(str,lstr)
                ierr=nf_copy_att(ncsrc, tvar_in, str(1:lstr),
     & nctarg, tvar_out)
                if (ierr == nf_noerr) then
                  write(*,*) 'copied attribute ''',str(1:lstr),''''
                else
                  write(*,*) '### ERROR 10: copy_att' ; stop
                endif
              else
                write(*,*) '### ERROR 9: inq_attname' ; stop
              endif
            enddo
          else
            write(*,*) '### ERROR 8: inq_varnatts'
          endif
        else
          write(*,*) '### ERROR 7: inq_varid'
        endif
      else
        write(*,*) '### ERROR 6: inq_varid'
      endif
      if (ierr /= nf_noerr) stop

      ierr=nf_enddef(nctarg) !<-- set to input mode
      if (ierr == nf_noerr) then
        ierr=nf_get_vara_double(ncsrc, tvar_in, rec,1, cff1)
        if (ierr == nf_noerr) then
          ierr=nf_put_var_double(nctarg, tvar_out, cff1)
          if (ierr == nf_noerr) then
            write(*,*) 'wrote ''', time_var_name, ''''
          else
            write(*,*) '### ERROR 13: put_var'
          endif
        else
          write(*,*) '### ERROR 12: get_vara'
        endif
      else
        write(*,*) '### ERROR 11: nf_enddef'
      endif
      if (ierr /= nf_noerr) stop


! Process variables
      call get_patch_by_name_real(ncsrc, prnt_data, 'zeta',
     & iwestpg,jsouthpg, nx,ny,0,rec, wrk1)
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1)
      call etch_into_land_thread(nx,ny, mskp, wrk1)
C$OMP END PARALLEL
      call spln2d_real(nx,ny, wrk1, srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, wrk1, srX,srY,sXY, ncx,ncy,
C$OMP& ip,jp, xi,eta, wrk2)
      call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY,
     & ncx,ncy, ip,jp, xi,eta, wrk2)
C$OMP END PARALLEL
      call apply_mask(ncx,ncy,1, mask,wrk2)
      call put_rec_by_name_real(nctarg,roms_init, 'zeta',
     & ncx,ncy,0,1, wrk2)
# 720 "r2r_init.F"
      call get_patch_by_name_real(ncsrc, prnt_data, 'u',iwestpg,
     & jsouthpg,nx-1,ny,Np, rec, wrk1)
      do k=1,Np
        isrc=1+(k-1)*(nx-1)*ny ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, umsp, wrk1, isrc)
        call etch_into_land_thread(nx-1,ny, umsp, wrk1(isrc))
C$OMP END PARALLEL
        call spln2d_real(nx-1,ny, wrk1(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc, wrk1, srX,srY,sXY,
C$OMP& ncx,ncy, itrg, ipu,jpu, xiu,etau, wrmx)
        call spln2d_interp_real(nx-1,ny, wrk1(isrc), srX,srY,sXY,
     & ncx,ncy, ipu,jpu, xiu,etau, wrmx(itrg))
C$OMP END PARALLEL
        write(6,'(A)',advance='no') '.' ; flush(6)
      enddo
      write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask, Np,wrmx, N,kprnt, wrk1)
      call r2r_vrtint_thread(ncx,ncy, 0,mask, btm_slp, Np,wrmx,
     & N,kprnt, wrk1)
C$OMP END PARALLEL
# 751 "r2r_init.F"
      call get_patch_by_name_real(ncsrc, prnt_data, 'v',iwestpg,
     & jsouthpg,nx,ny-1,Np, rec, wrk2)
      do k=1,Np
        isrc=1+(k-1)*nx*(ny-1) ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, vmsp, wrk2, isrc)
        call etch_into_land_thread(nx,ny-1, vmsp, wrk2(isrc))
C$OMP END PARALLEL
        call spln2d_real(nx,ny-1, wrk2(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc, wrk2, srX,srY,sXY,
C$OMP& ncx,ncy, itrg, ipv,jpv, xiv,etav, wrmx)
        call spln2d_interp_real(nx,ny-1, wrk2(isrc), srX,srY,sXY,
     & ncx,ncy, ipv,jpv, xiv,etav, wrmx(itrg))
C$OMP END PARALLEL
        write(6,'(A)',advance='no') '.' ; flush(6)
      enddo
      write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask, Np,wrmx, N,kprnt, wrk2)
      call r2r_vrtint_thread(ncx,ncy, 0,mask, btm_slp, Np,wrmx,
     & N,kprnt, wrk2)
C$OMP END PARALLEL
# 783 "r2r_init.F"
C$OMP PARALLEL SHARED(ncx,ncy, N, csA,snA, wrk1,wrk2, wrk3,wrmx)
      call r2r_rotate_shift_thread(ncx,ncy,N, csA,snA, wrk1,wrk2,
     & wrk3,wrmx)
C$OMP END PARALLEL
      call apply_mask(ncx-1,ncy,N, umask, wrk3)
      call apply_mask(ncx,ncy-1,N, vmask, wrmx)

      call put_rec_by_name_real(nctarg, roms_init, 'u',
     & ncx-1,ncy,N,1, wrk3)
      call put_rec_by_name_real(nctarg, roms_init, 'v',
     & ncx,ncy-1,N,1, wrmx)


C$OMP PARALLEL SHARED(ncx,ncy, N, hc,Cs_w, h, wrk3,wrmx, wrk1,wrk2)
      if ( hmin<=0) then
        h=h+hmin-0.2
      endif
      call compute_uvbar_thread(ncx,ncy, N, hc,Cs_w, h, wrk3,wrmx,
     & wrk1,wrk2)
C$OMP END PARALLEL
      call apply_mask(ncx-1,ncy,1, umask, wrk1)
      call apply_mask(ncx,ncy-1,1, vmask, wrk2)

      call put_rec_by_name_real(nctarg, roms_init, 'ubar',
     & ncx-1,ncy,0,1, wrk1)
      call put_rec_by_name_real(nctarg, roms_init, 'vbar',
     & ncx,ncy-1,0,1, wrk2)
# 819 "r2r_init.F"
      do itrc=1,ntrc/ipt_char_len
        call get_patch_by_name_real(ncsrc, prnt_data,trim(tracer(itrc)),
     & iwestpg,jsouthpg, nx,ny,Np,rec, wrk1)
        do k=1,Np
          isrc=1+(k-1)*nx*ny ; itrg=1+(k-1)*ncx*ncy
C$OMP PARALLEL SHARED(nx,ny, wrk1, isrc)
          call etch_into_land_thread(nx,ny, mskp, wrk1(isrc))
C$OMP END PARALLEL
          call spln2d_real(nx,ny, wrk1(isrc), srX,srY,sXY,sYX)
C$OMP PARALLEL SHARED(nx,ny, isrc, wrk1, srX,srY,sXY,
C$OMP& ncx,ncy, itrg, ip,jp, xi,eta, wrmx)
          call spln2d_interp_real(nx,ny, wrk1(isrc), srX,srY,sXY,
     & ncx,ncy, ip,jp, xi,eta, wrmx(itrg))
C$OMP END PARALLEL
          write(6,'(A)',advance='no') '.' ; flush(6)
        enddo
        write(6,*)
C$OMP PARALLEL SHARED(ncx,ncy, mask, Np,wrmx, N,kprnt,wrk1)
      call r2r_vrtint_thread(ncx,ncy, 1,mask, btm_trc, Np,wrmx,
     & N,kprnt,wrk1)
C$OMP END PARALLEL
        call apply_mask(ncx,ncy,N, mask,wrk1)
        call put_rec_by_name_real(nctarg,roms_init,trim(tracer(itrc)),
     & ncx,ncy,N,1, wrk1)
# 852 "r2r_init.F"
      enddo !<--itrc

      ierr=nf_close(nctarg)
# 873 "r2r_init.F"
      call lenstr(roms_init,lstr)
      write(*,'(/1x,5A/)') 'Files ''', roms_init(1:lstr),
     & ''' and ''', 'croco_init_diag.nc', ''' are ready.'
      end subroutine r2r_init
# 61 "tools_fort.F" 2
# 1 "r2r_bry.F" 1
      subroutine r2r_bry (ntrc,nhists,chldgrd,theta_s,theta_b,hc,N, WESN
     & ,prntgrd,prnt_data,tracer)

! Creates open boundary forcing file by extracting and interpolating
! relevant data from parent-grid solution. Should be used as
!
! r2r_bry child_grid.nc theta_s theta_b hc N WESN parent_grid.nc
! parent_his1.nc parent_his2.nc ...
!
! where theta_s theta_b hc N specify vertical coordinate of child-grid
! configuration; argument WESN (stands for West, East, North, South)
! specifies which sides are open; arguments beyond the parent grid file
! name are parent-grid history files, may be one or more, may be named
! individually or via wildcards.

! CPP-switch make this program read only relevant
! subdomain within the parent grid as opposite to reading the whole
! data. Refer to "r2r_init.F" for pros and cons of selecting it.

! Method: 2D horizontal interpolation is by bi-cubic pseudo-splines
! via Hermite basis functions followed by vertical splines by spline
! inversion (z-levels of child grid are translated into continuous
! k-index space of parent grid in such a way spline interpolation of
! the resultant kp=kp(z_r) values back into z-space of child grid
! yields the initial z_r=z_r(k) exactly. In doing so verticall spline
! derivatives for both z_r of parent grid and for the 3D fields to be
! interpolated are constructed in k-index space of horizontally-child
! but still vertically-parent grid). All horizontal interpolations are
! direct from source to destination locations on C-grids:
! RHO_parent --> RHO_child, U-->U, V-->V,
! V-->U, V-->U, RHO-->U, and RHO-->V
! overall a total of 7 permutations. Cross-interpolations, V-->U and
! V-->U, are needed for the rotation of vector components; RHO-->U and
! RHO-->V are is for computing barotropic mode (vertical averaging).
! Before horizontal interpolation data is extended to land-masked areas
! by etching.

! Parameters setting "btm_slp" and "btm_trc" control type of boundary
! condition at bottom for the purpose of computing spline derivatives:
! btm_slp=-1 is no-slip, meaning that u,v=0 exactly at the bottom, i.e.
! half-grid-interval below the k=1 value; This applies to velocities
! only; for tracers the meaningful b.c. are natural on Neumann (0 and
! +1 respectively). Top b.c. is always natural.

! Created and maintained by Alexander Shchepetkin, old_galaxy@yahoo.com


c--#define VERBOSE

      use mod_io_size_acct
      implicit none
      character(len=80) :: chldgrd, prntgrd
      integer :: nhists,ipt_char_len
      character(len=160), dimension(nhists) :: prnt_data
      character(len=160) :: prt_tmp
      character(len=16) :: roms_bry
      character(len=10) :: time_var_name
      character(len=16) :: VertCoordType, vname, trgname
      character(len=180) :: time_units, str
      integer :: ntrc,ipt_trc_len
      character(len=4) :: WESN
      character(len=20), dimension(ntrc):: tracer
      integer :: btm_slp ! no-slip=-1; +1 free; 0 natural
     & , btm_trc ! +1 Neumann; 0 natural;
      logical :: OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH
      real(kind=8) :: theta_s,theta_b, hc,hcp, xcmin,xcmax, time,
     & cff, cs,csP,sn,snP
      real(kind=8), allocatable, dimension(:,:) :: csA,snA, h,hprnt,
     & hp, xc,yc, xp,yp, xpu,ypu, xpv,ypv
      integer(kind=2), allocatable, dimension(:,:) :: mskp, umsp,vmsp
      real(kind=8), allocatable, dimension(:) :: Cs_w,Cs_r,Csp_w,Csp_r,
     & xr_west,yr_west, xu_west,yu_west, xiur_west,etaur_west,
     & xiu_west,etau_west, xiuv_west,etauv_west,
     & xi_west,eta_west, xv_west,yv_west, xivr_west,etavr_west,
     & xiv_west,etav_west, xivu_west,etavu_west,
     & csAu_west,snAu_west, csAv_west,snAv_west,
     & h_west,hp_west, hu_west,hpu_west, hv_west,hpv_west,
     & kpr_west, kpu_west, kpv_west,

     & xr_east,yr_east, xu_east,yu_east, xiur_east,etaur_east,
     & xiu_east,etau_east, xiuv_east,etauv_east,
     & xi_east,eta_east, xv_east,yv_east, xivr_east,etavr_east,
     & xiv_east,etav_east, xivu_east,etavu_east,
     & csAu_east,snAu_east, csAv_east,snAv_east,
     & h_east,hp_east, hu_east,hpu_east, hv_east,hpv_east,
     & kpr_east, kpu_east, kpv_east

       real(kind=8), allocatable, dimension(:) ::
     & xr_south,yr_south, xu_south,yu_south, xiur_south,etaur_south,
     & xiu_south,etau_south, xiuv_south,etauv_south,
     & xi_south,eta_south, xv_south,yv_south, xivr_south,etavr_south,
     & xiv_south,etav_south, xivu_south,etavu_south,
     & csAu_south,snAu_south, csAv_south,snAv_south,
     & h_south,hp_south, hu_south,hpu_south, hv_south,hpv_south,
     & kpr_south, kpu_south, kpv_south,

     & xr_north,yr_north, xu_north,yu_north, xiur_north,etaur_north,
     & xiu_north,etau_north, xiuv_north,etauv_north,
     & xi_north,eta_north, xv_north,yv_north, xivr_north,etavr_north,
     & xiv_north,etav_north, xivu_north,etavu_north,
     & csAu_north,snAu_north, csAv_north,snAv_north,
     & h_north,hp_north, hu_north,hpu_north, hv_north,hpv_north,
     & kpr_north, kpu_north, kpv_north

      integer(kind=4), allocatable, dimension(:) ::
     & ir_west,jr_west, iu_west,ju_west, iv_west,jv_west,
     & iur_west,jur_west, ivr_west,jvr_west,
     & iuv_west,juv_west, ivu_west,jvu_west,

     & ir_east,jr_east, iu_east,ju_east, iv_east,jv_east,
     & iur_east,jur_east, ivr_east,jvr_east,
     & iuv_east,juv_east, ivu_east,jvu_east,

     & ir_south,jr_south, iu_south,ju_south, iv_south,jv_south,
     & iur_south,jur_south, ivr_south,jvr_south,
     & iuv_south,juv_south, ivu_south,jvu_south,

     & ir_north,jr_north, iu_north,ju_north, iv_north,jv_north,
     & iur_north,jur_north, ivr_north,jvr_north,
     & iuv_north,juv_north, ivu_north,jvu_north
      integer(kind=2), allocatable, dimension(:) ::
     & mskr_west, mskr_east, mskr_south, mskr_north,
     & msku_west, msku_east, msku_south, msku_north,
     & mskv_west, mskv_east, mskv_south, mskv_north

      real(kind=8), allocatable, dimension(:) :: srX, srY, sXY, sYX
      real(kind=4), allocatable, dimension(:) :: wrk1,wrk2, wrk3,wrk4,
     & wrk5,wrk6, wrk7,wrk8
      integer :: net_alloc_mem, nargs, nx,ny,Np, ncx,ncy,N, ihis,
     & nccgrd, ncpgrd, ncsrc, nctarg, nrecs,rec, recout, varid,
     & vtype, ndims, vdimids(4), natts, tvar_in,tvar_out, ierr,
     & itrc,size,isrc, i,j,k, lstr,lvar,ltgv,lcgrd,lpgrd,lprnt
      real(kind=8) :: hmin
      character*30 :: orig_date

      integer :: iwestpg,jsouthpg, imin,imax,jmin,jmax
      integer :: margn




      real(kind=4) tstart,tend
      integer nclk, iclk(2), clk_rate, clk_max
      integer(kind=8) :: inc_clk, net_clk=0

      include "phys_const.h"
      include "netcdf.inc"
Cf2py intent(in) ntrc,nhists,chldgrd,theta_s,theta_b,hc,N, WESN , prntgrd, prnt_data,tracer

      call cpu_time(tstart)
      nclk=1
      call system_clock(iclk(nclk), clk_rate, clk_max)




      btm_slp=-1 ; btm_trc=+1;margn=8 ! parameters
      roms_bry='croco_chd_bry.nc'
      time_var_name='scrum_time'
      ipt_char_len=160;ipt_trc_len=20
      ihis=1 ; rec=1; recout=0 ; !time_units='                  '
      OBC_WEST=.false. ; OBC_SOUTH=.false.
      OBC_EAST=.false. ; OBC_NORTH=.false.

      call lenstr(chldgrd,lcgrd)
      ierr=nf_open(chldgrd(1:lcgrd), nf_nowrite, nccgrd)
      if (ierr /= nf_noerr) then
          write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #1 :: Cannot ',
     & 'open ''', chldgrd(1:lcgrd), '''.', nf_strerror(ierr)
          stop
      endif

      call lenstr(WESN,lvar)
      do i=1,lvar
        if (WESN(i:i) == 'W' .or. WESN(i:i) == 'w') then
          OBC_WEST=.true.
        elseif (WESN(i:i) == 'E' .or. WESN(i:i) == 'e') then
          OBC_EAST=.true.
        elseif (WESN(i:i) == 'S' .or. WESN(i:i) == 's') then
          OBC_SOUTH=.true.
        elseif (WESN(i:i) == 'N' .or. WESN(i:i) == 'n') then
          OBC_NORTH=.true.
        else
          write(*,'(/1x,6A/12x,A/)') '### ERROR: Wrong ',
     & 'argument # 6 ''',WESN(1:lvar),''': letter ''',
     & WESN(i:i), ''' should not be present.',
     & 'Only "W", "E", "S", and "N" are allowed.'
          stop
        endif
      enddo

      call lenstr(prntgrd,lpgrd)
      ierr=nf_open(prntgrd(1:lpgrd), nf_nowrite, ncpgrd)
      if (ierr /= nf_noerr) then
          write(*,'(/1x,4A/22x,A/)') '### ERROR: arg #6 ::',
     & ' Cannot open ''', prntgrd(1:lpgrd), '''.',
     & nf_strerror(ierr)
          stop
      endif


      do i=1,nhists/ipt_char_len ! Open and close
         prt_tmp=trim(prnt_data(i))
         call lenstr(prt_tmp,lprnt)
         if (i == 1) then
          ierr=nf_open(prt_tmp(1:lprnt),
     & nf_nowrite, ncsrc)
         else ! using "j"
           ierr=nf_open(prt_tmp(1:lprnt), ! as netCDF
     & nf_nowrite, j) ! file ID
         endif ! here is
         if (ierr == nf_noerr) then ! just for
           if (i > 1) ierr=nf_close(j) ! testing
         else
           write(*,'(/1x,A,I3,1x,3A/24x,A/)')
     & '### ERROR: arg #', i+7, ':: Cannot open ''',
     & prt_tmp(1:lprnt),'''.', nf_strerror(ierr)
           stop
         endif

      enddo

      call roms_find_dims(nccgrd, chldgrd, i,j,k)
      ncx=i+2 ; ncy=j+2
      call roms_find_dims(ncsrc, prnt_data(1), i,j,k)
      call roms_check_dims(ncpgrd, prntgrd, i,j,0)
      nx=i+2 ; ny=j+2 ; Np=k
      allocate(Cs_r(N), Cs_w(0:N), Csp_w(0:Np),Csp_r(Np))
      call set_scoord(theta_s,theta_b, N, Cs_r,Cs_w)
      call read_scoord(ncsrc, Np, Csp_r,Csp_w,hcp, VertCoordType)
      net_alloc_mem=2*(2*N+1)+ 2*(2*Np+1) ; ierr=0

      if (OBC_WEST) then
        allocate( kpr_west(ncy*N), kpu_west(ncy*N), kpv_west(ncy*N-N),
     & xr_west(ncy),yr_west(ncy), xu_west(ncy),yu_west(ncy),
     & csAu_west(ncy),snAu_west(ncy),
     & h_west(ncy),hp_west(ncy), hu_west(ncy),hpu_west(ncy),
     & xv_west(ncy-1),yv_west(ncy-1),
     & csAv_west(ncy-1),snAv_west(ncy-1),
     & hv_west(ncy-1),hpv_west(ncy-1),

     & xi_west(ncy),eta_west(ncy), xiu_west(ncy),etau_west(ncy),
     & xiuv_west(ncy),etauv_west(ncy),
     & xiur_west(ncy),etaur_west(ncy),
     & xivr_west(ncy-1),etavr_west(ncy-1),
     & xivu_west(ncy-1),etavu_west(ncy-1),
     & xiv_west(ncy-1),etav_west(ncy-1),

     & ir_west(ncy),jr_west(ncy), iu_west(ncy),ju_west(ncy),
     & iuv_west(ncy),juv_west(ncy),
     & iur_west(ncy),jur_west(ncy),
     & ivr_west(ncy-1),jvr_west(ncy-1),
     & ivu_west(ncy-1),jvu_west(ncy-1),
     & iv_west(ncy-1),jv_west(ncy-1),
     & mskr_west(ncy), msku_west(ncy), mskv_west(ncy-1),
     & stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncy*N + 2*(N*ncy-1)
     & +18*2*ncy +12*2*(ncy-1) +8*ncy +6*(ncy-1) +(3*ncy-1)/2
        write(*,*) 'allocated western  side coordinate arrays'
      endif

      if (OBC_EAST .and. ierr == 0 ) then
        allocate(kpr_east(ncy*N), kpu_east(ncy*N), kpv_east(ncy*N-N),
     & xr_east(ncy),yr_east(ncy), xu_east(ncy),yu_east(ncy),
     & csAu_east(ncy),snAu_east(ncy),
     & h_east(ncy),hp_east(ncy), hu_east(ncy),hpu_east(ncy),
     & xv_east(ncy-1),yv_east(ncy-1),
     & csAv_east(ncy-1),snAv_east(ncy-1),
     & hv_east(ncy-1),hpv_east(ncy-1),

     & xi_east(ncy),eta_east(ncy), xiu_east(ncy),etau_east(ncy),
     & xiuv_east(ncy),etauv_east(ncy),
     & xiur_east(ncy),etaur_east(ncy),
     & xivr_east(ncy-1),etavr_east(ncy-1),
     & xivu_east(ncy-1),etavu_east(ncy-1),
     & xiv_east(ncy-1),etav_east(ncy-1),

     & ir_east(ncy),jr_east(ncy), iu_east(ncy),ju_east(ncy),
     & iuv_east(ncy),juv_east(ncy),
     & iur_east(ncy),jur_east(ncy),
     & ivr_east(ncy-1),jvr_east(ncy-1),
     & ivu_east(ncy-1),jvu_east(ncy-1),
     & iv_east(ncy-1),jv_east(ncy-1),
     & mskr_east(ncy), msku_east(ncy), mskv_east(ncy-1),
     & stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncy*N + 2*(N*ncy-1)
     & +18*2*ncy +12*2*(ncy-1) +8*ncy +6*(ncy-1) +(3*ncy-1)/2
        write(*,*) 'allocated eastern  side coordinate arrays'
      endif

      if (OBC_SOUTH .and. ierr == 0 ) then
        allocate(kpr_south(ncx*N),kpu_south(ncx*N-N),kpv_south(ncx*N),
     & xr_south(ncx),yr_south(ncx), xv_south(ncx),yv_south(ncx),
     & csAv_south(ncx),snAv_south(ncx),
     & h_south(ncx),hp_south(ncx), hv_south(ncx),hpv_south(ncx),
     & xu_south(ncx-1),yu_south(ncx-1),
     & csAu_south(ncx-1),snAu_south(ncx-1),
     & hu_south(ncx-1),hpu_south(ncx-1),
     & xi_south(ncx),eta_south(ncx),
     & xiv_south(ncx),etav_south(ncx),
     & xivu_south(ncx),etavu_south(ncx),
     & xivr_south(ncx),etavr_south(ncx),
     & xiur_south(ncx-1),etaur_south(ncx-1),
     & xiuv_south(ncx-1),etauv_south(ncx-1),
     & xiu_south(ncx-1),etau_south(ncx-1),
     & ir_south(ncx),jr_south(ncx),
     & iv_south(ncx),jv_south(ncx),
     & ivu_south(ncx),jvu_south(ncx),
     & ivr_south(ncx),jvr_south(ncx),
     & iur_south(ncx-1),jur_south(ncx-1),
     & iuv_south(ncx-1),juv_south(ncx-1),
     & iu_south(ncx-1),ju_south(ncx-1),
     & mskr_south(ncx), msku_south(ncx-1), mskv_south(ncx),
     & stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncx*N + 2*(N*ncx-1)
     & +18*2*ncx +12*2*(ncx-1) +8*ncx +6*(ncx-1)+(3*ncx-1)/2
        write(*,*) 'allocated southern side coordinate arrays'
      endif

      if (OBC_NORTH .and. ierr == 0 ) then
        allocate(kpr_north(ncx*N),kpu_north(ncx*N-N),kpv_north(ncx*N),
     & xr_north(ncx),yr_north(ncx), xv_north(ncx),yv_north(ncx),
     & csAv_north(ncx),snAv_north(ncx),
     & h_north(ncx),hp_north(ncx), hv_north(ncx),hpv_north(ncx),
     & xu_north(ncx-1),yu_north(ncx-1),
     & csAu_north(ncx-1),snAu_north(ncx-1),
     & hu_north(ncx-1),hpu_north(ncx-1),
     & xi_north(ncx),eta_north(ncx),
     & xiv_north(ncx),etav_north(ncx),
     & xivu_north(ncx),etavu_north(ncx),
     & xivr_north(ncx),etavr_north(ncx),
     & xiur_north(ncx-1),etaur_north(ncx-1),
     & xiuv_north(ncx-1),etauv_north(ncx-1),
     & xiu_north(ncx-1),etau_north(ncx-1),
     & ir_north(ncx),jr_north(ncx),
     & iv_north(ncx),jv_north(ncx),
     & ivu_north(ncx),jvu_north(ncx),
     & ivr_north(ncx),jvr_north(ncx),
     & iur_north(ncx-1),jur_north(ncx-1),
     & iuv_north(ncx-1),juv_north(ncx-1),
     & iu_north(ncx-1),ju_north(ncx-1),
     & mskr_north(ncx), msku_north(ncx-1), mskv_north(ncx),
     & stat=ierr )
        net_alloc_mem=net_alloc_mem +2*2*N*ncx*N + 2*(N*ncx-1)
     & +18*2*ncx +12*2*(ncx-1) +8*ncx +6*(ncx-1) +(3*ncx-1)/2
        write(*,*) 'allocated northern side coordinate arrays'
      endif

      write(*,'(8x,A,F9.2,1x,A)') 'allocated', dble(net_alloc_mem)
     & /262144.D0,'MB total'

      allocate( csA(ncx,ncy),snA(ncx,ncy), h(ncx,ncy), hprnt(ncx,ncy),
     & stat=ierr)
      if (ierr == 0 ) then
        net_alloc_mem=net_alloc_mem +4*2*ncx*ncy
        allocate(xc(ncx,ncy),yc(ncx,ncy), stat=ierr)
        if (ierr == 0 ) then
          net_alloc_mem=net_alloc_mem +2*2*ncx*ncy

          allocate(xp(nx,ny),yp(nx,ny), stat=ierr)
          if (ierr == 0) then
            net_alloc_mem=net_alloc_mem +2*2*nx*ny
          endif

        endif
      endif

      if (ierr == 0) then
        write(*,'(8x,A,F9.2,1x,A)') 'allocated', dble(net_alloc_mem)
     & /262144.D0,'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 1.'
        stop
      endif

! Setup child-grid land mask: ! This program does not use 2D arrays
!------ ---------- ---- ----- ! for child-grid mask, but rather sets up
                              ! individual 1D arrays for RHO- U- and V-
      if (OBC_WEST) then ! points at each open boundary. Because
        mskr_west(:)=1 ! land mask may or may not be present for
        msku_west(:)=1 ! for the particular configuration, set
        mskv_west(:)=1 ! all masks to "water" everywhere status
      endif ! by assigning mask=1 (on the left), then
      if (OBC_EAST) then ! using native netCDF functions check
        mskr_east(:)=1 ! whether variable "mask_rho" is present
        msku_east(:)=1 ! in netCDF file: if yes, proceed reading
        mskv_east(:)=1 ! and set mask=0 where it is land; if no
      endif ! then all-water status on the left will
      if (OBC_SOUTH) then ! remain unchanged.
        mskr_south(:)=1
        mskv_south(:)=1 ! Temporarily use array "h" to read mask
        msku_south(:)=1 ! from the file as it is expected to be
      endif ! of real type there -- thereafter this
      if (OBC_NORTH) then ! program uses only 2-byte integer
        mskr_north(:)=1 ! version of mask.
        mskv_north(:)=1
        msku_north(:)=1
      endif ! Note that code

                                                    ! below relies on
      ierr=nf_inq_varid(nccgrd, 'mask_rho', varid) ! having all mask
      if (ierr == nf_noerr) then ! arrays mskX_side
        ierr=nf_get_var_double(nccgrd, varid, h) ! initialized to 1
        if (ierr == nf_noerr) then ! everywhere.

          if (OBC_WEST) then
            do j=1,ncy
              if (h(1,j) < 0.5D0) mskr_west(j)=0
              if (h(1,j) < 0.5D0 .and. h(2,j) < 0.5D0) msku_west(j)=0
            enddo
            do j=1,ncy-1
              if (h(1,j) < 0.5D0 .and. h(1,j+1) < 0.5D0) mskv_west(j)=0
            enddo
          endif
          if (OBC_EAST) then
            do j=1,ncy
              if (h(ncx,j) < 0.5D0) mskr_east(j)=0
              if (h(ncx,j) < 0.5D0 .and. h(ncx-1,j) < 0.5D0)
     & msku_east(j)=0
            enddo
            do j=1,ncy-1
              if (h(ncx,j) < 0.5D0 .and. h(ncx,j+1) < 0.5D0)
     & mskv_east(j)=0
            enddo
          endif
          if (OBC_SOUTH) then
            do i=1,ncx
              if (h(i,1) < 0.5D0) mskr_south(i)=0
              if (h(i,1) < 0.5D0 .and. h(i,2) < 0.5D0) mskv_south(i)=0
            enddo
            do i=1,ncx-1
              if (h(i,1) < 0.5D0 .and. h(i+1,1)< 0.5D0) msku_south(i)=0
            enddo
          endif
          if (OBC_NORTH) then
            do i=1,ncx
              if (h(i,ncy) < 0.5D0) mskr_north(i)=0
              if (h(i,ncy) < 0.5D0 .and. h(i,ncy-1) < 0.5D0)
     & mskv_north(i)=0
            enddo
            do i=1,ncx-1
              if (h(i,ncy) < 0.5D0 .and. h(i+1,ncy) < 0.5D0)
     & msku_north(i)=0
            enddo
          endif

        else
          write(*,'(/1x,5A/)') '### ERROR: Cannot read ''mask_rho'' ',
     & 'from ''', chldgrd(1:lcgrd), ''', ', nf_strerror(ierr)
          stop
        endif
      else
        write(*,'(9x,4A)') 'Variable ''mask_rho'' does not exist ',
     & 'in ''', chldgrd(1:lcgrd), ''', assuming mask=1 everywhere.'
      endif


! Read child-grid horizontal coordinates, extract and save them along
!------ ---------- --------- ----------- the open boundary lines. Only
! "r"-versions, xr_side,yr_side, are used initially for the parent grid
! index search in the following segment below. The other two, "u"- and
! "v"-types will used during the final search and initialization of the
! horizontal interpolation after identifying the relevant subdomain
! within the parent-grid.

      call get_var_by_name_double(nccgrd, 'lon_rho', xc)
      call get_var_by_name_double(nccgrd, 'lat_rho', yc)

      if (OBC_WEST) then
        do j=1,ncy
          xr_west(j)=xc(1,j); yr_west(j)=yc(1,j)
          xu_west(j)=0.5D0*(xc(1,j)+xc(2,j))
          yu_west(j)=0.5D0*(yc(1,j)+yc(2,j))
        enddo
        do j=1,ncy-1
          xv_west(j)=0.5D0*(xc(1,j)+xc(1,j+1))
          yv_west(j)=0.5D0*(yc(1,j)+yc(1,j+1))
        enddo
      endif

      if (OBC_EAST) then
        do j=1,ncy
          xr_east(j)=xc(ncx,j); yr_east(j)=yc(ncx,j)
          xu_east(j)=0.5D0*(xc(ncx,j)+xc(ncx-1,j))
          yu_east(j)=0.5D0*(yc(ncx,j)+yc(ncx-1,j))
        enddo
        do j=1,ncy-1
          xv_east(j)=0.5D0*(xc(ncx,j)+xc(ncx,j+1))
          yv_east(j)=0.5D0*(yc(ncx,j)+yc(ncx,j+1))
        enddo
      endif

      if (OBC_SOUTH) then
        do i=1,ncx
          xr_south(i)=xc(i,1); yr_south(i)=yc(i,1)
          xv_south(i)=0.5D0*(xc(i,1)+xc(i,2))
          yv_south(i)=0.5D0*(yc(i,1)+yc(i,2))
        enddo
        do i=1,ncx-1
          xu_south(i)=0.5D0*(xc(i,1)+xc(i+1,1))
          yu_south(i)=0.5D0*(yc(i,1)+yc(i+1,1))
        enddo
      endif

      if (OBC_NORTH) then
        do i=1,ncx
          xr_north(i)=xc(i,ncy) ; yr_north(i)=yc(i,ncy)
          xv_north(i)=0.5D0*(xc(i,ncy)+xc(i,ncy-1))
          yv_north(i)=0.5D0*(yc(i,ncy)+yc(i,ncy-1))
        enddo
        do i=1,ncx-1
          xu_north(i)=0.5D0*(xc(i,ncy)+xc(i+1,ncy))
          yu_north(i)=0.5D0*(yc(i,ncy)+yc(i+1,ncy))
        enddo
      endif

      write(*,'(1x,A)',advance='no') 'child grid longitude '
      call compute_min_max(ncx,ncy, xc, xcmin,xcmax)



! Preliminary step: read horizontal coordinates for the parent
!------------ ----- grid. then pretend initializing parent --> child
! interpolation, but actually all what we need is (ip,jp)-indices: use
! them to find the smallest logically rectangular subdomain within the
! parent grid which encloses the child (to be more specifically the
! unmasked portion of it). This region is characterized by bounds
! imin,imax,jmin,jmax defined within the parent grid index space.
! Expand the region in all four directions by making "margn"-wide halo
! to allow more points to reduce the influence of artificial boundary
! conditions associated with spline interpolations, then save
! coordinates of south-west corner of the expanded region and redefine
! nx,ny consistently with its size. Thereafter deallocate coordinate
! arrays for both parent and child -- basically all what this step is
! needed for is just 4 integer numbers: iwestpg,jsouthpg,nx,ny.

      call get_var_by_name_double(ncpgrd, 'lon_rho', xp)
      call get_var_by_name_double(ncpgrd, 'lat_rho', yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)

      imin=nx+1 ; imax=-1 !<-- initialize to unrealistic
      jmin=ny+1 ; jmax=-1 !<-- values outside the range

      if (OBC_WEST) then
        call r2r_bry_search( nx,ny, xp,yp, ncy, xr_west,yr_west,
     & ir_west,jr_west)
        call r2r_bry_index_bounds(ncy, ir_west,jr_west, mskr_west,
     & imin,imax,jmin,jmax)
      endif

      if (OBC_EAST) then
        call r2r_bry_search( nx,ny, xp,yp, ncy, xr_east,yr_east,
     & ir_east,jr_east)
        call r2r_bry_index_bounds(ncy, ir_east,jr_east, mskr_east,
     & imin,imax,jmin,jmax)
      endif

      if (OBC_SOUTH) then
        call r2r_bry_search( nx,ny, xp,yp, ncx, xr_south,yr_south,
     & ir_south,jr_south)
        call r2r_bry_index_bounds(ncx, ir_south,jr_south, mskr_south,
     & imin,imax,jmin,jmax)
      endif

      if (OBC_NORTH) then
        call r2r_bry_search( nx,ny, xp,yp, ncx, xr_north,yr_north,
     & ir_north,jr_north)
        call r2r_bry_index_bounds(ncx, ir_north,jr_north, mskr_north,
     & imin,imax,jmin,jmax)
      endif

      imax=imax+1 ; jmax=jmax+1
      write(*,'(1x,2A/2(4x,A,2I5,1x,A,I5))') 'minimal parent-grid ',
     & 'index bounds to accommodate child grid:',
     & 'imin,imax =', imin,imax, 'of nx =', nx,
     & 'jmin,jmax =', jmin,jmax, 'of ny =', ny
      imin=max(imin-margn,1) ; jmin=max(jmin-margn,1)
      imax=min(imax+margn,nx) ; jmax=min(jmax+margn,ny)
      write(*,'(1x,A/4x,A,2I5,17x,A,2I5)') 'adjusted to',
     & 'imin,imax =',imin,imax, 'jmin,jmax =',jmin,jmax

      iwestpg=imin ; nx=imax-imin+1 ; jsouthpg=jmin ; ny=jmax-jmin+1
      write(*,'(/2(2x,A,I5)/)') 'setting subdomain sizes  nx =', nx,
     & 'ny =', ny
      deallocate(xp,yp)
      net_alloc_mem=net_alloc_mem -2*2*nx*ny
      write(*,'(6x,A,F9.2,1x,A)') 'deallocated',
     & dble(2*2*nx*ny)/262144.D0,'MB'


! Note that the only outcome of the segment above is 4 integer numbers,
! iwestpg,jsouthpg, nx,ny
! which are the indices of south-west corner of the subdomain within
! the parent grid, and the sizes of the subdomain. Everything else is
! discarded. The next step is to allocate subdomain-sized 2D arrays
! for parent-grid variables: note that xp(nx,ny),yp(nx,ny) were
! deallocated above and now are allocated again with a different
! (expected to be smaller) size).

      allocate( srX(nx*ny), srY(nx*ny), sXY(nx*ny), sYX(nx*ny),
     & stat=ierr )
      if (ierr == 0) then
        allocate(hp(nx,ny), mskp(nx,ny), umsp(nx-1,ny), vmsp(nx,ny-1),
     & xp(nx,ny), yp(nx,ny), xpu(nx-1,ny), xpv(nx,ny-1),
     & ypu(nx-1,ny), ypv(nx,ny-1),
     & stat=ierr )
      endif
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem +(4+7)*2*nx*ny +3*(nx*ny)/2
        write(*,'(10x,A,F9.2,1x,A)') 'allocated',dble(net_alloc_mem)
     & /262144.D0, 'MB total'
      else
        write(*,'(/1x,A/)') '### ERROR: Memory allocation failure 2.'
        stop
      endif

! Re-read (read first time if CPP is undefined)
! horizontal coordinates for parent grid. however this time only within
! the subdomain of parent grid. Add/subtract 360 degrees to/from
! longitude if necessary to be consistent with child-grid xcmin,xcmax
! range determined above.

      call get_patch_by_name_double(ncpgrd, prntgrd, 'lon_rho',
     & iwestpg,jsouthpg, nx,ny,0,0, xp)
      call get_patch_by_name_double(ncpgrd, prntgrd, 'lat_rho',
     & iwestpg,jsouthpg, nx,ny,0,0, yp)
      write(*,'(1x,A)',advance='no') 'parent grid longitude'
      call adjust_lon_into_range(nx,ny, xp, xcmin,xcmax)
      do j=1,ny
        do i=1,nx-1
          xpu(i,j)=0.5D0*(xp(i,j)+xp(i+1,j))
          ypu(i,j)=0.5D0*(yp(i,j)+yp(i+1,j))
        enddo
      enddo
      do j=1,ny-1
        do i=1,nx
          xpv(i,j)=0.5D0*(xp(i,j)+xp(i,j+1))
          ypv(i,j)=0.5D0*(yp(i,j)+yp(i,j+1))
        enddo
      enddo

! Initialize horizontal interpolation: The purpose of this stage is
!----------- ---------- -------------- to find (i) parent grid indices
! ip=ip(ic),jp=jp(ic) enclosing child-grid point "ic" into the parent
! grid element defined by 4 vertices, (ip,jp), (ip+1,jp), (ip,jp+1),
! (ip,jp+1)), and (ip+1,jp+1), and (ii) fractions xi,eta such that
! bi-linear interpolation of coordinates of the 4 vertices of parent
! grid yields coordinates (xc,yc) of the child grid.

! The semantic rules for variable names are:
!
! role placement side
! [i,j,xi,eta][r,u,v.ur,vr,uv.vu]_[west,east,south,north]
!
! allearing with all possible permutations (there are 4*7*4=102 arrays
! computed as sets of 4 in 28 calls to "bry_init_interp").
! Single-letter placement "r" indicates RHO-points of child grid
! interpolated from RHO-points of parent; similarly for "u" and "v".
! Dual-letter placements use first letter to identify location on
! the child grid and while second letter is type of source on parent.


      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(/F8.2,1x,A)') dble(net_clk)/dble(clk_rate),



     & 'initializing horizontal interpolation...'


      imin=nx+1 ; imax=-1 !<-- initialize to unrealistic
      jmin=ny+1 ; jmax=-1 !<-- values outside the range


      if (OBC_WEST) then
        call bry_init_interp(nx,ny, xp,yp, ncy, xr_west,yr_west,
     & ir_west,jr_west, xi_west,eta_west)

        call r2r_bry_index_bounds( ncy, ir_west,jr_west, mskr_west,
     & imin,imax,jmin,jmax)

        call bry_init_interp(nx,ny, xp,yp, ncy, xu_west,yu_west,
     & iur_west,jur_west, xiur_west,etaur_west)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncy, xu_west,yu_west,
     & iu_west,ju_west, xiu_west,etau_west)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncy, xu_west,yu_west,
     & iuv_west,juv_west, xiuv_west,etauv_west)

        call bry_init_interp(nx,ny, xp,yp, ncy-1, xv_west,yv_west,
     & ivr_west,jvr_west, xivr_west,etavr_west)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncy-1,xv_west,yv_west,
     & ivu_west,jvu_west, xivu_west,etavu_west)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncy-1,xv_west,yv_west,
     & iv_west,jv_west, xiv_west,etav_west)
      endif

      if (OBC_EAST) then
        call bry_init_interp(nx,ny, xp,yp, ncy, xr_east,yr_east,
     & ir_east,jr_east, xi_east,eta_east)

        call r2r_bry_index_bounds( ncy, ir_east,jr_east, mskr_east,
     & imin,imax,jmin,jmax)

        call bry_init_interp(nx,ny, xp,yp, ncy, xu_east,yu_east,
     & iur_east,jur_east, xiur_east,etaur_east)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncy, xu_east,yu_east,
     & iu_east,ju_east, xiu_east,etau_east)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncy, xu_east,yu_east,
     & iuv_east,juv_east, xiuv_east,etauv_east)

        call bry_init_interp(nx,ny, xp,yp, ncy-1, xv_east,yv_east,
     & ivr_east,jvr_east, xivr_east,etavr_east)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncy-1,xv_east,yv_east,
     & ivu_east,jvu_east, xivu_east,etavu_east)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncy-1,xv_east,yv_east,
     & iv_east,jv_east, xiv_east,etav_east)
      endif

      if (OBC_SOUTH) then
        call bry_init_interp(nx,ny, xp,yp, ncx, xr_south,yr_south,
     & ir_south,jr_south, xi_south,eta_south)

        call r2r_bry_index_bounds( ncx, ir_south,jr_south, mskr_south,
     & imin,imax,jmin,jmax)

        call bry_init_interp(nx,ny, xp,yp, ncx-1, xu_south,yu_south,
     & iur_south,jur_south, xiur_south,etaur_south)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncx-1,xu_south,yu_south,
     & iu_south,ju_south, xiu_south,etau_south)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncx-1,xu_south,yu_south,
     & iuv_south,juv_south, xiuv_south,etauv_south)

        call bry_init_interp(nx,ny, xp,yp, ncx, xv_south,yv_south,
     & ivr_south,jvr_south, xivr_south,etavr_south)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncx, xv_south,yv_south,
     & ivu_south,jvu_south, xivu_south,etavu_south)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncx,xv_south,yv_south,
     & iv_south,jv_south, xiv_south,etav_south)
      endif

      if (OBC_NORTH) then
        call bry_init_interp(nx,ny, xp,yp, ncx, xr_north,yr_north,
     & ir_north,jr_north, xi_north,eta_north)

        call r2r_bry_index_bounds( ncx, ir_north,jr_north, mskr_north,
     & imin,imax,jmin,jmax)

        call bry_init_interp(nx,ny, xp,yp, ncx-1, xu_north,yu_north,
     & iur_north,jur_north, xiur_north,etaur_north)

        call bry_init_interp(nx-1,ny,xpu,ypu, ncx-1,xu_north,yu_north,
     & iu_north,ju_north, xiu_north,etau_north)

        call bry_init_interp(nx,ny-1,xpv,ypv, ncx-1,xu_north,yu_north,
     & iuv_north,juv_north, xiuv_north,etauv_north)

        call bry_init_interp(nx,ny, xp,yp, ncx, xv_north,yv_north,
     & ivr_north,jvr_north, xivr_north,etavr_north)

        call bry_init_interp(nx-1,ny, xpu,ypu, ncx, xv_north,yv_north,
     & ivu_north,jvu_north, xivu_north,etavu_north)

        call bry_init_interp(nx,ny-1, xpv,ypv, ncx, xv_north,yv_north,
     & iv_north,jv_north, xiv_north,etav_north)
      endif


      imax=imax+1 ; jmax=jmax+1
      write(*,'(1x,2A/2(4x,A,2I5,1x,A,I5))') 're-checking parent-',
     & 'grid index bounds to accommodate child grid:',
     & 'imin,imax =', imin,imax, 'of nx =', nx,
     & 'jmin,jmax =', jmin,jmax, 'of ny =', ny
      imin=max(imin-margn,1) ; jmin=max(jmin-margn,1)
      imax=min(imax+margn,nx) ; jmax=min(jmax+margn,ny)
      write(*,'(1x,A/4x,A,2I5,17x,A,2I5)') 'adjusted to',
     & 'imin,imax =',imin,imax, 'jmin,jmax =',jmin,jmax
      if ( imin == 1 .and. imax == nx .and.
     & jmin == 1 .and. jmax == ny ) then
        write(*,*) 'parent-child bounding check passed'
      else
        write(*,*) '### ERROR: Algorithm failure.' ; stop
      endif


      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     & 'horizontal initialization complete in',
     & dble(inc_clk)/dble(clk_rate), 'sec'


! Setup (cosA,sinA) arrays to rotate velocity components at each open
!------ ----------- ----- boundary. Variable "angle" stored in netCDF
! grid files the angle between true East and local direction of ROMS
! XI-coordinate of the grid. Read it then compute csA=cos(alpha) and
! snA=sin(alpha), first for child grid then for the parent. For the
! latter temporarily place the outcome into arrays xp,yp -- after all
! horizontal interpolations have been initialized above the content
! arrays is no longer needed. Then interpolate cos(A) and sin(A) of
! the parent into child grid and compute cos and sin of the child-
! parent difference of angles; these will be used to rotate velocity
! vector components.

      call read_angle(nccgrd, chldgrd, 1,1, ncx,ncy, csA,snA)
      call read_angle(ncpgrd, prntgrd, iwestpg,jsouthpg, nx,ny, xp,yp)

      call spln2d_double(nx,ny, xp, srX,srY,sXY,sYX)
      if (OBC_WEST) then
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncy,1,
     & iur_west,jur_west, xiur_west,etaur_west, csAu_west)
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncy-1,1,
     & ivr_west,jvr_west, xivr_west,etavr_west, csAv_west)
      endif
      if (OBC_EAST) then
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncy,1,
     & iur_east,jur_east, xiur_east,etaur_east, csAu_east)
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncy-1,1,
     & ivr_east,jvr_east, xivr_east,etavr_east, csAv_east)
      endif
      if (OBC_SOUTH) then
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncx-1,1,
     & iur_south,jur_south, xiur_south,etaur_south, csAu_south)
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncx,1,
     & ivr_south,jvr_south, xivr_south,etavr_south, csAv_south)
      endif
      if (OBC_NORTH) then
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncx-1,1,
     & iur_north,jur_north, xiur_north,etaur_north, csAu_north)
        call spln2d_interp_double(nx,ny, xp, srX,srY,sXY, ncx,1,
     & ivr_north,jvr_north, xivr_north,etavr_north, csAv_north)
      endif

      call spln2d_double(nx,ny, yp, srX,srY,sXY,sYX)
      if (OBC_WEST) then
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncy,1,
     & iur_west,jur_west, xiur_west,etaur_west, snAu_west)
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncy-1,1,
     & ivr_west,jvr_west, xivr_west,etavr_west, snAv_west)
      endif
      if (OBC_EAST) then
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncy,1,
     & iur_east,jur_east, xiur_east,etaur_east, snAu_east)
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncy-1,1,
     & ivr_east,jvr_east, xivr_east,etavr_east, snAv_east)
      endif
      if (OBC_SOUTH) then
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncx-1,1,
     & iur_south,jur_south, xiur_south,etaur_south, snAu_south)
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncx,1,
     & ivr_south,jvr_south, xivr_south,etavr_south, snAv_south)
      endif
      if (OBC_NORTH) then
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncx-1,1,
     & iur_north,jur_north, xiur_north,etaur_north, snAu_north)
        call spln2d_interp_double(nx,ny, yp, srX,srY,sXY, ncx,1,
     & ivr_north,jvr_north, xivr_north,etavr_north, snAv_north)
      endif
                                                    ! Thus far arrays
      if (OBC_WEST) then ! csAu_west and
        do j=1,ncy ! snAu_west contain
          csP=csAu_west(j) ; cs=csA(1,j)+csA(2,j) ! cos(parent) and
          snP=snAu_west(j) ; sn=snA(1,j)+snA(2,j)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_west(j)=cff*(cs*csP+sn*snP)
          snAu_west(j)=cff*(sn*csP-cs*snP) ! sin(parent)
        enddo ! interpolated
        do j=1,ncy-1 ! from parent
          csP=csAv_west(j) ; cs=csA(1,j)+csA(1,j+1) ! to child grid.
          snP=snAv_west(j) ; sn=snA(1,j)+snA(1,j+1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_west(j)=cff*(cs*csP+sn*snP)
          snAv_west(j)=cff*(sn*csP-cs*snP) ! Convert them into
        enddo ! cos(child-parent)
      endif ! and
                                                    ! sin(child-parent)
      if (OBC_EAST) then
        do j=1,ncy
          csP=csAu_east(j) ; cs=csA(ncx,j)+csA(ncx-1,j)
          snP=snAu_east(j) ; sn=snA(ncx,j)+snA(ncx-1,j)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_east(j)=cff*(cs*csP+sn*snP)
          snAu_east(j)=cff*(sn*csP-cs*snP)
        enddo
        do j=1,ncy-1
          csP=csAv_east(j) ; cs=csA(ncx,j)+csA(ncx,j+1)
          snP=snAv_east(j) ; sn=snA(ncx,j)+snA(ncx,j+1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_east(j)=cff*(cs*csP+sn*snP)
          snAv_east(j)=cff*(sn*csP-cs*snP)
        enddo
      endif

      if (OBC_SOUTH) then
        do i=1,ncx
          csP=csAv_south(i) ; cs=csA(i,1)+csA(i,2)
          snP=snAv_south(i) ; sn=snA(i,1)+snA(i,2)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_south(i)=cff*(cs*csP+sn*snP)
          snAv_south(i)=cff*(sn*csP-cs*snP)
        enddo
        do i=1,ncx-1
          csP=csAu_south(i) ; cs=csA(i,1)+csA(i+1,1)
          snP=snAu_south(i) ; sn=snA(i,1)+snA(i+1,1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_south(i)=cff*(cs*csP+sn*snP)
          snAu_south(i)=cff*(sn*csP-cs*snP)
        enddo
      endif

      if (OBC_NORTH) then
        do i=1,ncx
          csP=csAv_north(i) ; cs=csA(i,ncy)+csA(i,ncy-1)
          snP=snAv_north(i) ; sn=snA(i,ncy)+snA(i,ncy-1)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAv_north(i)=cff*(cs*csP+sn*snP)
          snAv_north(i)=cff*(sn*csP-cs*snP)
        enddo
        do i=1,ncx-1
          csP=csAu_north(i) ; cs=csA(i,ncy)+csA(i+1,ncy)
          snP=snAu_north(i) ; sn=snA(i,ncy)+snA(i+1,ncy)
          cff=1.D0/sqrt((cs*cs+sn*sn)*(csP*csP+snP*snP))
          csAu_north(i)=cff*(cs*csP+sn*snP)
          snAu_north(i)=cff*(sn*csP-cs*snP)
        enddo
      endif

! Read parent-grid land mask. Similarly to the above: temporarily
! use array "hp" to read "mask_rho" from the file, thereafter this
! program uses only integer(kind=2) version of mask.

      ierr=nf_inq_varid(ncpgrd, 'mask_rho', varid)
      if (ierr == nf_noerr) then
        call get_patch_by_name_double(ncpgrd, prntgrd, 'mask_rho',
     & iwestpg,jsouthpg, nx,ny,0,0, hp)
        call set_mask(nx,ny, hp, mskp,umsp,vmsp)
      else
        mskp=1 ; umsp=1 ; vmsp=1
        write(*,'(9x,4A)') 'No land mask ''mask_rho'' is present ',
     & 'in ''', prntgrd(1:lpgrd), ''', assuming mask=1 everywhere.'
      endif

! Read parent-grid topography and interpolate it to child grid, then
! initialize vertical interpolation for each line along each boundary
! (there are 4 boundaries, each having 3 lines: RHO-, U-, and V-points
! all of which have distinct type of vertical interpolations, so the
! ultimate goal is to initialize 12 arrays, kpX_size, where X={r,y,v}
! and _side={west,east,north,south}). Note that it is possible that
! some portion of child grid is outside the parent, so in that areas
! horizontal interpolation results to zero values. Even though in any
! meaningful model configuration these areas are under land mask, and
! therefore make no effect, leaving zero-valued "hp" there confuses
! vertical search algorithm inside "bry_init_vertinterp" resulting
! in out-of bound array index. For performance reason it is better
! not to add extra protective logic there, but rather fill-in these
! areas with child-grid topography.

      call get_var_by_name_double(nccgrd, 'h', h)
      call get_patch_by_name_double(ncpgrd, prntgrd, 'h',
     & iwestpg,jsouthpg, nx,ny,0,0, hp)

      if ( minval(h) <=0 ) then ! To prevent errors from
          hmin=minval(h) ! wet and drying, move
          hp=hp-hmin+0.2 ! topography so it is positive
          h=h-hmin+0.2 ! everywhere
      endif

      call spln2d_double(nx,ny, hp, srX,srY,sXY,sYX)

      if (OBC_WEST) then
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncy,1,
     & ir_west,jr_west, xi_west,eta_west, hp_west)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncy,1,
     & iur_west,jur_west, xiur_west,etaur_west, hpu_west)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncy-1,1,
     & ivr_west,jvr_west, xivr_west,etavr_west, hpv_west)
        do j=1,ncy
          h_west(j)=h(1,j)
          if (hp_west(j) < 0.0001D0) hp_west(j)=h_west(j)
          hu_west(j)=0.5D0*(h(1,j)+h(2,j))
          if (hpu_west(j) < 0.0001D0) hpu_west(j)=hu_west(j)
        enddo
        do j=1,ncy-1
          hv_west(j)=0.5D0*(h(1,j)+h(1,j+1))
          if (hpv_west(j) < 0.0001D0) hpv_west(j)=hv_west(j)
        enddo
        call bry_init_vertinterp(ncy, hp_west, Np,hcp,Csp_r,
     & h_west, N,hc,Cs_r, kpr_west)
        call bry_init_vertinterp(ncy, hpu_west, Np,hcp,Csp_r,
     & hu_west, N,hc,Cs_r, kpu_west)
        call bry_init_vertinterp(ncy-1, hpv_west, Np,hcp,Csp_r,
     & hv_west, N,hc,Cs_r, kpv_west)

        if ( hmin <=0 ) then
          hp_west=hp_west+hmin-0.2
          hu_west=hp_west+hmin-0.2
          hv_west=hp_west+hmin-0.2
        endif

      endif

      if (OBC_EAST) then
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncy,1,
     & ir_east,jr_east, xi_east,eta_east, hp_east)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncy,1,
     & iur_east,jur_east, xiur_east,etaur_east, hpu_east)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncy-1,1,
     & ivr_east,jvr_east, xivr_east,etavr_east, hpv_east)
        do j=1,ncy
          h_east(j)=h(ncx,j)
          hu_east(j)=0.5D0*( h(ncx,j)+h(ncx-1,j))
          if (hp_east(j) < 0.0001D0) hp_east(j)=h_east(j)
          hu_east(j)=0.5D0*( h(ncx,j)+h(ncx-1,j))
          if (hpu_east(j) < 0.0001D0) hpu_east(j)=hu_east(j)
        enddo
        do j=1,ncy-1
          hv_east(j)=0.5D0*(h(ncx,j)+h(ncx,j+1))
          if (hpv_east(j) < 0.0001D0) hpv_east(j)=hv_east(j)
        enddo
        call bry_init_vertinterp(ncy, hp_east, Np,hcp,Csp_r,
     & h_east, N,hc,Cs_r, kpr_east)
        call bry_init_vertinterp(ncy, hpu_east, Np,hcp,Csp_r,
     & hu_east, N,hc,Cs_r, kpu_east)
        call bry_init_vertinterp(ncy-1, hpv_east, Np,hcp,Csp_r,
     & hv_east, N,hc,Cs_r, kpv_east)
        if ( hmin <=0 ) then
          hp_east=hp_east+hmin-0.2
          hu_east=hp_east+hmin-0.2
          hv_east=hp_east+hmin-0.2
        endif

      endif

      if (OBC_SOUTH) then
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncx,1,
     & ir_south,jr_south, xi_south,eta_south, hp_south)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncx-1,1,
     & iur_south,jur_south, xiur_south,etaur_south, hpu_south)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncx,1,
     & ivr_south,jvr_south, xivr_south,etavr_south, hpv_south)
        do i=1,ncx
          h_south(i)=h(i,1)
          if (hp_south(i) < 0.0001D0) hp_south(i)=h_south(i)
          hv_south(i)=0.5D0*(h(i,1)+h(i,2))
          if (hpv_south(i) < 0.0001D0) hpv_south(i)=hv_south(i)
        enddo
        do i=1,ncx-1
          hu_south(i)=0.5D0*(h(i,1)+h(i+1,1))
          if (hpu_south(i) < 0.0001D0) hpu_south(i)=hu_south(i)
        enddo
        call bry_init_vertinterp(ncx, hp_south, Np,hcp,Csp_r,
     & h_south, N,hc,Cs_r, kpr_south)
        call bry_init_vertinterp(ncx-1, hpu_south, Np,hcp,Csp_r,
     & hu_south, N,hc,Cs_r, kpu_south)
        call bry_init_vertinterp(ncx, hpv_south, Np,hcp,Csp_r,
     & hv_south, N,hc,Cs_r, kpv_south)

        if ( hmin <=0 ) then
          hp_south=hp_south+hmin-0.2
          hu_south=hp_south+hmin-0.2
          hv_south=hp_south+hmin-0.2
        endif
      endif

      if (OBC_NORTH) then
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncx,1,
     & ir_north,jr_north, xi_north,eta_north, hp_north)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncx-1,1,
     & iur_north,jur_north, xiur_north,etaur_north, hpu_north)
        call spln2d_interp_double(nx,ny, hp, srX,srY,sXY, ncx,1,
     & ivr_north,jvr_north, xivr_north,etavr_north, hpv_north)
        do i=1,ncx
          h_north(i)=h(i,ncy)
          if (hp_north(i) < 0.0001D0) hp_north(i)=h_north(i)
          hv_north(i)=0.5D0*(h(i,ncy) +h(i,ncy-1))
          if (hpv_north(i) < 0.0001D0) hpv_north(i)=hv_north(i)
        enddo
        do i=1,ncx-1
          hu_north(i)=0.5D0*(h(i,ncy)+h(i+1,ncy))
          if (hpu_north(i) < 0.0001D0) hpu_north(i)=hu_north(i)
        enddo
        call bry_init_vertinterp(ncx, hp_north, Np,hcp,Csp_r,
     & h_north, N,hc,Cs_r, kpr_north)
        call bry_init_vertinterp(ncx-1, hpu_north, Np,hcp,Csp_r,
     & hu_north, N,hc,Cs_r, kpu_north)
        call bry_init_vertinterp(ncx, hpv_north, Np,hcp,Csp_r,
     & hv_north, N,hc,Cs_r, kpv_north)
        if ( hmin <=0 ) then
          hp_north=hp_north+hmin-0.2
          hu_north=hp_north+hmin-0.2
          hv_north=hp_north+hmin-0.2
        endif

      endif

      write(*,'(/1x,A)') 'Checking...'
      if (OBC_WEST) then
        write(*,'(2x,A)') 'western'
        call bry_check_init_vertinterp(ncy, hp_west, Np,hcp,Csp_r,
     & h_west, N,hc,Cs_r, kpr_west)
        call bry_check_init_vertinterp(ncy, hpu_west, Np,hcp,Csp_r,
     & hu_west, N,hc,Cs_r, kpu_west)
        call bry_check_init_vertinterp(ncy-1, hpv_west, Np,hcp,Csp_r,
     & hv_west, N,hc,Cs_r, kpv_west)
      endif
      if (OBC_EAST) then
        write(*,'(2x,A)') 'eastern'
        call bry_check_init_vertinterp(ncy, hp_east, Np,hcp,Csp_r,
     & h_east, N,hc,Cs_r, kpr_east)
        call bry_check_init_vertinterp(ncy, hpu_east, Np,hcp,Csp_r,
     & hu_east, N,hc,Cs_r, kpu_east)
        call bry_check_init_vertinterp(ncy-1, hpv_east, Np,hcp,Csp_r,
     & hv_east, N,hc,Cs_r, kpv_east)
      endif
      if (OBC_SOUTH) then
        write(*,'(2x,A)') 'southern'
        call bry_check_init_vertinterp(ncx, hp_south, Np,hcp,Csp_r,
     & h_south, N,hc,Cs_r, kpr_south)
        call bry_check_init_vertinterp(ncx-1, hpu_south, Np,hcp,Csp_r,
     & hu_south, N,hc,Cs_r, kpu_south)
        call bry_check_init_vertinterp(ncx, hpv_south, Np,hcp,Csp_r,
     & hv_south, N,hc,Cs_r, kpv_south)
      endif
      if (OBC_NORTH) then
        write(*,'(2x,A)') 'northern'
        call bry_check_init_vertinterp(ncx, hp_north, Np,hcp,Csp_r,
     & h_north, N,hc,Cs_r, kpr_north)
        call bry_check_init_vertinterp(ncx-1, hpu_north, Np,hcp,Csp_r,
     & hu_north, N,hc,Cs_r, kpu_north)
        call bry_check_init_vertinterp(ncx, hpv_north, Np,hcp,Csp_r,
     & hv_north, N,hc,Cs_r, kpv_north)
      endif


      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      write(*,'(/F8.2,1x,A,F8.2,1x,A/)') dble(net_clk)/dble(clk_rate),
     & 'vertical initialization complete in',
     & dble(inc_clk)/dble(clk_rate), 'sec'


! Create boundary forcing file.

      ierr=nf_create(roms_bry, nf_netcdf4, nctarg)
      if (ierr == nf_noerr) then
        call def_bry_file(roms_bry, nctarg, ncx,ncy,N,
     & OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH,
     & theta_s, theta_b, hc, Cs_w,Cs_r,ntrc,tracer,ncsrc)
        ierr=nf_inq_varid(nctarg, 'bry_time', tvar_out)
        if (ierr /= nf_noerr) then
          write(*,*) '### ERROR 2: nf_inq_varid(. tvar_out .)'
        endif

! Copy all attributes for time variable while the target file is still
! in redefinition mode, then switch it into input and copy time itself.

        ierr=nf_inq_varid(nctarg, 'bry_time', tvar_out)
        if (ierr == nf_noerr) then
          ierr=nf_inq_varid(ncsrc, time_var_name, tvar_in)
          if (ierr == nf_noerr) then
            ierr=nf_inq_varnatts(ncsrc, tvar_in, natts)
            if (ierr == nf_noerr) then
              do i=1,natts
                ierr=nf_inq_attname(ncsrc, tvar_in, i, str)
                if (ierr == nf_noerr) then
                  call lenstr(str,lstr)
! If units change 'seconds since...' to 'days since...'
                  if (str(1:lstr) .eq. 'units') then
                    ierr=nf_get_att_text(ncsrc,tvar_in, 'units',str)
                    if (ierr == nf_noerr) then
                      call lenstr(str,lstr)
                      write(orig_date,'(2A)') 'days',str(8:lstr)
                      call lenstr(orig_date,lstr)
                      ierr= nf_put_att_text(nctarg, tvar_out,'units',
     & lstr,orig_date(1:lstr))
                    else
                      write(*,'(/1x,7A/12x,A/)') '### ERROR: ',
     & 'Found, but cannot read attribute ''units'' ',
     & 'for variable in ''',
     & prt_tmp(1:lprnt), '''.',
     & nf_strerror(ierr)
                      stop
                    endif
                  else
                    ierr=nf_copy_att(ncsrc, tvar_in, str(1:lstr),
     & nctarg, tvar_out)
                  endif
                  if (ierr == nf_noerr) then
                    write(*,*) 'copied attribute ''',str(1:lstr),''''
                  else
                    write(*,*) '### ERROR 10: copy_att' ; stop
                  endif
                else
                  write(*,*) '### ERROR 9: inq_attname' ; stop
                endif
              enddo
            else
              write(*,*) '### ERROR 8: inq_varnatts'
            endif
          else
            write(*,*) '### ERROR 7: inq_varid'
          endif
        else
          write(*,*) '### ERROR 6: inq_varid'
        endif
        if (ierr /= nf_noerr) stop

      else
        write(*,*) '### ERROR 1: nf_create'
      endif
      if (ierr /= nf_noerr) stop
      ierr=nf_enddef(nctarg) !<-- set to input mode
      ierr=nf_sync(nctarg)


! Finally clean up memory which is no longer needed and allocate
! big buffer arrays to read the actual data.


      deallocate(hp, xp,yp, xpu,ypu, xpv,ypv, xc,yc)
      net_alloc_mem=net_alloc_mem -(6*2*nx*ny+2*2*ncx*ncy)
      write(*,'(8x,A,F9.2,1x,A)') 'deallocated',
     & dble(6*2*nx*ny+2*2*ncx*ncy)/262144.D0,'MB'

      size=max(ncx,ncy)*max(Np,N)
      allocate(wrk3(size), wrk4(size), wrk5(size), wrk6(size),
     & wrk7(size), wrk8(size), stat=ierr)
      if (ierr == 0) then
        net_alloc_mem=net_alloc_mem + 6*size
        write(*,'(8x,A,F10.2,1x,A)') 'allocated wrk3,3.4, reaching',
     & dble(net_alloc_mem)/262144.D0, 'MB total'
        size=nx*ny*(Np+1)
        allocate(wrk1(size), wrk2(size), stat=ierr)
        if (ierr == 0) then
          net_alloc_mem=net_alloc_mem + 2*size
          write(*,'(8x,A,F10.2,1x,A)') 'allocated wrk1,2, reaching',
     & dble(net_alloc_mem)/262144.D0, 'MB total'
        endif
      endif
      write(*,'(2/7(/4x,A)/)')

     &'      *****    *********    ******   *******    *********  ',
     &'    ***   ***  *  ***  *   **  ***   ***   ***  *  ***  *  ',
     &'    ***           ***     **   ***   ***   ***     ***     ',
     &'      *****       ***    ***   ***   ***   **      ***     ',
     &'          ***     ***    *********   ******        ***     ',
     &'    ***   ***     ***    ***   ***   ***  **       ***     ',
     &'      *****       ***    ***   ***   ***   ***     ***     '


! At this moment the first history file is already open, while all
! others need to be open when their time comes. Determine the number
! of records in the file, netCDF ID and units for timing variable (time
! in seconds needs to be converted ito days).

      do ihis=1,nhists/ipt_char_len
        if (ihis > 1) then
          prt_tmp=trim(prnt_data(ihis))
          call lenstr(prt_tmp,lprnt)
          ierr=nf_open(prt_tmp(1:lprnt), nf_nowrite, ncsrc)
        endif
        ierr=nf_inq_varid(ncsrc, time_var_name, tvar_in)
        if (ierr == nf_noerr) then
          ierr=nf_inq_var(ncsrc, tvar_in, vname, vtype,
     & ndims, vdimids, natts)
          if (ierr == nf_noerr) then
            call lenstr(vname,lvar)
            if (ndims == 1) then
              ierr=nf_inq_dimlen(ncsrc, vdimids, nrecs)
              if (ierr == nf_noerr) then





                do i=1,natts
                  ierr=nf_inq_attname(ncsrc, tvar_in, i, str)
                  if (ierr == nf_noerr) then
                    call lenstr(str,lstr)
                    if (str(1:lstr) == 'units') then
                      ierr=nf_get_att_text(ncsrc,tvar_in, 'units',str)
                      if (ierr == nf_noerr) then
                        call lenstr(str,lstr) ; time_units=str(1:lstr)
                      else
                        write(*,'(/1x,7A/12x,A/)') '### ERROR: ',
     & 'Found, but cannot read attribute ''units'' ',
     & 'for variable ''', vname(1:lvar), ''' in ''',
     & prt_tmp(1:lprnt), '''.',
     & nf_strerror(ierr)
                        stop
                      endif
                    endif
                  else
                    write(*,'(/1x,2A,I3,1x,5A/12x,A/)') '### ERROR: ',
     & 'Cannot make inquiry for attribute #', i,
     & 'for variable ''', vname(1:lvar), ''' in ''',
     & prt_tmp(1:lprnt),'''.', nf_strerror(ierr)
                    stop
                  endif
                enddo
              else
                write(*,'(1x,4A/12x,A/)') '### ERROR: Cannot ',
     & 'determine number of time records in ''',
     & prt_tmp(1:lprnt),'''.', nf_strerror(ierr)
              endif
            elseif (ndims == 0) then
              write(*,'(2(/1x,3A)/)') 'Time variable ''',
     & vname(1:lvar), ''' does not have dimensions.',
     & 'Presuming that ''', prt_tmp(1:lprnt),
     & ''' is a single-record file.'
              nrecs=1
            else
              write(*,'(/1x,5A/)') '### ERROR: Time variable ''',
     & vname(1:lvar), ''' in ''', prt_tmp(1:lprnt),
     & ''' has more than one dimension.'
            endif
          else
            write(*,'(/1x,6A/12x,A/)') '### ERROR: Cannot make ',
     & 'inquiry for variable ''', vname(1:lvar), ''' in ''',
     & prt_tmp(1:lprnt), '''.', nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,6A/12x,A/)') '### ERROR: Cannot find ',
     & 'variable ''', vname(1:lvar), ''' in ''',
     & prt_tmp(1:lprnt), '''.', nf_strerror(ierr)
        endif
        if (ierr == nf_noerr) then
          call lenstr(time_units,lstr)
          if (lstr > 0) then
            write(*,'(/1x,A,I5,1x,6A)') 'Found', nrecs, 'time ',
     & 'records in ''', prt_tmp(1:lprnt), ''', ',
     & 'time units = ', time_units(1:lstr)
          else
            write(*,'(/1x,5A/)') '### ERROR: Time variable ''',
     & vname(1:lvar), ''' in ''', prt_tmp(1:lprnt),
     & ''' does not have attribute ''units''.'
            stop
          endif
        else
          stop
        endif

        do rec=1,nrecs
          recout=recout+1

          nclk=3-nclk
          call system_clock(iclk(nclk), clk_rate,clk_max)
          inc_clk=iclk(nclk)-iclk(3-nclk)
          if (inc_clk < 0) inc_clk=inc_clk+clk_max
          net_clk=net_clk+inc_clk
          write(*,'(/F10.2,1x,A,I6/)') dble(net_clk)/dble(clk_rate),
     & 'sec  Processing record ', recout


! Time
          ierr=nf_get_vara_double(ncsrc, tvar_in, rec,1, time)
          if (ierr == nf_noerr) then
            if (time_units(1:6) == 'second') time=time/86400.D0
            ierr=nf_put_vara_double(nctarg, tvar_out, recout,1, time)
            if (ierr == nf_noerr) then
              write(*,*) '    bry_time =', time, ' days'
            else
              write(*,'(/1x,3A/)') '### ERROR: Cannot write ''',
     & 'bry_time'', ', nf_strerror(ierr)
            endif
          else
            write(*,'(/1x,5A/)') '### ERROR: Cannot read ''',
     & time_var_name,''' from ''',
     & prt_tmp(1:lprnt), ''' ,', nf_strerror(ierr)
          endif

! Free surface ...

          call get_patch_by_name_real(ncsrc, prt_tmp,'zeta',
     & iwestpg,jsouthpg, nx,ny,0,rec, wrk1)
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1)
          call etch_into_land_thread(nx,ny, mskp, wrk1)
C$OMP END PARALLEL
          call spln2d_real(nx,ny, wrk1, srX,srY,sXY,sYX)

          if (OBC_WEST) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncy,1,
     & ir_west,jr_west, xi_west,eta_west, wrk2)
            call bry_apply_mask(ncy,1, mskr_west, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_west',
     & ncy,0,0,recout, wrk2)
          endif
          if (OBC_EAST) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncy,1,
     & ir_east,jr_east, xi_east,eta_east, wrk2)
            call bry_apply_mask(ncy,1, mskr_east, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_east',
     & ncy,0,0,recout, wrk2)
          endif
          if (OBC_SOUTH) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncx,1,
     & ir_south,jr_south, xi_south,eta_south, wrk2)
            call bry_apply_mask(ncx,1, mskr_south, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_south',
     & ncx,0,0,recout, wrk2)
          endif
          if (OBC_NORTH) then
            call spln2d_interp_real(nx,ny, wrk1, srX,srY,sXY, ncx,1,
     & ir_north,jr_north, xi_north,eta_north, wrk2)
            call bry_apply_mask(ncx,1, mskr_north, wrk2)
            call put_rec_by_name_real(nctarg, roms_bry, 'zeta_north',
     & ncx,0,0,recout, wrk2)
          endif

! Horizontal velocities, both u,v and barotropic...

          call get_patch_by_name_real(ncsrc, prt_tmp, 'u',
     & iwestpg,jsouthpg, nx-1,ny,Np,rec, wrk1)
          call get_patch_by_name_real(ncsrc, prt_tmp, 'v',
     & iwestpg,jsouthpg, nx,ny-1,Np,rec, wrk2)
          do k=1,Np
            isrc=1+(k-1)*(nx-1)*ny
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1, isrc)
            call etch_into_land_thread(nx-1,ny, umsp, wrk1(isrc))
C$OMP END PARALLEL
          enddo
          do k=1,Np
            isrc=1+(k-1)*nx*(ny-1)
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1, isrc)
            call etch_into_land_thread(nx,ny-1, vmsp, wrk2(isrc))
C$OMP END PARALLEL
          enddo

          if (OBC_WEST) then
            call bry_interp(nx-1,ny,Np, wrk1, ncy, iu_west,ju_west,
     & xiu_west,etau_west, msku_west, wrk3)
            call bry_interp(nx,ny-1,Np, wrk2, ncy, iuv_west,juv_west,
     & xiuv_west,etauv_west, msku_west, wrk4)
            call bry_rotate_u_in_place(ncy,Np, csAu_west,snAu_west,
     & wrk3,wrk4)
            call bry_vertinterp(ncy, btm_slp, Np, wrk3, N, kpu_west,
     & wrk5)
            call bry_vert_average(ncy, hu_west, N,hc,Cs_w, wrk5,wrk7)
            call bry_apply_mask(ncy,N, msku_west, wrk5)
            call bry_apply_mask(ncy,1, msku_west, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_west',
     & ncy,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry, 'ubar_west',
     & ncy,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np, wrk2, ncy-1, iv_west,jv_west,
     & xiv_west,etav_west, mskv_west, wrk4)
            call bry_interp(nx-1,ny,Np, wrk1, ncy-1, ivu_west,jvu_west,
     & xivu_west,etavu_west, mskv_west, wrk3)
            call bry_rotate_v_in_place(ncy-1,Np, csAv_west,snAv_west,
     & wrk3,wrk4)
            call bry_vertinterp(ncy-1,btm_slp, Np, wrk4, N, kpv_west,
     & wrk6)
            call bry_vert_average(ncy-1,hv_west, N,hc,Cs_w, wrk6,wrk8)
            call bry_apply_mask(ncy-1,N, mskv_west, wrk6)
            call bry_apply_mask(ncy-1,1, mskv_west, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_west',
     & ncy-1,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry, 'vbar_west',
     & ncy-1,0,0,recout, wrk8)
          endif

          if (OBC_EAST) then
            call bry_interp(nx-1,ny,Np, wrk1, ncy, iu_east,ju_east,
     & xiu_east,etau_east, msku_east, wrk3)
            call bry_interp(nx,ny-1,Np, wrk2, ncy, iuv_east,juv_east,
     & xiuv_east,etauv_east, msku_east, wrk4)
            call bry_rotate_u_in_place(ncy,Np, csAu_east,snAu_east,
     & wrk3,wrk4)
            call bry_vertinterp(ncy, btm_slp, Np,wrk3, N,kpu_east,
     & wrk5)
            call bry_vert_average(ncy, hu_east, N,hc,Cs_w, wrk5,wrk7)
            call bry_apply_mask(ncy,N, msku_east, wrk5)
            call bry_apply_mask(ncy,1, msku_east, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_east',
     & ncy,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry, 'ubar_east',
     & ncy,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np,wrk2, ncy-1, iv_east,jv_east,
     & xiv_east,etav_east, mskv_east, wrk4)
            call bry_interp(nx-1,ny,Np,wrk1, ncy-1, ivu_east,jvu_east,
     & xivu_east,etavu_east, mskv_east, wrk3)
            call bry_rotate_v_in_place(ncy-1,Np, csAv_east,snAv_east,
     & wrk3,wrk4)
            call bry_vertinterp(ncy-1, btm_slp, Np,wrk4, N,kpv_east,
     & wrk6)
            call bry_vert_average(ncy-1,hv_east, N,hc,Cs_w,wrk6,wrk8)
            call bry_apply_mask(ncy-1,N, mskv_east, wrk6)
            call bry_apply_mask(ncy-1,1, mskv_east, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_east',
     & ncy-1,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry, 'vbar_east',
     & ncy-1,0,0,recout, wrk8)
          endif

          if (OBC_SOUTH) then
            call bry_interp(nx-1,ny,Np,wrk1, ncx-1,iu_south,ju_south,
     & xiu_south,etau_south, msku_south, wrk3)
            call bry_interp(nx,ny-1,Np,wrk2, ncx-1,iuv_south,juv_south,
     & xiuv_south,etauv_south, msku_south, wrk4)
            call bry_rotate_u_in_place(ncx-1,Np, csAu_south,snAu_south,
     & wrk3,wrk4)
            call bry_vertinterp(ncx-1,btm_slp, Np,wrk3, N,kpu_south,
     & wrk5)
            call bry_vert_average(ncx-1,hu_south, N,hc,Cs_w,wrk5,wrk7)
            call bry_apply_mask(ncx-1,N, msku_south, wrk5)
            call bry_apply_mask(ncx-1,1, msku_south, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_south',
     & ncx-1,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry,'ubar_south',
     & ncx-1,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np, wrk2, ncx, iv_south,jv_south,
     & xiv_south,etav_south, mskv_south, wrk4)
            call bry_interp(nx-1,ny,Np, wrk1, ncx, ivu_south,jvu_south,
     & xivu_south,etavu_south, mskv_south, wrk3)
            call bry_rotate_v_in_place(ncx,Np, csAv_south,snAv_south,
     & wrk3,wrk4)
            call bry_vertinterp(ncx, btm_slp, Np,wrk4, N,kpv_south,
     & wrk6)
            call bry_vert_average(ncx, hv_south, N,hc,Cs_w, wrk6,wrk8)
            call bry_apply_mask(ncx,N, mskv_south, wrk6)
            call bry_apply_mask(ncx,1, mskv_south, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_south',
     & ncx,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry,'vbar_south',
     & ncx,0,0,recout, wrk8)
          endif

          if (OBC_NORTH) then
            call bry_interp(nx-1,ny,Np,wrk1, ncx-1, iu_north,ju_north,
     & xiu_north,etau_north, msku_north, wrk3)
            call bry_interp(nx,ny-1,Np,wrk2, ncx-1,iuv_north,juv_north,
     & xiuv_north,etauv_north, msku_north, wrk4)
            call bry_rotate_u_in_place(ncx-1,Np, csAu_north,snAu_north,
     & wrk3,wrk4)
            call bry_vertinterp(ncx-1, btm_slp, Np, wrk3, N,kpu_north,
     & wrk5)
            call bry_vert_average(ncx-1,hu_north, N,hc,Cs_w,wrk5,wrk7)
            call bry_apply_mask(ncx-1,N, msku_north, wrk5)
            call bry_apply_mask(ncx-1,1, msku_north, wrk7)
            call put_rec_by_name_real(nctarg, roms_bry, 'u_north',
     & ncx-1,N,0,recout, wrk5)
            call put_rec_by_name_real(nctarg, roms_bry, 'ubar_north',
     & ncx-1,0,0,recout, wrk7)


            call bry_interp(nx,ny-1,Np, wrk2, ncx, iv_north,jv_north,
     & xiv_north,etav_north, mskv_north, wrk4)
            call bry_interp(nx-1,ny,Np, wrk1, ncx, ivu_north,jvu_north,
     & xivu_north,etavu_north, mskv_north, wrk3)
            call bry_rotate_v_in_place(ncx,Np, csAv_north,snAv_north,
     & wrk3,wrk4)
            call bry_vertinterp(ncx, btm_slp, Np,wrk4, N,kpv_north,
     & wrk6)
            call bry_vert_average(ncx, hv_north, N,hc,Cs_w, wrk6,wrk8)
            call bry_apply_mask(ncx,N, mskv_north, wrk6)
            call bry_apply_mask(ncx,1, mskv_north, wrk8)
            call put_rec_by_name_real(nctarg, roms_bry, 'v_north',
     & ncx,N,0,recout, wrk6)
            call put_rec_by_name_real(nctarg, roms_bry,'vbar_north',
     & ncx,0,0,recout, wrk8)
          endif

! Tracers...

          do itrc=1,ntrc/ipt_trc_len
            vname=trim(tracer(itrc)) ; call lenstr(vname,lvar)
            call get_patch_by_name_real(ncsrc, prt_tmp,vname,
     & iwestpg,jsouthpg, nx,ny,Np,rec, wrk1)
            do k=1,Np
              isrc=1+(k-1)*nx*ny
C$OMP PARALLEL SHARED(nx,ny, mskp, wrk1, isrc)
              call etch_into_land_thread(nx,ny, mskp, wrk1(isrc))
C$OMP END PARALLEL
            enddo

            if (OBC_WEST) then
              trgname=vname(1:lvar)/ /'_west'; call lenstr(trgname,ltgv)
              call bry_interp(nx,ny,Np, wrk1, ncy, ir_west,jr_west,
     & xi_west,eta_west, mskr_west, wrk3)
              call bry_vertinterp(ncy, btm_trc, Np,wrk3, N,kpr_west,
     & wrk5)
              call bry_apply_mask(ncy,N, mskr_west, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     & 0,ncy,N,recout, wrk5)
            endif

            if (OBC_EAST) then
              trgname=vname(1:lvar)/ /'_east'; call lenstr(trgname,ltgv)
              call bry_interp(nx,ny,Np, wrk1, ncy, ir_east,jr_east,
     & xi_east,eta_east, mskr_east, wrk3)
              call bry_vertinterp(ncy, btm_trc, Np,wrk3, N,kpr_east,
     & wrk5)
              call bry_apply_mask(ncy,N, mskr_east, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     & 0,ncy,N,recout, wrk5)
            endif

            if (OBC_SOUTH) then
              trgname=vname(1:lvar)/ /'_south';call lenstr(trgname,ltgv)
              call bry_interp(nx,ny,Np, wrk1, ncx, ir_south,jr_south,
     & xi_south,eta_south, mskr_south, wrk3)
              call bry_vertinterp(ncx, btm_trc, Np,wrk3, N,kpr_south,
     & wrk5)
              call bry_apply_mask(ncx,N, mskr_south, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     & ncx,0,N,recout, wrk5)
            endif

            if (OBC_NORTH) then
              trgname=vname(1:lvar)/ /'_north';call lenstr(trgname,ltgv)
              call bry_interp(nx,ny,Np, wrk1, ncx, ir_north,jr_north,
     & xi_north,eta_north, mskr_north, wrk3)
              call bry_vertinterp(ncx, btm_trc, Np,wrk3, N,kpr_north,
     & wrk5)
              call bry_apply_mask(ncx,N, mskr_north, wrk5)
              call put_rec_by_name_real(nctarg, roms_bry, trgname,
     & ncx,0,N,recout, wrk5)
            endif
          enddo !<-- itrc ! forcefully sync target file once
        enddo !<-- rec ! in a while to make it readable if
        ierr=nf_close(ncsrc)
        if (mod(recout,32) == 0) ierr=nf_sync(nctarg)
      enddo !<-- ihis
      ierr=nf_close(nctarg) ! the program is interrupted.

      write(*,'(/4x,A,F12.3,1x,A,F10.1,1x,A)') 'total data read ',
     & dble(sz_read_acc)/dble(1024**2), 'MBytes in',
     & dble(read_clk)/dble(clk_rate), 'sec'
      write(*,'(1x,A,F12.3,1x,A,F10.1,1x,A/)') 'total data written ',
     & dble(sz_write_acc)/dble(1024**2), 'MBytes in',
     & dble(write_clk)/dble(clk_rate), 'sec'
      nclk=3-nclk
      call system_clock(iclk(nclk), clk_rate,clk_max)
      inc_clk=iclk(nclk)-iclk(3-nclk)
      if (inc_clk < 0) inc_clk=inc_clk+clk_max
      net_clk=net_clk+inc_clk
      call cpu_time(tend)
      write(*,'(/1x,A,F8.2,1x,A,4x,A,F8.2,1x,A,F8.1,1x,A)')
     & 'Wall Clock time:', dble(net_clk)/dble(clk_rate), 'sec',
     & 'CPU time:', tend-tstart, 'sec',
     & (tend-tstart)*dble(clk_rate)/dble(net_clk)*100.D0,'% CPU'

      end subroutine r2r_bry
# 62 "tools_fort.F" 2
# 1 "etch_into_land.F" 1
! Content of this package is user-callable "etch_into_land" and/or
!-------- -- ---- ------- "etch_into_land_thread",

      subroutine etch_into_land(nx,ny, mask, qfld)
      implicit none
      integer nx,ny
      integer(kind=2) mask(nx,ny)
      real(kind=4) qfld(nx,ny)
C$OMP PARALLEL SHARED(nx,ny, mask, qfld)
      call etch_into_land_thread(nx,ny, mask, qfld)
C$OMP END PARALLEL
      end

! which progressively fill in land masked and/or special/missing value
! areas of a given 2D field "qfld" in order to subsequently interpolate
! it by an algorithm not designed to handle missing values. All other
! routines in this package are for its own internal use and not meant
! to be user callable. This package replaces "extend_on_land.F" by
! with the same functionality at a significantly reduced computational
! cost. While mathematically equivalent, this package uses radically
! different approach to code organization (and optimization) by
! replacing 2D-i,j-index sweeps with progressively reduction of a
! precomputed list of indices special-valued points:

! (1) make list of indices of all points with special-valued
! (2) select the ones from list (1) above which have immediate
! neighbors with valid values - these are the points which
! potentially can be filled in;
!
! (3) attempt to fill in all points from the list (2) - in order
! to be filled in a point must pass certain threshold of having
! enough valid neighbors - not all from list (2) are eventually
! be selected; compute values for those which will, an place
! these values into a special array (not directly into "qfld" at
! this moment) in order to avoid interference); mark positions
! of these finally selected points within list (1) by zeroing
! out the first index;
!
! (4) once step (3) is complete for the entire array and by all
! the threads, fill in selected values into "qfld";
!
! (5) shorten the list (1) by excluding zero-index points which
! were marked (3) and filled in by step (4);
!
! keep repeating steps (2)-(5) until nothing left in the list.
!

! The roles of the individual routines in the above algorithm are:
!---- -----
! call copy_extend_tile ! preliminary step
! call set_qext_bc_tile
! call init_ijmsk_tile !<-- list (1) of all sp.val. points
! do while(mskd_pts>0)
! call select_coast_pts !<-- list (2) as subset of (1)
! call comp_patch_pts !<-- list (3) as subset of (2)
! call apply_patch_pts !<-- (4)
! call set_qext_bc_tile
! call shortlist !<--(5)
! enddo
! call strip_halo_tile

! Meaning of the variables: [dimensional variables "nmsk", "ncst", and
!-------- -- --- ---------- "npths" indicate meaningful portions of the
! arrays as rather than allocated sizes;
! qext(0:nx+1,0:ny+1) working array with one row of ghost points all
! around, initially copy of "qfld" with Neumann b.c. applied;
! ijmsk(2,nmsk) list of i,j-indices of all special-valued points of
! array "qext", this is list (1) above;
! nijmsk(ncst) subset of indices within ijmsk(2,*) identifying points
! which can be potentially filled during the current
! iteration (list (2) above;
! ijptch(2:npths) list of i,j-indices of "qext" which are about to be
! filled; this is list (3) which is a subset nijmsk;
! ptch(npths) a set of fill values computed by step (3) and applied by
! step(4);

! Parallelization: Because of the mathematical nature of the problem,
!----------------- the standard tiling approach is would not be most
! efficient because the work done by each thread depends on the number
! of special valued points rather than total number of points within a
! subdomain, and therefore not predictable. This may cause significant
! load miss-balance. Note that the algorithm within each thread is
! fundamentally sequential as it involves counting to form index lists.
! So the approach below is to one-dimensionally divide all the domain
! into set of narrow stripes (parameter "jsize" below) and make them
! interleaved among the threads to, hopefully, make all of them get
! approximately the same amount of land. All the lists of indices are
! strictly PRIVATE in the code below, so each thread is responsible for
! processing only points within its own pre-determined set of stripes.
! Because the first list "ijmsk" is formed by a "_tile" (i.e. striped)
! routine, this establishes AFFINITY between the stripes and the
! subsequent going-by-the-list loops in sense that all the indices
! from the lists are guaranteed to be within the set of stripes of
! the thread. This eliminates the necessity to have barrier copying
! fill values into "qext" by a going-along-index-list routine and
! applying Neumann b.c. by tiled. Thus, the only operations needed
! to be barrier-separated are the ones where "qext" is input (forming
! lists and calculating fill values) and where "qext" is output
! (applying fill values into "qext").
# 109 "etch_into_land.F"
      module etch_on_land_vars
        real(kind=4), allocatable, dimension(:) :: qext
! integer, save :: allc_ext_size=0, trd_count=0, mskd_pts=0

! integer, save :: alloc_msk_size=0
C$OMP THREADPRIVATE(alloc_msk_size)
        real(kind=4), allocatable, dimension(:) :: ptch
C$OMP THREADPRIVATE(ptch)
        integer(kind=2), allocatable, dimension(:,:) :: ijmsk
C$OMP THREADPRIVATE(ijmsk)
        integer(kind=2), allocatable, dimension(:,:) :: ijptch
C$OMP THREADPRIVATE(ijptch)
        integer(kind=4), allocatable, dimension(:) :: nijmsk
C$OMP THREADPRIVATE(nijmsk)
      end module etch_on_land_vars




      subroutine etch_into_land_thread(nx,ny, mask, qfld)
      use etch_on_land_vars
      implicit none
      integer nx,ny
      integer(kind=2) mask(nx,ny)
      real(kind=4) qfld(nx,ny)
!>
      integer ntrds,trd, nmsk,nmsk_new, ncst, npths,
     & istr,iend,jstr,jend, j,j0, jskip, ierr
      integer, parameter :: jsize=3
      integer :: allc_ext_size, trd_count, mskd_pts
      integer :: alloc_msk_size

C$ integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$ ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      allc_ext_size=0; trd_count=0; mskd_pts=0;alloc_msk_size=0
      j0=trd*jsize ; jskip=ntrds*jsize ; istr=1 ; iend=nx


C$OMP MASTER
      mskd_pts=1!<-- set to start while loop
      if ((nx+2)*(ny+2) > allc_ext_size) then
        allc_ext_size=(nx+2)*(ny+2)
        if (allocated(qext)) deallocate(qext)
        allocate(qext(allc_ext_size))
! write(*,'(1x,2A,F16.8,1x,A)') 'etch_into_land_thread :: ',
! & 'allocated', dble((nx+2)*(ny+2))/dble(262144),
! & 'MB shared workspace array'
      endif
C$OMP END MASTER
C$OMP BARRIER

      nmsk=0
      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call copy_extend_tile(istr,iend,jstr,jend, nx,ny, mask,
     & qfld,qext, nmsk)
         call set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
      enddo
C$OMP BARRIER

      if (nmsk > alloc_msk_size) then
        alloc_msk_size=nmsk
        if (allocated(ptch)) then
            deallocate(ptch,ijptch, nijmsk,ijmsk)
        endif
        allocate( ijmsk(2,alloc_msk_size), nijmsk(alloc_msk_size),
     & ijptch(2,alloc_msk_size), ptch(alloc_msk_size),
     & stat=ierr )
! if (ierr == 0) then
C$OMP CRITICAL(etch_cr_rgn)
! write(*,'(1x,2A,F16.8,1x,A,I3)') 'etch_into_land_thread :: ',
! & 'allocated', dble(4*alloc_msk_size)/dble(262144),
! & 'MB private workspace, trd =', trd
C$OMP END CRITICAL(etch_cr_rgn)
! else
         if (ierr /= 0) then
          write(*,*) '### ERROR: etch_into_land_thread :: ',
     & 'memory allocation error.'
        endif
      endif


      nmsk=0 !<-- number of masked points
      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call init_ijmsk_tile(istr,iend,jstr,jend, nx,ny, qext,
     & alloc_msk_size, ijmsk, nmsk)
      enddo
C$OMP BARRIER

      do while (mskd_pts > 0)
        ncst=0 !<-- number of coastal points
        call select_coast_pts(nx,ny, qext, nmsk,ijmsk, ncst,nijmsk)

        npths=0 !<-- number of points to be patched
        call comp_patch_pts(nx,ny, qext, nmsk,ijmsk, ncst,nijmsk,
     & npths, ijptch,ptch)
C$OMP BARRIER
        call apply_patch_pts(npths, ijptch,ptch, nx,ny, qext)
        do j=j0,ny,jskip
          jstr=max(1,j) ; jend=min(j+jsize-1,ny)
          call set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
        enddo
        nmsk_new=0
        call shortlist(nmsk,ijmsk, nmsk_new)
        nmsk=nmsk_new
C$OMP CRITICAL(etch_cr_rgn)
        if (trd_count == 0) then
            mskd_pts=0
        endif
        trd_count=trd_count+1
        mskd_pts=mskd_pts+nmsk
        if (trd_count == ntrds) then
          trd_count=0



        endif
C$OMP END CRITICAL(etch_cr_rgn)
C$OMP BARRIER
      enddo !<-- while

      do j=j0,ny,jskip
        jstr=max(1,j) ; jend=min(j+jsize-1,ny)
        call strip_halo_tile(istr,iend,jstr,jend, nx,ny, qext,qfld)
      enddo
C$OMP BARRIER
      end

      subroutine copy_extend_tile(istr,iend,jstr,jend, nx,ny, mask,
     & qsrc,qext, nmsk)

! Copy array "qsrc" into "qext" with has one row of ghost points all
! around, while setting land points to special values. The logic here
! is designed to work both ways: either there a non-trivial land mask
! array, or masked data is already set to some special value, while
! mask(i,j) is trivially set to all-water mask(:,:)=1 status and makes
! no effect. The secondary goal is to determine the total number of
! special-value points encountered by this thread so an appropiate
! sized arrays can be allocated to hold list of indices.

      implicit none
      integer istr,iend,jstr,jend, nx,ny, nmsk, i,j
      integer(kind=2) mask(nx,ny)
      real(kind=4) qsrc(nx,ny), qext(0:nx+1,0:ny+1)
!> write(*,*) 'enter copy_extend_tile'
      do j=jstr,jend
        do i=istr,iend
          if (mask(i,j) > 0 .and. abs(qsrc(i,j)) < abs(-9.9E+9)) then
            qext(i,j)=qsrc(i,j)
          else
            qext(i,j)=-1.D+10 ; nmsk=nmsk+1
          endif
        enddo
      enddo
      end

      subroutine strip_halo_tile(istr,iend,jstr,jend, nx,ny, qext,qsrc)
      implicit none
      integer istr,iend,jstr,jend, nx,ny, i,j
      real(kind=4) qext(0:nx+1,0:ny+1), qsrc(nx,ny)
!> write(*,*) 'enter strip_halo_tile'
      do j=jstr,jend
        do i=istr,iend
          qsrc(i,j)=qext(i,j)
        enddo
      enddo
      end

      subroutine set_qext_bc_tile(istr,iend,jstr,jend, nx,ny, qext)
      implicit none
      integer istr,iend,jstr,jend, nx,ny, i,j
      real(kind=4) qext(0:nx+1,0:ny+1)
!> write(*,*) 'enter set_qext_bc_tile'
      if (istr==1) then
        do j=jstr,jend
          qext(istr-1,j)=qext(istr,j)
        enddo
      endif
      if (iend==nx) then
        do j=jstr,jend
          qext(iend+1,j)=qext(iend,j)
        enddo
      endif
      if (jstr==1) then
        do i=istr,iend
          qext(i,jstr-1)=qext(i,jstr)
        enddo
      endif
      if (jend==ny) then
        do i=istr,iend
          qext(i,jend+1)=qext(i,jend)
        enddo
      endif
      if (istr==1 .and. jstr==1) then
        qext(istr-1,jstr-1)=qext(istr,jstr)
      endif
      if (istr==1 .and. jend==ny) then
        qext(istr-1,jend+1)=qext(istr,jend)
      endif
      if (iend==nx .and. jstr==1) then
        qext(iend+1,jstr-1)=qext(iend,jstr)
      endif
      if (iend==nx .and. jend==ny) then
        qext(iend+1,jend+1)=qext(iend,jend)
      endif
      end

      subroutine init_ijmsk_tile(istr,iend,jstr,jend, nx,ny, src,
     & max_pts, ijmsk, nmsk)
      implicit none
      integer istr,iend,jstr,jend, nx,ny, max_pts, nmsk, i,j
      real(kind=4) src(0:nx+1,0:ny+1)
      integer(kind=2) ijmsk(2,max_pts)
!> write(*,*) 'enter init_ijmsk_tile'
      do j=jstr,jend ! Form list of indices of
        do i=istr,iend ! points with special value.
          if (src(i,j) < -9.9E+9) then
            nmsk=nmsk+1
            ijmsk(1,nmsk)=i
            ijmsk(2,nmsk)=j
          endif
        enddo
      enddo
      end

      subroutine select_coast_pts(nx,ny, src, nmsk,ijmsk, ncst,nijmsk)
      implicit none
      integer nx,ny, nmsk, ncst, i,j,n ! Take list of previously
      real(kind=4) src(0:nx+1,0:ny+1) ! identified masked points
      integer(kind=2) ijmsk(2,nmsk) ! and select the ones among
      integer(kind=4) nijmsk(nmsk) ! them which have at least
!> write(*,*) 'enter select_coast_pts' ! one water neighbor.
      do n=1,nmsk
        i=ijmsk(1,n) ; j=ijmsk(2,n)
        if (src(i,j) < -9.9E+9) then
          if (src(i+1,j) > -9.9E+9 .or. src(i,j+1) > -9.9E+9 .or.
     & src(i-1,j) > -9.9E+9 .or. src(i,j-1) > -9.9E+9 ) then
            ncst=ncst+1 ; nijmsk(ncst)=n
          endif
        endif
      enddo
      end

      subroutine comp_patch_pts(nx,ny, src, nmsk,ijmsk, ncst,nijmsk,
     & npths, ijptch, ptch)
      implicit none
      integer nx,ny, nmsk, ncst, npths, i,j,n
      real(kind=4) src(0:nx+1,0:ny+1), ptch(ncst), wgt,vlu
      integer(kind=2) ijmsk(2,nmsk), ijptch(2,ncst)
      integer(kind=4) nijmsk(nmsk)

      real(kind=4), parameter :: grad=1./3., corn=0.707106781186547,
     & corngrad=0.5*corn*grad, threshold=2.4

!> write(*,*) 'enter comp_patch_pts, ncst =',ncst

      do n=1,ncst ! check surrounding points:
        i=ijmsk(1,nijmsk(n)) ; wgt=0. ! counterclockwise direction
        j=ijmsk(2,nijmsk(n)) ; vlu=0. ! starting from the east.

        if (src(i+1,j) > -9.9E+9) then
          wgt=wgt+1. ; vlu=vlu + src(i+1,j)
          if (i < nx) then
            if (src(i+2,j) > -9.9E+9) then
                vlu=vlu+grad*(src(i+1,j)-src(i+2,j))
            endif
          endif
        endif


        if (src(i+1,j+1) > -9.9E+9) then
          wgt=wgt+corn ; vlu=vlu + corn*src(i+1,j+1)
          if (i < nx) then
            if (src(i+2,j+1) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i+1,j+1)-src(i+2,j+1))
            endif
          endif
          if (j < ny) then
            if (src(i+1,j+2) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i+1,j+1)-src(i+1,j+2))
            endif
          endif
        endif


        if (src(i,j+1) > -9.9E+9) then
          wgt=wgt+1. ; vlu=vlu + src(i,j+1)
          if (j < ny) then
            if (src(i,j+2) > -9.9E+9) then
                vlu=vlu+grad*(src(i,j+1)-src(i,j+2))
            endif
          endif
        endif


        if (src(i-1,j+1) > -9.9E+9) then
          wgt=wgt+corn ; vlu=vlu + corn*src(i-1,j+1)
          if (j < ny) then
            if (src(i-1,j+2) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i-1,j+1)-src(i-1,j+2))
            endif
          endif
          if (i > 1) then
            if (src(i-2,j+1) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i-1,j+1)-src(i-2,j+1))
            endif
          endif
        endif


        if (src(i-1,j) > -9.9E+9) then
          wgt=wgt+1. ; vlu=vlu + src(i-1,j)
          if (i > 1) then
            if (src(i-2,j) > -9.9E+9) then
                vlu=vlu+grad*(src(i-1,j)-src(i-2,j))
            endif
          endif
        endif


        if (src(i-1,j-1) > -9.9E+9) then
          wgt=wgt+corn ; vlu=vlu + corn*src(i-1,j-1)
          if (i > 1) then
            if (src(i-2,j-1) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i-1,j-1) -src(i-2,j-1))
            endif
          endif
          if (j > 1) then
            if (src(i-1,j-2) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i-1,j-1)-src(i-1,j-2))
            endif
          endif
        endif


        if (src(i,j-1) > -9.9E+9) then
          wgt=wgt+1. ; vlu=vlu + src(i,j-1)
          if (j > 1) then
            if (src(i,j-2) > -9.9E+9) then
                vlu=vlu+grad*(src(i,j-1)-src(i,j-2))
            endif
          endif
        endif


        if (src(i+1,j-1) > -9.9E+9) then
          wgt=wgt+corn ; vlu=vlu + corn*src(i+1,j-1)
          if (j > 1) then
            if (src(i+1,j-2) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i+1,j-1)-src(i+1,j-2))
            endif
          endif
          if (i < nx) then
            if (src(i+2,j-1) > -9.9E+9) then
                vlu=vlu + corngrad*(src(i+1,j-1)-src(i+2,j-1))
            endif
          endif
        endif

        if (wgt > threshold) then
          npths=npths+1 ! At the end set "ijmsk" i-index
          ptch(npths)=vlu/wgt ! to zero to signal that the point
          ijptch(1,npths)=i ! is no longer a special value.
          ijptch(2,npths)=j
          ijmsk(1,nijmsk(n))=0
        endif
      enddo
      end

      subroutine apply_patch_pts(npths, ijptch,ptch, nx,ny, src)
      implicit none
      integer npths, nx,ny, i,j,n
      integer(kind=2) ijptch(2,npths)
      real(kind=4) ptch(npths), src(0:nx+1,0:ny+1)
!> write(*,*) 'enter apply_patch_pts'
      do n=1,npths
        i=ijptch(1,n) ; j=ijptch(2,n)
        src(i,j)=ptch(n)
      enddo
      end


      subroutine shortlist(nmsk,ijmsk, nmsk_new)
      implicit none
      integer nmsk,nmsk_new,n
      integer(kind=2) ijmsk(2,nmsk)
      do n=1,nmsk
        if (ijmsk(1,n) > 0) then
          nmsk_new=nmsk_new+1
          ijmsk(1,nmsk_new)=ijmsk(1,n)
          ijmsk(2,nmsk_new)=ijmsk(2,n)
        endif
      enddo
!> write(*,*) 'shortlist, nmsk=', nmsk, ' -->', nmsk_new
      end
# 63 "tools_fort.F" 2
# 1 "set_depth.F" 1
! The content of this file is a collection of routines associated
! with vertical coordinate transformation, and sould matches exactly
! the actual transform in ROMS code.

c--#define VERBOSE


      subroutine set_scoord(theta_s,theta_b, N, Cs_r,Cs_w)
      implicit none
      integer N, k ! Setup S-coordinate system:
      real(kind=8) theta_s, theta_b, ! compute vertical stretching
     & Cs_w(0:N), Cs_r(N), ! curves Cs_w,Cs_r at W- and
     & ds, sc, CSF ! RHO-points.
      ds=1.D0/dble(N) ! input: theta_s, theta_b, N
      Cs_w(N)=0.D0 ! output: Cs_w, Cs_r
      do k=N-1,1,-1 ! -1 < Cs_r,Cs_w < 0
        sc=ds*dble(k-N)
        Cs_w(k)=CSF(sc, theta_s,theta_b)
      enddo
      Cs_w(0)=-1.D0
      do k=1,N
        sc=ds*(dble(k-N)-0.5D0)
        Cs_r(k)=CSF(sc, theta_s,theta_b)
      enddo
# 34 "set_depth.F"
      end

      subroutine print_scoord(N, Cs_w, hc,hmin,hmax)
      implicit none
      integer N, k
      real(kind=8) Cs_w(0:N), hc,hmin,hmax, ds, sc, z1,zhc,z2,z3
      ds=1.D0/dble(N)
      write(*,'(/1x,A/,/2x,A,7x,A/)')
     & 'Vertical S-coordinate system (z at W-points):',
     & 'level   S-coord    Cs-curve    Z at hmin',
     & 'at hc    half way     at hmax'
      do k=N,0,-1
        sc=ds*dble(k-N)

        z1=hmin*(hc*sc + hmin*Cs_w(k))/(hc+hmin)
        zhc=0.5*hc*(sc + Cs_w(k))
        z2=0.5*hmax*(hc*sc + 0.5*hmax*Cs_w(k))/(hc+0.5*hmax)
        z3=hmax*(hc*sc + hmax*Cs_w(k))/(hc+hmax)






        if (hc < 1.E+4) then
          write(*,'(I7,F11.6,F12.7,4F12.3)')
     & k, ds*(k-N),Cs_w(k), z1,zhc,z2,z3
        else
          write(*,'(I7,F11.6,F12.7,F12.3,12x,2F12.3)')
     & k, ds*(k-N),Cs_w(k), z1, z2,z3
        endif
      enddo
      end



! Vertical stretching functions: In principle stretching curves can be
! selected independently from the vertical transformation type, which
! would require a separate name for CPP-switch here. This however has
! little practical value (merely to reproduce my own Pacific solutions
! of 2005 KPP presentation where the SM09 coordinate transformation
! formula was used in combination with legacy SH94 stretching curves),
! so for simplicity selection of SM09 stretching curve is tied with
! SM09 coordinate transform.




                                              ! Note that mathematical
      function CSF(sc, theta_s,theta_b) ! limits of csrf,CSF for
      implicit none ! theta_s, theta_b --> 0
      real*8 CSF, sc, theta_s,theta_b,csrf ! match that under "else"
                                              ! logical branches.
      if (theta_s > 0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2
      endif
      if (theta_b > 0.D0) then
        CSF=(exp(theta_b*csrf)-1.D0)/(1.D0-exp(-theta_b))
      else
        CSF=csrf ! Reference: This form of
      endif ! CSF exactly corresponds
      end ! to Eq.(2.4) from SM2009
# 165 "set_depth.F"
! The following routine just retrieves Cs_r,Cs_w, hc, VertCoordType
! from a previously opened necCDF file. Note that it is preferable to
! read Cs-curves from the file (here saved as attributes) rather than
! try to recompute them from transform parameters "theta_s, theta_b"
! because the specific details of the transformation may be changed
! at any time.


      subroutine read_scoord(ncid, N, Cs_r, Cs_w, hc, VertCoordType)
      implicit none
      integer ncid, N, ierr
      real(kind=8) Cs_r(N), Cs_w(0:N), hc
      character(len=*) VertCoordType
      include "netcdf.inc"

      VertCoordType=' '
      ierr=nf_get_att_text (ncid, nf_global, 'VertCoordType',
     & VertCoordType)
      if (ierr == nf_noerr) then
        write(*,*) 'VertCoordType=', VertCoordType
        if (VertCoordType=='NEW') VertCoordType='SM09'
      else
        write(*,*) 'Global attribute VertCoordType not found.'
        ierr=nf_get_att_text (ncid, nf_global, 'sc_r', Cs_r)
        if (ierr == nf_noerr) then
          write(*,*) 'Found  ''sc_r'':==>  old s-coord.'
        else
          ierr=nf_get_att_text (ncid, nf_global, 's_rho', Cs_r)
          if (ierr == nf_noerr) then
            write(*,*) 'Found  ''s_rho'':==>  old s-coord.'
          else
            write(*,*) '''s_rho'' not found ==> VertCoordType=''SM09''.'
            VertCoordType='SM09'
          endif
        endif
      endif

      ierr=nf_get_att_double (ncid, nf_global, 'Cs_r', Cs_r)
      if (ierr /= nf_noerr) then
        ierr=nf_get_att_double (ncid, nf_global, 'Cs_rho', Cs_r)
      endif

      if (ierr == nf_noerr) then
        ierr=nf_get_att_double (ncid, nf_global, 'Cs_w', Cs_w)
        if (ierr == nf_noerr) then
          ierr=nf_get_att_double (ncid, nf_global, 'hc', hc)
          if (ierr == nf_noerr) then
            write(*,*) 'read Cs_r,Cs_w,hc from netCDF file.'
          return !--> successful return

          else
            write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot get ',
     & 'global attribute ''hc''', nf_strerror(ierr)
          endif
        else
          write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot get global ',
     & 'attribute ''Cs_w''', nf_strerror(ierr)
        endif
      else
          write(*,'(/1x,2A/12x,A/)') '### ERROR: Cannot get global ',
     & 'attribute ''Cs_r'' or ''Cs_rho''', nf_strerror(ierr)
      endif
      stop
      end


      subroutine set_depth(Lm,Mm,N, Cs_r,Cs_w,hc, VertCoordType,
     & h, z_r,z_w)
      implicit none
      integer Lm,Mm,N, i,j,k, ierr
      real(kind=8) h(0:Lm+1,0:Mm+1), hinv(0:Lm+1), hc, ds,
     & z_r(0:Lm+1,0:Mm+1,N), Cs_r(N), Cs_r_k, cff_r,
     & z_w(0:Lm+1,0:Mm+1,0:N), Cs_w(0:N), Cs_w_k, cff_w
      character(len=*) VertCoordType
# 248 "set_depth.F"
      if (VertCoordType == 'SM09') then
        write(*,*) 'new VertCoordType = SM09'
      else
        write(*,*) 'defaulting to SH94 S-coordinate'
      endif

      ds=1.D0/dble(N)
      do j=0,Mm+1
        if (VertCoordType == 'SM09') then
          Cs_r_k=Cs_r(N) ; cff_r=-0.5D0*hc*ds
          do i=0,Lm+1
            z_w(i,j,N)=0.D0
            hinv(i)=h(i,j)/(hc+h(i,j))
            z_r(i,j,N)=hinv(i)*( cff_r + Cs_r_k*h(i,j) )
            z_w(i,j,0)=-h(i,j)
          enddo
          do k=N-1,1,-1
            Cs_w_k=Cs_w(k) ; Cs_r_k=Cs_r(k)
            cff_w=hc*ds*dble(k-N) ; cff_r=hc*ds*(dble(k-N)-0.5D0)
            do i=0,Lm+1
              z_w(i,j,k)=hinv(i)*( cff_w + Cs_w_k*h(i,j) )
              z_r(i,j,k)=hinv(i)*( cff_r + Cs_r_k*h(i,j) )
            enddo
          enddo
        else
          ierr=0
          Cs_r_k=Cs_r(N) ; cff_r=hc*(-0.5D0*ds -Cs_r_k)
          do i=0,Lm+1
            if (hc > h(i,j)) ierr=1
            z_w(i,j,N)=0.D0
            z_r(i,j,N)=cff_r + Cs_r_k*h(i,j)
            z_w(i,j,0)=-h(i,j)
          enddo
          if (ierr /= 0) then
            write(*,'(/1x,2A/)') '### ERROR: set_depth: hc > hmin ',
     & 'occurrence while selecting legacy s-coordinate.'
            stop
          endif
          do k=N-1,1,-1
            Cs_w_k=Cs_w(k) ; cff_w=hc*(ds* dble(k-N) -Cs_w_k)
            Cs_r_k=Cs_r(k) ; cff_r=hc*(ds*(dble(k-N)-0.5D0) -Cs_r_k)
            do i=0,Lm+1
              z_w(i,j,k)=cff_w + Cs_w_k*h(i,j)
              z_r(i,j,k)=cff_r + Cs_r_k*h(i,j)
            enddo
          enddo
        endif
      enddo !<-- j



      end



      subroutine compute_ubar(Lm,Mm,N, z_w, u,ubar)
      implicit none
      integer Lm,Mm,N, i,j,k
      real(kind=8) z_w(0:Lm+1,0:Mm+1,0:N)
      real(kind=4) u(1:Lm+1,0:Mm+1,N), ubar(1:Lm+1,0:Mm+1)
      do j=0,Mm+1
        do i=1,Lm+1
          ubar(i,j)=0.0
        enddo
        do k=N,1,-1
          do i=1,Lm+1
            ubar(i,j)=ubar(i,j) + u(i,j,k)*( z_w(i,j,k)+z_w(i-1,j,k)
     & -z_w(i,j,k-1)-z_w(i-1,j,k-1))
          enddo
        enddo
        do i=1,Lm+1
          ubar(i,j)=ubar(i,j)/( z_w(i,j,N)+z_w(i-1,j,N)
     & -z_w(i,j,0)-z_w(i-1,j,0))
        enddo
      enddo
      end

      subroutine compute_vbar (Lm,Mm,N, z_w, v,vbar)
      implicit none
      integer Lm,Mm,N, i,j,k
      real(kind=8) z_w(0:Lm+1,0:Mm+1,0:N)
      real(kind=4) v(0:Lm+1,1:Mm+1,N), vbar(0:Lm+1,1:Mm+1)
      do j=1,Mm+1
        do i=0,Lm+1
          vbar(i,j)=0.0
        enddo
        do k=N,1,-1
          do i=0,Lm+1
            vbar(i,j)=vbar(i,j) + v(i,j,k)*( z_w(i,j,k)+z_w(i,j-1,k)
     & -z_w(i,j,k-1)-z_w(i,j-1,k-1))
          enddo
        enddo
        do i=0,Lm+1
          vbar(i,j)=vbar(i,j)/( z_w(i,j,N)+z_w(i,j-1,N)
     & -z_w(i,j,0)-z_w(i,j-1,0))
        enddo
      enddo
      end
# 64 "tools_fort.F" 2
# 1 "def_roms_file.F" 1
      subroutine def_roms_file(ntrc,ncid, fname, xi_rho,eta_rho,s_rho,
     & theta_s,theta_b, hc, Cs_w,Cs_r,tracer,ncsrc)

! Set up netCDF structure for a file suitable to be ROMS initial or
! climatological input file. Note that the file itself is expected to
! be created externally - incoming argument "ncid" is a valid netCDF ID
! of a writable file, while this routine only takes care about creating
! dimensions, variables, and writing attributes. Similarly, it does
! not call nf_enddef to finish the definitions and switch into input
! mode leaving the calling program some room to customize the
! definitions without resorting to nf_redef (which may require moving
! the data).

      implicit none
      character(len=*) fname
      integer ncid, xi_rho,xi_u, eta_rho,eta_v, s_rho, ierr,
     & old_fill_mode, varid, lfnm,itrc
      real(kind=8) theta_s,theta_b, hc, Cs_w(0:s_rho),Cs_r(s_rho)
      integer :: ntrc,ipt_trc_len,ncsrc,trc_in ,natts,lstr,i
      character(len=20), dimension(ntrc) :: tracer
      character(len=16) :: str




      integer, parameter :: n2d=2, n3d=3

      integer r2dgrd(n2d), u2dgrd(n2d), v2dgrd(n2d),
     & r3dgrd(n3d), u3dgrd(n3d), v3dgrd(n3d)
      real(kind=8), parameter :: cycle_length=360.D0
      include "spval.h"
      include "netcdf.inc"

      call lenstr(fname,lfnm) ; write(*,'(1x,3A)', advance='no')
     & 'creating netCDF file ''', fname(1:lfnm), ''' ...'

      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)
      if (ierr /= nf_noerr) then
        write(*,'(/1x,A/)') '### WARNING: Cannot set nofill mode.'
      endif
      ipt_trc_len=20
! Define dimensions

      xi_u=xi_rho-1 ; eta_v=eta_rho-1
      ierr=nf_def_dim(ncid, 'xi_rho', xi_rho, r2dgrd(1))
      ierr=nf_def_dim(ncid, 'xi_u', xi_u, u2dgrd(1))
      ierr=nf_def_dim(ncid, 'eta_rho', eta_rho, r2dgrd(2))
      ierr=nf_def_dim(ncid, 'eta_v', eta_v, v2dgrd(2))
      ierr=nf_def_dim(ncid, 's_rho', s_rho, r3dgrd(3))




      v2dgrd(1)=r2dgrd(1) ; u2dgrd(2)=r2dgrd(2)

      r3dgrd(1)=r2dgrd(1) ; r3dgrd(2)=r2dgrd(2)
      u3dgrd(1)=u2dgrd(1) ; u3dgrd(2)=u2dgrd(2)
      v3dgrd(1)=v2dgrd(1) ; v3dgrd(2)=v2dgrd(2)

      u3dgrd(3)=r3dgrd(3) ; v3dgrd(3)=r3dgrd(3)






! Recoord S-coordinate control parameters "theta_s", "theta_b", "hc",
! and stretching curves "Cs_w", "Cs_r" at vertical W- and RHO-points.

      ierr=nf_put_att_text(ncid,nf_global,'VertCoordType',4,'SM09')

      ierr=nf_put_att_double(ncid, nf_global,'theta_s', nf_double,
     & 1, theta_s)
      ierr=nf_put_att_double(ncid, nf_global,'theta_b', nf_double,
     & 1, theta_b)
      ierr=nf_put_att_double(ncid, nf_global, 'hc',nf_double,1,hc)

      ierr=nf_put_att_double(ncid, nf_global, 'Cs_w', nf_double,
     & s_rho+1, Cs_w)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_r', nf_double,
     & s_rho, Cs_r)

! Time.
# 116 "def_roms_file.F"
      ierr=nf_def_var (ncid, 'scrum_time', nf_double, n2d-2,
     & r3dgrd(n3d), varid)
      ierr=nf_put_att_text (ncid, varid, 'long_name', 22,
     & 'time since initialization')
      ierr=nf_put_att_text (ncid, varid, 'units', 4, 'days')


! Free-surface.

      ierr=nf_def_var(ncid, 'zeta', nf_real, n2d, r2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     & 'free-surface elevation')
      ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)


! 2D momenta in XI- and ETA-directions.

      ierr=nf_def_var(ncid, 'ubar', nf_real, n2d, u2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     & 'barotropic XI-velocity')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)


      ierr=nf_def_var(ncid, 'vbar', nf_real, n2d, v2dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 23,
     & 'barotropic ETA-velocity')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

! 3D momenta in XI- and ETA-directions.

      ierr=nf_def_var(ncid, 'u', nf_real, n3d, u3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 21,
     & 'XI-velocity component')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)


      ierr=nf_def_var(ncid, 'v', nf_real, n3d, v3dgrd, varid)
      ierr=nf_put_att_text(ncid, varid, 'long_name', 22,
     & 'ETA-velocity component')
      ierr=nf_put_att_text(ncid, varid, 'units',12,'meter/second')
      ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

! Tracer variables.
      do itrc=1,ntrc/ipt_trc_len
          ! Create tracer var
          ierr=nf_def_var(ncid, trim(tracer(itrc)), nf_real, n3d, r3dgrd
     & , varid)
          ! Check if data exists in prt
          ierr=nf_inq_varid(ncsrc, trim(tracer(itrc)), trc_in)
          if (ierr == nf_noerr) then
            ierr=nf_inq_varnatts(ncsrc, trc_in, natts)
            if (ierr == nf_noerr) then
              do i=1,natts
                ierr=nf_inq_attname(ncsrc, trc_in, i, str)
                if (ierr == nf_noerr) then
                  call lenstr(str,lstr)
                  ierr=nf_copy_att(ncsrc, trc_in, str(1:lstr),
     & ncid, varid)
                endif
              enddo
            endif
          endif
          ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     & 1, FillValue)
      enddo
      end
# 65 "tools_fort.F" 2
# 1 "def_bry_file.F" 1
      subroutine def_bry_file(fname, ncid, xi_rho,eta_rho,N,
     & OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH,
     & theta_s, theta_b, hc, Cs_w,Cs_r,ntrc,tracer,ncsrc)
      implicit none
      character(len=*) fname
      integer ncid, xi_rho,eta_rho,N
      logical OBC_WEST, OBC_EAST, OBC_SOUTH, OBC_NORTH
      real(kind=8) theta_s,theta_b,hc, Cs_w(0:N),Cs_r(N)
      integer r2d_EW(2), u2d_EW(2), r3d_EW(3), u3d_EW(3),
     & r2d_NS(2), v2d_NS(2), r3d_NS(3), v3d_NS(3),
     & vert_dim, time_dim, varid, ierr,
     & old_fill_mode,lfname,itrc
      integer :: ntrc,ipt_trc_len,ncsrc,trc_in ,natts,lstr,i
      character(len=20), dimension(ntrc) :: tracer
      character(len=16) :: str
      include "spval.h"
      include "netcdf.inc"
      ipt_trc_len=20
      call lenstr(fname,lfname)
      write(*,'(1x,3A)',advance='no') 'def_bry_file :: creating ''',
     & fname(1:lfname), ''' ...'

      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)

! Dimensions

      ierr=nf_def_dim(ncid, 'xi_rho', xi_rho, r2d_EW(1))
      ierr=nf_def_dim(ncid, 'xi_u', xi_rho-1, u2d_EW(1))
      ierr=nf_def_dim(ncid, 'eta_rho', eta_rho, r2d_NS(1))
      ierr=nf_def_dim(ncid, 'eta_v', eta_rho-1, v2d_NS(1))
      ierr=nf_def_dim(ncid, 's_rho', N, vert_dim)
      ierr=nf_def_dim(ncid, 'bry_time',NF_UNLIMITED, time_dim)

      r3d_EW(1)=r2d_EW(1) ; r2d_EW(2)=time_dim
      u3d_EW(1)=u2d_EW(1) ; u2d_EW(2)=time_dim
      r3d_NS(1)=r2d_NS(1) ; r2d_NS(2)=time_dim
      v3d_NS(1)=v2d_NS(1) ; v2d_NS(2)=time_dim

      r3d_EW(2)=vert_dim ; r3d_EW(3)=time_dim
      u3d_EW(2)=vert_dim ; u3d_EW(3)=time_dim
      r3d_NS(2)=vert_dim ; r3d_NS(3)=time_dim
      v3d_NS(2)=vert_dim ; v3d_NS(3)=time_dim


! Recoord S-coordinate control parameters "theta_s", "theta_b", "hc",
! and stretching curves "Cs_w", "Cs_r" at vertical W- and RHO-points.

      ierr=nf_put_att_text(ncid, nf_global, 'VertCoordType',4,'SM09')

      ierr=nf_put_att_double(ncid, nf_global,'theta_s', nf_double, 1,
     & theta_s)
      ierr=nf_put_att_double(ncid, nf_global,'theta_b', nf_double, 1,
     & theta_b)
      ierr=nf_put_att_double(ncid, nf_global, 'hc', nf_double, 1, hc)

      ierr=nf_put_att_double(ncid, nf_global, 'Cs_w', nf_double, N+1,
     & Cs_w)
      ierr=nf_put_att_double(ncid, nf_global, 'Cs_r', nf_double, N,
     & Cs_r)

! Time

      ierr=nf_def_var (ncid, 'bry_time', nf_double, 1, time_dim, varid)
      ierr=nf_put_att_text (ncid, varid, 'long_name', 22,
     & 'time since initialization')
! ierr=nf_put_att_text (ncid, varid, 'units', 4,'days')

! Side boundary forcing variables

      if (OBC_WEST) then
        ierr=nf_def_var(ncid, 'zeta_west',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'ubar_west',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'vbar_west',nf_real, 2,v2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'u_west',nf_real, 3, r3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'v_west',nf_real, 3, v3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)
        do itrc=1,ntrc/ipt_trc_len
          ! Create tracer var
          ierr=nf_def_var(ncid, trim(tracer(itrc))/ /'_west',nf_real, 3,
     & r3d_NS ,varid)
          ! Check if data exists in prt
          ierr=nf_inq_varid(ncsrc, trim(tracer(itrc)), trc_in)
          if (ierr == nf_noerr) then
            ierr=nf_inq_varnatts(ncsrc, trc_in, natts)
            if (ierr == nf_noerr) then
              do i=1,natts
                ierr=nf_inq_attname(ncsrc, trc_in, i, str)
                if (ierr == nf_noerr) then
                  call lenstr(str,lstr)
                  ierr=nf_copy_att(ncsrc, trc_in, str(1:lstr),
     & ncid, varid)
                endif
              enddo
            endif
          endif
          ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     & 1, FillValue)
        enddo
      endif

      if (OBC_EAST) then
        ierr=nf_def_var(ncid, 'zeta_east',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'ubar_east',nf_real, 2,r2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'vbar_east',nf_real, 2,v2d_NS,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'u_east',nf_real, 3, r3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'v_east',nf_real, 3, v3d_NS, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)
        do itrc=1,ntrc/ipt_trc_len
          ! Create tracer var
          ierr=nf_def_var(ncid, trim(tracer(itrc))/ /'_east',
     & nf_real,3,r3d_NS ,varid)
          ! Check if data exists in prt
          ierr=nf_inq_varid(ncsrc, trim(tracer(itrc)), trc_in)
          if (ierr == nf_noerr) then
            ierr=nf_inq_varnatts(ncsrc, trc_in, natts)
            if (ierr == nf_noerr) then
              do i=1,natts
                ierr=nf_inq_attname(ncsrc, trc_in, i, str)
                if (ierr == nf_noerr) then
                  call lenstr(str,lstr)
                  ierr=nf_copy_att(ncsrc, trc_in, str(1:lstr),
     & ncid, varid)
                endif
              enddo
            endif
          endif
          ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     & 1, FillValue)
        enddo
      endif

      if (OBC_SOUTH) then
        ierr=nf_def_var(ncid, 'zeta_south',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'ubar_south',nf_real,2,u2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'vbar_south',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'u_south',nf_real, 3, u3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'v_south',nf_real, 3, r3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)
        do itrc=1,ntrc/ipt_trc_len
          ! Create tracer var
          ierr=nf_def_var(ncid, trim(tracer(itrc))/ /'_south', nf_real,
     & 3,r3d_EW ,varid)
          ! Check if data exists in prt
          ierr=nf_inq_varid(ncsrc, trim(tracer(itrc)), trc_in)
          if (ierr == nf_noerr) then
            ierr=nf_inq_varnatts(ncsrc, trc_in, natts)
            if (ierr == nf_noerr) then
              do i=1,natts
                ierr=nf_inq_attname(ncsrc, trc_in, i, str)
                if (ierr == nf_noerr) then
                  call lenstr(str,lstr)
                  ierr=nf_copy_att(ncsrc, trc_in, str(1:lstr),
     & ncid, varid)
                endif
              enddo
            endif
          endif
          ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     & 1, FillValue)
        enddo
      endif

      if (OBC_NORTH) then
        ierr=nf_def_var(ncid, 'zeta_north',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid, varid, 'units', 5, 'meter')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'ubar_north',nf_real,2,u2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'vbar_north',nf_real,2,r2d_EW,varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'u_north',nf_real, 3, u3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)

        ierr=nf_def_var(ncid, 'v_north',nf_real, 3, r3d_EW, varid)
        ierr=nf_put_att_text(ncid,varid, 'units', 12,'meter/second')
        ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real, 1,
     & FillValue)
        do itrc=1,ntrc/ipt_trc_len
          ! Create tracer var
          ierr=nf_def_var(ncid, trim(tracer(itrc))/ /'_north', nf_real,
     & 3,r3d_EW ,varid)
          ! Check if data exists in prt
          ierr=nf_inq_varid(ncsrc, trim(tracer(itrc)), trc_in)
          if (ierr == nf_noerr) then
            ierr=nf_inq_varnatts(ncsrc, trc_in, natts)
            if (ierr == nf_noerr) then
              do i=1,natts
                ierr=nf_inq_attname(ncsrc, trc_in, i, str)
                if (ierr == nf_noerr) then
                  call lenstr(str,lstr)
                  ierr=nf_copy_att(ncsrc, trc_in, str(1:lstr),
     & ncid, varid)
                endif
              enddo
            endif
          endif
          ierr=nf_put_att_real(ncid, varid, '_FillValue', nf_real,
     & 1, FillValue)
        enddo
      endif
      write(*,*) '...done.'
      end
# 66 "tools_fort.F" 2
# 1 "r2r_subs.F" 1


      subroutine set_mask(nx,ny, rmsk, mask,umask,vmask)
      implicit none
      integer nx,ny, i,j
      real(kind=8) rmsk(nx,ny)
      integer(kind=2) mask(nx,ny), umask(2:nx,ny), vmask(nx,2:ny)
      do j=1,ny
        do i=1,nx ! Convert RHO-point mask,
          if (rmsk(i,j) < 0.5D0) then ! real(kind=8) array into
            mask(i,j)=0 ! short 2-byte integer.
          else
            mask(i,j)=1 ! Note that here U- and V-
          endif ! masks are set to zero only
        enddo ! for U- and V-point points
        do i=2,nx
          if (rmsk(i,j) < 0.5D0 .and. rmsk(i-1,j) < 0.5D0) then
            umask(i,j)=0
          else ! which are FULLY INSIDE land
            umask(i,j)=1 ! (hence ".and." logic instead
          endif ! of ".or." in the actual ROMS
        enddo ! model), while points of
        if (j > 1) then ! normal velocity components
          do i=1,nx
            if (rmsk(i,j) < 0.5D0 .and. rmsk(i,j-1) < 0.5D0) then
              vmask(i,j)=0
            else ! on the coast line - the
              vmask(i,j)=1 ! ones where normal velocities
            endif ! must be set to zero due to
          enddo ! no-flux boundary conditions
        endif ! - are considered to be part
      enddo ! of the the solution and
      end ! therefore left unmasked.


      subroutine set_contour_mask(nx,ny, rmsk, mask,umask,vmask)
      implicit none
      integer nx,ny, i,j, iw2,iw,ie,js2,js,jn
      real(kind=8) rmsk(nx,ny)
      integer(kind=2) mask(nx,ny), umask(2:nx,ny), vmask(nx,2:ny)
      do j=1,ny
        js2=max(1,j-2) ; js=max(1,j-1) ; jn=min(j+1,ny)
        do i=1,nx
          iw=max(1,i-1) ; ie=min(i+1,nx)
          if (rmsk(i,j) < 0.5D0 .and. ( rmsk(ie,j) > 0.5D0
     & .or. rmsk(i,jn) > 0.5D0 .or. rmsk(iw,j) > 0.5D0
     & .or. rmsk(i,js) > 0.5D0 )) then
            mask(i,j)=0
          else
            mask(i,j)=1
          endif
        enddo
        do i=2,nx
          iw2=max(1,i-2) ; iw=max(1,i-1) ; ie=min(i+1,nx)
          if ( rmsk(i-1,j) < 0.5D0 .and. rmsk(i,j) < 0.5D0 .and.
     & ( rmsk(i-1,jn) > 0.5D0 .or. rmsk(i,jn) > 0.5D0 .or.
     & rmsk(iw2,j) > 0.5D0 .or. rmsk(ie,j) > 0.5D0 .or.
     & rmsk(i-1,js) > 0.5D0 .or. rmsk(i,js) > 0.5D0 )) then
            umask(i,j)=0
          else
            umask(i,j)=1
          endif
        enddo
        if (j > 1) then
          do i=1,nx
            iw=max(1,i-1) ; ie=min(i+1,nx)
            if ( rmsk(i,j-1) < 0.5D0 .and. rmsk(i,j) < 0.5D0 .and.
     & ( rmsk(i,jn) > 0.5D0 .or.
     & rmsk(iw,j ) > 0.5D0 .or. rmsk(ie,j ) > 0.5D0 .or.
     & rmsk(iw,j-1) > 0.5D0 .or. rmsk(ie,j-1) > 0.5D0 .or.
     & rmsk(i,js2) > 0.5D0 )) then
              vmask(i,j)=0
            else
              vmask(i,j)=1
            endif
          enddo
        endif
      enddo
      end

      subroutine apply_mask(nx,ny,N, mask,q)
      integer nx,ny,N, i,j,k
      integer(kind=2) mask(nx,ny)
      real (kind=4) q(nx,ny,N)





      if (N > 1) then
        do j=1,ny
          do i=1,nx
            if (mask(i,j) == 0) then
              do k=1,N
                q(i,j,k)=0.
              enddo
            endif
          enddo
        enddo
      else
        k=1
        do j=1,ny
          do i=1,nx
            if (mask(i,j) == 0) q(i,j,k)=0.
          enddo
        enddo
      endif
      end

      subroutine bry_apply_mask(nx,N, mask,q)
      integer nx,N, i,k
      integer(kind=2) mask(nx) ! Apply land mask
      real (kind=4) q(nx,N) ! to horizontally



      if (N > 1) then ! one-dimensional
        do i=1,nx ! array. Note that
          if (mask(i) == 0) then ! mask array itself
            do k=1,N ! is horizontally
              q(i,k)=0. ! one-dimensional
            enddo ! too.
          endif
        enddo
      else
        k=1
        do i=1,nx
          if (mask(i) == 0) q(i,k)=0.
        enddo
      endif
      end
# 228 "r2r_subs.F"
      subroutine compute_min_max(nx,ny, x, xmin,xmax)
      implicit none
      integer nx,ny, i,j
      real(kind=8) x(nx,ny), xmin,xmax
      xmin=x(1,1) ; xmax=x(1,1)
      do j=1,ny ! Compute min,max values of
        do i=1,nx ! array "x" irrespective of
          if (x(i,j) < xmin) then ! what is its meaning.
            xmin=x(i,j)
          elseif (x(i,j) > xmax) then
            xmax=x(i,j)
          endif
        enddo
      enddo
      write(*,*) 'min,max =', xmin,xmax
      end

      subroutine adjust_lon_into_range(nx,ny, x, west,east)
      implicit none
      integer nx,ny, i,j
      real(kind=8) x(nx,ny), west,east, xmin,xmax

      xmin=x(1,1) ; xmax=x(1,1) ! Compute min,max values of
      do j=1,ny ! longitude (array "x") whose
        do i=1,nx ! definition is presumed to be
          if (x(i,j) < xmin) then ! unknown in advance, e.g., it
            xmin=x(i,j) ! may be from 0 to 360 or from
          elseif (x(i,j) > xmax) then ! -180 to +180 or whatever --
            xmax=x(i,j) ! it is however presumed to be
          endif ! continuous within the array.
        enddo ! Check whether there is an
      enddo ! overlap with the externally
      write(*,*) 'min,max =', xmin,xmax ! specified west,east range.
      if (xmin > east) then ! Adjust by adding/subtracting
        do j=1,ny ! 360 degrees is necessary to
          do i=1,nx ! make it overlap.
            x(i,j)=x(i,j)-360.D0
          enddo
        enddo
        write(*,*) 'adjusted to', xmin-360.D0,xmax-360.D0
      elseif (xmax < west) then
        do j=1,ny
          do i=1,nx
            x(i,j)=x(i,j)+360.D0
          enddo
        enddo
        write(*,*) 'adjusted to', xmin+360.D0, xmax+360.D0
      else
        write(*,*) 'no need to adjust parent-grid longitude'
      endif
      end

! The following routines find bounding indices on parent grid which
! define minimal logically rectangular patch fully containing unmasked
! portion of child grid. It is almost the same finding min,max for two
! integer arrays, expect that there are two caveats:
!
! (i) some of indices ip,jp may be set to non-positive values to
! indicate that that portion of child grid is outside the parent,
! so interpolation is not possible there -- it is still OK as
! long as these areas are under land mask -- these points are
! ignored for the purpose of search below. For this reason the
! search is done by full 2D-sweep rather than moving along the
! perimeter of child grid (which may not be even closed because
! of having special-valued (ip,jp)-s; another good reason for
! needing 2D-sweep is land mask; and
!
! (ii) indices [ip(i,j),jp(i,j)] mean that child-grid point i,j is
! located somewhere within the parent-grid area bounded by 4
! vertices
!
! (xp(ip,jp+1),yp(ip,jp+1)) --- (xp(ip+1,jp+1),yp(ip+1,jp+1))
! | |
! | |
! (xp(ip,jp),yp(ip,jp)) --- (xp(ip+1,jp),yp(ip+1,jp))
!
! so +1 is added to both illegal and illegal at the very end.


      subroutine compute_index_bounds(ncx,ncy, ip,jp, mask, imin,imax,
     & jmin,jmax)
      implicit none
      integer ncx,ncy,ip(ncx,ncy),jp(ncx,ncy), imin,imax,jmin,jmax,i,j
      integer(kind=2) mask(ncx,ncy)
      imin=10000000 ; imax=-1 !<-- initialize to unrealistic
      jmin=10000000 ; jmax=-1 !<-- values outside the range
      do j=1,ncy
        do i=1,ncx
          if (mask(i,j) > 0 .and. ip(i,j) > 0 .and. jp(i,j) > 0) then
            imin=min(imin,ip(i,j)) ; jmin=min(jmin,jp(i,j))
            imax=max(imax,ip(i,j)) ; jmax=max(jmax,jp(i,j))
          endif
        enddo
      enddo
      imax=imax+1 ; jmax=jmax+1 !<-- because of (ii) above
      end


! Same as above, but to be applied along a 1D line. Note that this
! time imin,imax,jmin,jmax are expected to be initialized externally
! and there is no adding +1 at the end: this is because consecutive
! calls to this routine applied to the different sides of perimeter
! of the grid (e.g., open boundaries) contribute to finding the same
! bounds and it needs to be summarized externally after all done.

      subroutine r2r_bry_index_bounds(ncx, ip,jp, mask, imin,imax,
     & jmin,jmax)
      implicit none
      integer ncx, ip(ncx),jp(ncx), imin,imax,jmin,jmax, ic
      integer(kind=2) mask(ncx)
      do ic=1,ncx
        if (mask(ic) > 0 .and. ip(ic) > 0 .and. jp(ic) > 0) then
          imin=min(imin,ip(ic)) ; jmin=min(jmin,jp(ic))
          imax=max(imax,ip(ic)) ; jmax=max(jmax,jp(ic))
        endif
      enddo
      end
# 411 "r2r_subs.F"
! The following routine is called from "r2r_init" and writes out all
! the incoming 2D arrays into a special-purpose netCDF file it creates
! internally. This is needed only for diagnostic purposes only and does
! not affect outcome of any computations. Consequently the routine is
! written in a relaxed way without performing checks for netCDF errors.


      subroutine r2r_init_diag_file(ncx,ncy, N, ip,jp,xi,eta,
     & ipu,jpu,xiu,etau, ipv,jpv,xiv,etav,
     & csA,snA, h, hprnt, kprnt)
      implicit none
      integer ncx,ncy, N
      integer(kind=4), dimension(ncx,ncy) :: ip,jp, ipu,jpu, ipv,jpv
      real(kind=8), dimension(ncx,ncy) :: csA,snA, h,hprnt, xi,eta,
     & xiu,etau, xiv,etav
      real(kind=8), dimension(ncx,ncy,N) :: kprnt
      integer ncid, ierr,i,j,dimids(4), ipvar, jpvar, xivar, etavar,
     & old_fill_mode, csvar, ipuvar,jpuvar, xiuvar,etauvar,
     & snvar, ipvvar,jpvvar, xivvar,etavvar,
     & hvar, hpvar, dhvar, kpvar
      real(kind=8), allocatable, dimension(:,:) :: dh

      include "netcdf.inc"

      write(*,'(1x,2A)',advance='no') 'creating diagnostic file ',
     & '''croco_init_diag.nc''...'
      ierr=nf_create('croco_init_diag.nc', nf_netcdf4, ncid)
      ierr=nf_set_fill(ncid, nf_nofill, old_fill_mode)

      ierr=nf_def_dim(ncid, 'xi_rho', ncx, dimids(1))
      ierr=nf_def_dim(ncid, 'eta_rho', ncy, dimids(2))
      ierr=nf_def_dim(ncid, 's_rho', N, dimids(3))
!>
      ierr=nf_def_var(ncid, 'ip', nf_int, 2, dimids, ipvar)
      ierr=nf_put_att_text(ncid, ipvar, 'long_name', 37,
     & 'parent-to-child interpolation i-index')
      ierr=nf_put_att_int(ncid, ipvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'jp', nf_int, 2, dimids, jpvar)
      ierr=nf_put_att_text(ncid, jpvar, 'long_name', 37,
     & 'parent-to-child interpolation j-index')
      ierr=nf_put_att_int(ncid, jpvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'xi', nf_double, 2, dimids, xivar)
      ierr=nf_put_att_text(ncid, xivar, 'long_name', 38,
     & 'parent-to-child interpolation offset X')
      ierr=nf_put_att_double(ncid, xivar, '_FillValue', nf_double,
     & 1, -1.D0)

      ierr=nf_def_var(ncid, 'eta', nf_double, 2, dimids, etavar)
      ierr=nf_put_att_text(ncid, etavar, 'long_name', 38,
     & 'parent-to-child interpolation offset Y')
      ierr=nf_put_att_double(ncid, etavar,'_FillValue', nf_double,
     & 1, -1.D0)
!>
      ierr=nf_def_var(ncid, 'ipu', nf_int, 2, dimids, ipuvar)
      ierr=nf_put_att_text(ncid, ipuvar, 'long_name', 39,
     & 'parent-to-child U-interpolation i-index')
      ierr=nf_put_att_int(ncid, ipuvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'jpu', nf_int, 2, dimids, jpuvar)
      ierr=nf_put_att_text(ncid, jpuvar, 'long_name', 39,
     & 'parent-to-child U-interpolation j-index')
      ierr=nf_put_att_int(ncid, jpuvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'xiu', nf_double, 2, dimids, xiuvar)
      ierr=nf_put_att_text(ncid, xiuvar, 'long_name', 40,
     & 'parent-to-child U-interpolation offset X')
      ierr=nf_put_att_double(ncid, xiuvar, '_FillValue', nf_double,
     & 1, -1.D0)

      ierr=nf_def_var(ncid, 'etau', nf_double, 2, dimids, etauvar)
      ierr=nf_put_att_text(ncid, etauvar, 'long_name', 40,
     & 'parent-to-child U-interpolation offset Y')
      ierr=nf_put_att_double(ncid, etauvar,'_FillValue', nf_double,
     & 1, -1.D0)
!>
      ierr=nf_def_var(ncid, 'ipv', nf_int, 2, dimids, ipvvar)
      ierr=nf_put_att_text(ncid, ipvvar, 'long_name', 39,
     & 'parent-to-child V-interpolation i-index')
      ierr=nf_put_att_int(ncid, ipvvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'jpv', nf_int, 2, dimids, jpvvar)
      ierr=nf_put_att_text(ncid, jpvvar, 'long_name', 39,
     & 'parent-to-child V-interpolation j-index')
      ierr=nf_put_att_int(ncid, jpvvar, '_FillValue', nf_int, 1,-1)

      ierr=nf_def_var(ncid, 'xiv', nf_double, 2, dimids, xivvar)
      ierr=nf_put_att_text(ncid, xivvar, 'long_name', 40,
     & 'parent-to-child V-interpolation offset X')
      ierr=nf_put_att_double(ncid, xivvar, '_FillValue', nf_double,
     & 1, -1.D0)

      ierr=nf_def_var(ncid, 'etav', nf_double, 2, dimids, etavvar)
      ierr=nf_put_att_text(ncid, etauvar, 'long_name', 40,
     & 'parent-to-child V-interpolation offset Y')
      ierr=nf_put_att_double(ncid, etavvar,'_FillValue', nf_double,
     & 1, -1.D0)
!>
      ierr=nf_def_var(ncid, 'csA', nf_double, 2, dimids, csvar)
      ierr=nf_put_att_text(ncid, csvar, 'long_name', 44,
     & 'cosine of child-parent grid angle difference')

      ierr=nf_def_var(ncid, 'snA', nf_double, 2, dimids, snvar)
      ierr=nf_put_att_text(ncid, snvar, 'long_name', 42,
     & 'sine of child-parent grid angle difference')

      ierr=nf_def_var(ncid, 'h', nf_double, 2, dimids, hvar)
      ierr=nf_put_att_text(ncid, hvar, 'long_name', 35,
     & 'child grid native bottom topography')

      ierr=nf_def_var(ncid, 'hprnt', nf_double, 2, dimids, hpvar)
      ierr=nf_put_att_text(ncid, hpvar, 'long_name', 35,
     & 'interpolated parent grid topography')

      ierr=nf_def_var(ncid, 'dh', nf_double, 2, dimids, dhvar)
      ierr=nf_put_att_text(ncid, dhvar, 'long_name', 34,
     & 'child-parent topography difference')

      ierr=nf_def_var(ncid, 'kprnt', nf_real, 3, dimids, kpvar)
      ierr=nf_put_att_text(ncid, kpvar, 'long_name', 56,
     & 'index coordinate for parent-child vertical interpolation')

!>
      ierr=nf_enddef(ncid)

      ierr=nf_put_var_int(ncid, ipvar, ip)
      ierr=nf_put_var_int(ncid, jpvar, jp)
      ierr=nf_put_var_double(ncid, xivar, xi)
      ierr=nf_put_var_double(ncid, etavar, eta)

      ierr=nf_put_var_int(ncid, ipuvar, ipu)
      ierr=nf_put_var_int(ncid, jpuvar, jpu)
      ierr=nf_put_var_double(ncid, xiuvar, xiu)
      ierr=nf_put_var_double(ncid, etauvar, etau)

      ierr=nf_put_var_int(ncid, ipvvar, ipv)
      ierr=nf_put_var_int(ncid, jpvvar, jpv)
      ierr=nf_put_var_double(ncid, xivvar, xiv)
      ierr=nf_put_var_double(ncid, etavvar, etav)

      ierr=nf_put_var_double(ncid, csvar, csA)
      ierr=nf_put_var_double(ncid, snvar, snA)

      ierr=nf_put_var_double(ncid, hvar, h)
      ierr=nf_put_var_double(ncid, hpvar, hprnt)

      ierr=nf_put_var_double(ncid, kpvar, kprnt)

      allocate(dh(ncx,ncy))
      do j=1,ncy
        do i=1,ncx
          dh(i,j)=h(i,j)-hprnt(i,j)
        enddo
      enddo
      ierr=nf_put_var_double(ncid, dhvar, dh)
      deallocate(dh)
      ierr=nf_close(ncid)
      write(*,'(2x,A)') '...done'
      end
# 67 "tools_fort.F" 2
# 1 "r2r_vert_interp.F" 1
! Content of this package:
!------------------------- ! initialization of vertical interpolation
! r2r_init_vertint_thread ! by inverse mapping by cubic splines: given
! r2r_check_vertint_thread ! parent-grid zp_r(kp=1,Np) and child-grid
! r2r_init_vrtint_tile ! z_r(k=1,N) compute parent grid index-space
! r2r_check_vrtint_tile ! coordinate kprnt(k=1,N) such that
! comp_zspline_tile ! zspline[kp=kprnt(k)] = z_r(k)
!
! ! interpolation routines for fields:
! r2r_vrtint_thread ! compute vertical derivatives in parent
! r2r_vertint_tile ! grid index-coordinate space, then knowing
! r2r_vsplnint_tile ! kp=kprnt(k=1,N), interpolate parent-grid
! ! data onto child grid.
!
! r2r_set_depth_tile ! compute z_r=z_r[h(j,j), Cs(k)]
!
! compute_uvbar_thread ! vertical integration of u,v-velocity
! compute_uvbar_tile ! components to compute barotripic mode
# 29 "r2r_vert_interp.F"
      subroutine r2r_init_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer :: ncx,ncy, Np,N
      real(kind=8) :: hprnt(ncx,ncy),h(ncx,ncy),hcp,hc,Csp_r(Np),
     & Cs_r(N),kprnt(ncx,ncy,N)
      integer :: ntrds,trd, nsub_x,nsub_y, ntls, isize,
     & my_first,my_last,tile, istr,iend,jstr,jend
C$ integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$ ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zc_size < isize*(N+1) .or.
     & allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) then
            deallocate(zp_r,drv,zc)
        endif
        allocate( zc(alloc_zc_size), zp_r(allc_zpr_size),
     & drv(allc_zpr_size) )
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)') 'r2r_init_vertint_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     & 'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)

      elseif (trd == 0) then
        write(*,'(2x,A)',advance='no')
     & 'entering r2r_init_vertint_thread ...'

      endif

      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend )
        call r2r_init_vrtint_tile( istr,iend,jstr,jend, ncx,ncy,
     & hprnt,Np,hcp,Csp_r,zp_r,drv,
     & h, N,hc,Cs_r,zc, kprnt)
      enddo

      if (trd == 0) write(*,'(2x,A)')
     & 'leaving r2r_init_vertint_thread'

C$OMP BARRIER
      end

      subroutine r2r_check_vertint_thread(ncx,ncy, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer ncx,ncy, Np,N
      real(kind=8) hprnt(ncx,ncy),h(ncx,ncy), hcp,hc,Csp_r(Np),Cs_r(N),
     & kprnt(ncx,ncy,N),my_error
      integer ntrds,trd, nsub_x,nsub_y, ntls, isize,
     & my_first,my_last, tile, istr,iend,jstr,jend
C$ integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$ ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zc_size < isize*(N+1) .or.
     & allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) then
            deallocate(zp_r,drv,zc)
        endif
        allocate(zc(alloc_zc_size), zp_r(allc_zpr_size),
     & drv(allc_zpr_size))
      endif
      my_error=0.D0 !<-- initialize
      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend )
        call r2r_check_vrtint_tile(istr,iend,jstr,jend, ncx,ncy,
     & hprnt, Np,hcp,Csp_r,zp_r,drv,
     & h, N,hc,Cs_r,zc, kprnt, my_error)
      enddo
C$OMP CRITICAL(r2r_vert_crgn)
      if (trd_count == 0) vert_int_error=0.D0
      trd_count=trd_count+1
      vert_int_error=max(vert_int_error, my_error)
      if (trd_count == ntrds) then
        trd_count=0
        write(*,*)
        write(*,*) 'maximum vert_int_error =', vert_int_error
        write(*,*)
      endif
C$OMP END CRITICAL(r2r_vert_crgn)
      end



      subroutine r2r_init_vrtint_tile( istr,iend,jstr,jend, ncx,ncy,
     & hprnt,Np,hcp,Csp_r, zp_r,drv,
     & h,N, hc, Cs_r,z_r, kprnt)

! Given a set of parent-grid z-levels zp_r=z_r(kp), kp=1,..,Np defined
! at vertical RHO-points, as well as knowing that zp=zeta at surface
! kp=Np+1/2, (always presume zeta=0 for this purpose), and zp=-h_parent
! at bottom, kp=1/2, construct

! (i) a set of spine derivatives drv=drv(kp) = d zp_r/dkp (vertical
! derivative in index space of parent grid) defined at the same
! RHO-points, then

! (ii) construct an inverse mapping function
!
! kprnt=kprnt(z), 1/2 < kprnt < Np+1/2
!
! which maps child-grid z-levels z_r(k), k=1,..,N onto
! parent-grid vertical index kp in such a way that index-space
! spline-interpolated value
!
! z* = spline[p, zp_r(kp), zp_r(kp+1), drv(kp), drv(kp+1)]
!
! where index "kp" in the integer part of kprnt(z) and "p" is
! the fractional reminder, is equal to "z" itself.
!
! Once this mapping kp = kp(k) = kprnt[z_r(k)] is established, any
! field "qsrc" defined on the parent grid, hence qsrc(kp), kp=1,...,Np
! can be spline interpolated in index space,
!
! q(k)= spline[p, qsrc(kp), src(kp+1), drvq(kp), drvq(kp+1)]
!
! where drvq = d qsrc/d kp is vertical derivative of qsrc computed
! in index space of parent grid.

! Task (ii) above is achieved by first using linear interpolation to
! map child-grid z-levels z_r(k) onto parent-grid zp_r(kp). This leads
! to an an initial approximation for kprnt(k). Subsequently use Newton
! iterations to compute fractional offsets "p" in such a way that spline
! version of z_r(k) = spline[p, zp_r(kp),...] item (ii) above, holds.

! The linear stage constructs "continuous index" mapping coordinate
! kp=kprnt(k) such that interpolation of field "qsrc" defined on the
! parent grid z=zp_r(kp), kp=1,...,Np into child grid z=z_r(k),
! k=1,...,N is approximated as
!
! kp=int(kprnt(k)) ; p = kprnt(k)-float(kp)
! q(k)=p*qsrc(kp+1) + (1-p)*qsrc(kp)
!
! where, "kp" and "p" are the integer and fractional parts of "kprnt"
! defined in such a way that zp_r(kp), kp=1,...,Np interpolated the
! same way into the location z_r(k) should yield z_r(k) itself, that
! is
! z_r(k)=p*zp_r(kp+1) + (1-p)*zp_r(kp)
!
! where index "kp" is chosen to make z_r(k) bounded as
!
! zp_r(kp) <= z_r(k) < zp_r(kp+1)
!
! and the fractional part
!
! p = [z_r(k)-zp_r(kp)]/[zp_r(kp+1)-zp_r(kp)]
!
! The caveats are:
!
! (iii) it is possible that some child-grid "z_r" points may be above
! the uppermost available zp_r(Np); in this case kp=Np while
! fraction "p" is set using unperturbed free surface z=0 instead
! of nonexisting zp_r(kp+1); interpolation becomes extrapolation
! toward surface;
!
! (iv) for a similar reason some "z_r" may be below the lowest
! available zp_r(kp=1), and even worse, because of inconsistency
! between parent and child topography some "z_r" points may be
! even below bottom of the parent grid.

! To address both (or at least mitigate in the case of (iv)) spline
! interpolation algorithm for zp_r (see "comp_zspline_tile" below with
! CPP-switches STAGG_BCS and/or STAGG_NOT_A_KNOT_BCS activated) is
! designed to utilize boundary conditions
!
! zp(kprnt=Np+1/2) = zeta !<-- surface
! zp(kprnt=1/2) = -h_parent !<-- bottom
!
! in addition to the usual Np equations for spline derivatives: i.e.,
! the extremal kp=1 and kp=Np are no-longer treated as boundaries, but
! are formed from the same assumption of second-derivative continuity
! as for all the interior points. Furthermore, to simplify the searches
! and interpolations "comp_zspline_tile" produces ghost point values
!
! zp_r(kp=0), drv(0), zp_r(Np+1), drv(Np+1)
!
! which are at RHO-points half-grid-interval below the bottom and
! above the surface. These ghost points are constructed in such a way
! that spline polynomial of zp_r taken at kprnt=1/2 yields exactly
! -h_parent (bottom) and at kprnt=N+1/2 yields zeta (=0 in all cases
! here, surface). Then, finally,
!
! kprnt > Np+1/2 above surface, should never occur
! Np < kprnt < Np+1/2 extrapolation toward surface; accepted
! 1 < kprnt < Np interpolation within the vertical column
! 1/2 < kprnt < 1 extrapolation toward bottom;
! kprnt < 1/2 below bottom of parent grid; use filling


      implicit none
      integer :: istr,iend,jstr,jend, ncx,ncy, Np,N
      real(kind=8) :: hprnt(ncx,ncy), h(ncx,ncy), Csp_r(Np), Cs_r(N),
     & zp_r(istr:iend,0:Np+1), z_r(istr:iend,N), hcp,hc,
     & drv(istr:iend,0:Np+1), kprnt(ncx,ncy,N), p
      integer :: i,j,k, kp

      integer :: iter
      real(kind=8) :: q,pq, zerr, dZds

     & , dZs2, d2Zds2


      do j=jstr,jend
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     & h, N, hc,Cs_r, z_r)
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     & hprnt, Np,hcp,Csp_r, zp_r(istr,1) )
        do i=istr,iend
          zp_r(i,Np+1)=0.D0 !<-- surface
          zp_r(i,0)=-hprnt(i,j) !<-- bottom
        enddo

        call comp_zspline_tile(istr,iend, Np, zp_r,drv)

        do i=istr,iend !--> search loop
          kp=Np
          do k=N,1,-1 !--> recursive because of kp
            do while(z_r(i,k) < zp_r(i,kp) .and. kp > 0)
              kp=kp-1
            enddo
            kprnt(i,j,k)=dble(kp)+0.5D0 !<-- temporarily
          enddo
        enddo ! set initial approximation
        do k=1,N ! for fractional distance "p"
          do i=istr,iend ! by linear interpolation then
            kp=int(kprnt(i,j,k))
            p=(z_r(i,k)-zp_r(i,kp))/(zp_r(i,kp+1)-zp_r(i,kp))
            kprnt(i,j,k)=dble(kp)+p
          enddo ! apply Newton iterations

          do i=istr,iend
            if (kprnt(i,j,k) > 0.499999999999D0) then
              kp=int(kprnt(i,j,k)) ; p=kprnt(i,j,k)-dble(kp)
              do iter=1,8
                q=1.D0-p ; pq=p*q

                zerr=p*((p+2.D0*pq)*zp_r(i,kp+1) -pq*drv(i,kp+1))
     & +q*((q+2.D0*pq)*zp_r(i,kp ) +pq*drv(i,kp ))
     & -z_r(i,k)
                dZds=6.D0*pq*(zp_r(i,kp+1)-zp_r(i,kp))
     & +(p*p-2.D0*pq)*drv(i,kp+1)
     & +(q*q-2.D0*pq)*drv(i,kp )




                d2Zds2=(2.D0*p-q)*drv(i,kp+1) -(2.D0*q-p)*drv(i,kp)
     & -3.D0*(p-q)*(zp_r(i,kp+1)-zp_r(i,kp))

                dZs2=dZds*dZds
                p=p-zerr*(dZs2+d2Zds2*zerr)/(dZs2*dZds)

              enddo
              kprnt(i,j,k)=dble(kp)+p !<-- corrected value
            else
              kprnt(i,j,k)=0.D0
            endif
          enddo

        enddo !<-- k
      enddo !<-- j
      end

      subroutine r2r_check_vrtint_tile(istr,iend,jstr,jend, ncx,ncy,
     & hprnt, Np,hcp,Csp_r, zp_r,drv,
     & h, N,hc,Cs_r,z_r, kprnt, my_error)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, Np,N, i,j,k, kp
      real(kind=8) hprnt(ncx,ncy), h(ncx,ncy), Csp_r(Np), Cs_r(N),
     & zp_r(istr:iend,0:Np+1), z_r(istr:iend,N), hcp,hc,
     & drv(istr:iend,0:Np+1), kprnt(ncx,ncy,N), my_error
     & , p

     & , q, pq

      do j=jstr,jend
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     & h, N, hc,Cs_r, z_r)
        call r2r_set_depth_tile(istr,iend, j,j, ncx,ncy,
     & hprnt, Np,hcp,Csp_r, zp_r(istr,1) )
        do i=istr,iend
          zp_r(i,Np+1)=0.D0 !<-- surface
          zp_r(i,0)=-hprnt(i,j) !<-- bottom
        enddo

        call comp_zspline_tile(istr,iend, Np, zp_r,drv)

        do i=istr,iend
          do k=1,N
            kp=max(int(kprnt(i,j,k)), 0) ; p=kprnt(i,j,k)-float(kp)

            if (kprnt(i,j,k) > 0.499999999999D0) then
              q=1.D0-p ; pq=p*q
              my_error=max( my_error, abs(
     & p*((p+2.D0*pq)*zp_r(i,kp+1) -pq*drv(i,kp+1))
     & +q*((q+2.D0*pq)*zp_r(i,kp ) +pq*drv(i,kp ))
     & -z_r(i,k) ) )
            endif







          enddo
        enddo
      enddo !<-- j
      end







      subroutine comp_zspline_tile(istr,iend, Np, zp_r,drv)

! Compute spline derivatives of drv=d zp_r/d kp for the parent grid
! vertical coordinate. The incoming argument zp_r=zp_r(kp) contains
! z-levels at RHO-point for indices k=1,..,Np, AND, IN ADDITION TO
! THAT indices kp=0 and kp=Np+1 contain bottom zp_r(0)=-h_parent and
! surface zp_r(Np+1)=zeta value. On output zp_r(0) and zp_r(Np+1)
! are converted into EXTRAPOLATED ghost-point values, kind of half-
! grid interval below the bottom and above the surface. Derivatives
! drv(0) and drv(Np+1) are also at ghost points, while all others
! for kp=1,...,Np are for regular RHO-points.

      implicit none
      integer istr,iend, Np
      real(kind=8) zp_r(istr:iend,0:Np+1), drv(istr:iend,0:Np+1)
      integer i,k
      real(kind=8) CF(1:Np-1), cff,cff1


      CF(1)=1.D0/3.D0 ; cff=8.D0/9.D0
      do i=istr,iend
        drv(i,1)=cff*(zp_r(i,2) -zp_r(i,0))
      enddo
# 398 "r2r_vert_interp.F"
      do k=2,Np-1,+1 !--> forward elimination
        cff=1.D0/(4.D0-CF(k-1)) ; CF(k)=cff
        do i=istr,iend
          drv(i,k)=cff*( 3.D0*(zp_r(i,k+1)-zp_r(i,k-1)) -drv(i,k-1))
        enddo
      enddo


      cff=1.D0/(3.D0-CF(Np-1)) ; cff1=8.D0/3.D0
      do i=istr,iend
        drv(i,Np)=cff*( cff1*(zp_r(i,Np+1)-zp_r(i,Np-1)) -drv(i,Np-1) )
      enddo
# 423 "r2r_vert_interp.F"
      do k=Np-1,1,-1 !--> backsubstitution
        cff=CF(k)
        do i=istr,iend
          drv(i,k)=drv(i,k)-cff*drv(i,k+1)
        enddo
      enddo

      do i=istr,iend

        drv(i,0)=2.D0*drv(i,2)-drv(i,1) -8.D0*zp_r(i,0)+12.D0*zp_r(i,1)
     & -4.D0*zp_r(i,2)
        zp_r(i,0)=4.D0*zp_r(i,0)-4.D0*zp_r(i,1)+zp_r(i,2)
     & -0.5D0*(drv(i,2)-drv(i,1))

        drv(i,Np+1)=2.D0*drv(i,Np-1)-drv(i,Np) +8.D0*zp_r(i,Np+1)
     & -12.D0*zp_r(i,Np)+4.D0*zp_r(i,Np-1)
        zp_r(i,Np+1)=4.D0*zp_r(i,Np+1)-4.D0*zp_r(i,Np)+zp_r(i,Np-1)
     & -0.5D0*(drv(i,Np)-drv(i,Np-1))







      enddo
      end





      subroutine r2r_vrtint_thread(ncx,ncy, lmsk,mask, btm_bc,
     & Np,qsrc, N,kprnt,qtr)
      use r2r_vertint_vars
      implicit none
      integer :: ncx,ncy, lmsk, btm_bc, Np,N
      integer(kind=2) :: mask(ncx,ncy)
      real(kind=8) :: kprnt(ncx,ncy,N)
      real(kind=4) :: qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
!>
      integer :: ntrds,trd, nsub_x,nsub_y, ntls, my_first,my_last,
     & tile, istr,iend,jstr,jend
C$ integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$ ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      do tile=my_first,my_last,+1
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

        call r2r_vsplnint_tile( istr,iend,jstr,jend, ncx,ncy,
     & lmsk,mask, btm_bc, Np,qsrc,drv, N,kprnt,qtr)




      enddo
      end

      subroutine r2r_vretint_tile(istr,iend,jstr,jend, ncx,ncy,
     & lmsk,mask, Np,qsrc, N,kprnt,qtr)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, lmsk,Np,N, i,j,k,kp
      integer(kind=2) mask(ncx,ncy)
      real(kind=8) kprnt(ncx,ncy,N), p
      real(kind=4) qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)

       include "spval.h"

      do j=jstr,jend
        do k=1,N
          do i=istr,iend
            kp=int(kprnt(i,j,k)) ; p=kprnt(i,j,k)-float(kp)
            if (kp > Np-1) then
              qtr(i,j,k)=qsrc(i,j,Np)
            elseif (kp > 0) then
              qtr(i,j,k)=p*qsrc(i,j,kp+1)+(1.D0-p)*qsrc(i,j,kp)
            else
              qtr(i,j,k)=qsrc(i,j,1)
            endif
          enddo
        enddo
        if (lmsk == 1) then
          do i=istr,iend
            if (mask(i,j) == 0) then
              do k=1,N
                qtr(i,j,k)=0.
              enddo
            endif
          enddo
        endif
      enddo
      end



      subroutine r2r_vsplnint_tile(istr,iend,jstr,jend, ncx,ncy, lmsk,
     & mask, btm_bc,Np,qsrc,drv, N,kprnt,qtr)
      implicit none
      integer :: istr,iend,jstr,jend, ncx,ncy, lmsk,btm_bc,Np,N,
     & i,j,k,kp
      integer(kind=2) ::mask(ncx,ncy)
      real(kind=4) :: qsrc(ncx,ncy,Np), qtr(ncx,ncy,N)
      real(kind=8) :: kprnt(ncx,ncy,N), CF(Np), drv(istr:iend,Np),
     & cff, p,q,pq

      include "spval.h"

      do j=jstr,jend
        if (btm_bc == -1) then ! staggered no-slip b.c.
          CF(1)=1.D0/3.D0 ; cff=8.D0/9.D0 ! combined with not-a-knot
          do i=istr,iend ! at botton
            drv(i,1)=cff*qsrc(i,j,2)
          enddo
        elseif (btm_bc == +1) then ! staggered Neumann b.c.
          CF(1)=7.D0/15.D0 ; cff=1.2D0 ! combined with not-a-knot
          do i=istr,iend
            drv(i,1)=cff*(qsrc(i,j,2)-qsrc(i,j,1))
          enddo
        else
          CF(1)=0.5D0 ! "natural" lower b.c.
          do i=istr,iend
            drv(i,1)=1.5D0*(qsrc(i,j,2)-qsrc(i,j,1))
          enddo
        endif

        do k=2,Np-1,+1 !--> forward elimination
          CF(k)=1.D0/(4.D0-CF(k-1)) ; cff=CF(k)
          do i=istr,iend
            drv(i,k)=cff*( 3.D0*(qsrc(i,j,k+1)-qsrc(i,j,k-1))
     & -drv(i,k-1) )
          enddo
        enddo

        cff=1.D0/(2.D0-CF(Np-1)) !<-- upper b.c.
        do i=istr,iend
          drv(i,Np)=cff*( 3.D0*(qsrc(i,j,Np)-qsrc(i,j,Np-1))
     & -drv(i,Np-1) )
        enddo
        do k=Np-1,1,-1 !--> backsubstitution
          cff=CF(k)
          do i=istr,iend
            drv(i,k)=drv(i,k)-cff*drv(i,k+1)
          enddo
        enddo
        do k=1,N
          do i=istr,iend
            kp=int(kprnt(i,j,k))
            kp=max( 1, int(kprnt(i,j,k)) )
            p=max(-0.5D0, kprnt(i,j,k)-float(kp) )
            if (kp > Np-1) then ! extrapolate within
              qtr(i,j,k)=qsrc(i,j,Np)+p*drv(i,Np) ! upper half of upper
            elseif (kp > 0) then ! most grid box.
              q=1.D0-p ; pq=p*q
              qtr(i,j,k)=p*((p+2.D0*pq)*qsrc(i,j,kp+1) -pq*drv(i,kp+1))
     & +q*((q+2.D0*pq)*qsrc(i,j,kp ) +pq*drv(i,kp ))
c** else
c** qtr(i,j,k)=qsrc(i,j,1)
            endif
          enddo
        enddo
        if (lmsk == 1) then
          do i=istr,iend
            if (mask(i,j) == 0) then
              do k=1,N
                qtr(i,j,k)=0.
              enddo
            endif
          enddo
        endif
      enddo
      end



      subroutine r2r_set_depth_tile(istr,iend,jstr,jend, ncx,ncy,
     & h, N,hc,Cs, z)

! Compute
! hc*s(k) + Cs(k)*h(i,j)
! z(i,j,k) = h(i,j) * ------------------------
! hc + h(i,j)
!
! at vertical RHO- or W-points which is determined by looking at the
! first value of Cs(k): for for W-type -1 <= Cs <= 0 reaching both -1
! and 0 at the ends; for RHO-type it stays fully inside the interval.
! Place the outcome into private array.

      implicit none
      integer :: istr,iend,jstr,jend, ncx,ncy,N, i,j,k
      real(kind=8) :: h(ncx,ncy), hc, Cs(N), z(istr:iend,jstr:jend,N),
     & hcds,bias,cf1,cf2

      if (Cs(1) > -0.999999D0) then
        hcds=hc/dble(N) ; bias=0.5D0 !<-- RHO-points
      else
        hcds=hc/dble(N-1) ; bias=0.D0 !<-- W-points
      endif
      do j=jstr,jend
        do i=istr,iend
          z(i,j,N)=h(i,j)/(hc+h(i,j)) !<-- temporarily
        enddo
        do k=1,N
          cf1=hcds*(dble(k-N)-bias) ; cf2=Cs(k)
          do i=istr,iend
            z(i,j,k)= z(i,j,N)*(cf1 + cf2*h(i,j))
          enddo
        enddo
      enddo
      end


      subroutine compute_uvbar_thread(ncx,ncy, N, hc,Cs_w, h, u,v,
     & ubar,vbar)

! Compute barotropic velocities "ubar" and "vbar" from "u" and "v" --
! basically it is just vertical integration, except that the use of 3D
! arrays z_w Hz is avoided for memory reason: using vertical slice of
! "z_w" instead, resulting in somewhat awkward-looking code.

      use comp_uvbar_vars
      implicit none
      integer ncx,ncy, N
      real(kind=8) hc, Cs_w(0:N), h(ncx,ncy)
      real(kind=4) u(2:ncx,ncy,N), ubar(2:ncx,ncy),
     & v(ncx,2:ncy,N), vbar(ncx,2:ncy)
      integer ntrds,trd, nsub_x,nsub_y, tile,my_first,my_last,
     & ntls, isize, istr,iend,jstr,jend, istrU,jstrV
C$ integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$ ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)
      isize=(ncx+nsub_x-1)/nsub_x
      if (alloc_zw_size < (isize+2)*(N+1)) then
        alloc_zw_size=(isize+2)*(N+1)
        if (allocated(z_w)) deallocate(Hz,z_w,D)
        allocate(D(isize+2), z_w(alloc_zw_size), Hz(2*alloc_zw_size))
C$OMP CRITICAL(comp_uvbar_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)') 'compute_uvbar_thread: ',
     & 'allocated', (2*(isize+1)+6*alloc_zw_size)/float(262144),
     & 'MB private workspace, trd =', trd
C$OMP END CRITICAL(comp_uvbar_crgn)
      endif

      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
        istrU=max(istr,2) ; jstrV=max(jstr,2)
        call compute_uvbar_tile(istr,iend,jstr,jend, istrU,jstrV,
     & ncx,ncy, N, hc,Cs_w, h, u,v, ubar,vbar, D,z_w,Hz)
      enddo
      end

! Note that istrU,jstrV are computed above by the driver rather than
! inside the working routine below because the first dimension of "z_w"
! must be set from "istrU-1" rather than based on "istr" (say "istr-1")
! because private array dimensions inside "r2r_set_depth_tile" are tied
! to its loop ranges, leaving no other choice but to have horizontal
! dimensions of "z_w" be the same in both "compute_uvbar_tile" and
! "r2r_set_depth_tile".

      subroutine compute_uvbar_tile(istr,iend,jstr,jend, istrU,jstrV,
     & ncx,ncy, N, hc,Cs_w, h, u,v, ubar,vbar, D,z_w,Hz)
      implicit none
      integer istr,iend,jstr,jend, istrU,jstrV, ncx,ncy, N, i,j,k,j1
      real(kind=8) hc, Cs_w(0:N), h(ncx,ncy), cff, D(istrU-1:iend),
     & z_w(istrU-1:iend,0:N), Hz(istrU-1:iend,2,N)
      real(kind=4) u(2:ncx,ncy,N), ubar(2:ncx,ncy),
     & v(ncx,2:ncy,N), vbar(ncx,2:ncy)

      j1=2
      do j=jstrV-1,jend,+1 !--> recursive because of Hz.
        call r2r_set_depth_tile(istrU-1,iend, j,j, ncx,ncy,
     & h, N+1,hc,Cs_w, z_w)
        j1=3-j1 !<-- rotate index
        do k=1,N
          do i=istrU-1,iend
            Hz(i,j1,k)=z_w(i,k)-z_w(i,k-1)
          enddo
        enddo
        if (j > jstr-1) then
          do i=istrU,iend
            ubar(i,j)=0. ; D(i)=0.D0
          enddo
          do k=N,1,-1
            do i=istrU,iend
              cff=Hz(i,j1,k)+Hz(i-1,j1,k)
              D(i)=D(i)+cff
              ubar(i,j)=ubar(i,j)+cff*u(i,j,k)
            enddo
          enddo
          do i=istrU,iend
            ubar(i,j)=ubar(i,j)/D(i)
          enddo
        endif
        if (j > jstrV-1) then
          do i=istr,iend
            vbar(i,j)=0. ; D(i)=0.D0
          enddo
          do k=N,1,-1 ! Here Hz(i,3-j1,k) is
            do i=istr,iend ! the j-slice computed
              cff=Hz(i,j1,k)+Hz(i,3-j1,k) ! during the previous j.
              D(i)=D(i)+cff
              vbar(i,j)=vbar(i,j)+cff*v(i,j,k)
            enddo
          enddo
          do i=istr,iend
            vbar(i,j)=vbar(i,j)/D(i)
          enddo
        endif
      enddo !<-- j
      end
# 68 "tools_fort.F" 2
# 1 "r2r_rotate.F" 1
! module r2r_rotate_scratch
! integer, save :: alloc_FX_size=0
!C$OMP THREADPRIVATE(alloc_FX_size)
! real(kind=4), allocatable, dimension(:):: FX,FE
!C$OMP THREADPRIVATE(FX,FE)
! end module r2r_rotate_scratch


      subroutine r2r_rotate_shift_thread(nx,ny,N, csA,snA, ur,vr, u,v)

! Rotate vector components (ur,vr) colocated ar RHO-points while
! simultaneously shifting them into U- and V-locations on C-grid.
! inputs csA,snA, ur,vr (all remain unchanged); outpits u,v

! use r2r_rotate_scratch
      implicit none
      integer nx,ny,N
      real(kind=8) csA(nx,ny),snA(nx,ny)
      real(kind=4) ur(nx,ny,N),vr(nx,ny,N), u(nx-1,ny,N),v(nx,ny-1,N)
!>
      integer ntrds,trd, nsub_x,nsub_y, ntls,isize,jsize,
     & my_first,my_last, tile, istr,iend,jstr,jend
      integer :: alloc_FX_size
      real(kind=4), allocatable, dimension(:):: FX,FE
C$ integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$ ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      alloc_FX_size=0
      call set_tiles(nx,ny, nsub_x,nsub_y)

      isize=(nx+nsub_x-1)/nsub_x ; jsize=(ny+nsub_y-1)/nsub_y
      if (alloc_FX_size < (isize+4)*(jsize+4)) then
        alloc_FX_size=(isize+4)*(jsize+4)
        if (allocated(FX)) deallocate(FX,FE)
        allocate(FX(alloc_FX_size),FE(alloc_FX_size))
C$OMP CRITICAL(ext_wr_rgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)') 'r2r_rotate_shift_thread: ',
     & 'allocated', float(2*alloc_FX_size)/float(262144),
     & 'MB of private workspace, trd =', trd
!#if > 1
! write(*,'(4x,2(A,I5,1x),2(1x,A,I4),1x,2(1x,A,I4))')
! & 'nx =',nx, 'ny =',ny, 'nsub_x =',nsub_x, 'nsub_y =',nsub_y,
! & 'isize =',isize, 'jsize =',jsize
!#endif
C$OMP END CRITICAL(ext_wr_rgn)

      elseif (trd == 0) then
        write(*,'(1x,A)',advance='no')
     & 'enter r2r_rotate_shift_thread...'

      endif

      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)

      do tile=my_first,my_last,+1
        call comp_tile_bounds( tile, nx,ny, nsub_x,nsub_y,
     & istr,iend,jstr,jend )
        call r2r_rotate_shift_tile( istr,iend,jstr,jend,
     & nx,ny,N, csA,snA, ur,vr, u,v, FX,FE )
      enddo

      if (trd == 0) write(*,'(2x,A)') 'leaving r2r_rotate_shift_thread'

      end


      subroutine r2r_rotate_shift_tile(istr,iend,jstr,jend, nx,ny,N,
     & csA,snA, ur,vr, u,v, FX,FE)
      implicit none
      integer istr,iend,jstr,jend, nx,ny,N, istrU,jstrV, i,j,k
      real(kind=8) csA(nx,ny),snA(nx,ny)
      real(kind=4) ur(nx,ny,N), u(2:nx,ny,N),
     & vr(nx,ny,N), v(nx,2:ny,N),
     & FX(istr-1:iend,jstr-1:jend),
     & FE(istr-1:iend,jstr-1:jend)
      istrU=max(istr,2) ; jstrV=max(jstr,2)
      do k=1,N
        do j=jstrV-1,jend
          do i=istrU-1,iend
            FX(i,j)=ur(i,j,k)*csA(i,j)+vr(i,j,k)*snA(i,j)
            FE(i,j)=vr(i,j,k)*csA(i,j)-ur(i,j,k)*snA(i,j)
          enddo
        enddo
        do j=jstr,jend
          do i=istrU,iend
            u(i,j,k)=0.5*(FX(i,j)+FX(i-1,j))
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            v(i,j,k)=0.5*(FE(i,j)+FE(i,j-1))
          enddo
        enddo
      enddo !<-- k
      end



      subroutine r2r_rotate_in_place_thread(nx,ny,N, csA,snA, ur,vr)

! Rotate vector components (ur,vr) colocated ar RHO-points and place
! the outcome back into the same arrays.

      implicit none
      integer nx,ny,N
      real(kind=8) csA(nx,ny), snA(nx,ny)
      real(kind=4) ur(nx,ny,N), vr(nx,ny,N)
!>
      real(kind=4) cff1,cff2
      integer ntrds,trd, nsub_x,nsub_y, istr,iend,jstr,jend,
     & ntls,tile, my_first,my_last, i,j,k
C$ integer omp_get_thread_num, omp_get_num_threads
      ntrds=1 ; trd=0
C$ ntrds=omp_get_num_threads() ; trd=omp_get_thread_num()
      if (trd == 0) write(*,'(2x,A)',advance='no')
     & 'entering r2r_rotate_r_thread...'

      call set_tiles(nx,ny, nsub_x,nsub_y)

      ntls=(nsub_x*nsub_y+ntrds-1)/ntrds
      my_first=trd*ntls -(ntls*ntrds-nsub_x*nsub_y)/2
      my_last=min(my_first+ntls-1, nsub_x*nsub_y-1)
      my_first=max(my_first, 0)

      do tile=my_first,my_last,+1
        call comp_tile_bounds(tile, nx,ny, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
        do j=jstr,jend
          do k=1,N
            do i=istr,iend
              cff1=ur(i,j,k)*csA(i,j)+vr(i,j,k)*snA(i,j)
              cff2=vr(i,j,k)*csA(i,j)-ur(i,j,k)*snA(i,j)
              ur(i,j,k)=cff1 ; vr(i,j,k)=cff2
            enddo
          enddo
        enddo
      enddo
      if (trd == 0) write(*,'(2x,A)') '...done'
      end
# 69 "tools_fort.F" 2
# 1 "r2r_bry_rotate.F" 1
! This is a rather silly situation where the velocity components u,v
! need to be rotated, but only one is needed as the output. So the two
! versions are generated by CPP and are different by having _u_ or _v_
! in the middle of the subroutine name, and by selecting one line or
! the other. (Note that both lines cannot coexist because the first
! one overwrites u).
# 18 "r2r_bry_rotate.F"
      subroutine bry_rotate_u_in_place(ncx,N, csA,snA, u,v)
      implicit none
      integer ncx,N
      real(kind=8) csA(ncx),snA(ncx)
      real(kind=4) u(ncx,N),v(ncx,N)
C$OMP PARALLEL SHARED(ncx,N, csA,snA, u,v)
      call bry_rotate_u_inp_thread(ncx,N, csA,snA, u,v)
C$OMP END PARALLEL
      end

      subroutine bry_rotate_u_inp_thread(ncx,N, csA,snA, u,v)
      implicit none
      integer ncx,N
      real(kind=8) csA(ncx),snA(ncx), cosA,sinA
      real(kind=4) u(ncx,N),v(ncx,N)
      integer icmin,icmax,isize, istr,iend,tile, i,k
      integer numthreads, trd, chunk_size
C$ integer omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        do i=istr,iend
          cosA=csA(i) ; sinA=snA(i) ! <-- invariant for k-index
          do k=1,N

            u(i,k)=u(i,k)*cosA +v(i,k)*sinA



          enddo
        enddo
      enddo
      end





# 1 "r2r_bry_rotate.F" 1
! This is a rather silly situation where the velocity components u,v
! need to be rotated, but only one is needed as the output. So the two
! versions are generated by CPP and are different by having _u_ or _v_
! in the middle of the subroutine name, and by selecting one line or
! the other. (Note that both lines cannot coexist because the first
! one overwrites u).
# 18 "r2r_bry_rotate.F"
      subroutine bry_rotate_v_in_place(ncx,N, csA,snA, u,v)
      implicit none
      integer ncx,N
      real(kind=8) csA(ncx),snA(ncx)
      real(kind=4) u(ncx,N),v(ncx,N)
C$OMP PARALLEL SHARED(ncx,N, csA,snA, u,v)
      call bry_rotate_v_inp_thread(ncx,N, csA,snA, u,v)
C$OMP END PARALLEL
      end

      subroutine bry_rotate_v_inp_thread(ncx,N, csA,snA, u,v)
      implicit none
      integer ncx,N
      real(kind=8) csA(ncx),snA(ncx), cosA,sinA
      real(kind=4) u(ncx,N),v(ncx,N)
      integer icmin,icmax,isize, istr,iend,tile, i,k
      integer numthreads, trd, chunk_size
C$ integer omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        do i=istr,iend
          cosA=csA(i) ; sinA=snA(i) ! <-- invariant for k-index
          do k=1,N



            v(i,k)=v(i,k)*cosA -u(i,k)*sinA

          enddo
        enddo
      enddo
      end
# 61 "r2r_bry_rotate.F" 2
# 70 "tools_fort.F" 2

!!! tools used for smoothing !!!
# 1 "smooth_utils.F" 1
      subroutine rx_diag_tile(istr,iend,jstr,jend, Lm,Mm, h,msk,
     & iter,ntiles,ifrst_call)

      implicit none

      real(kind=8), dimension(0:Lm+1,0:Mm+1) :: h
      real(kind=8) cff, my_rx, my_ry
      real(kind=8) :: rx_max=0.D0, ry_max=0.D0

      integer istr,iend,jstr,jend, Lm,Mm, iter,ntiles
      integer(kind=2), dimension(0:Lm+1,0:Mm+1) :: msk
      integer i,j, my_i_rx, my_i_ry, my_j_rx, my_j_ry
      integer :: trd_count=0, ifrst_call,
     & irx_max=0, jrx_max=0, iry_max=0, jry_max=0

                                          ! This part is purely for
      my_rx=0.D0 ; my_i_rx=0 ; my_j_rx=0 ! diagnostics and causes no
      my_ry=0.D0 ; my_i_ry=0 ; my_j_ry=0 ! effect onto the resultant
                                          ! topography.
      do j=jstr,jend+1
        do i=istr,iend+1
          if (msk(i,j)> 0 .and. msk(i-1,j) > 0) then
            cff=abs(h(i,j)-h(i-1,j))/(h(i,j)+h(i-1,j))
            if (cff > my_rx) then
              my_rx=cff ; my_i_rx=i ; my_j_rx=j
            endif
          endif
          if (msk(i,j) > 0 .and. msk(i,j-1) > 0) then
            cff=abs(h(i,j)-h(i,j-1))/(h(i,j)+h(i,j-1))
            if (cff > my_ry) then
              my_ry=cff ; my_i_ry=i ; my_j_ry=j
            endif
          endif
        enddo
      enddo
! write(*,*) 'my_rx:',my_rx,'/ my_ry:',my_ry
C$OMP CRITICAL(rxdiag_cr_rgn)
      if (ifrst_call == 1) then
        ifrst_call=0
        write(*,'(/1x,A,2x,A,1x,A,14x,A,1x,A/1x,2A)') 'iter',
     & 'i,jrx_max', 'rx_max', 'i,jry_max', 'ry_max',
     & '---------------------------------',
     & '---------------------------------'
      endif
      trd_count=trd_count+1
      if (my_rx > rx_max) then
        rx_max=my_rx
        irx_max=my_i_rx ; jrx_max=my_j_rx
      endif
      if (my_ry > ry_max) then
        ry_max=my_ry
        iry_max=my_i_ry ; jry_max=my_j_ry
      endif
      if (trd_count == ntiles) then
        trd_count=0
        write(*,'(I5,2(2x,I4,I5,1x,A,F10.7))') iter, irx_max,jrx_max,
     & 'rx_max =',rx_max, iry_max,jry_max, 'ry_max =',ry_max
        rx_max=0.D0
        ry_max=0.D0
      endif
C$OMP END CRITICAL(rxdiag_cr_rgn)
      end



      subroutine neumann_bc_tile(istr,iend,jstr,jend, Lm,Mm,A)
      implicit none
      integer istr,iend,jstr,jend, Lm,Mm, i,j
      real(kind=8) A(0:Lm+1,0:Mm+1)
      if (istr == 1) then
        do j=jstr,jend
          A(istr-1,j)=A(istr,j)
        enddo
      endif
      if (iend == Lm) then
        do j=jstr,jend
          A(iend+1,j)=A(iend,j)
        enddo
      endif
      if (jstr == 1) then
        do i=istr,iend
          A(i,jstr-1)=A(i,jstr)
        enddo
      endif
      if (jend == Mm) then
        do i=istr,iend
          A(i,jend+1)=A(i,jend)
        enddo
      endif
      if (istr == 1 .and. jstr == 1) then
        A(istr-1,jstr-1)=A(istr,jstr)
      endif
      if (istr == 1 .and. jend == Mm) then
        A(istr-1,jend+1)=A(istr,jend)
      endif
      if (iend == Lm .and. jstr == 1) then
        A(iend+1,jstr-1)=A(iend,jstr)
      endif
      if (iend == Lm .and. jend == Mm) then
        A(iend+1,jend+1)=A(iend,jend)
      endif
      end
# 73 "tools_fort.F" 2
# 1 "tiling.F" 1
! A set of generic tiling tools for 2D subdomain decomposition needed
! for a coarse-grained shared-memory OpenMP parallelization approach
! similar to that of the actual ROMS code with the exception that now
! it is adapted for dynamically changing dimensions, so the numbers of
! partitions are selected automatically from the grid array dimensions
! and number of threads, with an effort to minimize the miss-balance
! (if the dimensions are not evenly divisible by the desired number of
! tiles), as well as to bring subdomain sizes close to the optimum for
! cache performance, subject to judicial compromise in the case of
! conflict between the two goals.
!
! Content:
!
! set_tiles(nx,ny, nsub_x,nsub_y)
!
! input: nx,ny -- grid dimensions;
! number of threads (determined internally);
!
! output: nsub_x -- numbers of tiles in each direction (their
! nsub_y product is guaranteed to be divisible by
! the number of threads;
!
! comp_tile_bounds(tile, nx,ny,nsub_x,nsub_y, istr,iend,jstr,jend)
!
! input: tile -- tile number ranging from 0 to nsub_x*nsub_y-1;
! nx,ny, nsub_x,nsub_y -- same as above;
!
! output: istr,iend,jstr,jend -- starting and ending indices
! of subdomain "tile";

c--#define TEST
# 45 "tiling.F"
      subroutine set_tiles(nx,ny, nsub_x,nsub_y)
      implicit none
      integer nx,ny, nsub_x,nsub_y
      integer ntrds, ntx,nty, nsb, i,size,excess,max_exc
C$ integer omp_get_num_threads
      integer,parameter :: min_vect_lenght=96, targ_length=128

      ntrds=1
C$ ntrds=omp_get_num_threads()

      if (nx>8*min_vect_lenght-1 .and. mod(ntrds,8)==0) then
        ntx=8
      elseif (nx>6*min_vect_lenght-1 .and. mod(ntrds,6)==0) then
        ntx=6
      elseif (nx>4*min_vect_lenght-1 .and. mod(ntrds,4)==0) then
        ntx=4
      elseif (nx>3*min_vect_lenght-1 .and. mod(ntrds,3)==0) then
        ntx=3
      elseif (nx>2*min_vect_lenght-1 .and. mod(ntrds,2)==0) then
        ntx=2
      else ! Scan possible tile sizes within
        ntx=1 ! the range of "i" index and select
      endif ! one which yields the least number
      nty=ntrds/ntx ! of excess points, subject to the
                                  ! constraint that "nsb" can be
      max_exc=ny ! evenly divided by "nty"...
      do i=9,25
        nsb=(ny+i-1)/i !<-- prospective tile size
        nsb=nsb-mod(nsb,nty) !<-- make "nsb" divisible by "nty"
        if (nsb == 0) nsb=nty
        size=(ny+nsb-1)/nsb !--> actual tile size
        excess=nsb*size-ny
        if (excess<max_exc) then
          max_exc=excess
          nsub_y=nsb
        endif



      enddo

      if (mod(nsub_y,ntrds) > 0) then
        nsb=nsub_y/nty
        if (mod(nsb,3)==0 .and. mod(ntx,3)==0) then
          nsb=nsb/3
          ntx=ntx/3
        endif
        do while(mod(nsb,2)==0 .and. mod(ntx,2)==0)
          nsb=nsb/2
          ntx=ntx/2
        enddo
      else
        ntx=1
      endif

      nsub_x=(nx+targ_length-1)/targ_length
      if (ntx > 1) then
        i=mod(nsub_x,ntx) !--> make sure that "nsub_x"
        if (i > ntx/2) then ! is divisible by "ntx"
          nsub_x=nsub_x+ntx-i
        else
          nsub_x=nsub_x-i
        endif
      endif
# 123 "tiling.F"
C$OMP CRITICAL (tiling_cr_rgn)
      write(*,'(2(1x,A,I4),3(1x,A,I3))') 'set_tiles: nx =', nx,
     & 'ny =',ny, 'numthreads =', ntrds, 'nsub_x =', nsub_x,
     & 'nsub_y =', nsub_y
C$OMP END CRITICAL (tiling_cr_rgn)

      end

      subroutine comp_tile_bounds(tile, nx,ny,nsub_x,nsub_y,
     & istr,iend,jstr,jend)
      implicit none
      integer tile, nx,ny,nsub_x,nsub_y, istr,iend,jstr,jend,
     & i,j,size
      j=tile/nsub_x
      i=tile-nsub_x*j ! This tiling algorithm
                                         ! generally follows ROMS,
      size=(nx+nsub_x-1)/nsub_x ! where tile sizes in each
      istr=1+i*size -(nsub_x*size-nx)/2 ! direction are computed by
      iend=min(istr+size-1,nx) ! integer division rounded
      istr=max(istr,1) ! upward, while the excess
                                         ! (in the case when the
      size=(ny+nsub_y-1)/nsub_y ! domain cannot be divided
      jstr=1+j*size -(nsub_y*size-ny)/2 ! evenly) is half-and-half
      jend=min(jstr+size-1,ny) ! split between the ends.
      jstr=max(jstr,1) ! As the result tiles the
                                         ! first and the last in each
      end ! direction may be smaller.
# 74 "tools_fort.F" 2
# 1 "lenstr.F" 1
      subroutine lenstr(str,lstr)

! Find the position of the last non-blank character in input string
! after removing all leading blanks, if any. At first, find the length
! of input string using intrinsic function "len" and search for the
! last and the first non-blank character, "ie" and "is". Move the whole
! string to the beginning if there are leading blanks (is>1). Returned
! value "lenstr" is the position of the last non-blanc character of the
! modified string.

! WARNING: if there are leading blank characters, user must ensure
! that the string is "writable", i.e., there is a character variable
! in the calling program which holds the string: otherwise call to
! lenstr results in segmentation fault, i.e. passing directly typed
! argument like
! lstr=lenstr(' x...')
!
! is not allowed, however
!
! lstr=lenstr('x...')
!
! is OK because lenstr makes no attempt to shift the string.

! implicit none ! In the code below there
      character(len=*), intent(inout) :: str ! are two possible outcomes
      integer, intent(out) :: lstr
      integer :: is,ie ! of the search for the first
      ie=len(str)
      do while(ie > 1 .and. str(ie:ie) == ' ')
        ie=ie-1
      enddo ! non-blank character "is":
      is=1
      do while(is < ie .and. str(is:is) == ' ')
        is=is+1
      enddo ! it either finds one, or
      if (str(is:is) /= ' ') then ! the search is terminated
        if (is > 1) str=str(is:ie) ! by reaching the condition
        lstr=ie-is+1 ! (is == ie), while the
      else ! character is still blank,
        lstr=0 ! which means that the
      endif ! whole string consists of
      end subroutine lenstr ! blank characters only.
# 75 "tools_fort.F" 2
# 1 "mrg_zone_subs.F" 1





      subroutine etch_mgz_weights_thread(ncx,ncy, mask, mgz,wrk, wdth)

! Assuming that "mgz" is initialized as mgz=1 just on the row of
! perimeter points along the unmasked part of open boundary etch it
! into the interior of the domain by increasing "mgz" values at
! unmasked points by 1 at every step in the location where the value
! iself or at least one of its 4 immediate neighbors is already
! positive. After completion this procedure creates a constant-slope
! shape function function connected to the boundary by water points.
! Argument "mask" is input-only land mask; "mgz" is input-output;
! "wrk" is used as work array;

      implicit none
      integer ncx,ncy, wdth
      integer(kind=2), dimension(ncx,ncy) :: mask, mgz,wrk
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     & range, istr,iend,jstr,jend, iter, max_iters, i,j
C$ integer omp_get_num_threads, omp_get_thread_num
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c*** nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

! Note that "mgz" and "wrk" switch roles in the two code segments
! below: output of the first becomes input of the second, and vice
! versa. Other than that and tile reversal the two are identical.

      max_iters=wdth/2
      if (wdth > 2*max_iters) max_iters=max_iters+1
      do iter=1,max_iters
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          call etch_weights_tile(istr,iend,jstr,jend, ncx,ncy,
     & mask, mgz,wrk)
        enddo
C$OMP BARRIER
        do tile=my_first,my_last+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)

          if (iter < max_iters .or. wdth == 2*max_iters) then
            call etch_weights_tile(istr,iend,jstr,jend, ncx,ncy,
     & mask, wrk,mgz)
          else
            do j=jstr,jend ! Copy the outcome from the previous
              do i=istr,iend ! call to etch_weights_tile back into
                mgz(i,j)=wrk(i,j) ! array "mgz" in the case when "wdth"
              enddo ! is an odd number. No need to do so
            enddo ! for even "wdth" as the final state
          endif ! is naturally there.
        enddo
C$OMP BARRIER
      enddo
      end


      subroutine etch_weights_tile(istr,iend,jstr,jend, ncx,ncy,
     & mask, ms1,ms2)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, i,j, iw,ie,js,jn
      integer(kind=2), dimension(ncx,ncy) :: mask, ms1, ms2
      do j=jstr,jend
        js=max(j-1,1) ; jn=min(j+1,ncy)
        do i=istr,iend
          iw=max(i-1,1) ; ie=min(i+1,ncx)
          if ( mask(i,j) > 0 .and. ( ms1(i,j) > 0 .or.
     & ms1(iw,j) > 0 .or. ms1(ie,j) > 0 .or.
     & ms1(i,js) > 0 .or. ms1(i,jn) > 0 )) then
            ms2(i,j)=ms1(i,j)+1
          else
            ms2(i,j)=0
          endif
        enddo
      enddo
      end



      subroutine mrg_zone_cont_thread(ncx,ncy, mgz,ms1,ms2)

! Enforce the property that every water point within the merging
! zone identified as mgz=1 can be reached from the boundary by water
! (unreachable points mgz=1 are reset to mgz=0). The procedure is
! essentially the same as in "single_connect" with the exception that
! the initial points are set on the perimeter, and, in principle,
! there is a possibility of multiple unconnected merging zones
! separated from each other by land.

      implicit none
      integer ncx,ncy
      integer(kind=2), dimension(ncx,ncy) :: mgz, ms1, ms2
      integer trd_count, npts, npts_bak
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     & range, istr,iend,jstr,jend, iter, my_sum, i,j

      trd_count=0; npts=0; npts_bak=-1
C$ integer omp_get_num_threads, omp_get_thread_num
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c*** nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      if (trd == 0) then
        write(*,*) 'Enforcing continuity of merging zone, ',
     & 'numthreads =', numthreads
        write(*,*) 'total number of points in grid',
     & (ncx-2)*(ncy-2), ' excluding perimeter rows'
        trd_count=0 ; npts=0 ; npts_bak=-1 !<-- initialize
      endif

      do tile=my_first,my_last
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
        do j=jstr,jend
          do i=iend,istr
            ms1(i,j)=0 ; ms2(i,j)=0 !<-- initialize
          enddo
        enddo
        if (istr==1) then ! initialize search
          do j=jstr,jend ! by setting merging
            if (mgz(istr,j) > 0) then ! zone masks along
              ms1(istr,j)=1 ; ms2(istr,j)=1 ! the perimeter
            endif
          enddo
        endif
        if (iend==ncx) then
          do j=jstr,jend
            if (mgz(iend,j) > 0) then
              ms1(iend,j)=1 ; ms2(iend,j)=1
            endif
          enddo
        endif
        if (jstr==1) then
          do i=istr,iend
            if (mgz(i,jstr) > 0) then
              ms1(i,jstr)=1 ; ms2(i,jstr)=1
            endif
          enddo
        endif
        if (jend==ncy) then
          do i=istr,iend
            if (mgz(i,jend) > 0) then
              ms1(i,jend)=1 ; ms2(i,jend)=1
            endif
          enddo
        endif
      enddo !<-- tile
C$OMP BARRIER

      iter=0
      do while(npts /= npts_bak)
        my_sum=0 ; iter=iter+1
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          if (istr==1) istr=istr+1
          if (iend==ncx) iend=iend-1
          if (jstr==1) jstr=jstr+1
          if (jend==ncy) jend=jend-1

          do j=jstr,jend
            do i=istr,iend
              if ( mgz(i,j) > 0 .and. ( ms1(i,j ) > 0 .or.
     & ms1(i-1,j) > 0 .or. ms1(i+1,j) > 0 .or.
     & ms1(i,j-1) > 0 .or. ms1(i,j+1) > 0 )) then
                ms2(i,j)=1
                my_sum=my_sum+1
              endif
            enddo
          enddo
        enddo !<-- tile
C$OMP CRITICAL(cr_region)
        if (trd_count == 0) then
          npts_bak=npts ; npts=0
        endif
        npts=npts+my_sum
        trd_count=trd_count+1
        if (trd_count == numthreads) then
          trd_count=0
          if (mod(iter,20) == 0 .or. npts == npts_bak) then
            write(*,'(8x,A,I7,2(2x,A,I10))') 'iter =', iter,
     & 'npts =', npts, 'changes =', npts-npts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER

        my_sum=0 ; iter=iter+1
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          if (istr==1) istr=istr+1
          if (iend==ncx) iend=iend-1
          if (jstr==1) jstr=jstr+1
          if (jend==ncy) jend=jend-1

          do j=jstr,jend
            do i=istr,iend
              if ( mgz(i,j) > 0 .and. ( ms2(i,j ) > 0 .or.
     & ms2(i-1,j) > 0 .or. ms2(i+1,j) > 0 .or.
     & ms2(i,j-1) > 0 .or. ms2(i,j+1) > 0 )) then
                ms1(i,j)=1
                my_sum=my_sum+1
              endif
            enddo
          enddo
        enddo
C$OMP CRITICAL(cr_region)
        if (trd_count == 0) then
          npts_bak=npts ; npts=0
        endif
        npts=npts+my_sum
        trd_count=trd_count+1
        if (trd_count == numthreads) then
          trd_count=0
          if (mod(iter,20) == 0 .or. npts == npts_bak) then
            write(*,'(8x,A,I7,2(2x,A,I10))') 'iter =', iter,
     & 'npts =', npts, 'changes =', npts-npts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo !<-- while

      do tile=my_last,my_first,-1
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
        do j=jstr,jend
          do i=istr,iend ! cancel points
             if (ms1(i,j) == 0) mgz(i,j)=0 ! which cannot be
          enddo ! reached by water
        enddo
      enddo
C$OMP BARRIER
      end


      subroutine etch_mgz_into_land_thread(ncx,ncy, mask, mgz,ms1,ms2)
      implicit none
      integer ncx,ncy
      integer(kind=2), dimension(ncx,ncy) :: mask, mgz,ms1,ms2
      integer :: trd_count, npts, npts_bak
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     & range, istr,iend,jstr,jend, iter, my_sum, i,j, iw,ie,js,jn

      trd_count=0
      npts=0
      npts_bak=-1
C$ integer omp_get_num_threads, omp_get_thread_num
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c** nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first+range-1, nsub_x*nsub_y-1)

      if (trd == 0) then
        write(*,*) 'Etching merging zone area into land...'
        trd_count=0 ; npts=0 ; npts_bak=-1 !<-- initialize
      endif

      do tile=my_first,my_last
        call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
        do j=jstr,jend
          do i=istr,iend
            if (mask(i,j) > 0) then ! initialize etching
              if (mgz(i,j) == 0) mgz(i,j)=-1 ! procedure: -1
            endif ! -1 = water
            ms1(i,j)=mgz(i,j) ! 0 = land
            ms2(i,j)=mgz(i,j) ! +1 = merging zone
          enddo
        enddo
      enddo
C$OMP BARRIER

      iter=0
      do while(npts /= npts_bak)
        iter=iter+1 ; my_sum=0
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              ms2(i,j)=ms1(i,j)
              if (ms1(i,j) == 0) then
                if (ms1(iw,j)+ms1(ie,j)+ms1(i,js)+ms1(i,jn) < 0) then
                  if ( ms1(iw,j) <= 0 .and. ms1(ie,j) <= 0 .and.
     & ms1(i,js) <= 0 .and. ms1(i,jn) <= 0) then
                    ms2(i,j)=-1
                    my_sum=my_sum+1
                  endif
                elseif (ms1(iw,j)+ms1(ie,j)+ms1(i,js)+ms1(i,jn) > 0
     & ) then
                  if ( ms1(iw,j) >= 0 .and. ms1(ie,j) >= 0 .and.
     & ms1(i,js) >= 0 .and. ms1(i,jn) >= 0) then
                    ms2(i,j)=+1
                    my_sum=my_sum+1
                  endif
                endif
              endif
            enddo
          enddo
        enddo !<-- tile
C$OMP BARRIER

        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              ms1(i,j)=ms2(i,j)
              if (ms2(i,j) == 0) then
                if (ms2(iw,j)+ms2(ie,j)+ms2(i,js)+ms2(i,jn) < 0) then
                  if ( ms2(iw,j) <= 0 .and. ms2(ie,j) <= 0 .and.
     & ms2(i,js) <= 0 .and. ms2(i,jn) <= 0) then
                    ms1(i,j)=-1
                    my_sum=my_sum+1
                  endif
                elseif (ms2(iw,j)+ms2(ie,j)+ms2(i,js)+ms2(i,jn) > 0
     & ) then
                  if ( ms2(iw,j) >= 0 .and. ms2(ie,j) >= 0 .and.
     & ms2(i,js) >= 0 .and. ms2(i,jn) >= 0) then
                    ms1(i,j)=+1
                    my_sum=my_sum+1
                  endif
                endif
              endif
            enddo
          enddo
        enddo !<-- tile

C$OMP CRITICAL(cr_region)
        if (trd_count == 0) then
          npts_bak=npts ; npts=my_sum
        else
          npts=npts+my_sum
        endif
        trd_count=trd_count+1
        if (trd_count == numthreads) then
          trd_count=0
          if (mod(iter,20) == 0 .or. npts == npts_bak) then
            write(*,'(8x,A,I7,2(2x,A,I10))') 'iter =', iter,
     & 'npts =', npts, 'changes =', npts-npts_bak
          endif
        endif
C$OMP END CRITICAL(cr_region)
C$OMP BARRIER
      enddo !<-- while
      end




      subroutine etch_weights_into_land_thread(wdth, ncx,ncy, mgz,
     & ms1,ms2, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
      implicit none
      integer ncx,ncy, wdth
      integer(kind=2), dimension(ncx,ncy) :: mgz, ms1,ms2
      logical OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH
      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last, tile,
     & range, istr,iend,jstr,jend, iter, max_iters
C$ integer omp_get_num_threads, omp_get_thread_num
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)
c*** nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only
      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)
      if (trd == 0) then
        write(*,*) 'Etching merging weight function into land...'
      endif
      max_iters=(wdth+1)/2
      do iter=1,max_iters
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          call mgz_step_into_land_tile( istr,iend,jstr,jend,
     & ncx,ncy, mgz, ms1,ms2,
     & wdth, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
        enddo
C$OMP BARRIER
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          call mgz_step_into_land_tile( istr,iend,jstr,jend,
     & ncx,ncy, mgz, ms2,ms1,
     & wdth, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)
        enddo
C$OMP BARRIER
        if (trd == 0 .and. mod(iter,20) == 0) then
          write(*,'(8x,A,I7)') 'iter =', iter
        endif
      enddo
      end

      subroutine mgz_step_into_land_tile(istr,iend,jstr,jend, ncx,ncy,
     & mgz, ms1,ms2, wdth, OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH)

! Perform one step of etching "ms2" into allowed land area using "ms1"
! as input. Note that in the calling sequence above "ms1" and "ms2"
! alternate, so the algorithm here is quasi time stepping.
! At the beginning of this procedure "ms1" and "ms2" are initialized
! as integer constant-slope functions: starting with the maximum value
! of "width" at the open boundary and decreasing by 1 for every row of
! points as proceeding into the interior of the domain. However, this
! applies only to water points; on all the land points both "ms1" and
! "ms2" are set to zero; Condition "mgz > 0" specifies merging zone
! which already extended into some land (etched), so now "ms1" and
! "ms2" will be allowed to intrude into land areas where "mgz > 0" in
! such a way that their values will be decreasing by 1 every step
! moving away from the perimeter of the grid and from the coastline
! into the land, however, the expansion is not be allowed to go
! beyond "mgz > 0".

      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, wdth
      integer(kind=2), dimension(ncx,ncy) :: mgz, ms1,ms2
      logical OBC_WEST,OBC_EAST,OBC_SOUTH,OBC_NORTH
      integer i,j, iw,ie,js,jn, imgz,jmgz

      if (OBC_WEST) then
        imgz=min(iend, wdth)
        if (istr <= imgz) then
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,imgz
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max( ms2(i,j), ms1(i,j), ms1(ie,j),
     & max(ms1(iw,j), ms1(i,jn), ms1(i,js))-1
     & )
              endif
            enddo
          enddo
        endif
      endif
      if (OBC_SOUTH) then
        jmgz=min(jend, wdth)
        if (jstr <= jmgz) then
          do j=jstr,jmgz
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max( ms2(i,j), ms1(i,j), ms1(i,jn),
     & max(ms1(i,js), ms1(ie,j), ms1(iw,j))-1
     & )
              endif
            enddo
          enddo
        endif
      endif
      if (OBC_EAST) then
        imgz=max(istr, ncx-wdth+1)
        if (imgz <= iend) then
          do j=jstr,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=imgz,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max( ms2(i,j), ms1(i,j), ms1(iw,j),
     & max(ms1(ie,j), ms1(i,jn), ms1(i,js))-1
     & )
              endif
            enddo
          enddo
        endif
      endif
      if (OBC_NORTH) then
        jmgz=max(jstr, ncy-wdth+1)
        if (jmgz <= jend) then
          do j=jmgz,jend
            js=max(1,j-1) ; jn=min(j+1,ncy)
            do i=istr,iend
              iw=max(1,i-1) ; ie=min(i+1,ncx)
              if (mgz(i,j) > 0) then
                ms2(i,j)=max( ms2(i,j), ms1(i,j), ms1(i,js),
     & max(ms1(i,jn), ms1(ie,j), ms1(iw,j))-1
     & )
              endif
            enddo
          enddo
        endif
      endif
      end


      subroutine smooth_wgt_thread(max_iters, ncx,ncy, mgz,wgt,hwg)
      implicit none
      integer max_iters, ncx,ncy
      integer(kind=2) mgz(ncx,ncy)
      real(kind=8) wgt(ncx,ncy), hwg (ncx,ncy)

      integer numthreads,trd, nsub_x,nsub_y, my_first,my_last,
     & tile, range, istr,iend,jstr,jend, iter
C$ integer omp_get_num_threads, omp_get_thread_num
      numthreads=1 ; trd=0
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      call set_tiles(ncx,ncy, nsub_x,nsub_y)

c*** nsub_x=1 ; nsub_y=1 !<-- for testing parallel correctness only

      range=(nsub_x*nsub_y +numthreads-1)/numthreads
      my_first=trd*range
      my_last=min(my_first + range-1, nsub_x*nsub_y-1)

      if (trd == 0) then
        write(*,*) 'Smoothing weight function, iters =',
     & max_iters, '  numthreads =', numthreads
      endif

      do iter=1,max_iters
        do tile=my_first,my_last,+1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          call gauss_sm_wgt_tile(istr,iend,jstr,jend,
     & ncx,ncy, mgz, wgt,hwg)
        enddo
C$OMP BARRIER
        do tile=my_last,my_first,-1
          call comp_tile_bounds(tile, ncx,ncy, nsub_x,nsub_y,
     & istr,iend,jstr,jend)
          call smooth_wgt_tile(istr,iend,jstr,jend,
     & ncx,ncy, mgz, hwg,wgt)
        enddo
C$OMP BARRIER
        if (mod(iter,20) == 0 .and. trd == 0) then
          write(*,'(8x,A,I7)') 'iter =', iter
        endif
      enddo
      end

      subroutine smooth_wgt_tile(istr,iend,jstr,jend,
     & ncx,ncy, mgz, src,targ)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, i,j
      integer(kind=2) mgz(ncx,ncy)
      real(kind=8) src(ncx,ncy), targ(ncx,ncy)

      if (istr==1) istr=istr+1 ! isotropic smoothing operator
      if (iend==ncx) iend=iend-1 !
      if (jstr==1) jstr=jstr+1 ! 1/32 1/8 1/32
      if (jend==ncy) jend=jend-1 !
                                      ! 1/8 3/8 1/8
      do j=jstr,jend !
        do i=istr,iend ! 1/32 1/8 1/32
          if (mgz(i,j) > 0) then
            targ(i,j)=0.125D0*( 3.D0*src(i,j)+src(i-1,j)+src(i,j-1)
     & +src(i+1,j)+src(i,j+1)
     & +0.25D0*( src(i-1,j+1)+src(i+1,j+1)
     & +src(i-1,j-1)+src(i+1,j-1)
     & ) )
          else
            targ(i,j)=0.D0
          endif
        enddo
      enddo
      end


      subroutine gauss_sm_wgt_tile(istr,iend,jstr,jend,
     & ncx,ncy, mgz, src,targ)
      implicit none
      integer istr,iend,jstr,jend, ncx,ncy, i,j
      integer(kind=2) mgz(ncx,ncy)
      real(kind=8) src(ncx,ncy), targ(ncx,ncy)

      if (istr==1) istr=istr+1 ; if (jstr==1) jstr=jstr+1
      if (iend==ncx) iend=iend-1 ; if (jend==ncy) jend=jend-1

      do j=jstr,jend
        do i=istr,iend
          if (mgz(i,j) > 0) then
            targ(i,j)=0.2D0*( src(i-1,j)+src(i,j-1)
     & +src(i+1,j)+src(i,j+1)
     & +0.25D0*( src(i-1,j+1)+src(i+1,j+1)
     & +src(i-1,j-1)+src(i+1,j-1)
     & ) )
          else
            targ(i,j)=0.D0
          endif
        enddo
      enddo
      end
# 76 "tools_fort.F" 2
# 1 "r2r_interp_init.F" 1
      subroutine r2r_interp_init_thread( nx,ny,x,y, ncx,ncy, xc,yc,
     & ip,jp, xi,eta)

! Initialize interpolation between two arbitrary oriented curvilinear
! grids: given arrays of coordinates [x(nx,ny),y(nx,ny)] for points of
! a non-Cartesian "source" grid, and similar [xc(ncx,ncy),yc(ncx,ncy)]
! for the "target" grid (hereafter "c" stands for "child") find arrays
! of indices [ip(ncx,ncy),jp(ncx,ncy)] and arrays fractional distances
! [xi(ncx,ncy),eta(ncx,ncy)] such that each point of the target grid
! (xc,yc) is surrounded by 4 points of the source,
!
! [x,y](ip,jp+1) --__
! / --__ [x,y](ip+1,jp+1)
! / [xc,yc] /
! / /
! [x,y](ip,jp) -- __ /
! --__ [x,y](ip+1,jp)
!
! such that bi-linear interpolation of [x,y] into location of [xc,yc]
! yields [xc,yc] themselves,
!
! xc = (1-xi)*(1-eta)*x(i,j) + xi*(1-eta)*x(i+1,j)
! +(1-xi)* eta *x(i,j+1) + xi* eta *x(i+1,j+1)
!
! yc = (1-xi)*(1-eta)*y(i,j) + xi*(1-eta)*y(i+1,j)
! +(1-xi)* eta *y(i,j+1) + xi* eta *y(i+1,j+1)
!
! in other words, (ip+xi,jp+eta) are coordinates of point (xc,yc) in
! "continuous index" space (i,j) of the source grid.

! Normally the target grid is expected to be entirely within the area
! covered by the source, however it is not strictly required by the
! algorithm below: if some portions of the target grid are outside the
! source (hence interpolation is impossible for these locations) then
! ip=-1 and jp=-1 are set as special values.

      implicit none
      integer nx,ny, ncx,ncy
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer icmin,icmax,jcmin,jcmax
C$ integer numthreads, trd, chunk_size
C$ integer omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncy+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx ; jcmin=1 ; jcmax=ncy
C$ jcmin=1+trd*chunk_size ; jcmax=min(jcmin+chunk_size-1,ncy)

      call search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     & ncx,ncy, xc,yc, ip,jp)

      call comp_offsets_tile( icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     & ncx,ncy, xc,yc, ip,jp, xi,eta)
      end



      subroutine r2r_interp_search_thread( nx,ny,x,y, ncx,ncy, xc,yc,
     & ip,jp)
      implicit none
      integer nx,ny, ncx,ncy ! Search for parent
      real(kind=8), dimension(nx,ny) :: x,y ! grid indices only,
      real(kind=8), dimension(ncx,ncy) :: xc,yc ! do not compute
      integer(kind=4), dimension(ncx,ncy) :: ip,jp ! xi,eta offsets.
      integer icmin,icmax,jcmin,jcmax
C$ integer numthreads, trd, chunk_size
C$ integer omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncy+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx ; jcmin=1 ; jcmax=ncy
C$ jcmin=1+trd*chunk_size ; jcmax=min(jcmin+chunk_size-1,ncy)

      call search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny, x,y,
     & ncx,ncy, xc,yc, ip,jp)
      end



      subroutine search_indices_tile(icmin,icmax,jcmin,jcmax, nx,ny,
     & x,y, ncx,ncy, xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc ! 4 <---- 3
      integer(kind=4), dimension(ncx,ncy) :: ip,jp ! ! ^
      integer i,j, ic,jc, inew,jnew, iinc ! ! !
      real(kind=8) dx4,dy4, dx3,dy3, r41,r34, ! v !
     & dx1,dy1, dx2,dy2, r12,r23 ! 1 ----> 2

      write(*,*) 'enter search_tile, jcmin =',jcmin, ' jcmax =', jcmax


! The search algorithm is organized as follows:
!
! initialization: select starting point ic,jc
! and compute cyclic vector products r12...r41
!
! do while(.true.)
! do while( r12...r41 >= 0 -- meaning that point
! ic,jc is inside of [i:i+1]x[j:j+1] )
!
! record ip,jp for this point, and proceed to the next point
! of target grid using by incrementing only ONE index, either
! "ic" or "jc" -- if "ic" already reached its bound, while
! simultaneously reversing the direction of "ic" sweep -- the
! reversing is done to make sure that the next target point is
! always nearby, to maximize the probability that it is still
! inside of [i:i+1]x[j:j+1] cell, so the inner while-loop would
! not break off too often. Recompute r12...r41 for the next
! ic,jc point to check logical condition during next iteration
! of while loop; [Note that there is no attempt to change
! source-grid indices i,j inside this loop.]
!
! enddo
!
! once the while(r12...r41 >= 0) loop breaks off, it means that
! the point ic,jc is no longer inside [i:i+1]x[j:j+1], so i and/or
! j index of source grid must be incremented, depending which of
! the r12...r41 is negative. The indices i,j are incremented by
! one and are restricted to be within the range of source grid.
! This leads to two possibilities:
! either
! (i) at least one of them, i,j is incremented, then keep
! ic,jc the same and proceed with recomputing r12...r41
! to resume while(r12...r41 >= 0) loop;
! or
! (ii) both increments of i and j are canceled by the
! restriction, so i,j get "stuck" which is detected
! by inew,jnew both having the same values. This means
! that this ic,jc point cannot be bounded properly
! because its location is outside the source grid, so
! mark it by special values of ip,jp (essentially skip
! this point) and proceed to the next one by
! incrementing ip,jp;
!
! enddo !<-- while(.true.)
!
! Note (a) that not having special care for (i) leads to an infinite
! loop in the algorithm if some of the target points cannot be
! bounded; and
! (b) assuming that the target grid has finer resolution than the
! source, it is expected that most of the computing time is
! spent inside while(r12...r41 >= 0) loop, hence its body is
! minimized at the expense of somewhat awkwardness of what is
! around it.

      i=nx/2 ; j=ny/2 ; ic=icmin ; jc=jcmin ; iinc=+1

      dx1=x(i+1,j )-xc(ic,jc) ; dy1=y(i+1,j )-yc(ic,jc) ! 3 <-- 2
      dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc) ! ! !
      dx3=x(i ,j+1)-xc(ic,jc) ; dy3=y(i ,j+1)-yc(ic,jc) ! ! !
      dx4=x(i ,j )-xc(ic,jc) ; dy4=y(i ,j )-yc(ic,jc) ! 4 --> 1

      r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
      r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      do while(.true.)
        do while( r12 >= 0.D0 .and. r23 >= 0.D0
     & .and. r34 >= 0.D0 .and. r41 >= 0.D0 )

          ip(ic,jc)=i ; jp(ic,jc)=j !<-- record bounding indices

          if (iinc > 0 .and. ic < icmax) then !--> proceed to the
            ic=ic+1 ! next target point
          elseif (iinc < 0 .and. ic > icmin) then
            ic=ic-1
          else
            jc=jc+1 ; iinc=-iinc ; if (jc > jcmax) return
!#ifdef
! write(*,*) 'jc =', jc
!#endif
          endif

          dx1=x(i+1,j )-xc(ic,jc) ; dy1=y(i+1,j )-yc(ic,jc)
          dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
          dx3=x(i ,j+1)-xc(ic,jc) ; dy3=y(i ,j+1)-yc(ic,jc)
          dx4=x(i ,j )-xc(ic,jc) ; dy4=y(i ,j )-yc(ic,jc)

          r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
          r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4
        enddo

        inew=i ; jnew=j
        if (r12 < 0.D0) inew=min(inew+1, nx-1)
        if (r23 < 0.D0) jnew=min(jnew+1, ny-1)
        if (r34 < 0.D0) inew=max(inew-1, 1)
        if (r41 < 0.D0) jnew=max(jnew-1, 1)

        if (inew == i .and. jnew == j) then
          ip(ic,jc)=-1 ; jp(ic,jc)=-1 !<-- cannot be bounded

          if (iinc > 0 .and. ic < icmax) then !--> proceed to the
            ic=ic+1 ! next target point
          elseif (iinc < 0 .and. ic > icmin) then
            ic=ic-1
          else
            jc=jc+1 ; iinc=-iinc ; if (jc > jcmax) return
!#ifdef
! write(*,*) 'jc =', jc
!#endif
          endif
        else
          i=inew ; j=jnew !<-- accept move of source location.
        endif

        dx1=x(i+1,j )-xc(ic,jc) ; dy1=y(i+1,j )-yc(ic,jc)
        dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
        dx3=x(i ,j+1)-xc(ic,jc) ; dy3=y(i ,j+1)-yc(ic,jc)
        dx4=x(i ,j )-xc(ic,jc) ; dy4=y(i ,j )-yc(ic,jc)

        r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
        r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      enddo
      end



      subroutine comp_offsets_tile(icmin,icmax,jcmin,jcmax, nx,ny,
     & x,y, ncx,ncy, xc,yc, ip,jp, xi,eta)
      implicit none
      integer nx,ny, ncx,ncy, icmin,icmax,jcmin,jcmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer ic,jc, i,j, iter
      real(kind=8) a11,a12,a21,a22, det, dX,dY, p,p1,q,q1

      do jc=jcmin,jcmax
        do ic=icmin,icmax
          i=ip(ic,jc) ; j=jp(ic,jc)
          if (i > 0 .and. j > 0) then
            a11=0.5D0*(x(i+1,j+1)-x(i,j+1) +x(i+1,j)-x(i,j))
            a12=0.5D0*(x(i+1,j+1)-x(i+1,j) +x(i,j+1)-x(i,j))
            a21=0.5D0*(y(i+1,j+1)-y(i,j+1) +y(i+1,j)-y(i,j))
            a22=0.5D0*(y(i+1,j+1)-y(i+1,j) +y(i,j+1)-y(i,j))

            dX=xc(ic,jc)-0.25D0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))
            dY=yc(ic,jc)-0.25D0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))

            det=1.D0/(a11*a22-a12*a21)
            xi(ic,jc) =0.5D0 + det*(a22*dX-a12*dY)
            eta(ic,jc)=0.5D0 + det*(a11*dY-a21*dX)
          else
            xi(ic,jc)=-1.D0 ; eta(ic,jc)=-1.D0 !<-- special values
          endif
        enddo
        do iter=1,10
          do ic=icmin,icmax
            i=ip(ic,jc) ; j=jp(ic,jc)
            if (i > 0 .and. j > 0) then
              p=xi(ic,jc) ; p1=1.D0-p
              q=eta(ic,jc) ; q1=1.D0-q

              a11=q*(x(i+1,j+1)-x(i,j+1)) +q1*(x(i+1,j)-x(i,j))
              a12=p*(x(i+1,j+1)-x(i+1,j)) +p1*(x(i,j+1)-x(i,j))
              a21=q*(y(i+1,j+1)-y(i,j+1)) +q1*(y(i+1,j)-y(i,j))
              a22=p*(y(i+1,j+1)-y(i+1,j)) +p1*(y(i,j+1)-y(i,j))

              dX=xc(ic,jc) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     & -p*q1*x(i+1,j) -p1*q1*x(i,j)
              dY=yc(ic,jc) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     & -p*q1*y(i+1,j) -p1*q1*y(i,j)

              det=1.D0/(a11*a22-a12*a21)

              xi(ic,jc) =p + det*(a22*dX-a12*dY)
              eta(ic,jc)=q + det*(a11*dY-a21*dX)
            endif
          enddo
        enddo
      enddo
      end


      subroutine check_search_indices(nx,ny,x,y, ncx,ncy, xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx,ncy
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer ic,jc, i,j
      real(kind=8) dx1,dy1, dx2,dy2, dx3,dy3, dx4,dy4, r12,r23,r34,r41

      write(*,'(2x,A)',advance='no') 'check_search_indices ...'
      do jc=1,ncy
        do ic=1,ncx
          i=ip(ic,jc) ; j=jp(ic,jc)
          if (i > 0 .and. j > 0) then
            dx1=x(i+1,j )-xc(ic,jc) ; dy1=y(i+1,j )-yc(ic,jc)
            dx2=x(i+1,j+1)-xc(ic,jc) ; dy2=y(i+1,j+1)-yc(ic,jc)
            dx3=x(i ,j+1)-xc(ic,jc) ; dy3=y(i ,j+1)-yc(ic,jc)
            dx4=x(i ,j )-xc(ic,jc) ; dy4=y(i ,j )-yc(ic,jc)

            r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
            r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

            if (r12<0.D0 .or. r23<0.D0 .or. r34<0.D0 .or. r41<0.D0)
     & write(*,*) '### ERROR: Search algorithm failure at ',
     & 'ic =', ic, ' jc =', jc
          endif
        enddo
      enddo
      write(*,'(2x,A)') '...done'
      end



      subroutine check_offsets(nx,ny, x,y, ncx,ncy, xc,yc,
     & ip,jp, xi,eta)
      implicit none
      integer nx,ny, ncx,ncy
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx,ncy) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx,ncy) :: ip,jp
      integer ic,jc, i,j
      real(kind=8) dX,dY, p,p1,q,q1, errX,errY

      write(*,'(2x,A)',advance='no') 'checking offsets...'
      errX=0.D0 ; errY=0.D0
      do jc=1,ncy
        do ic=1,ncx ! Notice simple semantic rule
          i=ip(ic,jc) ; j=jp(ic,jc) ! here: one and only one "1" is
          if (i > 0 .and. j > 0) then ! always present: either as +1
            p=xi(ic,jc) ; p1=1.D0-p ! in index, or as p1,q1 in the
            q=eta(ic,jc) ; q1=1.D0-q ! corresponding coefficient.

            dX=xc(ic,jc) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     & -p*q1*x(i+1,j) -p1*q1*x(i,j)
            dY=yc(ic,jc) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     & -p*q1*y(i+1,j) -p1*q1*y(i,j)

            if (abs(dX) > errX) errX=abs(dX)
            if (abs(dY) > errY) errY=abs(dY)

c** xi(ic,jc)=dX ; eta(ic,jc)=dY
          endif
        enddo
      enddo
      write(*,'(2x,A,2ES22.15)') 'max errors =', errX,errY
      end
# 77 "tools_fort.F" 2
# 1 "r2r_bry_interp.F" 1
      subroutine bry_init_interp(nx,ny, x,y, ncx, xc,yc, ip,jp, xi,eta)
      implicit none
      integer :: nx,ny, ncx
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
C$OMP PARALLEL SHARED(nx,ny,x,y, ncx, xc,yc, ip,jp, xi,eta)
      call bry_init_interp_thread(nx,ny,x,y, ncx,xc,yc,
     & ip,jp, xi,eta)
C$OMP END PARALLEL
      call check_search_line(nx,ny, x,y, ncx, xc,yc, ip,jp)
      call check_offsts_line(nx,ny, x,y, ncx, xc,yc, ip,jp,
     & xi,eta)
      end

      subroutine bry_init_interp_thread(nx,ny,x,y, ncx, xc,yc, ip,jp,
     & xi,eta)
      implicit none
      integer :: nx,ny, ncx, icmin,icmax
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
      integer numthreads,trd, chunk_size
C$ integer omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      call bry_search_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc, ip,jp)
      call compute_offsts_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc,
     & ip,jp, xi,eta)
      end

      subroutine r2r_bry_search(nx,ny,x,y, ncx,xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc
      integer(kind=4), dimension(ncx) :: ip,jp
C$OMP PARALLEL SHARED(nx,ny,x,y, ncx,xc,yc, ip,jp)
      call r2r_bry_search_thread(nx,ny,x,y, ncx,xc,yc, ip,jp)
C$OMP END PARALLEL
      call check_search_line( nx,ny,x,y, ncx,xc,yc, ip,jp)
      end

      subroutine r2r_bry_search_thread(nx,ny,x,y, ncx,xc,yc, ip,jp)
      implicit none
      integer :: nx,ny, ncx, icmin,icmax ! Search for parent
      real(kind=8), dimension(nx,ny) :: x,y ! grid indices only,
      real(kind=8), dimension(ncx) :: xc,yc ! do not compute
      integer(kind=4), dimension(ncx) :: ip,jp ! xi,eta offsets.
      integer :: numthreads, chunk_size,trd
C$ integer :: omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      call bry_search_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc, ip,jp)
      end



      subroutine bry_search_line(icmin,icmax, nx,ny, x,y, ncx, xc,yc,
     & ip,jp)
      implicit none
      integer icmin,icmax, nx,ny, ncx, ic,i,j,inew,jnew
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) dx4,dy4, dx3,dy3, r41,r34,
     & dx1,dy1, dx2,dy2, r12,r23

      write(*,*) 'enter search_line, icmin =',icmin, ' icmax =', icmax

      i=nx/2 ; j=ny/2 ; ic=icmin ! 3 <---- 2
      dx1=x(i+1,j )-xc(ic) ; dy1=y(i+1,j )-yc(ic) ! ! ^
      dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic) ! ! !
      dx3=x(i ,j+1)-xc(ic) ; dy3=y(i ,j+1)-yc(ic) ! v !
      dx4=x(i ,j )-xc(ic) ; dy4=y(i ,j )-yc(ic) ! 4 ----> 1

      r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
      r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      do while(.true.)
        do while( r12 >= 0.D0 .and. r23 >= 0.D0
     & .and. r34 >= 0.D0 .and. r41 >= 0.D0 )

          write(*,*) 'ic =', ic

          ip(ic)=i ; jp(ic)=j !<-- record bounding indices
          ic=ic+1 !--> proceed to the next target point
          if (ic > icmax) return

          dx1=x(i+1,j )-xc(ic) ; dy1=y(i+1,j )-yc(ic)
          dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic)
          dx3=x(i ,j+1)-xc(ic) ; dy3=y(i ,j+1)-yc(ic)
          dx4=x(i ,j )-xc(ic) ; dy4=y(i ,j )-yc(ic)

          r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
          r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4
        enddo

        inew=i ; jnew=j
        if (r12 < 0.D0) inew=min(inew+1, nx-1)
        if (r23 < 0.D0) jnew=min(jnew+1, ny-1)
        if (r34 < 0.D0) inew=max(inew-1, 1)
        if (r41 < 0.D0) jnew=max(jnew-1, 1)

        if (inew == i .and. jnew == j) then

          write(*,*) 'ic =', ic

          ip(ic)=-1 ; jp(ic)=-1 !<-- cannot be bounded
          ic=ic+1
          if (ic > icmax) return
        else
          i=inew ; j=jnew !--> accept move of source location.
        endif

        dx1=x(i+1,j )-xc(ic) ; dy1=y(i+1,j )-yc(ic)
        dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic)
        dx3=x(i ,j+1)-xc(ic) ; dy3=y(i ,j+1)-yc(ic)
        dx4=x(i ,j )-xc(ic) ; dy4=y(i ,j )-yc(ic)

        r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
        r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

      enddo !<-- while(.true.)
      end

      subroutine compute_offsts_line(icmin,icmax, nx,ny, x,y, ncx,
     & xc,yc,ip,jp, xi,eta)
      implicit none
      integer icmin,icmax, nx,ny, ncx, ic, i,j, iter
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) a11,a12,a21,a22, det, dX,dY, p,p1,q,q1

      do ic=icmin,icmax
        i=ip(ic) ; j=jp(ic)
        if (i > 0 .and. j > 0) then
          a11=0.5D0*(x(i+1,j+1)-x(i,j+1) +x(i+1,j)-x(i,j))
          a12=0.5D0*(x(i+1,j+1)-x(i+1,j) +x(i,j+1)-x(i,j))
          a21=0.5D0*(y(i+1,j+1)-y(i,j+1) +y(i+1,j)-y(i,j))
          a22=0.5D0*(y(i+1,j+1)-y(i+1,j) +y(i,j+1)-y(i,j))

          dX=xc(ic)-0.25D0*(x(i+1,j+1)+x(i,j+1)+x(i+1,j)+x(i,j))
          dY=yc(ic)-0.25D0*(y(i+1,j+1)+y(i,j+1)+y(i+1,j)+y(i,j))

          det=1.D0/(a11*a22-a12*a21)
          xi(ic) =0.5D0 + det*(a22*dX-a12*dY)
          eta(ic)=0.5D0 + det*(a11*dY-a21*dX)
        else
          xi(ic)=-1.D0 ; eta(ic)=-1.D0 !<-- special values
        endif
      enddo
      do iter=1,10
        do ic=icmin,icmax
          i=ip(ic) ; j=jp(ic)
          if (i > 0 .and. j > 0) then
            p=xi(ic) ; p1=1.D0-p
            q=eta(ic) ; q1=1.D0-q

            a11=q*(x(i+1,j+1)-x(i,j+1)) +q1*(x(i+1,j)-x(i,j))
            a12=p*(x(i+1,j+1)-x(i+1,j)) +p1*(x(i,j+1)-x(i,j))
            a21=q*(y(i+1,j+1)-y(i,j+1)) +q1*(y(i+1,j)-y(i,j))
            a22=p*(y(i+1,j+1)-y(i+1,j)) +p1*(y(i,j+1)-y(i,j))

            dX=xc(ic) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     & -p*q1*x(i+1,j) -p1*q1*x(i,j)
            dY=yc(ic) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     & -p*q1*y(i+1,j) -p1*q1*y(i,j)

            det=1.D0/(a11*a22-a12*a21)

            xi(ic) =p + det*(a22*dX-a12*dY)
            eta(ic)=q + det*(a11*dY-a21*dX)
          endif
        enddo
      enddo
      end

      subroutine check_search_line(nx,ny,x,y, ncx, xc,yc, ip,jp)
      implicit none
      integer nx,ny, ncx, ic, i,j
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) dx4,dy4, dx3,dy3, r41,r34,
     & dx1,dy1, dx2,dy2, r12,r23

      write(*,'(2x,A)',advance='no') 'check_search_line ...'
      do ic=1,ncx
        i=ip(ic) ; j=jp(ic)
        if (i > 0 .and. j > 0) then
          dx1=x(i+1,j )-xc(ic) ; dy1=y(i+1,j )-yc(ic)
          dx2=x(i+1,j+1)-xc(ic) ; dy2=y(i+1,j+1)-yc(ic)
          dx3=x(i ,j+1)-xc(ic) ; dy3=y(i ,j+1)-yc(ic)
          dx4=x(i ,j )-xc(ic) ; dy4=y(i ,j )-yc(ic)

          r12=dx1*dy2-dx2*dy1 ; r23=dx2*dy3-dx3*dy2
          r34=dx3*dy4-dx4*dy3 ; r41=dx4*dy1-dx1*dy4

          if (r12<0.D0 .or. r23<0.D0 .or. r34<0.D0 .or. r41<0.D0)
     & write(*,*) '### ERROR: Search algorithm failure at ',
     & 'ic =', ic
        endif
      enddo
      write(*,'(2x,A)') '...done'
      end

      subroutine check_offsts_line(nx,ny, x,y, ncx, xc,yc, ip,jp,
     & xi,eta)
      implicit none
      integer nx,ny, ncx, ic, i,j
      real(kind=8), dimension(nx,ny) :: x,y
      real(kind=8), dimension(ncx) :: xc,yc, xi,eta
      integer(kind=4), dimension(ncx) :: ip,jp
      real(kind=8) dX,dY, p,p1,q,q1, errX,errY

      write(*,'(2x,A)',advance='no') 'checking offsets...'
      errX=0.D0 ; errY=0.D0
      do ic=1,ncx ! Notice simple semantic rule
        i=ip(ic) ; j=jp(ic) ! here: one and only one "1" is
        if (i > 0 .and. j > 0) then ! always present: either as +1
          p=xi(ic) ; p1=1.D0-p ! in index, or as p1,q1 in the
          q=eta(ic) ; q1=1.D0-q ! corresponding coefficient.

          dX=xc(ic) -p*q*x(i+1,j+1) -p1*q*x(i,j+1)
     & -p*q1*x(i+1,j) -p1*q1*x(i,j)
          dY=yc(ic) -p*q*y(i+1,j+1) -p1*q*y(i,j+1)
     & -p*q1*y(i+1,j) -p1*q1*y(i,j)

          if (abs(dX) > errX) errX=abs(dX)
          if (abs(dY) > errY) errY=abs(dY)
        endif
      enddo
      write(*,'(2x,A,2ES22.15)') 'max errors =', errX,errY
      end

      subroutine bry_interp(nx,ny,N, src, ncx, ip,jp,xi,eta, msk, targ)
      implicit none
      integer nx,ny,N, ncx
      integer(kind=4) ip(ncx),jp(ncx)
      real(kind=8) xi(ncx), eta(ncx)
      integer(kind=2) msk(ncx)
      real(kind=4) src(nx,ny,N), targ(ncx,N)
C$OMP PARALLEL SHARED(nx,ny,N, src, ncx, ip,jp,xi,eta, msk, targ)
      call bry_interp_thread(nx,ny,N, src, ncx, ip,jp,xi,eta,msk, targ)
C$OMP END PARALLEL
      end

      subroutine bry_interp_thread(nx,ny,N, src, ncx, ip,jp,xi,eta,
     & msk,targ)
      implicit none
      integer :: nx,ny,N, ncx
      integer(kind=4) :: ip(ncx), jp(ncx)
      real(kind=8) :: xi(ncx), eta(ncx)
      integer(kind=2) :: msk(ncx)
      real(kind=4) :: src(nx,ny,N), targ(ncx,N)
      integer :: icmin,icmax,isize, istr,iend,tile
      integer :: numthreads, trd, chunk_size
C$ integer :: omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call bry_interp_line( istr,iend, nx,ny,N, src, ncx, ip,jp,
     & xi,eta, msk, targ)
      enddo
      end

! subroutine bry_interp_line(istr,iend, nx,ny,N, src, ncx, ip,jp,
! & xi,eta, msk, targ)
! implicit none
! integer istr,iend, nx,ny,N, ncx, ic, i,j,k
! integer(kind=4) ip(ncx), jp(ncx)
! real(kind=8) xi(ncx), eta(ncx), p,p1,q,q1
! integer(kind=2) msk(ncx)
! real(kind=4) src(nx,ny,N), targ(ncx,N)
!
! do ic=istr,iend
! if (ip(ic) > 0) then
! i=ip(ic) ; p=xi(ic) ; p1=1.D0-p
! j=jp(ic) ; q=eta(ic) ; q1=1.D0-q
! do k=1,N
! targ(ic,k)=p*q*src(i+1,j+1,k) +p1*q*src(i,j+1,k)
! & +p*q1*src(i+1,j,k) +p1*q1*src(i,j,k)
! enddo
! else
! do k=1,N
! targ(ic,k)=0.
! enddo
! endif
! enddo
! end


      subroutine bry_interp_line(istr,iend, nx,ny,N, src, ncx, ip,jp,
     & xi,eta, msk, targ)
      implicit none
      integer istr,iend, nx,ny,N, ncx, ic, i,j,k
      integer(kind=4) ip(ncx), jp(ncx)
      real(kind=8) xi(ncx), eta(ncx)
      integer(kind=2) msk(ncx)
      real(kind=4) src(nx,ny,N), targ(ncx,N)
      real(kind=8), parameter :: TwoThird=2.D0/3.D0,
     & FourNineth=4.D0/9.D0
      real(kind=8) px,qx,pqx, HxL,HxR,GxL,GxR, FxLL,FxRL,FxLR,FxRR,
     & py,qy,pqy, HyL,HyR,GyL,GyR, FyLL,FyRL,FyLR,FyRR,
     & FxyLL,FxyRL,FxyLR,FxyRR

      do ic=istr,iend
        i=ip(ic) ; j=jp(ic) ;
        if ( msk(ic) == 1 .and. 2 < i .and. i < nx-2 .and.
     & 2 < j .and. j < ny-2 ) then
          px=xi(ic) ; qx=1.D0-px; pqx=px*qx
          GxR=-px*pqx ; HxR=px*(px+2.D0*pqx)
          GxL= qx*pqx ; HxL=qx*(qx+2.D0*pqx)

          py=eta(ic); qy=1.D0-py; pqy=py*qy
          GyR=-py*pqy ; HyR=py*(py+2.D0*pqy)
          GyL= qy*pqy ; HyL=qy*(qy+2.D0*pqy)

          do k=1,N
            FxLL=TwoThird*( src(i+1,j,k)-src(i-1,j,k)
     & -0.125D0*(src(i+2,j,k)-src(i-2,j,k)))
            FxRL=TwoThird*( src(i+2,j,k)-src(i ,j,k)
     & -0.125D0*(src(i+3,j,k)-src(i-1,j,k)))
            FxLR=TwoThird*( src(i+1,j+1,k)-src(i-1,j+1,k)
     & -0.125D0*(src(i+2,j+1,k)-src(i-2,j+1,k)))
            FxRR=TwoThird*( src(i+2,j+1,k)-src(i ,j+1,k)
     & -0.125D0*(src(i+3,j+1,k)-src(i-1,j+1,k)))

            FyLL=TwoThird*( src(i,j+1,k)-src(i,j-1,k)
     & -0.125D0*(src(i,j+2,k)-src(i,j-2,k)))
            FyRL=TwoThird*( src(i+1,j+1,k)-src(i+1,j-1,k)
     & -0.125D0*(src(i+1,j+2,k)-src(i+1,j-2,k)))
            FyLR=TwoThird*( src(i,j+2,k)-src(i,j ,k)
     & -0.125D0*(src(i,j+3,k)-src(i,j-1,k)))
            FyRR=TwoThird*( src(i+1,j+2,k)-src(i+1,j ,k)
     & -0.125D0*(src(i+1,j+3,k)-src(i+1,j-1,k)))

            FxyLL=FourNineth*( -src(i-1,j+1,k) +src(i+1,j+1,k)
     & +src(i-1,j-1,k) -src(i+1,j-1,k)
     & +0.125D0*( src(i-1,j+2,k) -src(i+1,j+2,k)
     & +src(i-2,j+1,k) -src(i+2,j+1,k)
     & -src(i-2,j-1,k) +src(i+2,j-1,k)
     & -src(i-1,j-2,k) +src(i+1,j-2,k)
     & +0.125D0*( -src(i-2,j+2,k) +src(i+2,j+2,k)
     & +src(i-2,j-2,k) -src(i+2,j-2,k)
     & )))
            FxyRL=FourNineth*( -src(i ,j+1,k) +src(i+2,j+1,k)
     & +src(i ,j-1,k) -src(i+2,j-1,k)
     & +0.125D0*( src(i ,j+2,k) -src(i+2,j+2,k)
     & +src(i-1,j+1,k) -src(i+3,j+1,k)
     & -src(i-1,j-1,k) +src(i+3,j-1,k)
     & -src(i ,j-2,k) +src(i+2,j-2,k)
     & +0.125D0*( -src(i-1,j+2,k) +src(i+3,j+2,k)
     & +src(i-1,j-2,k) -src(i+3,j-2,k)
     & )))
            FxyLR=FourNineth*( -src(i-1,j+2,k) +src(i+1,j+2,k)
     & +src(i-1,j ,k) -src(i+1,j ,k)
     & +0.125D0*( src(i-1,j+3,k) -src(i+1,j+3,k)
     & +src(i-2,j+2,k) -src(i+2,j+2,k)
     & -src(i-2,j ,k) +src(i+2,j ,k)
     & -src(i-1,j-1,k) +src(i+1,j-1,k)
     & +0.125D0*( -src(i-2,j+3,k) +src(i+2,j+3,k)
     & +src(i-2,j-1,k) -src(i+2,j-1,k)
     & )))
            FxyRR=FourNineth*( -src(i ,j+2,k) +src(i+2,j+2,k)
     & +src(i ,j ,k) -src(i+2,j ,k)
     & +0.125D0*( src(i ,j+3,k) -src(i+2,j+3,k)
     & +src(i-1,j+2,k) -src(i+3,j+2,k)
     & -src(i-1,j ,k) +src(i+3,j ,k)
     & -src(i ,j-1,k) +src(i+2,j-1,k)
     & +0.125D0*( -src(i-1,j+3,k) +src(i+3,j+3,k)
     & +src(i-1,j-1,k) -src(i+3,j-1,k)
     & )))

            targ(ic,k) = HxL*HyL*src(i,j ,k) + HxR*HyL*src(i+1,j ,k)
     & + HxL*HyR*src(i,j+1,k) + HxR*HyR*src(i+1,j+1,k)

     & + GxL*HyL*FxLL + GxR*HyL*FxRL
     & + GxL*HyR*FxLR + GxR*HyR*FxRR

     & + HxL*GyL*FyLL + HxR*GyL*FyRL
     & + HxL*GyR*FyLR + HxR*GyR*FyRR

     & + GxL*GyL*FxyLL + GxR*GyL*FxyRL
     & + GxL*GyR*FxyLR + GxR*GyR*FxyRR
          enddo
        else
          do k=1,N
            targ(ic,k)=0. !<-- no special value masking here yet.
          enddo
        endif
      enddo
      end



! An alternative driver for "r2r_init_vrtint_tile" (r2r_vert_interp.F)
! designed for 1D line as needed by r2r_bry. The differences are:
! (i) the second horizontal dimension no longer exists; (ii) instead of
! using 2D tiling, it cuts "ncx" into a set of chunks according to the
! number of theads; (iii) external layer with parallel region is added;
! (iv) there is no option mode.


      subroutine bry_init_vertinterp(ncx, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
      implicit none
      integer ncx, Np,N
      real(kind=8) hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     & kprnt(ncx,N)
C$OMP PARALLEL SHARED(ncx,hprnt, Np,hcp,Csp_r,h,N,hc,Cs_r,kprnt)
      call bry_init_vinterp_thread(ncx, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
C$OMP END PARALLEL
      end

      subroutine bry_init_vinterp_thread(ncx, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer :: ncx, Np,N
      real(kind=8) :: hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     & kprnt(ncx,N)
      integer :: icmin,icmax,isize, istr,iend,tile
      integer :: numthreads, chunk_size,trd
      integer :: omp_get_num_threads, omp_get_thread_num
      numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
      chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
      icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      if (alloc_zc_size < isize*(N+1) .or.
     & allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) deallocate(zp_r,drv,zc)
        allocate( zc(alloc_zc_size), zp_r(allc_zpr_size),
     & drv(allc_zpr_size) )
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)') 'bry_init_vinterp_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     & 'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
      endif
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call r2r_init_vrtint_tile( istr,iend, 1, 1, ncx,1,
     & hprnt,Np,hcp,Csp_r, zp_r,drv,
     & h,N, hc, Cs_r, zc, kprnt)
      enddo
      end

! Same, but for checking routine.

      subroutine bry_check_init_vertinterp(ncx, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
      implicit none
      integer ncx, Np,N
      real(kind=8) hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     & kprnt(ncx,N)
C$OMP PARALLEL SHARED(ncx,hprnt, Np,hcp,Csp_r,h,N,hc,Cs_r,kprnt)
      call bry_check_vrtint_thread(ncx, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
C$OMP END PARALLEL
      end

      subroutine bry_check_vrtint_thread(ncx, hprnt, Np,hcp,Csp_r,
     & h, N,hc,Cs_r, kprnt)
      use r2r_vertint_vars
      implicit none
      integer :: ncx, Np,N
      real(kind=8) :: hprnt(ncx),h(ncx), hcp,hc, Csp_r(Np),Cs_r(N),
     & kprnt(ncx,N), my_error
      integer :: icmin,icmax,isize, istr,iend,tile
      integer :: numthreads, trd, chunk_size
C$ integer :: omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      if (alloc_zc_size < isize*(N+1) .or.
     & allc_zpr_size < isize*(Np+2)) then
        alloc_zc_size=isize*(N+1); allc_zpr_size=isize*(Np+2)
        if (allocated(zp_r)) then
            deallocate(zp_r,drv,zc)
        endif
        allocate( zc(alloc_zc_size), zp_r(allc_zpr_size),
     & drv(allc_zpr_size) )
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)') 'bry_check_vrtint_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     & 'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
      endif
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call r2r_check_vrtint_tile( istr,iend, 1,1, ncx,1,
     & hprnt, Np,hcp,Csp_r,zp_r,drv,
     & h, N,hc,Cs_r,zc, kprnt, my_error)
      enddo
C$OMP CRITICAL(r2r_vert_crgn)
      if (trd_count == 0) then
          vert_int_error=0.D0
      endif
      trd_count=trd_count+1
      vert_int_error=max(vert_int_error, my_error)
C$ if (trd_count == numthreads) then
        trd_count=0
        write(*,*) '          maximum vert_int_error =', vert_int_error
C$ endif
C$OMP END CRITICAL(r2r_vert_crgn)
      end






      subroutine bry_vertinterp(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
      implicit none
      integer :: ncx, btm_bc, Np,N
      real(kind=4) :: qsrc(ncx,Np), qtr(ncx,N)
      real(kind=8) :: kprnt(ncx,N)
C$OMP PARALLEL SHARED(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
      call bry_vinterp_thread(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
C$OMP END PARALLEL
      end

      subroutine bry_vinterp_thread(ncx, btm_bc, Np,qsrc, N,kprnt,qtr)
      use r2r_vertint_vars !<-- needed for "drv"
      implicit none
      integer :: ncx, btm_bc, Np,N
      real(kind=4) :: qsrc(ncx,Np), qtr(ncx,N)
      real(kind=8) :: kprnt(ncx,N)
      integer, parameter :: lmsk=0 !<-- suppress land masking
      integer(kind=2) :: mask(8) !<-- for compatibility; not used
      integer :: icmin,icmax,isize, istr,iend,tile
      integer :: numthreads, trd, chunk_size
C$ integer :: omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1,icmax)

        call r2r_vsplnint_tile( istr,iend, 1,1, ncx,1,
     & lmsk,mask, btm_bc, Np,qsrc,drv, N,kprnt,qtr)




      enddo
      end



      subroutine bry_vert_average(ncx, h, N,hc,Cs_w, qsrc,qbar)
      implicit none
      integer :: ncx,N
      real(kind=8) :: h(ncx), hc, Cs_w(0:N)
      real(kind=4) :: qsrc(ncx,N),qbar(ncx)
C$OMP PARALLEL SHARED(ncx,N, h, hc,Cs_w, qsrc,qbar)
      call bry_vert_average_thread(ncx, h, N, hc,Cs_w, qsrc,qbar)
C$OMP END PARALLEL
      end

      subroutine bry_vert_average_thread(ncx, h, N,hc,Cs_w, qsrc,qbar)
      use r2r_vertint_vars
      implicit none
      integer :: ncx,N
      real(kind=8) :: h(ncx), hc, Cs_w(0:N)
      real(kind=4) :: qsrc(ncx,N), qbar(ncx)
      integer :: icmin,icmax,isize, istr,iend,tile
      integer :: numthreads, trd, chunk_size
C$ integer :: omp_get_num_threads, omp_get_thread_num
C$ numthreads=omp_get_num_threads() ; trd=omp_get_thread_num()
C$ chunk_size=(ncx+numthreads-1)/numthreads
      icmin=1 ; icmax=ncx
C$ icmin=1+trd*chunk_size ; icmax=min(icmin+chunk_size-1,ncx)
      isize=(icmax-icmin+2)/2
      if (alloc_zc_size < isize*(N+1)) then
        alloc_zc_size=isize*(N+1)
        if (allocated(zc)) then
            deallocate(zc)
        endif
        allocate(zc(alloc_zc_size))
C$OMP CRITICAL(r2r_vert_crgn)
        write(*,'(1x,2A,F8.4,1x,A,I3)') 'bry_vert_average_thread: ',
     & 'allocated', dble(2*allc_zpr_size+alloc_zc_size)/dble(262144),
     & 'MB private workspace, trd =', trd
C$OMP END CRITICAL(r2r_vert_crgn)
      endif
      do tile=0,1
        istr=icmin+tile*isize ; iend=min(istr+isize-1, icmax)
        call bry_vert_average_line(istr,iend, ncx,h, N,hc,Cs_w,
     & qsrc,qbar, zc)
      enddo
      end

      subroutine bry_vert_average_line(istr,iend, ncx, h, N,hc,Cs_w,
     & qsrc,qbar, z_w)
      implicit none
      integer :: istr,iend, ncx,N, i,k
      real(kind=8) :: h(ncx), hc, Cs_w(0:N), z_w(istr:iend,0:N)
      real(kind=4) :: qsrc(ncx,N),qbar(ncx)

      call r2r_set_depth_tile(istr,iend,1,1, ncx,1,h, N+1,hc,Cs_w,z_w)
      do i=istr,iend
        qbar(i)=qsrc(i,N)*(z_w(i,N)-z_w(i,N-1))
      enddo
      do k=N-1,1,-1
        do i=istr,iend
          qbar(i)=qbar(i)+qsrc(i,k)*(z_w(i,k)-z_w(i,k-1))
        enddo
      enddo
      do i=istr,iend
        qbar(i)=qbar(i)/(z_w(i,N)-z_w(i,0))
      enddo
      end
# 78 "tools_fort.F" 2
!!!!!!!!!!
