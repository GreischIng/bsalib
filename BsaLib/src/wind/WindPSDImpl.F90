!! This file is part of BsaLib.
!! Copyright (C) 2024  Michele Esposito Marzino 
!!
!! BsaLib is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! BsaLib is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with BsaLib.  If not, see <https://www.gnu.org/licenses/>.
submodule(BsaLib_WindData) BsaLib_WindPSDImpl

#ifndef BSA_DEBUG
# define __use_concurrent_loops__
#endif

   implicit none (type, external)

   ! type :: arr_proc_pointer_t
   !    procedure(PSDfunc), pointer, nopass :: ptr => null()
   ! end type arr_proc_pointer_t

   ! ! TODO: might be statically initialised (parameter)
   ! type(arr_proc_pointer_t), dimension(5) :: psd_funcs

   procedure(PSDfunc), pointer  :: psd_func_internal_ => null()

contains


   module subroutine SetPSDType(this, ipsd)
      class(WindData_t), intent(inout) :: this
      integer, value :: ipsd
      integer :: istat
      character(len = 256) :: emsg

      if (ipsd < 1 .or. ipsd > 5) then
         print '(1x, 2a, i0, a)', &
            WARNMSG, 'Invalid psd ID ', ipsd, '. Setting default (5).'
         ipsd = 5
      endif

      ! NOTE: need to store it internally..
      this%i_psd_type_ = ipsd

      allocate(this%psd_, stat=istat, errmsg=emsg)
      if (istat /= 0) call allocKOMsg('this % psd_', istat, emsg)

      select case (ipsd)
      case (1)
         psd_func_internal_ => vonKarmanPSD_
      case (2)
         psd_func_internal_ => kaimalPSD_
      case (3)
         psd_func_internal_ => davenportPSD_Greisch_
      case (5)
         psd_func_internal_ => davenportPSD_Uliege_
      end select

      call this%psd_%SetPSDFunction(psd_func_internal_)

! #ifdef BSA_DEBUG
!       print *, INFOMSG, '@WindImpl::SetPSDType() : PSD type set to ', this%i_psd_type_
! #endif
   end subroutine





   module function getFullNodalPSD(this, innl, nodesl, PSDvec, f, idir) result(PSDmat)
      use BsaLib_Utility,  only: util_getCorrVectIndex
      use BsaLib_Data,     only: struct_data
      class(WindData_t),  intent(in) :: this
      integer, intent(in) :: innl, idir
      integer, intent(in) :: nodesl(innl)
      real,   intent(in) :: PSDvec(innl)
      real,   intent(in) :: f
      real :: PSDmat(innl, innl)
      real :: absf
      integer   :: i, j, ni, nj, id

      absf = abs(f)

      do i = 1, innl
         PSDmat(:, i) = sqrt(PSDvec * PSDvec(i))
         ni = nodesl(i)
         do j = 1, innl
            nj           = nodesl(j)
            id           = util_getCorrVectIndex(nj, ni, struct_data%nn_)
            PSDmat(j, i) = PSDmat(j, i) * (this%nod_corr_(id, idir)**(absf))
         enddo
      enddo
   end function




   module subroutine SetPSDFunction(this, func)
      class(psd_t) :: this
      procedure(PSDfunc), intent(in), pointer :: func

      this%psd_fct_ptr => func
   end subroutine




   module function evalPSD_(this, nf, f, innl, nnl, idir, itc) result(PSD)
      use BsaLib_Data, only: settings
      class(WindData_t),  intent(in) :: this
      integer, intent(in) :: nf, innl, idir, itc
      real,   intent(in) :: f(:)
      integer, intent(in) :: nnl(:)
      real :: PSD(nf, innl)

      if (idir /= 1) then
         print '(/ 1x, 2a, i1, a)', &
            WARNMSG, 'IDIR=  ', idir, ', when usually SHOULD be 1 (X wind direction).'
         print '(1x, 2a/, 2a)', &
            MSGCONT, '(It is uncommon to compute PSDs of wind turbulence assuming', &
            MSGCONT, ' that vortices do not move along X (idir=1), direction of mean wind)'
      endif

      ! invoking internal function pointer
      PSD = this%psd_%psd_fct_ptr(this, nf, f, innl, nnl, idir, itc)

      ! NOTE: PSD scaling (convention based)
      if (settings%i_def_scaling_ == 1) PSD = PSD / CST_PIt4
   end function




   function vonKarmanPSD_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t),  intent(in) :: wd
      integer, intent(in) :: nf       ! n. frequencies
      integer, intent(in) :: innl     ! n. actual nodes loaded
      integer, intent(in) :: idir     ! wind direction
      integer, intent(in) :: itc      ! turb. component id
      real,   intent(in) :: freqs(:) ! frequencies
      integer, intent(in) :: nnl(:)   ! list of actual loaded nodes
      real :: PSD(nf, innl)
      ! local
      real, dimension(1, innl) :: L
      real, allocatable :: rtmp1(:, :)
      integer :: i


      ! NOTE: how was programmed in FINELG
      ! BUG: check.
      if (itc == 1) then

         L(1, :) = wd%turb_scales_wz_(1, idir, wd%wz_node_(nnl))

         ! L/U
         rtmp1 = L / reshape(wd%u_node_(nnl), [1, innl])
         PSD   = matmul(reshape(abs(freqs), [nf, 1]), rtmp1)
         PSD   = PSD * PSD ! square
         rtmp1 = rtmp1 * reshape(wd%sigmaUVW_wz_(itc, wd%wz_node_(nnl))**2, [1, innl])
         PSD   = (1 + 70.7 * PSD)**(5./6.)


         do i = 1, innl
            PSD(:, i) = 4. * rtmp1(1, i) / PSD(:, i)
         enddo

      else

         block
            real :: dnlsu(nf, innl), rtmp2(nf, innl), rtmp3(nf, innl)

            dnlsu = 2. * &
               matmul(reshape(freqs, [nf, 1]), &
                  reshape(wd%turb_scales_wz_(itc, idir, wd%wz_node_(nnl)), [1, innl]) / &
                  reshape(wd%u_node_(nnl), [1, innl]) &
               )

            dnlsu = dnlsu*dnlsu

            rtmp1 = 1. + 70.7 * dnlsu
            rtmp2 = rtmp1 ** (11. / 6.)
            do i = 1, innl
               rtmp2(:, i) = rtmp2(:, i) * wd%u_node_(nnl(i))
               rtmp3(:, i) = wd%turb_scales_wz_(itc, idir, wd%wz_node_(nnl(i))) * &
                  (1. + 188.4 * dnlsu(:, i)) / &
                  rtmp2(:, i) * wd%sigmaUVW_wz_(itc, wd%wz_node_(nnl(i)))**2
            enddo
            rtmp1 = rtmp3 + rtmp3
            PSD   = rtmp1 + rtmp1
         end block
      endif

      if (allocated(rtmp1)) deallocate(rtmp1)
   end function vonKarmanPSD_




   function kaimalPSD_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t),  intent(in) :: wd
      integer, intent(in) :: nf            ! n. frequencies
      integer, intent(in) :: innl          ! n. actual nodes loaded
      integer, intent(in) :: idir          ! wind direction
      integer, intent(in) :: itc           ! 
      integer, intent(in) :: nnl(:)     ! list of actual loaded nodes
      real,   intent(in) :: freqs(:)     ! frequencies
      real :: PSD(nf, innl)


      PSD = 0.
   end function kaimalPSD_





   function davenportPSD_Greisch_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t),  intent(in) :: wd
      integer, intent(in) :: nf          ! n. frequencies
      integer, intent(in) :: innl        ! n. actual nodes loaded
      integer, intent(in) :: idir        ! wind direction
      integer, intent(in) :: itc         ! 
      integer, intent(in) :: nnl(:)   ! list of actual loaded nodes
      real,   intent(in) :: freqs(:)   ! frequencies
      real :: PSD(nf, innl)

      real, parameter :: cst1 = 0.65 * 1200.
      integer   :: i, n

      do i = 1, innl
         n         = nnl(i)
         PSD(:, i) = &
            cst1 * (wd%sigmaUVW_wz_(itc, wd%wz_node_(n)))**2 / wd%u_mean_ref_wz_(wd%wz_node_(n)) / &
            (1 + (cst1 * freqs / wd%u_mean_ref_wz_(wd%wz_node_(n)))**2.)**(5. / 6.)
      enddo
   end function davenportPSD_Greisch_





   function davenportPSD_Uliege_(wd, nf, freqs, innl, nnl, idir, itc) result(PSD)
      class(WindData_t),  intent(in) :: wd
      integer, intent(in) :: nf                    ! n. frequencies
      integer, intent(in) :: innl                  ! n. actual nodes loaded
      integer, intent(in) :: idir                  ! wind direction
      integer, intent(in) :: itc                   ! 
      integer, intent(in) :: nnl(:)     ! list of actual loaded nodes
      real,   intent(in) :: freqs(:)   ! frequencies
      real :: PSD(nf, innl)

      real :: cstL_U(1, innl), cstFL_U(nf, innl)
      integer   :: i, n


      cstL_U(1, :) = wd%turb_scales_wz_(itc, idir, wd%wz_node_(nnl)) / wd%u_node_(nnl)
      cstFL_U      = matmul(reshape(abs(freqs), [nf, 1]), cstL_U)

      do i = 1, innl
         n = nnl(i)

         PSD(:, i) = &
            (CST_2d3 * cstFL_U(:, i) * cstL_U(1, i) * wd%sigmaUVW_wz_(itc, wd%wz_node_(n))**2) / &
            ((1.d0 + (cstFL_U(:, i)**2))**(4.d0 / 3.d0))
      enddo
   end function davenportPSD_Uliege_




end submodule
