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
module BsaLib_Functions

   use BsaLib_CONSTANTS
   use BsaLib_Data, only: wd, struct_data, settings, dimM_bisp_, dimM_psd_
   implicit none (type, external)
   public
   private :: wd, struct_data, settings, dimM_bisp_, dimM_psd_

   ! make a local internal copy
   integer :: NFREQS, NNODES, NNODESL, NLIBS, NLIBSL
   integer :: NMODES, NMODES_EFF
   integer :: NPSDEL, NTCOMPS, NDIRS = 1
   integer, allocatable :: MODES(:)
   integer, allocatable :: TCOMPS(:), DIRS(:)

   integer              :: MSHR_SVD_LWORK = - 1
   integer, allocatable  :: MSHR_SVD_INFO        
   real, allocatable :: MSHR_SVD_WORK(:)


   interface

      module subroutine setBsaFunctionLocalVars()
      end subroutine


      module subroutine prefetchSVDWorkDim_()
      end subroutine

      module subroutine cleanSVDWorkInfo_()
      end subroutine



      module subroutine getFM_full_tnm_scalar_msh_(bfm, fi, fj)
         real, intent(inout), contiguous :: bfm(:, :)
         real, intent(in), contiguous :: fi(:), fj(:)
      end subroutine


      module subroutine getFM_full_tm_scalar_msh_POD_(bfm, fi, fj)
         real, intent(inout), contiguous :: bfm(:, :)
         real, intent(in), contiguous :: fi(:), fj(:)
      end subroutine


      module subroutine getRM_full_scalar_msh_(brm, fi, fj, bfm)
         real, intent(inout), contiguous :: brm(:, :)
         real, intent(in), contiguous :: fi(:), fj(:)
         real, intent(in), contiguous :: bfm(:, :)
      end subroutine



      module subroutine getFM_diag_tnm_scalar_msh_(bfm, fi, fj)
         real, intent(inout), contiguous :: bfm(:, :)
         real, intent(in), contiguous :: fi(:), fj(:)
      end subroutine


      module subroutine getRM_diag_scalar_msh_(brm, fi, fj, bfm)
         real, intent(inout), contiguous :: brm(:, :)
         real, intent(in), contiguous :: fi(:), fj(:)
         real, intent(in), contiguous :: bfm(:, :)
      end subroutine









      module subroutine getFM_full_tnm_vect_cls_(f, Suvw, psd, bisp)
         real, intent(in) :: f(NFREQS)
         real, intent(in) :: Suvw(NFREQS, NPSDEL)
         real, allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine



      module subroutine getRM_full_vect_cls_(f, psd, bisp)
         real, intent(in)                 :: f(NFREQS)
         real, allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine



      module subroutine getFM_diag_tnlm_vect_cls_(f, Suvw, psd, bisp)
         real, intent(in) :: f(NFREQS)
         real, intent(in) :: Suvw(NFREQS, NPSDEL)
         real, intent(inout), allocatable :: psd(:, :), bisp(:, :, :)
      end subroutine   



      module subroutine getRM_diag_vect_cls_(f, psd, bisp)
         real, intent(in)                 :: f(NFREQS)
         real, allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine







      pure module subroutine getFM_full_tnm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         integer, intent(in)  :: ii, ij
         real, intent(in)    :: fi, fj
         real, intent(in)    :: Suvw(NFREQS, NPSDEL)
         real, intent(in)    :: Suvw_pad(NPSDEL)
         real, intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine



      module subroutine getRM_full_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         integer, intent(in) :: ii, ij
         real, intent(in)   :: fi, fj
         real, intent(in)   :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real, intent(out)  :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine


      !> BUG: this routine is adapted to the case where we use
      !>      convention on PULSATION.
      !>      Please, adapt it to the case of convention over FREQUENCIES.
      pure module subroutine getFM_diag_tnlm_scalar_cls_(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         integer, intent(in)  :: ii, ij
         real, intent(in)    :: fi, fj
         real, intent(in)    :: Suvw(NFREQS, NPSDEL)
         real, intent(in)    :: Suvw_pad(NPSDEL)
         real, intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine



      module subroutine getRM_diag_scalar_cls_(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         integer, intent(in) :: ii, ij  ! freqs indexes
         real, intent(in)   :: fi, fj
         real, intent(in)   :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real, intent(out)  :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine




      pure module subroutine getBR_SFm_val_(nm, Suvw, fnat, im, m, psd)
         integer, intent(in)  :: im, m, nm
         real, intent(in)    :: Suvw(nm, NPSDEL), fnat
         real, intent(inout) :: psd
      end subroutine

   end interface

end module BsaLib_Functions
