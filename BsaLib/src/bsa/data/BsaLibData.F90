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
module BsaLib_Data

   use BsaLib_Timing
   use BsaLib_CONSTANTS
   use BsaLib_Settings
   use BsaLib_Structure
   use BsaLib_WindData
   !$ use omp_lib
#ifdef BSA_USE_GPU
   use BsaCL
#endif
   implicit none (type, external)
   public


   interface
      module function evaluatePSD(f, nf, itc) result(psd)
         integer, intent(in) :: nf, itc
         real, intent(in)   :: f(nf)
         real, allocatable, target :: psd(:, :)
      end function

      module subroutine cleanBSAData_()
      end subroutine

      module subroutine bsa_Abort(emsg)
         character(len = *), intent(in), optional :: emsg
      end subroutine
   end interface


! **********************************************************************
!   GLOBAL control data
! **********************************************************************

   type(settings_t),      allocatable, target :: settings
   type(WindData_t),      allocatable, target :: wd
   type(StructureData_t), allocatable, target :: struct_data


   !> If <code>.true.</code> reads BSA input parameters from
   !> <code>BsaLib</code> formatted file <code>bsa.bsadata</code>.
   !> This can be enabled by calling <code>bsa_readInputParamsFromBSAFile()</code>.
   !> If the file does not exist, an error is raised.
   logical :: do_read_input_internal_ = .false.

   !> Tracks if we have cleaned BSA data or not.
   logical :: is_data_cleaned_   = .false.

   !> Controls whether we need to close debug unit or not.
   !> This is for those cases where the debug unit is provided
   !> by the hosting unit, and so we do not own it.
   logical :: close_deb_unit_    = .true.

   !> Controls if we need to validate modal data or not.
   logical :: do_validate_modal_ = .true.

   integer :: dimNf_psd_ = 0, dimNf_bisp_ = 0
   integer :: dimNr_psd_ = 0, dimNr_bisp_ = 0
   integer :: dimM_psd_  = 0, dimM_bisp_  = 0

   real, allocatable, target :: PHItimesC_local_(:, :, :)

   !> If true, in Run, only generates compatible files for BSA executable
   logical :: do_gen_bsa_input_files_ = .false.

   !> If false, does not run BsaLib.
   logical :: do_run_bsalib_ = .true.

   !> Keeps track of max n. of OMP threads. 1 (main) thread by default.
   integer :: max_num_omp_threads_ = 1


! **********************************************************************
!   Exporting control data
! **********************************************************************

   !> If true, generates the "modal.txt" formatted file
   !> containing modal data
   logical :: do_export_modal_data_ = .false.

   ! BUG: unused!
   ! TODO: if visual for modal bisp, we don't need to compute the full set.
   !> Modal dimensions in Post phase.
   integer :: dimM_psd_post_, dimM_bisp_post_

   !> If true, post-mesh for generating Bispectrums.
   logical :: is_visual_ = .false.

   !> Stores visual indexes (both modal and nodal)
   integer, target :: visual_indexes_(3) = 0

   ! BUG: this only valid for nodal visual mode (make it more general)
   !> If visual, allows to export Nodal Bispectrum instead of modal
   logical :: is_brn_export_ = .false.

   !> Tracks I/O handle to which exporting bispectrum data
   integer :: bisp_export_iun_internal_ = 0

   !> Final offset indexing value in visual mode
   integer, target :: visual_idx_ = 1

   !> If true, exports all BRM data to file (NOT visual!)
   logical :: do_export_brm_ = .false.

   integer :: i_brmexport_mode_ = BSA_EXPORT_BRM_MODE_BASE

   !> Base exporting data type
   type, public :: BsaExportBaseData_t

      integer :: i_doNotPrintGenHeader_ = 0   ! == 0  means DO PRINT !!
      integer :: ncomb_  = 0
      integer :: ispsym_ = 0
      integer :: nzones_ = 0
      integer :: nm_     = 0
      integer, pointer :: modes_(:) => null()

      integer :: idZone_ = 0
      integer :: nI_     = 0
      integer :: nJ_     = 0
   end type

   !> Function pointer to Bispectrum exporting procedure
   procedure(exportInterf_vect_), pointer :: write_brm_fptr_  => null()
   type(BsaExportBaseData_t), target      :: export_data_base_

   !> If true, means we need to export, using base internal version.
   logical :: do_export_base_ = .false.

   procedure(getBRN_), pointer :: getBRN => null()
   abstract interface
      function getBRN_(brm)
         real, intent(in) :: brm(:)
         real :: getBRN_
      end function
   end interface



! **********************************************************************
!   GPU control data
! **********************************************************************

   logical :: is_gpu_enabled_ = .false.
#ifdef BSA_USE_GPU
   integer, target :: ierr_cl_
#endif



! **********************************************************************
!   CLASSIC control data
! **********************************************************************

   logical :: force_cls_execution_ = .false.
   integer(int64), parameter :: MAX_VECT_ALLOC_ELEMS = 1000000000 ! 1B -> almost 8Gb
   integer :: ifr = 0, jfr = 0
   ! real(RDP), pointer :: m2mf_cls_ptr_(:), m2mr_cls_ptr_(:)   ! 2nd order moments
   ! real(RDP), pointer :: m3mf_cls_ptr_(:), m3mr_cls_ptr_(:)   ! 3rd order moments

   procedure(getBFMClsVect), pointer :: getBFM_vect_cls => null()
   procedure(getBRMClsVect), pointer :: getBRM_vect_cls => null()
   abstract interface
      subroutine getBFMClsVect(f, Suvw, psd, bisp)
         import :: settings, struct_data, wd
         real, intent(in) :: f(settings%nfreqs_)
         real, intent(in) :: Suvw(settings%nfreqs_, struct_data%nn_load_ * wd%i_ndirs_ * wd%i_ntc_)
         real, allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine

      subroutine getBRMClsVect(f, psd, bisp)
         import :: settings
         real, intent(in)                 :: f(settings%nfreqs_)
         real, allocatable, intent(inout) :: psd(:, :), bisp(:, :, :)
      end subroutine
   end interface


   procedure(getBFMClsScalar), pointer :: getBFM_scalar_cls => null()
   procedure(getBRMClsScalar), pointer :: getBRM_scalar_cls => null()
   abstract interface
      pure subroutine getBFMClsScalar(ii, ij, fi, fj, Suvw, Suvw_pad, psd, bisp)
         import :: dimM_psd_, dimM_bisp_
         import :: settings, struct_data, wd
         integer, intent(in)  :: ii, ij
         real, intent(in)    :: fi, fj
         real, intent(in)    :: Suvw(settings%nfreqs_, struct_data%nn_load_ * wd%i_ndirs_ * wd%i_ntc_)
         real, intent(in)    :: Suvw_pad(struct_data%nn_load_ * wd%i_ndirs_ * wd%i_ntc_)
         real, intent(inout) :: psd(dimM_psd_), bisp(dimM_bisp_)
      end subroutine

      subroutine getBRMClsScalar(ii, ij, fi, fj, psdin, psdout, bispin, bispout)
         import :: dimM_psd_, dimM_bisp_
         integer, intent(in) :: ii, ij
         real, intent(in)   :: fi, fj
         real, intent(in)   :: psdin(dimM_psd_), bispin(dimM_bisp_)
         real, intent(out)  :: psdout(dimM_psd_), bispout(dimM_bisp_)
      end subroutine
   end interface




! **********************************************************************
!   MESHER control data
! **********************************************************************

   real, pointer :: m3mf_msh_ptr_(:) => null(), m3mr_msh_ptr_(:) => null()

   !> If true, stops at pre-mesh phase.
   logical :: is_only_premesh_ = .false.

   !> Premeshing scheme type
   integer, public :: ipre_mesh_type = BSA_PREMESH_TYPE_DIAG_CREST_NO

   !> Premeshing scheme mode
   integer, public :: ipre_mesh_mode = BSA_PREMESH_MODE_ZONE_REFINED

   integer, public :: msh_iZone

   !> Tracks zone with max N. of points
   integer, public :: msh_max_zone_NPts = 0

   !> Controls if checking zone's deltas or not.
   logical :: do_validate_deltas_ = .true.

   !> If true, restricts back-peak computation to only first column
   logical :: do_restrict_bkgpeak_ = .false.

   !> Stores width of background peak
   real :: bkg_peakw_ = 0.

   !> Array of resonant peak extensions (widths)
   real, allocatable :: peak_exts_(:)

   ! Total pre-mesh/post-mesh phase points
   integer, public :: msh_bfmpts_pre_
   integer, public :: msh_bfmpts_post_
   integer, public :: msh_brmpts_post_

   !> Code "-1" means that interest modes have to be read from the NEXT (right)
   !> limit only, since this is the first limit in the list
   integer, public, parameter :: CODE_PRE_PEAK_OK = -1

   !> Code "-2" means that interest modes have to be read from 
   !> NEXT (right) limit only, since this is the first limit in the list.
   !> However, info is MISSING from the left side (BKG) peak
   !> in which some modes fall within, but we cannot know which
   !> one is close to this zone's limit, in order to determine if
   !> to be added to its interest modes.
   !> Hence, it is a sort of "ALARM". Info will not be accurate in this case.
   integer, public, parameter :: CODE_PRE_PEAK_KO = -2

   !> Limit zones interest modes indexes
   integer, public, allocatable :: msh_ZoneLimsInterestModes(:)

   !> Tot n. of zones counter.
   integer, public, target :: msh_NZones = 0

   !> Flag to signal whether POD CACHING has been used
#ifdef BSA_USE_POD_DATA_CACHING
   integer, parameter :: POD_CACHING_FLAG = 1
#else
   integer, parameter :: POD_CACHING_FLAG = 0
#endif

   !> Controls whether employing new BFM MLR method or not
   logical :: test_no_bfm_mlr_ = .false.

   !> Controls whether to perform modal truncation or not
   logical          :: do_trunc_POD_  = .false.

   !> POD truncation limit
   real :: POD_trunc_lim_ = 0.

   !> If true, keeps specified n. of POD modes
   logical          :: nPODmodes_set_ = .false.

   !> N. of POD modes to be kept
   integer :: nmodes_POD_  = 0_int32

   !> If .true., exports POD truncation info to a file.
   logical :: do_export_POD_info_ = .false.
   logical, allocatable :: do_export_POD_trunc_(:)    !<-- BUG: this is because of OMP.
   integer,     parameter :: iun_POD_trunc_ = 659_int32
   character(len = *), parameter :: iun_POD_trunc_fname_ = 'POD_trunc_info.txt'

   ! Mesher function pointer (pre/post meshing)
   procedure(getMshBFM), pointer :: getBFM_msh => null()
   procedure(getMshBRM), pointer :: getBRM_msh => null()
   abstract interface
      subroutine getMshBFM(bfm, fi, fj)
         real, intent(inout), contiguous :: bfm(:, :)
         real, intent(in), contiguous    :: fi(:), fj(:)
      end subroutine

      subroutine getMshBRM(brm, fi, fj, bfm)
         real, intent(inout), contiguous :: brm(:, :)
         real, intent(in), contiguous :: fi(:), fj(:)
         real, intent(in), contiguous :: bfm(:, :)
      end subroutine
   end interface

end module BsaLib_Data
