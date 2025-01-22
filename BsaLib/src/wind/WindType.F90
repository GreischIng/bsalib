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
module BsaLib_WindData

   use BsaLib_CONSTANTS, only: INFOMSG, WARNMSG, ERRMSG, MSGCONT, DBGMSG
   use BsaLib_Settings
   use BsaLib_IO
   implicit none (type, external)
   private


   !**************************************************************************************
   !   WIND LABELs
   !**************************************************************************************
   character(len = *), dimension(6), public, parameter :: CST_PSD_TYPES = [&
        'VON_KARMAN        ' &
      , 'KAIMAL            ' &
      , 'DAVENPORT         ' &
      , 'EUROCODE          ' &
      , 'DAVENPORT_REFORM  ' &
      , 'ORNSTEIN_UHLENBECK' &
   &]



   !> NOTE: some wind profile definitions come from real measurement campaigns 
   !>       conducted by Greisch Design Office (Liège, Belgium) for the 
   !>       design of the Millau Bridge (Millau, France).
   character(len = *), dimension(5), public, parameter :: CST_WIND_V_PROFILES = [&
        'POWER     ' &
      , 'LOGARITHM ' &
      , 'MILLAU    ' &
      , 'MILLAU MAQ' &
      , 'GLOB. POW ' &
   &]


   !> Vertical wind direction axis labels.
   character(len = *), dimension(3), public, parameter :: CST_WIND_VERT_DIRS = ['X', 'Y', 'Z']


   ! ----- private variables
   real, public              :: air_dens_ = 1.225
   real, public, allocatable :: rot_W2G_(:, :)  ! rotation from GRS to WRS (global)


   type, public :: psd_t
      procedure(PSDfunc), nopass, private, pointer :: psd_fct_ptr => null()
   contains
      procedure, public, pass :: SetPSDFunction
   end type psd_t



   type, public :: WindData_t

      !> N. of wind zones
      integer :: nz_ = 0


      !> PSD type (VonKarman, Davenport, Kaimal, etc..).
      integer :: i_psd_type_ = 5

      !> If true, using same nodal wind speed everywhere.
      integer :: i_eq_nod_wind_speed_ = 0

      !> Wind vertical profile (Power, Log).
      integer :: i_wind_prof_ = 2

      !> Specifies which of the three is the main
      !> structural direction (for wind altidude!!)
      integer :: i_vert_ = 3




      !> Number of considered turbulent components
      !> By default ==1, i.e. only u(t) considered.
      integer :: i_ntc_ = 1

      !> List of considered wind turbulent components
      !> Range from 1 (u) to 3 (w).
      integer, allocatable :: tc_(:)

      !> How many spatial directions are considered
      !> By default ==1, only direction parallel to mean wind.
      integer :: i_ndirs_ = 1

      !> List of considered spatial directions.
      !> Range from 1 (u) to 3 (w).
      integer, allocatable :: dirs_(:)






      !> Actual nodal wind speeds (i.e. dipending on altitude, and wind zone)
      real, pointer :: u_node_(:) => null()

      !> Wind loaded nodes' wind zone
      integer, pointer :: wz_node_(:) => null()

      !> Nodal Wind Altitudes (needed for some PSDs computation)
      real, pointer :: wAlt_node_(:) => null()

      !> Nodal spatial coherence 
      !!
      !! NOTE: dimensions [(NN**2 + NN)/2, 3-dirs])
      !!       Hence, reconstruct 2D sym matrix if needed.
      real, pointer :: nod_corr_(:, :) => null()






      !> Wind forces coefficients (nlib_l, ndegw+3, nnodes_l)
      real, pointer :: wfc_(:, :, :) => null()

      !> Wind forces coefficients projected onto modal base (pre-computed)
      !> NOTE: dimensions (NM, NNL, NDEGW)
      real, pointer :: phi_times_A_ndegw_(:, :, :) => null()






      !> u(t), v(t), w(t) std deviations [m] (per wind zone)
      real, pointer :: sigmaUVW_wz_(:, :) => null()

      !> Lu_xyz, Lv_xyz, Lw_xyz wind turbulence scales [m] (per wind zone)
      real, pointer :: turb_scales_wz_(:, :, :) => null()

      !> Spatial turbulence correlation coefficients
      real, pointer :: corrCoeffs_wz_(:, :, :) => null()

      !> Spatial turbulence correlation exponents
      real, pointer :: corrExp_wz_(:, :, :) => null()

      !> reference altitude [m] (per wind zone)
      real, pointer :: Zref_wz_(:) => null()

      !> Mean Wind Speed at Zref (ref altitude) [m/s] (per wind zone)
      real, pointer :: u_mean_ref_wz_(:) => null()

      !> Rotation matrix from LWRS (Local Element) to GRS
      real, pointer :: rot_LW2G_wz_(:, :, :) => null()

      !> Wind zones limits 
      !> (might be along X, Y, Z, depending on relative orientation)
      real, pointer :: limits_wz_(:) => null()

      !> Mean wind incidence angle (per wind zone)
      real, pointer :: incAng_wz_(:) => null()


      !> internal PSD type variable
      type(psd_t), allocatable, private :: psd_

   contains

      procedure, public, nopass :: SetAirDensity
      procedure, public, nopass :: SetGlobalW2G
      procedure, public, pass :: SetWindVertProf
      procedure, public, pass :: SetPSDType
      procedure, public, pass :: SetMainVertDir
      procedure, public, pass :: SetNodalVel
      procedure, public, pass :: SetNodalWindZones
      procedure, public, pass :: SetSpatialNodalCorr
      procedure, public, pass :: getFull2DNodCorrMat
      procedure, public, pass :: SetWindFCoeffs
      procedure, public, pass :: SetPhitimesC
      procedure, public, pass :: SetNodalWindAltitudes
      procedure, public, pass :: SetWindZoneLimits
      procedure, public, pass :: SetWZMeanWindVel
      procedure, public, pass :: SetWZRefAlt
      procedure, public, pass :: SetTurbWindScales
      procedure, public, pass :: SetTurbWindSDT
      procedure, public, pass :: SetWindCorrCoeffs
      procedure, public, pass :: SetWindCorrExpnts
      procedure, public, pass :: SetLocalRotMatW2G
      procedure, public, pass :: SetIncidenceAngles

      procedure, public, pass :: setTurbComps
      procedure, public, pass :: setWindDirections
      procedure, public, pass :: SetTurbCompsAndDirsDefault

      procedure, public, pass :: evalPSD => evalPSD_

      procedure, public, pass :: getFullNodalPSD

      procedure, public, pass :: clean
   end type WindData_t




   abstract interface
      function PSDfunc(wd, nf, freqs, innl, nodes_loaded, idir, itc) result(PSD)
         import :: WindData_t
         class(WindData_t),  intent(in) :: wd
         integer, intent(in) :: nf          ! n. frequencies
         integer, intent(in) :: innl        ! n. actual nodes loaded
         integer, intent(in) :: idir        ! wind direction
         integer, intent(in) :: itc         ! turb component id
         real,   intent(in) :: freqs(:)           ! frequencies
         integer, intent(in) :: nodes_loaded(:)     ! list of actual loaded nodes
         real, dimension(nf, innl) :: PSD
      end function
   end interface



   interface

      module subroutine SetAirDensity(aird)
         real, value :: aird
      end subroutine


      module subroutine SetGlobalW2G(mat)
         real, intent(in) :: mat(3, 3)
      end subroutine


      module subroutine SetWindVertProf(this, iwprof)
         class(WindData_t), intent(inout) :: this
         integer, value        :: iwprof
      end subroutine


      module subroutine SetPSDType(this, ipsd)
         class(WindData_t), intent(inout) :: this
         integer, value        :: ipsd
      end subroutine


      module subroutine SetMainvertDir(this, ivert)
         class(WindData_t), intent(inout) :: this
         integer, value        :: ivert
      end subroutine


      module subroutine SetWindZoneLimits(this, lim, ilim)
         class(WindData_t), intent(inout) :: this
#if  ((defined(__INTEL_COMPILER_BUILD_DATE)) && (__INTEL_COMPILER_BUILD_DATE >= 20221019))
         real, intent(in), target     :: lim(..)
         integer, intent(in), optional :: ilim(..)
#else
         real, intent(in), target     :: lim(:)
         integer, intent(in), optional :: ilim(:)   ! limits' index passed
#endif
      end subroutine


      module subroutine SetWZMeanWindVel(this, UBref)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: UBref(this%nz_)
      end subroutine


      module subroutine SetWZRefAlt(this, Zref)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: Zref(this%nz_)
      end subroutine


      module subroutine SetTurbWindScales(this, L)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: L(3, 3, this%nz_)
      end subroutine


      module subroutine SetTurbWindSDT(this, wtstd)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: wtstd(3, this%nz_)
      end subroutine

      module subroutine SetWindCorrCoeffs(this, corrcoeff)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: corrcoeff(3, 3, this%nz_)
      end subroutine

      module subroutine SetWindCorrExpnts(this, correxpn)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: correxpn(3, 3, this%nz_)
      end subroutine

      module subroutine SetIncidenceAngles(this, incang)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: incang(this%nz_)
      end subroutine


      module subroutine SetLocalRotMatW2G(this, rotW2G_L)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: rotW2G_L(3, 3, this%nz_)
      end subroutine



      module subroutine setWindDirections(this, dirs, ndirs)
         class(WindData_t) :: this
         integer, intent(in) :: dirs(:)
         integer, value, optional :: ndirs
      end subroutine


      module subroutine setTurbComps(this, tc, ntc)
         class(WindData_t) :: this
         integer, intent(in)      :: tc(:)
         integer, value, optional :: ntc
      end subroutine



      module subroutine SetTurbCompsAndDirsDefault(this)
         class(WindData_t), intent(inout) :: this
      end subroutine


      module subroutine SetNodalVel(this, Unod)
         class(WindData_t), intent(inout)     :: this
         real, target, intent(in) :: Unod(:)
      end subroutine


      module subroutine SetNodalWindZones(this, NodWZ)
         class(WindData_t), intent(inout) :: this
         integer, target, intent(in) :: NodWZ(:)
      end subroutine


      module subroutine SetNodalWindAltitudes(this, WnodAlt)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: WnodAlt(:)
      end subroutine 


      module subroutine SetSpatialNodalCorr(this, nodCorr)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: nodCorr(:, :)
      end subroutine


      module subroutine getFull2DNodCorrMat(this, nn, nodcorr2d)
         class(WindData_t), intent(in) :: this
         integer, value     :: nn
         real, allocatable, intent(inout) :: nodcorr2d(:, :)
      end subroutine


      module subroutine SetWindFCoeffs(this, wfc)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: wfc(:, :, :)
      end subroutine


      module subroutine SetPhitimesC(this, phiTc)
         class(WindData_t), intent(inout) :: this
         real, target, intent(in) :: phiTc(:, :, :)
      end subroutine




      module subroutine SetPSDFunction(this, func)
         class(psd_t) :: this
         procedure(PSDfunc), intent(in), pointer :: func
      end subroutine





      module function evalPSD_(this, nf, f, innl, nnl, idir, itc) result(PSD)
         class(WindData_t),  intent(in) :: this
         integer, intent(in) :: nf, innl, idir, itc
         integer, intent(in) :: nnl(:)
         real,   intent(in) :: f(:)
         real :: PSD(nf, innl)
      end function



      module function getFullNodalPSD(this, innl, nodesl, PSDvec, f, idir) result(PSDmat)
         class(WindData_t),  intent(in) :: this
         integer, intent(in) :: innl, idir
         integer, intent(in) :: nodesl(innl)
         real,   intent(in) :: PSDvec(innl)
         real,   intent(in) :: f
         real :: PSDmat(innl, innl)
      end function



      module subroutine clean(this)
         class(WindData_t) :: this
      end subroutine

   end interface


end module
