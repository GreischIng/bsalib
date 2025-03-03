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

!**************************************************************************************
!   BSA  integer/real types selected kinds
!**************************************************************************************

   integer, parameter :: bsa_int_t  = int32
   integer, parameter :: bsa_real_t = real64
   
!**************************************************************************************
!   BSA  GENERICS
!**************************************************************************************

   integer, parameter :: BSA_SPATIAL_SYM_NONE = 0
   integer, parameter :: BSA_SPATIAL_SYM_HALF = 2
   integer, parameter :: BSA_SPATIAL_SYM_FOUR = 4

   integer, parameter :: BSA_PREMESH_TYPE_DIAG_CREST_NO  = 0
   integer, parameter :: BSA_PREMESH_TYPE_DIAG_CREST_YES = 1

   integer, parameter :: BSA_PREMESH_MODE_BASE         = 0
   integer, parameter :: BSA_PREMESH_MODE_ZONE_REFINED = 1

   integer, parameter :: BSA_CLASSIC_MODE_VECTOR = 0
   integer, parameter :: BSA_CLASSIC_MODE_SCALAR = 1

   integer, parameter :: BSA_PSD_CONVENTION_FREQ = 0
   integer, parameter :: BSA_PSD_CONVENTION_PULS = 1

!**************************************************************************************
!   BSA  I/O  DEFAULTS
!**************************************************************************************

   character(len = *), parameter :: BSA_OUT_DIRNAME_DEFAULT          = '.\bsaresults\'
   character(len = *), parameter :: BSA_OUT_FILENAME_PREFIX_DEFAULT_ = 'bsaout_def'

   character(len = *), parameter :: BSA_DATA_FNAME = "bsa.bsadata"
   character(len = *), parameter :: EXT_DATA_FNAME = "bsa.extdata"

#ifdef BSA_USE_GPU
# ifdef BSA_USE_CUDA
   character(len = *), parameter :: BSA_FILE_NAME_CL_SUFFIX = '_CUDA'
# else
   character(len = *), parameter :: BSA_FILE_NAME_CL_SUFFIX = '_CL'
# endif
#else
   character(len = 0), parameter :: BSA_FILE_NAME_CL_SUFFIX = ''
#endif



!**************************************************************************************
!   WIND
!**************************************************************************************

   integer, parameter :: BSA_WIND_VERT_PROFILE_POWER = 1
   integer, parameter :: BSA_WIND_VERT_PROFILE_LOG   = 2

   integer, parameter :: BSA_WIND_PSD_VONKARMAN = 1
   integer, parameter :: BSA_WIND_PSD_KAIMAL    = 2
   integer, parameter :: BSA_WIND_PSD_DAVENPORT = 5
   ! integer, parameter :: BSA_WIND_PSD_DAVENPORT_FIXED = 3


!**************************************************************************************
!   LOGGING TAGs
!**************************************************************************************

   character(len = *), parameter :: MSGCONT = '             '
   character(len = *), parameter :: INFOMSG = '  --[info]   '
   character(len = *), parameter :: NOTEMSG = '  --[note]   '
   character(len = *), parameter :: WARNMSG = '  --[warn]   '
   character(len = *), parameter :: ERRMSG  = '  --[error]  '
   character(len = *), parameter :: DBGMSG  = '  --[debug]  '

   character(len = *), parameter :: CONSOLE_CR_SEQ = ' '


!**************************************************************************
!  I/O CONSTANTs
!**************************************************************************

   character(len = *), parameter :: IO_ACCESS_DIRECT = 'DIRECT'
   character(len = *), parameter :: IO_ACCESS_SEQUEN = 'SEQUENTIAL'
   character(len = *), parameter :: IO_ACCESS_STREAM = 'STREAM'
   character(len = *), parameter :: IO_ACCESS_APPEND = 'APPEND'

   character(len = *), parameter :: IO_ACTION_WRITE     = 'WRITE'
   character(len = *), parameter :: IO_ACTION_READ      = 'READ'
   character(len = *), parameter :: IO_ACTION_READWRITE = 'READWRITE'

   character(len = *), parameter :: IO_ASYNC_YES = 'YES'
   character(len = *), parameter :: IO_ASYNC_NO  = 'NO'

   character(len = *), parameter :: IO_FORM_FORMATTED   = 'FORMATTED'
   character(len = *), parameter :: IO_FORM_UNFORMATTED = 'UNFORMATTED'
   character(len = *), parameter :: IO_FORM_BINARY      = 'BINARY'

   character(len = *), parameter :: IO_POSITION_ASIS   = 'ASIS'
   character(len = *), parameter :: IO_POSITION_REWIND = 'REWIND'
   character(len = *), parameter :: IO_POSITION_APPEND = 'APPEND'

   character(len = *), parameter :: IO_STATUS_OLD     = 'OLD'
   character(len = *), parameter :: IO_STATUS_NEW     = 'NEW'
   character(len = *), parameter :: IO_STATUS_SCRATCH = 'SCRATCH'
   character(len = *), parameter :: IO_STATUS_REPLACE = 'REPLACE'
   character(len = *), parameter :: IO_STATUS_UNKNOWN = 'UNKNOWN'



!**************************************************************************
!  EXPORT CONSTANTs
!**************************************************************************

   integer, parameter :: BSA_EXPORT_FORMAT_FORMATTED   = 0
   integer, parameter :: BSA_EXPORT_FORMAT_UNFORMATTED = 1
   integer, parameter :: BSA_EXPORT_MODE_REPLACE = 0
   integer, parameter :: BSA_EXPORT_MODE_APPEND  = 1

   integer, parameter :: BSA_EXPORT_BRM_MODE_NONE = 0
   integer, parameter :: BSA_EXPORT_BRM_MODE_BASE = 1
   integer, parameter :: BSA_EXPORT_BRM_MODE_USR  = 9


   character(len = *), parameter :: BSA_EXPORT_M2MF_CLS_FNAME   = "m2mf_cls"
   character(len = *), parameter :: BSA_EXPORT_M2MR_CLS_FNAME   = "m2mr_cls"
   character(len = *), parameter :: BSA_EXPORT_M2O2MR_CLS_FNAME = "m2o2mr_cls"
   character(len = *), parameter :: BSA_EXPORT_M3MF_CLS_FNAME   = "m3mf_cls"
   character(len = *), parameter :: BSA_EXPORT_M3MR_CLS_FNAME   = "m3mr_cls"
   character(len = *), parameter :: BSA_EXPORT_M2MF_MSH_FNAME   = "m2mf_msh"
   character(len = *), parameter :: BSA_EXPORT_M2MR_MSH_FNAME   = "m2mr_msh"
   character(len = *), parameter :: BSA_EXPORT_M2O2MR_MSH_FNAME = "m2o2mr_msh"
   character(len = *), parameter :: BSA_EXPORT_M3MF_MSH_FNAME   = "m3mf_msh"
   character(len = *), parameter :: BSA_EXPORT_M3MR_MSH_FNAME   = "m3mr_msh"

   character(len = *), parameter :: BRM_EXPORT_FNAME_CLS = 'bsaexport_cls.brm'
   character(len = *), parameter :: BRN_EXPORT_FNAME_CLS = 'bsaexport_cls.brn'
   character(len = *), parameter :: BRM_EXPORT_FNAME_MSH = 'bsaexport_msh.brm'
   character(len = *), parameter :: BRN_EXPORT_FNAME_MSH = 'bsaexport_msh.brn'


   abstract interface
      subroutine exportInterf_vect_(f1, f2, brm, pdata)
         real, intent(in)  :: f1(:)
            !! Array of frequencies along the X-axis
         real, intent(in)  :: f2(:)
            !! Array of frequencies along the Y-axis
         real, intent(in)  :: brm(:, :)
            !! Data of bispectra
         class(*), pointer, intent(in) :: pdata
            !! Unlimited polymorphic object allowinf the user to pass any kind of 
            !! object, holding the necessary data to be used when backfiring the 
            !! provided callback function when exporting bispectra information.
      end subroutine
   end interface




!**************************************************************************************
!   NUMERICs
!**************************************************************************************

   !> @note To avoid assertion errors in Mesher due to `sin()`/`cos()` intrinsic functions 
   !>  rounding errors.
   real, parameter :: MACHINE_PRECISION = 1e-12

   real, parameter :: CST_PIGREC = 4.0 * atan(1.0)

   real, parameter :: CST_PIt2   = CST_PIGREC * 2.0
   real, parameter :: CST_PIt4   = CST_PIGREC * 4.0

   real, parameter :: CST_PId2   = CST_PIGREC / 2.0
   real, parameter :: CST_PId4   = CST_PIGREC / 4.0


   real, parameter :: CST_2d3    = 2.0 / 3.0
   real, parameter :: CST_3d2    = 3.0 / 2.0
   real, parameter :: CST_PIt3d2 = CST_PIGREC * CST_3d2

