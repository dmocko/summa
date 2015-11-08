! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE data_struc
 ! used to define model data structures
 USE nrtype
 implicit none
 private
 ! ***********************************************************************************************************
 ! Define the model decisions
 ! ***********************************************************************************************************
 ! the model decision structure
 type,public  :: model_options
  character(len=64)                      :: cOption
  character(len=64)                      :: cDecision
  integer(i4b)                           :: iDecision
 end type model_options
 type(model_options),pointer,save,public :: model_decisions(:)      ! the decision structure
 ! ***********************************************************************************************************
 ! Define metadata for model forcing datafile
 ! ***********************************************************************************************************
 ! define a derived type for the data in the file
 !type,public  :: file_info
 ! character(len=256)                     :: filenmDesc='notPopulatedYet' ! name of file that describes the data
 ! character(len=256)                     :: filenmData='notPopulatedYet' ! name of data file
 ! integer(i4b)                           :: ncols                    ! number of columns in the file
 ! integer(i4b)                           :: ixFirstHRU               ! index of the first HRU to share the same data
 ! integer(i4b),pointer                   :: time_ix(:) => null()     ! column index for each time variable
 ! integer(i4b),pointer                   :: data_ix(:) => null()     ! column index for each forcing data variable
 !end type file_info
 ! and save all the data in a single data structure
 ! NOTE: vector (HRU dimension)
 !type(file_info),pointer,save,public     :: forcFileInfo(:) => null()   ! file info for model forcing data
 ! ***********************************************************************************************************
 ! Define metadata on model parameters
 ! ***********************************************************************************************************
 ! define a data type to store model parameter information
 type,public  :: par_info
  real(dp)                               :: default_val              ! default parameter value
  real(dp)                               :: lower_limit              ! lower bound
  real(dp)                               :: upper_limit              ! upper bound
 endtype par_info
 ! define a vector, with a separate element for each parameter (variable)
 type(par_info),pointer,save,public      :: localParFallback(:) => null() ! local column default parameters
 type(par_info),pointer,save,public      :: basinParFallback(:) => null() ! basin-average default parameters
 ! ***********************************************************************************************************
 ! Define variable metadata
 ! ***********************************************************************************************************
 ! define derived type for model variables, including name, decription, and units
 type,public :: var_info
  character(len=64)                      :: varname=''               ! variable name
  CHARACTER(len=128)                     :: vardesc=''               ! variable description
  character(len=64)                      :: varunit=''               ! variable units
  character(len=32)                      :: vartype=''               ! variable type (scalar, model layers, etc.)
  logical(lgt)                           :: v_write=.TRUE.          ! flag to write variable to the output file
 endtype var_info

 ! define arrays of metadata
 type(var_info),pointer,save,public      :: time_meta(:) => null()   ! model time information
 type(var_info),pointer,save,public      :: forc_meta(:) => null()   ! model forcing data
 type(var_info),pointer,save,public      :: attr_meta(:) => null()   ! local attributes
 type(var_info),pointer,save,public      :: type_meta(:) => null()   ! local classification of veg, soil, etc.
 type(var_info),pointer,save,public      :: mpar_meta(:) => null()   ! local model parameters for each HRU
 type(var_info),pointer,save,public      :: mvar_meta(:) => null()   ! local model variables for each HRU
 type(var_info),pointer,save,public      :: indx_meta(:) => null()   ! local model indices for each HRU
 type(var_info),pointer,save,public      :: bpar_meta(:) => null()   ! basin parameters for aggregated processes
 type(var_info),pointer,save,public      :: bvar_meta(:) => null()   ! basin parameters for aggregated processes
 ! ***********************************************************************************************************
 ! Define hierarchal derived data types
 ! ***********************************************************************************************************
 ! define named variables to describe the layer type
 integer(i4b),parameter,public      :: ix_soil=1001             ! named variable to denote a soil layer
 integer(i4b),parameter,public      :: ix_snow=1002             ! named variable to denote a snow layer
 integer(i4b),parameter,public      :: ix_mixd=1003             ! named variable to denote a mixed layer

 ! gru data structure
 type, public :: gru
  integer(i4b)                      :: gru_ix                   ! index of the gru
  integer(i4b)                      :: hruCount                 ! total number of hrus in the gru
  type(hru), pointer                :: hru(:)                   ! hrus within the gru
 endtype gru
 
 ! hru data structure
 type, public :: hru
  integer(i4b)                      :: hru_ix                   ! index of the hru in the entire domain
  integer(i4b)                      :: hru_id                   ! id (non-sequential number) of the hru
 endtype hru

 ! used to relate sequencial index and ids and others
 type, public :: mapping
  integer(i4b)                      :: gru_ix                   ! index of gru which the hru belongs to 
  integer(i4b)                      :: ihru                     ! index of a hru within a gru
 endtype mapping

 ! define derived types to hold multivariate data for a single variable (different variables have different length)
 ! NOTE: use derived types here to facilitate adding the "variable" dimension
 ! ** double precision type
 type, public :: dlength
  real(dp),pointer                       :: dat(:) => null()
 endtype dlength

 ! ** integer type
 type, public :: ilength
  integer(i4b),pointer                   :: dat(:) => null()
 endtype ilength

 ! define derived types to hold data for multiple variables
 ! NOTE: use derived types here to facilitate adding extra dimensions (e.g., spatial)
 ! ** double precision type of variable length
 type, public :: var_dlength
  type(dlength),pointer                  :: var(:) => null()
 endtype var_dlength

 ! ** integer type of variable length
 type, public :: var_ilength
  type(ilength),pointer                  :: var(:) => null()
 endtype var_ilength

 ! ** double precision type of fixed length
 type, public :: var_d
  real(dp),pointer                       :: var(:) => null()
 endtype var_d

 ! ** integer type of variable length
 type, public :: var_i
  integer(i4b),pointer                   :: var(:) => null()
 endtype var_i


 ! **hru with double precision type of variable length
 type, public :: hru_dlength
  type(var_dlength),pointer              :: hru(:) => null()
 endtype hru_dlength

 ! **hru with integer type of variable length
 type, public :: hru_ilength
  type(var_ilength),pointer              :: hru(:) => null()
 endtype hru_ilength

 ! **hru with double precision type of fixed length
 type, public :: hru_d
  type(var_d),pointer                    :: hru(:) => null()
 endtype hru_d

 ! ** hru with integer type of variable length
 type, public :: hru_i
  type(var_i),pointer                 :: hru(:) => null()
 endtype hru_i

 ! define top-level derived types
 ! NOTE: either allocate directly, or use to point to higher dimensional structures
 type(gru), pointer, save, public        :: gru_struc(:) => null()  ! gru-hru map     
 type(mapping),pointer,save,public       :: index_map(:) => null()   ! relating indexing
 type(hru_i),pointer,save,public         :: time_gru(:) => null()    ! model time data
 type(hru_d),pointer,save,public         :: forc_gru(:) => null()    ! model forcing data
 type(hru_d),pointer,save,public         :: attr_gru(:) => null()    ! local attributes for each HRU
 type(hru_i),pointer,save,public         :: type_gru(:) => null()    ! local classification of soil veg etc. for each HRU
 type(hru_d),pointer,save,public         :: mpar_gru(:) => null()    ! model parameters
 type(hru_dlength),pointer,save,public   :: mvar_gru(:) => null()    ! model variables
 type(hru_ilength),pointer,save,public   :: indx_gru(:) => null()    ! model indices
 !type(hru_d), pointer,save, public       :: dt_init(:)  => null()    ! used to initialize the length of the sub-step for each HRU
 !type(hru_d), pointer,save, public       :: upArea(:)   => null()    ! area upslope of each HRU
 ! define data types for individual HRUs, and for basin-average quantities
 type(var_i),pointer,save,public         :: time_data => null()      ! model time data
 type(var_d),pointer,save,public         :: forc_data => null()      ! model forcing data
 type(var_d),pointer,save,public         :: attr_data => null()      ! local attributes
 type(var_i),pointer,save,public         :: type_data => null()      ! local classification of veg, soil, etc.
 type(var_d),pointer,save,public         :: mpar_data => null()      ! local column model parameters
 type(var_dlength),pointer,save,public   :: mvar_data => null()      ! local column model variables
 type(var_ilength),pointer,save,public   :: indx_data => null()      ! local column model indices
 type(var_d),pointer,save,public         :: bpar_data => null()      ! basin-average model parameters
 type(var_dlength),pointer,save,public   :: bvar_data => null()      ! basin-average model variables
 ! ***********************************************************************************************************
 ! Define common variables
 ! ***********************************************************************************************************
 integer(i4b),save,public                :: nGRU                     ! number of total grus
 integer(i4b),save,public                :: nHRU                     ! number of total hrus in entire domain
 integer(i4b),save,public                :: nSnow                    ! number of snow layers
 integer(i4b),save,public                :: nSoil                    ! number of soil layers
 integer(i4b),save,public                :: nLayers                  ! total number of layers in the snow-soil system
 integer(i4b),save,public                :: numtim                   ! number of time steps
 integer(i4b),save,public                :: data_step                ! time step of the data in seconds
 integer(i4b),save,public                :: data_steps               ! total time steps in a forcing data file
 real(dp),save,public                    :: refJulday                ! reference time in fractional julian days
 real(dp),save,public                    :: fracJulday               ! fractional julian days since the start of year
 real(dp),save,public                    :: dJulianStart             ! julian day of start time of simulation
 real(dp),save,public                    :: dJulianFinsh             ! julian day of end time of simulation
 integer(i4b),save,public                :: yearLength               ! number of days in the current year
 integer(i4b),save,public                :: urbanVegCategory=1       ! vegetation category for urban areas
 logical(lgt),save,public                :: doJacobian=.false.       ! flag to compute the Jacobian
 logical(lgt),save,public                :: globalPrintFlag=.false.  ! flag to compute the Jacobian
 ! ***********************************************************************************************************
 ! Define ancillary data structures
 ! ***********************************************************************************************************
 type(var_i),pointer,save,public         :: refTime      => null()   ! reference time for the model simulation
 type(var_i),pointer,save,public         :: startTime    => null()   ! start time for the model simulation
 type(var_i),pointer,save,public         :: finshTime    => null()   ! end time for the model simulation
 ! ***********************************************************************************************************


END MODULE data_struc

