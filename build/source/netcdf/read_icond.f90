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

module read_icond_module
USE nrtype
USE netcdf
implicit none
private
public::read_icond
public::read_icond_nlayers
#if ( defined LIS_SUMMA_2_0 )
public::read_icond_lis
public::read_icond_nlayers_lis
#endif

! define single HRU restart file
integer(i4b), parameter :: singleHRU=1001
integer(i4b), parameter :: multiHRU=1002
integer(i4b), parameter :: restartFileType=multiHRU
contains

 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_nlayers(iconFile,nGRU,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookIndex                        ! variable lookup structure
 USE globalData,only:gru_struc                         ! gru-hru mapping structures
 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling
 USE data_types,only:gru_hru_intVec                    ! actual data
 USE data_types,only:var_info                          ! metadata
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)        ,intent(in)     :: iconFile       ! name of input (restart) file
 integer(i4b)        ,intent(in)     :: nGRU           ! total # of GRUs in run domain
 type(var_info)      ,intent(in)     :: indx_meta(:)   ! metadata
 integer(i4b)        ,intent(out)    :: err            ! error code
 character(*)        ,intent(out)    :: message        ! returned error message

 ! locals
 integer(i4b)             :: ncID                       ! netcdf file id
 integer(i4b)             :: dimID                      ! netcdf file dimension id
 integer(i4b)             :: fileHRU                    ! number of HRUs in netcdf file
 integer(i4b)             :: snowID, soilID             ! netcdf variable ids
 integer(i4b)             :: iGRU, iHRU                 ! loop indexes
 integer(i4b),allocatable :: snowData(:)                ! number of snow layers in all HRUs
 integer(i4b),allocatable :: soilData(:)                ! number of soil layers in all HRUs
 character(len=256)       :: cmessage                   ! downstream error message

 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_nlayers/'

 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncid,err,cmessage);
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimId);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimId,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! allocate sotrage for reading from file
 allocate(snowData(fileHRU))
 allocate(soilData(fileHRU))
 snowData = 0
 soilData = 0

 ! get variable ids
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookIndex%nSnow)%varName),snowid); call netcdf_err(err,message)
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookIndex%nSoil)%varName),soilid); call netcdf_err(err,message)

 ! get data
 err = nf90_get_var(ncid,snowid,snowData); call netcdf_err(err,message)
 err = nf90_get_var(ncid,soilid,soilData); call netcdf_err(err,message)
 !print*, 'snowData = ', snowData
 !print*, 'soilData = ', soilData

 ! assign to index structure - gru by hru
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   
   ! single HRU
   if(restartFileType==singleHRU)then   
    gru_struc(iGRU)%hruInfo(iHRU)%nSnow = snowData(1)
    gru_struc(iGRU)%hruInfo(iHRU)%nSoil = soilData(1)

   ! multi HRU
   else
    gru_struc(iGRU)%hruInfo(iHRU)%nSnow = snowData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)
    gru_struc(iGRU)%hruInfo(iHRU)%nSoil = soilData(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)
   endif 

  end do
 end do

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(snowData,soilData)

 end subroutine read_icond_nlayers


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(iconFile,                      & ! intent(in):    name of initial conditions file
                       nGRU,                          & ! intent(in):    number of GRUs
                       mparData,                      & ! intent(in):    model parameters
                       progData,                      & ! intent(inout): model prognostic variables
                       indxData,                      & ! intent(inout): model indices
                       err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookPROG                          ! variable lookup structure
 USE var_lookup,only:iLookPARAM                         ! variable lookup structure
 USE var_lookup,only:iLookINDEX                         ! variable lookup structure
 USE globalData,only:prog_meta                          ! metadata for prognostic variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globaldata,only:iname_soil,iname_snow              ! named variables to describe the type of layer
 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_doubleVec                  ! full double precision structure
 USE data_types,only:gru_hru_intVec                     ! full integer structure
 USE data_types,only:var_dlength                        ! double precision structure for a single HRU
 USE data_types,only:var_info                           ! metadata
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 USE updatState_module,only:updateSoil                  ! update soil states
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)     :: iconFile     ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)     :: nGRU         ! number of grouped response units in simulation domain
 type(gru_hru_doubleVec),intent(in)     :: mparData     ! model parameters
 type(gru_hru_doubleVec),intent(inout)  :: progData     ! model prognostic variables
 type(gru_hru_intVec)   ,intent(inout)  :: indxData     ! model indices
 integer(i4b)           ,intent(out)    :: err          ! error code
 character(*)           ,intent(out)    :: message      ! returned error message

 ! locals
 character(len=256)                     :: cmessage     ! downstream error message
 integer(i4b)                           :: fileHRU      ! number of HRUs in file
 integer(i4b)                           :: iVar         ! loop index
 integer(i4b)                           :: iGRU         ! loop index
 integer(i4b)                           :: iHRU         ! loop index
 integer(i4b)                           :: dimID        ! varible dimension ids
 integer(i4b)                           :: ncVarID      ! variable ID in netcdf file
 character(256)                         :: dimName      ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                           :: dimLen       ! data dimensions
 integer(i4b)                           :: ncID         ! netcdf file ID
 integer(i4b)                           :: ixFile       ! index in file
 real(dp),allocatable                   :: varData(:,:) ! variable data storage
 integer(i4b)                           :: nSoil, nSnow, nToto ! # layers
 integer(i4b)                           :: iLayer,jLayer ! layer indices
 integer(i4b),parameter                 :: nBand=2      ! number of spectral bands

 character(len=32),parameter            :: scalDimName   ='scalarv'  ! dimension name for scalar data
 character(len=32),parameter            :: midSoilDimName='midSoil'  ! dimension name for soil-only layers
 character(len=32),parameter            :: midTotoDimName='midToto'  ! dimension name for layered varaiables
 character(len=32),parameter            :: ifcTotoDimName='ifcToto'  ! dimension name for layered varaiables

 ! --------------------------------------------------------------------------------------------------------

 ! Start procedure here
 err=0; message="read_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimID,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! loop through prognostic variables
 do iVar = 1,size(prog_meta)

  ! skip variables that are computed later
  if(prog_meta(iVar)%varName=='scalarCanopyWat'           .or. &
     prog_meta(iVar)%varName=='spectralSnowAlbedoDiffuse' .or. &
     prog_meta(iVar)%varName=='scalarSurfaceTemp'         .or. &
     prog_meta(iVar)%varName=='mLayerVolFracWat'          .or. &
     prog_meta(iVar)%varName=='mLayerHeight'                   ) cycle

  ! get variable id
  err = nf90_inq_varid(ncID,trim(prog_meta(iVar)%varName),ncVarID); call netcdf_err(err,message)
  if(err/=0)then
   message=trim(message)//': problem with getting variable id, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get variable dimension IDs
  select case (prog_meta(iVar)%varType)
   case (iLookVarType%scalarv); err = nf90_inq_dimid(ncID,trim(scalDimName)   ,dimID); call netcdf_err(err,message)
   case (iLookVarType%midSoil); err = nf90_inq_dimid(ncID,trim(midSoilDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%midToto); err = nf90_inq_dimid(ncID,trim(midTotoDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcToto); err = nf90_inq_dimid(ncID,trim(ifcTotoDimName),dimID); call netcdf_err(err,message)
   case default
    message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
    err=20; return
  end select

  ! check errors
  if(err/=0)then
   message=trim(message)//': problem with dimension ids, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get the dimension length
  err = nf90_inquire_dimension(ncID,dimID,dimName,dimLen); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the dimension length'; return; endif

  ! iniitialize the variable data
  allocate(varData(fileHRU,dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating variable data'; return; endif

  ! get data
  err = nf90_get_var(ncID,ncVarID,varData); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

  ! store data in prognostics structure
  ! loop through GRUs
  do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount

    ! get the number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
    nToto = nSnow + nSoil

    ! get the index in the file: single HRU
    if(restartFileType==singleHRU)then   
     ixFile = 1  ! use for single HRU restart file

    ! get the index in the file: multi HRU
    else
     ixFile = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
    endif

    ! put the data into data structures and check that none of the values are set to nf90_fill_double
    select case (prog_meta(iVar)%varType)
     case (iLookVarType%scalarv)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1)       = varData(ixFile,1)
      if(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData))then; err=20; endif
     case (iLookVarType%midSoil)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nSoil) = varData(ixFile,1:nSoil)
      if(any(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nSoil) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
     case (iLookVarType%midToto)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nToto) = varData(ixFile,1:nToto)
      if(any(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nToto) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
     case (iLookVarType%ifcToto)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(0:nToto) = varData(ixFile,1:nToto+1)
      if(any(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(0:nToto) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
     case default
      message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
      err=20; return
    end select

    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(prog_meta(iVar)%varName)//"')"; return; endif

    ! fix the snow albedo
    if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < 0._dp)then
     progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMax)%dat(1)
    endif

    ! initialize the spectral albedo
    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nBand) = progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1)

   end do ! iHRU
  end do ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData, stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating variable data'; return; endif

 end do ! iVar

 ! --------------------------------------------------------------------------------------------------------
 ! (2) set number of layers
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount

   ! save the number of layers
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)   = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)   = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1) = gru_struc(iGRU)%hruInfo(iHRU)%nSnow + gru_struc(iGRU)%hruInfo(iHRU)%nSoil

   ! set layer type
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat(1:gru_struc(iGRU)%hruInfo(iHRU)%nSnow) = iname_snow
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat((gru_struc(iGRU)%hruInfo(iHRU)%nSnow+1):(gru_struc(iGRU)%hruInfo(iHRU)%nSnow+gru_struc(iGRU)%hruInfo(iHRU)%nSoil)) = iname_soil

  end do
 end do

 ! --------------------------------------------------------------------------------------------------------
 ! (3) update soil layers
 ! --------------------------------------------------------------------------------------------------------

 ! loop through GRUs and HRUs
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount

   ! loop through soil layers
   do iLayer = 1,indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)

    ! get layer in the total vector
    jLayer = iLayer+indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)

    ! update soil layers
    call updateSoil(&
                    ! input
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp          )%dat(jLayer),& ! intent(in): temperature vector (K)
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerMatricHead    )%dat(iLayer),& ! intent(in): matric head (m)
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_alpha          )%dat(iLayer),& ! intent(in): van Genutchen "alpha" parameter
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_n              )%dat(iLayer),& ! intent(in): van Genutchen "n" parameter
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_sat          )%dat(iLayer),& ! intent(in): soil porosity (-)
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_res          )%dat(iLayer),& ! intent(in): soil residual volumetric water content (-)
                    1._dp - 1._dp/mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_n)%dat(iLayer),& ! intent(in): van Genutchen "m" parameter (-)
                    ! output
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat    )%dat(jLayer),& ! intent(out): volumetric fraction of total water (-)
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq    )%dat(jLayer),& ! intent(out): volumetric fraction of liquid water (-)
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce    )%dat(jLayer),& ! intent(out): volumetric fraction of ice (-)
                    err,message)                                                                   ! intent(out): error control
    if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

   end do  ! looping through soil layers

  end do  ! looping throuygh HRUs
 end do  ! looping through GRUs

 end subroutine read_icond



#if ( defined LIS_SUMMA_2_0 )
 ! Added by Zhuo Wang on 12/17/2019 for LIS/SUMMA restart run mode
 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers_lis: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
  subroutine read_icond_nlayers_lis(iconFile,nGRU,typeData,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE globalData,only:gru_struc 

 USE var_lookup,only:iLookINDEX                        ! variable lookup structure
 USE var_lookup,only:iLookTYPE

 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling

 USE data_types,only:gru_hru_intVec                    ! actual data
 USE data_types,only:var_info
 USE data_types,only:gru_hru_int

 use LIS_coreMod,    only : LIS_rc, LIS_masterproc
 use LIS_historyMod, only : LIS_readvar_restart
 use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
                            LIS_getNextUnitNumber,   &
                            LIS_releaseUnitNumber,   &
                            LIS_verify

!USE summa_type, only:summa1_type_dec

 implicit none
   
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
!integer               ,intent(in)    :: n
!type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure

 character(*)        ,intent(in)     :: iconFile       ! name of input (restart) file
 integer(i4b)        ,intent(in)     :: nGRU           ! total # of GRUs in run domain
 type(gru_hru_int)   ,intent(in)     :: typeData
 integer(i4b)        ,intent(out)    :: err            ! error code
 character(*)        ,intent(out)    :: message        ! returned error message
  
 ! locals
!type(summa1_type_dec)    :: summa1_struc 
 integer(i4b)             :: ncID                       ! netcdf file id
 integer(i4b)             :: dimID                      ! netcdf file dimension id
 integer(i4b)             :: fileHRU                    ! number of HRUs in netcdf file
 integer(i4b)             :: snowID, soilID             ! netcdf variable ids
 integer(i4b)             :: iGRU, iHRU                 ! loop indexes

 real, allocatable        :: tmptilen(:)
 real, allocatable        :: tmptilen_1d(:)
 integer(i4b), allocatable  :: tmpint_1d(:)

 integer(i4b),allocatable :: nSnowData(:)                ! number of snow layers in all HRUs
 integer(i4b),allocatable :: nSoilData(:)                ! number of soil layers in all HRUs
 integer(i4b),allocatable :: hruIdData(:)

 character(len=256)       :: cmessage                   ! downstream error message

 integer             :: n
 integer             :: t, l
 integer             :: nc, nr, npatch
 integer             :: ftn
 integer             :: status

!real(dp) :: nSnow 
!real(dp) :: nSoil
 integer(i4b) :: hruId
 integer(i4b)        :: nGRU_lis
 integer(i4b)        :: nHRUrun
 
 integer :: nlayer_real             ! rean snow number + nsoil = real snow number + 8
 
 logical             :: file_exists
 character*20        :: wformat

 ! --------------------------------------------------------------------------------------------------------

 ! initialize error message
 err=0
 message = 'read_icond_nlayers_lis/'

 ! Revised on 02/24/2020
 do n=1, LIS_rc%nnest
 ! open netcdf file
!wformat = trim(summa1_struc(n)%rformat)
 wformat = "netcdf"

! if(trim(LIS_rc%startcode) == "restart") then
   ! open restart file
!  if(wformat .eq. "netcdf") then
!#if (defined USE_NETCDF3 || defined USE_NETCDF4)
     status = nf90_open(path=iconFile, mode=NF90_NOWRITE, ncid=ftn)
     call LIS_verify(status, "Error opening file "//iconFile)
!#endif
!  endif
 !---------------------------------------------------------------------------------------------

  nGRU_lis = LIS_rc%ntiles(n)
  nHRUrun = sum(gru_struc%hruCount)

  allocate(tmptilen(nHRUrun))
  allocate(tmptilen_1d(nHRUrun))
  allocate(tmpint_1d(nHRUrun))

  allocate(nSnowData(nHRUrun))
  allocate(nSoilData(nHRUrun))
  allocate(hruIdData(nHRUrun))

  nSnowData = 0
  nSoilData = 0
  hruIdData = 0

   tmpint_1d = 0
   call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmpint_1d, &
                            varname="nSoil", wformat=wformat)

   do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount

     gru_struc(iGRU)%hruInfo(iHRU)%nSoil = tmpint_1d(iGRU) 
    enddo
   enddo
   nSoilData = tmpint_1d 

   tmpint_1d = 0
   call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmpint_1d, &
                 varname="nSnow", wformat=wformat)

   do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount
     gru_struc(iGRU)%hruInfo(iHRU)%nSnow = tmpint_1d(iGRU) 
    enddo
   enddo
   nSnowData = tmpint_1d

   tmpint_1d = 0
   call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmpint_1d, &
                            varname="hruId", wformat=wformat)

   do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount
     hruIdData(iGRU) = tmpint_1d(iGRU)
    enddo
   enddo
   hruIdData = tmpint_1d

     status = nf90_close(ftn)
     call LIS_verify(status,"Error in nf90_close in summa2_readrst")

  deallocate(tmptilen,tmptilen_1d,tmpint_1d,nSnowData,nSoilData)

!  endif ! restart

 enddo  ! n

 end subroutine read_icond_nlayers_lis
#endif


#if ( defined LIS_SUMMA_2_0 )
! Added by Zhuo Wang on 12/18/2019 for LIS/SUMMA restart run mode
! ************************************************************************************************
! public subroutine read_icond_lis: read model initial conditions from LIS restart file
! ************************************************************************************************
 subroutine read_icond_lis(iconFile,                    & ! intent(in):    name of initial conditions file
                           nGRU,                          & ! intent(in):    number of GRUs
                           typeData,                      &
                           mparData,                      & ! intent(in):    model parameters
                           progData,                      & ! intent(inout): model prognostic variables
                           indxData,                      & ! intent(inout): model indices
                           err,message)                     ! intent(out):   error control
! --------------------------------------------------------------------------------------------------------
! modules
 USE nrtype
 USE var_lookup,only:iLookType
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookPROG                          ! variable lookup structure
 USE var_lookup,only:iLookPARAM                         ! variable lookup structure
 USE var_lookup,only:iLookINDEX                         ! variable lookup structure

 USE globalData,only:prog_meta                          ! metadata for prognostic variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globaldata,only:iname_soil,iname_snow              ! named variables to describe the type of layer

 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_doubleVec                  ! full double precision structure
 USE data_types,only:gru_hru_intVec                     ! full integer structure
 USE data_types,only:var_dlength                        ! double precision structure for a single HRU
 USE data_types,only:var_info                           ! metadata
 USE data_types,only:gru_hru_int

 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 USE updatState_module,only:updateSoil                  ! update soil states

 use LIS_coreMod,    only : LIS_rc, LIS_masterproc
 use LIS_historyMod, only : LIS_readvar_restart
 use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
                            LIS_getNextUnitNumber,   &
                            LIS_releaseUnitNumber,   &
                            LIS_verify
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)     :: iconFile     ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)     :: nGRU         ! number of grouped response units in simulation domain
 type(gru_hru_int)      ,intent(inout)     :: typeData
 type(gru_hru_doubleVec),intent(in)     :: mparData     ! model parameters
 type(gru_hru_doubleVec),intent(inout)  :: progData     ! model prognostic variables
 type(gru_hru_intVec)   ,intent(inout)  :: indxData     ! model indices
 integer(i4b)           ,intent(out)    :: err          ! error code
 character(*)           ,intent(out)    :: message      ! returned error message

 ! locals
 character(len=256)                     :: cmessage     ! downstream error message
 integer(i4b)                           :: fileHRU      ! number of HRUs in file
 integer(i4b)                           :: iVar         ! loop index
 integer(i4b)                           :: iGRU         ! loop index
 integer(i4b)                           :: iHRU         ! loop index
 integer(i4b)                           :: dimID        ! varible dimension ids
 integer(i4b)                           :: ncVarID      ! variable ID in netcdf file
 character(256)                         :: dimName      ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                           :: dimLen       ! data dimensions
 integer(i4b)                           :: ncID         ! netcdf file ID
 integer(i4b)                           :: ixFile       ! index in file
 real(dp),allocatable                   :: varData(:,:) ! variable data storage
 integer(i4b)                           :: nSoil, nSnow, nToto ! # layers
 integer(i4b)                           :: iLayer,jLayer ! layer indices
 integer(i4b),parameter                 :: nBand=2      ! number of spectral bands

 character(len=32),parameter            :: scalDimName   ='time'  ! dimension name for scalar data
 character(len=32),parameter            :: midSoilDimName='dim2'  ! dimension name for soil-only layers
 character(len=32),parameter            :: midTotoDimName='dim3'  ! dimension name for layered varaiables
 character(len=32),parameter            :: ifcTotoDimName='dim5'  ! dimension name for layered varaiables

 integer             :: n
 integer             :: t, l
 integer             :: nc, nr, npatch
 integer             :: ftn
 integer             :: status

 real, allocatable        :: tmptilen(:)
 real, allocatable        :: tmptilen_1d(:)
 integer(i4b),allocatable :: tmpint_1d(:)

 integer(i4b),allocatable :: nSnowData(:)
 integer(i4b),allocatable :: nSoilData(:)
 integer(i4b),allocatable :: hruIdData(:)

 integer                  :: midToto_new, ifcToto_new, soilLay_new

 integer(i4b) :: nlayer_real   ! rean snow number + nsoil = real snow number + 8
 integer(i4b) :: nHRUrun, num

 logical             :: file_exists
 character*20        :: wformat

! --------------------------------------------------------------------------------------------------------
! Start procedure here
 err=0; message="read_icond_lis/"

 midToto_new = 8     ! midToto_new = 13
 ifcToto_new = 9     ! ifcToto_new = 14
 soilLay_new = 3
! --------------------------------------------------------------------------------------------------------
! (1) read the file
! --------------------------------------------------------------------------------------------------------

 do n=1, LIS_rc%nnest

! wformat = trim(summa1_struc(n)%rformat)
  wformat = "netcdf"

!  if(trim(LIS_rc%startcode) == "restart") then

  nHRUrun = sum(gru_struc%hruCount)
  allocate(tmptilen(nHRUrun))
  allocate(tmptilen_1d(nHRUrun))
  allocate(tmpint_1d(nHRUrun))
  allocate(nSnowData(nHRUrun))
  allocate(nSoilData(nHRUrun))
  allocate(hruIdData(nHRUrun))

!#if (defined USE_NETCDF3 || defined USE_NETCDF4)
     status = nf90_open(path=iconFile, mode=NF90_NOWRITE, ncid=ftn) 
!#endif
!-------------------------------------------------------------------------------------------
!!! These data are read in LIS hurId order
  ! read: length of initial time step at start of next data interval 
  tmpint_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmpint_1d, &
                           varname="hruId", wformat=wformat)
  typeData%gru(:)%hru(1)%var(iLookTYPE%hruId) = tmpint_1d
  hruIdData = tmpint_1d

  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="dt_init", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%dt_init)%dat(1) = tmptilen_1d

  ! read: mass of ice on the vegetation canopy
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarCanopyIce", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarCanopyIce)%dat(1) = tmptilen_1d

  ! read: mass of liquid water on the vegetation canopy
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarCanopyLiq", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarCanopyLiq)%dat(1) = tmptilen_1d

  ! read: mass of total water on the vegetation canopy
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarCanopyWat", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarCanopyWat)%dat(1) = tmptilen_1d

  ! read: temperature of the canopy air space
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarCanairTemp", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarCanairTemp)%dat(1) = tmptilen_1d

  ! read: temperature of the vegetation canopy
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarCanopyTemp", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarCanopyTemp)%dat(1) = tmptilen_1d

  ! read: snow albedo for the entire spectral band
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarSnowAlbedo", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = tmptilen_1d

  ! read: total snow depth
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarSnowDepth", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarSnowDepth)%dat(1) = tmptilen_1d

  ! read: snow water equivalent
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarSWE", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarSWE)%dat(1) = tmptilen_1d

  ! read: ponded water caused by melt of the snow without a layer
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarSfcMeltPond", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarSfcMeltPond)%dat(1) = tmptilen_1d

  ! read: relative aquifer storage -- above bottom of the soil profile
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarAquiferStorage", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarAquiferStorage)%dat(1) = tmptilen_1d

  ! read: surface temperature
  tmptilen_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen_1d, &
                           varname="scalarSurfaceTemp", wformat=wformat)
  progData%gru(:)%hru(1)%var(iLookPROG%scalarSurfaceTemp)%dat(1) = tmptilen_1d

  ! read: number of snow layers
  tmpint_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmpint_1d, &
                           varname="nSnow", wformat=wformat)
  indxData%gru(:)%hru(1)%var(iLookINDEX%nSnow)%dat(1) = tmpint_1d

  ! read: number of soil layers
  tmpint_1d=0
  call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmpint_1d, &
                           varname="nSoil", wformat=wformat)

  indxData%gru(:)%hru(1)%var(iLookINDEX%nSoil)%dat(1) = tmpint_1d

  nSnowData = indxData%gru(:)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
  nSoilData = indxData%gru(:)%hru(1)%var(iLookINDEX%nSoil)%dat(1)

  ! Set number of layers
  do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
     nlayer_real = nSnowData(t) + nSoilData(t)

     indxData%gru(t)%hru(1)%var(iLookINDEX%nLayers)%dat(1) = nSnowData(t) + nSoilData(t)

     ! Set layer type
     indxData%gru(t)%hru(1)%var(iLookINDEX%layerType)%dat(1:nSnowData(t)) = iname_snow
     indxData%gru(t)%hru(1)%var(iLookINDEX%layerType)%dat((nSnowData(t) + 1):nlayer_real) = iname_soil
  enddo
!---------------------------------------------------------------------------------------------
  ! read: diffuse snow albedo for individual spectral bands 
  do l=1, 2 ! spectral 
   tmptilen=0
   call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="spectralSnowAlbedoDiffuse", &
                            dim=l, vlevels = 2, wformat=wformat)
   do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
    progData%gru(t)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(l) = tmptilen(t)

    ! fix the snow albedo
    if(progData%gru(t)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < 0._dp) then
      progData%gru(t)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(t)%hru(1)%var(iLookPARAM%albedoMax)%dat(1)
    endif
            
    ! initialize the spectral albedo
    progData%gru(t)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(l) = &
                        progData%gru(t)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1)
    enddo
   enddo
 
   ! read: temperature of each layer
!  do l=1, 13 ! midToto
   do l=1, 8  ! midToto
     tmptilen=0
     call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="mLayerTemp", &
                              dim=l, vlevels = midToto_new, wformat=wformat)
     do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
      nlayer_real = nSnowData(t) + nSoilData(t)
      if(l .le. nlayer_real) then
        progData%gru(t)%hru(1)%var(iLookPROG%mLayerTemp)%dat(l) = tmptilen(t)
      endif
     enddo
   enddo
 
   ! read: volumetric fraction of ice in each layer
!  do l=1, 13 ! midToto
   do l=1, 8  ! midToto
     tmptilen=0
     call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="mLayerVolFracIce", &
                              dim=l, vlevels = midToto_new, wformat=wformat)
     do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
      nlayer_real = nSnowData(t) + nSoilData(t)
      if(l .le. nlayer_real) then
        progData%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat(l) = tmptilen(t)
      endif
     enddo
    enddo

    ! read: volumetric fraction of liquid water in each layer
!   do l=1, 13 ! midToto
    do l=1, 8  ! midToto
       tmptilen=0
       call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="mLayerVolFracLiq", &
                                dim=l, vlevels = midToto_new, wformat=wformat)
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        nlayer_real = nSnowData(t) + nSoilData(t)
        if(l .le. nlayer_real) then
          progData%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat(l) = tmptilen(t)
        endif
       enddo
    enddo

    ! read: volumetric fraction of total water in each layer
!   do l=1, 13 ! midToto
    do l=1, 8  ! midToto
      tmptilen=0
      call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="mLayerVolFracWat", &
                               dim=l, vlevels = midToto_new, wformat=wformat)
      do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
       nlayer_real = nSnowData(t) + nSoilData(t)
       if(l .le. nlayer_real) then
         progData%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat(l) = tmptilen(t)
       endif
      enddo
    enddo

    ! read: matric head of water in the soil
    !!! Note: For 3-Layer soil, vlevels = soilLay_new=3 !!!!
!   do l=1, 8 ! midSoil
    do l=1, 3 ! midSoil
       tmptilen=0
       call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="mLayerMatricHead", &
                                dim=l, vlevels = soilLay_new, wformat=wformat)
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
          progData%gru(t)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat(l) = tmptilen(t)
       enddo
    enddo

    ! read: depth of each layer
!   do l=1, 13 ! midToto
    do l=1, 8  ! midToto
       tmptilen=0
       call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="mLayerDepth", &
                                dim=l, vlevels = midToto_new, wformat=wformat)
       do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        nlayer_real = nSnowData(t) + nSoilData(t)
        if(l .le. nlayer_real) then
          progData%gru(t)%hru(1)%var(iLookPROG%mLayerDepth)%dat(l) = tmptilen(t)
        endif
       enddo
    enddo

!   do l=1, 13 ! midToto
    do l=1, 8  ! midToto
      tmptilen=0
      call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="mLayerHeight", &
                               dim=l, vlevels = midToto_new, wformat=wformat)
      do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        nlayer_real = nSnowData(t) + nSoilData(t)
        if(l .le. nlayer_real) then
          progData%gru(t)%hru(1)%var(iLookPROG%mLayerHeight)%dat(l) = tmptilen(t)
        endif
      enddo
    enddo

    ! read: height of the layer interface; top of soil = 0
!   do l=0, 13 ! ifcToto -- 0:midToto
    do l=0, 8  ! ifcToto -- 0:midToto
      tmptilen=0
      call LIS_readvar_restart(ftn, n, LIS_rc%lsm_index, tmptilen, varname="iLayerHeight", &
                               dim=l+1, vlevels = ifcToto_new, wformat=wformat)
      do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index)
        nlayer_real = nSnowData(t) + nSoilData(t)
        if(l .le. nlayer_real) then
          progData%gru(t)%hru(1)%var(iLookPROG%iLayerHeight)%dat(l) = tmptilen(t)
        endif
      enddo
    enddo

! --------------------------------------------------------------------------------------------------------
! (3) update soil layers
! --------------------------------------------------------------------------------------------------------
 
  do t=1, LIS_rc%npatch(n, LIS_rc%lsm_index) 
    ! loop through soil layers
    do iLayer = 1,indxData%gru(t)%hru(1)%var(iLookINDEX%nSoil)%dat(1)

    ! get layer in the total vector
    jLayer = iLayer+indxData%gru(t)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
 
    ! update soil layers
    !!! progData are read in from LIS restart file and in LIS hruId order, but mparData is in stand-alone hruId order
    call updateSoil(&
                    ! input
                    progData%gru(t)%hru(1)%var(iLookPROG%mLayerTemp          )%dat(jLayer),& ! intent(in): temperature vector (K)
                    progData%gru(t)%hru(1)%var(iLookPROG%mLayerMatricHead    )%dat(iLayer),& ! intent(in): matric head (m)
                    mparData%gru(t)%hru(1)%var(iLookPARAM%vGn_alpha          )%dat(iLayer),& ! intent(in): van Genutchen "alpha" parameter
                    mparData%gru(t)%hru(1)%var(iLookPARAM%vGn_n              )%dat(iLayer),& ! intent(in): van Genutchen "n" parameter
                    mparData%gru(t)%hru(1)%var(iLookPARAM%theta_sat          )%dat(iLayer),& ! intent(in): soil porosity (-)
                    mparData%gru(t)%hru(1)%var(iLookPARAM%theta_res          )%dat(iLayer),& ! intent(in): soil residual volumetric water content (-)
                    1._dp - 1._dp/mparData%gru(t)%hru(1)%var(iLookPARAM%vGn_n)%dat(iLayer),& ! intent(in): van Genutchen "m" parameter (-)
                    ! output
                    progData%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracWat    )%dat(jLayer),& ! intent(out): volumetric fraction of total water (-)
                    progData%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracLiq    )%dat(jLayer),& ! intent(out): volumetric fraction of liquid water (-)
                    progData%gru(t)%hru(1)%var(iLookPROG%mLayerVolFracIce    )%dat(jLayer),& ! intent(out): volumetric fraction of ice (-)
                    err,message)                                                                   ! intent(out): error control
   if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
   end do  ! looping through soil layers
 end do

 ! close restart file
 if(wformat .eq. "binary") then
    call LIS_releaseUnitNumber(ftn)
 elseif(wformat .eq. "netcdf") then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    status = nf90_close(ftn)
    call LIS_verify(status, "Error in nf90_close in summa2_readrst")
#endif
  endif

   deallocate(tmptilen,tmptilen_1d,tmpint_1d,nSnowData,nSoilData)
! endif ! restart

 enddo ! n

  end subroutine read_icond_lis
 
#endif

end module read_icond_module
