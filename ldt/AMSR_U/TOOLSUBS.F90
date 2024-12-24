!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!=============================i==========================================
!  MODULE, TOOLSUBS P.W.LIU, 09/12/18
!  Contains subroutines that are necessry for downscaling program
!-----------------------------------------------------------------------
!  AVHRR_NDVI
!  NOTES : Modified P.W. Liu, 09/12/18
!=======================================================================

#include "LDT_misc.h"

MODULE TOOLSUBS
    USE FUNCTIONS
#if (defined USE_HDF5)
    USE HDF5
#endif
    IMPLICIT NONE

    CONTAINS

      ! Forked version of GetSMAP_L1B_NRT_subset, and Modified it for AMSR2 
      ! L1R fields.
      ! Ehsan Jalilvand.
      SUBROUTINE GetAMSR_L1R_NRT_subset(filename, tb_time_seconds, &
          tb_10v, tb_10h, tb_18v, tb_18h, &
          tb_23v, tb_23h, tb_36v, tb_36h, &
          tb_89v, tb_89h, lat, lon, &
          scan_qual_flag, pixel_qual_flag, Land_Ocean_flag, &
          n, m, ierr)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun ! EMK

      ! Arguments
      character(*), intent(in) :: filename
      real*8, allocatable, intent(out) :: tb_time_seconds(:,:)
      real*4, allocatable, intent(out) :: tb_10v(:,:), tb_10h(:,:), tb_18v(:,:), tb_18h(:,:), tb_23v(:,:), tb_23h(:,:), tb_36v(:,:),
      tb_36h(:,:), tb_89v(:,:), tb_89h(:,:)
      real*4, allocatable, intent(out) :: lat_89A(:,:), lon_89A(:,:)
      integer*4, allocatable, intent(out) :: scan_qual_flag(:,:)
      integer*4, allocatable, intent(out) :: pixel_qual_flag(:,:), Land_Ocean_flag(:,:)
      !real*4, allocatable, intent(out) :: sc_nadir_angle(:)
      !real*4, allocatable, intent(out) :: antenna_scan_angle(:,:)
      integer, intent(out) :: m, n, m89, n89 ! m89 & n89 are for the 89GHz band for which lat and lon are provided
      integer, intent(out) :: ierr

#if (defined USE_HDF5)

      ! Locals
      character(100) :: dataset
      integer(HID_T) :: file_id, dataset_id, dspace_id ! HID_T type is used for HDF5 obj
      integer(HSIZE_T) :: dims(2), maxdims(2) ! HSIZE_T for dims of HDF5 file
      integer :: rank
      integer :: hdferr
      logical :: exists, ishdf5, link_exists

      ierr = 0
      m = 0
      n = 0

      ! Make sure file exists
      inquire(file=trim(filename), exist=exists)
      if (.not. exists) then
         write(LDT_logunit,*)'[ERR] Cannot find file ', trim(filename)
         ierr = 1
         return
      end if

      ! Initialize HDF5
      call h5open_f(hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Cannot initialize HDF5 Fortran interface!'
         call h5close_f(hdferr)
         ierr = 1
         return
      end if

      ! Make sure the file is HDF5
      call h5fis_hdf5_f(trim(filename), ishdf5, hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Problem checking if ', trim(filename), &
              ' is HDF5'
         call h5close_f(hdferr)
         ierr = 1
         return
      end if
      if (.not. ishdf5) then
         write(LDT_logunit,*)'[ERR] File ', trim(filename), ' is not HDF5!'
         call h5close_f(hdferr)
         ierr = 1
         return
      end if

      ! Open the file
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, hdferr)
      if (hdferr == -1) then
         write(LDT_logunit,*)'[ERR] Cannot open ', trim(filename)
         call h5close_f(hdferr)
         ierr = 1
         return
      end if

      ! Get the data

      dataset = "Scan Time"
      call get_dataset_real8_2d(file_id, dataset, n, m, tb_time_seconds, &
           ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if    
      
      dataset = "Brightness Temperature (res10,10.7GHz,V)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_10v, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if
      
      dataset = "Brightness Temperature (res10,18.7GHz,H)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_18h,ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if
      
      dataset = "Brightness Temperature (res10,18.7GHz,V)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_18v, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Brightness Temperature (res10,23.8GHz,H)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_23h, ierr)     
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if
      dataset = "Brightness Temperature (res10,23.8GHz,V)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_23v, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Brightness Temperature (res10,36.5GHz,H)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_36h, ierr)    
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if
      dataset = "Brightness Temperature (res10,36.5GHz,V)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_36v, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Brightness Temperature (res10,89.0GHz,H)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_89h, ierr)   
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if
   
      dataset = "Brightness Temperature (res10,89.0GHz,V)"
      call get_dataset_real4_2d(file_id, dataset, n, m, tb_89v, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Latitude of Observation Point for 89A"
      call get_dataset_real4_2d(file_id, dataset, n89, m89,lat_89A, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Longitude of Observation Point for 89A"
      call get_dataset_real4_2d(file_id, dataset, n89, m89,lon_89A, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if     
      
      ! Resample lat/lon using zoom

      call zoom_2d(lat_89A, (/ n89, m89 /), lat, (/ n, m /))
      call zoom_2d(lon_89A, (/ n89, m89 /), lon, (/ n, m /))
      
      dataset = "Scan Data Quality"
      call get_dataset_integer2_2d(file_id, dataset, n, m, &
           scan_qual_flag, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Pixel Data Quality 6 to 36"
      call get_dataset_integer2_2d(file_id, dataset, n, m, &
           pixel_qual_flag, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Land_Ocean Flag 6 to 36"
      call get_dataset_integer2_2d(file_id, dataset, n, m, &
           Land_Ocean_flag, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if


      ! Clean up
      call h5fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

      ierr = 0

      return

    contains

      ! Internal subroutine
      subroutine get_dataset_integer2_2d(file_id, dataset, n, m, var2d, &
           ierr)

        ! Defaults
        implicit none

        ! Arguments
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: dataset
        integer, intent(out) :: n
        integer, intent(out) :: m
        integer*4, allocatable, intent(out) :: var2d(:,:)
        integer, intent(out) :: ierr

        ! Locals
        integer(HID_T) :: dataset_id, datatype_id
        logical :: link_exists
        integer :: hdferr
        integer(HID_T) :: dspace_id
        integer(HSIZE_T) :: dims(2), maxdims(2)
        integer :: rank
        integer :: class
        integer(SIZE_T) :: size
        integer :: sign

        ierr = 0

        ! See if the dataset is in the file
        call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Problem finding ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif
        if (.not. link_exists) then
           write(LDT_logunit,*)'[ERR] Nonexistent dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Get the dataset id
        call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot open dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the datatype id
        call h5dget_type_f(dataset_id, datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the datatype class
        call h5tget_class_f(datatype_id, class, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get class for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (class .ne. H5T_INTEGER_F) then
           write(LDT_logunit,*)'[ERR] Bad class for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected ', H5T_INTEGER_F, &
                ', found ', class
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the size of the datatype.
        call h5tget_size_f(datatype_id, size, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get size for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (size .ne. 2) then
           write(LDT_logunit,*)'[ERR] Wrong byte size found for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected 2, found ', size
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the sign type of the datatype.  Should be unsigned.
        call h5tget_sign_f(datatype_id, sign, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get sign type for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (sign .ne. H5T_SGN_NONE_F) then
           write(LDT_logunit,*)'[ERR] Wrong sign type found for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected ', H5T_SGN_NONE_F, &
                ', found ', sign
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close the datatype
        call h5tclose_f(datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot close datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dspace id for the variable dimensions
        call h5dget_space_f(dataset_id, dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot find dimensions for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the rank of the dataset in the file.
        call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get rank for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check that the rank is 2.
        if (rank .ne. 2) then
           write(LDT_logunit,*) &
                '[ERR] Wrong rank for ', trim(dataset), &
                ', expected 2, found ', rank
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dimensions
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get dimensions for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close access to the dataspace.
        call h5sclose_f(dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot close access to dataspace for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Allocate and initialize the array
        n = dims(1)
        m = dims(2)
        allocate(var2d(n,m))
        var2d = 0

        ! Read the dataset.  Fortran doesn't have unsigned integers,
        ! so we save the 16-bit unsigned integer in a 32-bit signed
        ! integer (should have the room).
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, var2d, dims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot read dataset ', trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Close access to the dataset
        call h5dclose_f(dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_Logunit,*)'[ERR] Problem closing dataset ', &
                trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        return
      end subroutine get_dataset_integer2_2d

      ! Internal subroutine
      subroutine get_dataset_real4_2d(file_id, dataset, n, m, var2d, ierr)

        ! Defaults
        implicit none

        ! Arguments
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: dataset
        integer, intent(out) :: n
        integer, intent(out) :: m
        real*4, allocatable, intent(out) :: var2d(:,:)
        integer, intent(out) :: ierr

        ! Locals
        integer(HID_T) :: dataset_id, datatype_id
        logical :: link_exists
        integer :: hdferr
        integer(HID_T) :: dspace_id
        integer(HSIZE_T) :: dims(2), maxdims(2)
        integer :: rank
        integer :: class
        integer(SIZE_T) :: size

        ierr = 0

        ! See if the dataset is in the file
        call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Problem finding ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif
        if (.not. link_exists) then
           write(LDT_logunit,*)'[ERR] Nonexistent dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Get the dataset id
        call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot open dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the datatype id
        call h5dget_type_f(dataset_id, datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the datatype class
        call h5tget_class_f(datatype_id, class, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get class for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (class .ne. H5T_FLOAT_F) then
           write(LDT_logunit,*)'[ERR] Bad class for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected ', H5T_FLOAT_F, &
                ', found ', class
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the size of the datatype
        call h5tget_size_f(datatype_id, size, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get size for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (size .ne. 4) then
           write(LDT_logunit,*)'[ERR] Wrong byte size found for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected 4, found ', size
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close the datatype
        call h5tclose_f(datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot close datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dspace id for the variable dimensions
        call h5dget_space_f(dataset_id, dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot find dimensions for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the rank of the dataset in the file.
        call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get rank for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check that the rank is 2.
        if (rank .ne. 2) then
           write(LDT_logunit,*) &
                '[ERR] Wrong rank for ', trim(dataset), &
                ', expected 2, found ', rank
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dimensions
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get dimensions for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close access to the dataspace.
        call h5sclose_f(dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot close access to dataspace for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Allocate and initialize the array
        n = dims(1)
        m = dims(2)
        allocate(var2d(n,m))
        var2d = 0

        ! Read the dataset
        call h5dread_f(dataset_id, H5T_NATIVE_REAL, var2d, dims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot read dataset ', trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Close access to the dataset
        call h5dclose_f(dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_Logunit,*)'[ERR] Problem closing dataset ', &
                trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        return
      end subroutine get_dataset_real4_2d

      ! Internal subroutine
      subroutine get_dataset_real8_2d(file_id, dataset, n, m, var2d, ierr)

        ! Defaults
        implicit none

        ! Arguments
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: dataset
        integer, intent(out) :: n
        integer, intent(out) :: m
        real*8, allocatable, intent(out) :: var2d(:,:)
        integer, intent(out) :: ierr

        ! Locals
        integer(HID_T) :: dataset_id, datatype_id
        logical :: link_exists
        integer :: hdferr
        integer(HID_T) :: dspace_id
        integer(HSIZE_T) :: dims(2), maxdims(2)
        integer :: rank
        integer :: class
        integer(SIZE_T) :: size

        ierr = 0

        ! See if the dataset is in the file
        call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Problem finding ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif
        if (.not. link_exists) then
           write(LDT_logunit,*)'[ERR] Nonexistent dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Get the dataset id
        call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot open dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the datatype id
        call h5dget_type_f(dataset_id, datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the datatype class
        call h5tget_class_f(datatype_id, class, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get class for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (class .ne. H5T_FLOAT_F) then
           write(LDT_logunit,*)'[ERR] Bad class for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected ', H5T_FLOAT_F, &
                ', found ', class
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the size of the datatype
        call h5tget_size_f(datatype_id, size, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get size for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (size .ne. 8) then
           write(LDT_logunit,*)'[ERR] Wrong byte size found for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected 8, found ', size
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close the datatype
        call h5tclose_f(datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot close datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dspace id for the variable dimensions
        call h5dget_space_f(dataset_id, dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot find dimensions for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the rank of the dataset in the file.
        call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get rank for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check that the rank is 2.
        if (rank .ne. 2) then
           write(LDT_logunit,*) &
                '[ERR] Wrong rank for ', trim(dataset), &
                ', expected 2, found ', rank
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dimensions
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get dimensions for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close access to the dataspace.
        call h5sclose_f(dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*) &
                '[ERR] Cannot close access to dataspace for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Allocate and initialize the array
        n = dims(1)
        m = dims(2)
        allocate(var2d(n,m))
        var2d = 0

        ! Read the dataset
        call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, var2d, dims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot read dataset ', trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Close access to the dataset
        call h5dclose_f(dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_Logunit,*)'[ERR] Problem closing dataset ', &
                trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        return
      end subroutine get_dataset_real8_2d

      ! Internal subroutine
      subroutine get_dataset_real4_1d(file_id, dataset, n, var1d, ierr)

        ! Defaults
        implicit none

        ! Arguments
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: dataset
        integer, intent(out) :: n
        real*4, allocatable, intent(out) :: var1d(:)
        integer, intent(out) :: ierr

        ! Locals
        integer(HID_T) :: dataset_id, datatype_id
        logical :: link_exists
        integer :: hdferr
        integer(HID_T) :: dspace_id
        integer(HSIZE_T) :: dims(1), maxdims(1)
        integer :: rank
        integer :: class
        integer(SIZE_T) :: size

        ierr = 0

        ! See if the dataset is in the file
        call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Problem finding ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif
        if (.not. link_exists) then
           write(LDT_logunit,*)'[ERR] Nonexistent dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Get the dataset id
        call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot open dataset ', trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the datatype id
        call h5dget_type_f(dataset_id, datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the datatype class
        call h5tget_class_f(datatype_id, class, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get class for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (class .ne. H5T_FLOAT_F) then
           write(LDT_logunit,*)'[ERR] Bad class for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected ', H5T_FLOAT_F, &
                ', found ', class
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check the size of the datatype
        call h5tget_size_f(datatype_id, size, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get size for ', &
                trim(dataset)
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if
        if (size .ne. 4) then
           write(LDT_logunit,*)'[ERR] Wrong byte size found for ', &
                trim(dataset)
           write(LDT_logunit,*)'[ERR] Expected 4, found ', size
           call h5tclose_f(datatype_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close the datatype
        call h5tclose_f(datatype_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot close datatype for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dspace id for the variable dimensions
        call h5dget_space_f(dataset_id, dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot find dimensions for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the rank of the dataset in the file.
        call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get rank for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Check that the rank is 1.
        if (rank .ne. 1) then
           write(LDT_logunit,*) &
                '[ERR] Wrong rank for ', trim(dataset), &
                ', expected 1, found ', rank
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Get the dimensions
        call h5sget_simple_extent_dims_f(dspace_id, dims, maxdims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot get dimensions for ', &
                trim(dataset)
           call h5sclose_f(dspace_id, hdferr)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Close access to the dataspace.
        call h5sclose_f(dspace_id, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot close access to dataspace for ', &
                trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        ! Allocate and initialize the array
        n = dims(1)
        allocate(var1d(n))
        var1d = 0

        ! Read the dataset
        call h5dread_f(dataset_id, H5T_NATIVE_REAL, var1d, dims, hdferr)
        if (hdferr == -1) then
           write(LDT_logunit,*)'[ERR] Cannot read dataset ', trim(dataset)
           call h5dclose_f(dataset_id, hdferr)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        endif

        ! Close access to the dataset
        call h5dclose_f(dataset_id, hdferr)
        if (hdferr == -1) then
           write(LDT_Logunit,*)'[ERR] Problem closing dataset ', &
                trim(dataset)
           call h5fclose_f(file_id, hdferr)
           call h5close_f(hdferr)
           ierr = 1
           return
        end if

        return
      end subroutine get_dataset_real4_1d

        ! Internal subroutine.  Warning -- deallocates memory in
        ! parent subroutine and resets two variables.  This is intended
        ! for gracefully handling errors returned from HDF5.
        subroutine freeall(ierr)
          implicit none
          integer, intent(out) :: ierr
          if (allocated(tb_10h)) deallocate(tb_10h)
          if (allocated(tb_10v)) deallocate(tb_10v)
          if (allocated(tb_18h)) deallocate(tb_18h)
          if (allocated(tb_18v)) deallocate(tb_18v)
          if (allocated(tb_23h)) deallocate(tb_23h)
          if (allocated(tb_23v)) deallocate(tb_23v)
          if (allocated(tb_36h)) deallocate(tb_36h)
          if (allocated(tb_36v)) deallocate(tb_36v)
          if (allocated(tb_89h)) deallocate(tb_89h)
          if (allocated(tb_89v)) deallocate(tb_89v)
          if (allocated(lat_89A)) deallocate(lat_89A)
          if (allocated(lon_89A)) deallocate(lon_89A)
          if (allocated(tb_time_seconds)) deallocate(tb_time_seconds)
          if (allocated(scan_qual_flag)) deallocate(scan_qual_flag)
          if (allocated(pixel_qual_flag)) deallocate(pixel_qual_flag)
          if (allocated(Land_Ocean_flag)) deallocate(Land_Ocean_flag)
          m = 0
          n = 0
          ierr = 1
          return
        end subroutine freeall

        subroutine zoom_2d(input, dims_in, output, dims_out) 
           implicit none
           integer(hsize_t), intent(in) :: dims_in(2), dims_out(2)
           real, intent(in) :: input(dims_in(1), dims_in(2))
           real, intent(out) :: output(dims_out(1), dims_out(2))
           
           real :: x_scale, y_scale, x, y
           integer :: i, j, x1, x2, y1, y2
           real :: dx, dy
           real :: c11, c12, c21, c22
           real :: f1, f2
           
           ! Compute scaling factors
           x_scale = real(dims_in(1) - 1) / real(dims_out(1) - 1)
           y_scale = real(dims_in(2) - 1) / real(dims_out(2) - 1)
           
           do j = 1, dims_out(2)
               do i = 1, dims_out(1)
                   ! Get input coordinates
                   x = 1.0 + (i-1) * x_scale
                   y = 1.0 + (j-1) * y_scale
                   
                   ! Get surrounding points
                   x1 = int(x)
                   x2 = min(x1 + 1, int(dims_in(1)))
                   y1 = int(y)
                   y2 = min(y1 + 1, int(dims_in(2)))
                   
                   ! Get interpolation weights
                   dx = x - x1
                   dy = y - y1
                   
                   ! Get corner values
                   c11 = input(x1, y1)
                   c12 = input(x1, y2)
                   c21 = input(x2, y1)
                   c22 = input(x2, y2)
                   
                   ! Bilinear interpolation
                   f1 = (1.0-dx)*c11 + dx*c21
                   f2 = (1.0-dx)*c12 + dx*c22
                   output(i,j) = (1.0-dy)*f1 + dy*f2
               end do
           end do
           
        end subroutine zoom_2d
        
        #else
        ! Dummy version if LDT was compiled w/o HDF5 support.
        write(LDT_logunit,*) &
             '[ERR] GetSMAP_L1B_NRT called without HDF5 support!'
        write(LDT_logunit,*) &
             '[ERR] Recompile LDT with HDF5 support and try again!'
        call LDE_endrun()
#endif
      end subroutine GetAMSR_L1R_NRT_subset 

END MODULE TOOLSUBS
