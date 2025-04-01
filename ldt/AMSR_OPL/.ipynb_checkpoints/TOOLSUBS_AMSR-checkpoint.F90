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
!  MODULE, TOOLSUBS_AMSR E.J for extracting data from AMSR L1R
!=======================================================================

#include "LDT_misc.h"

MODULE TOOLSUBS_AMSR
    USE FUNCTIONS
#if (defined USE_HDF5)
    USE HDF5
#endif
    IMPLICIT NONE

    CONTAINS

      ! Forked version of GetSMAP_L1B_NRT_subset, and Modified it for AMSR2 
      ! L1R fields.
      ! Ehsan Jalilvand.
      SUBROUTINE get_amsr_l1r(filename, tb_time_seconds, &
          tb_10v, tb_10h, tb_18v, tb_18h, &
          tb_23v, tb_23h, tb_36v, tb_36h, &
          tb_89v, tb_89h, &
          lat, lon, lat89, lon89, &
          land_water_frac, snow, precip, &
          pixel_qual_flag, &
          n, m, m89, n89, ierr)

      ! Imports
      use LDT_logMod, only: LDT_logunit, LDT_endrun ! EMK

      ! Arguments
      character(*), intent(in) :: filename
      real*8, allocatable, intent(out) :: tb_time_seconds(:,:)
      real*4, allocatable, intent(out) :: tb_10v(:,:), tb_10h(:,:), tb_18v(:,:), tb_18h(:,:), tb_23v(:,:), tb_23h(:,:), tb_36v(:,:), &
      tb_36h(:,:), tb_89v(:,:), tb_89h(:,:)
      integer*4, allocatable, intent(out) :: land_water_frac(:,:)
      real*4, allocatable, intent(out) :: lat89(:,:), lon89(:,:)
      real*4, allocatable, intent(out) :: lat(:,:), lon(:,:)
    ! integer*4, allocatable, intent(out) :: scan_qual_flag(:,:)
      integer*4, allocatable, intent(out) :: pixel_qual_flag(:,:)
    !  integer*4, allocatable, intent(out) :: land_ocean_flag(:,:)
      integer, intent(out) :: m, n, m89, n89 ! m89 & n89 are for the 89GHz band for which lat and lon are provided
      integer :: i, j 
      integer, intent(out) :: ierr
      
      real :: sil, tt18  ! For intermediate snow and precip qual flag calculations
      integer*4, allocatable :: snow(:,:), precip(:,:)
      
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
      call get_dataset_real4_2d(file_id, dataset, n89, m89,lat89, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

      dataset = "Longitude of Observation Point for 89A"
      call get_dataset_real4_2d(file_id, dataset, n89, m89,lon89, ierr)
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if     
      
      ! Resample lat/lon using zoom

      call zoom_2d(lat89, (/ n89, m89 /), lat, (/ n, m /))
      call zoom_2d(lon89, (/ n89, m89 /), lon, (/ n, m /))
     
      ! Reading quality flag data
      !dataset = "Scan Data Quality"
      !call get_dataset_integer2_2d(file_id, dataset, n, m, &
      !     scan_qual_flag, ierr)
      !if (ierr == 1) then
      !   call h5fclose_f(file_id, hdferr)
      !   call h5close_f(hdferr)
      !   call freeall(ierr)
      !   return
      !end if
      
! TODO E.J: read the pixel quality flag for the RFI in c-band both H and V polarization, this won't be used for filtering of the footprints before making a 2D grid rather we just want to make a 2D grid of RFI flag for user information [According to Rajat], keep in mind that this data has the dims of (#scans x 486) so we need the (lat89, lon89) and (m89,n89) for making a 2D grid from it!

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
      call get_dataset_integer1_3d(file_id, dataset, n, m, 4, land_water_frac,ierr) ! here we sliced the 4th layer or 36 GHz channel that has a more detailed land water fraction map
      if (ierr == 1) then
         call h5fclose_f(file_id, hdferr)
         call h5close_f(hdferr)
         call freeall(ierr)
         return
      end if

        ! Compute snow and precipitation flags
        ! tb_18h →  (18.7GHz,H)
        ! tb_18v →  (18.7GHz,V)
        ! tb_23v →  (23.8GHz,V)
        ! tb_36v →  (36.5GHz,V)
        ! tb_89v →  (89.0GHz,V)
        ! land_water_frac → land_water_frac 
        ! Precip and snow flag at footprint level
        do j = 1, int(m) !dims(2)
           do i = 1, int(n) ! dims(1)
              ! Check if it's land (using the Land_Ocean flag from your data(:,:,11))
              if (land_water_frac(i,j) >= 50) then
                 ! Calculate intermediate variables
                 sil = 451.88 - 0.44*tb_18v(i,j) - 1.775*tb_23v(i,j) + &
                       0.00574*tb_23v(i,j)**2 - tb_89v(i,j)
                 tt18 = tb_18v(i,j) - tb_18h(i,j)  ! tbv18 - tbh18
                 
                 if (sil > 10) then
                    if ((tb_23v(i,j) <= 264.0) .and. &
                        (tb_23v(i,j) <= (175.0 + 0.49*tb_89v(i,j)))) then
                       ! Snow branch
                       snow(i,j) = 1
                       if ((tt18 >= 18) .and. &
                           ((tb_18v(i,j) - tb_36v(i,j)) <= 10) .and. &
                           ((tb_36v(i,j) - tb_89v(i,j)) <= 10)) then
                          snow(i,j) = 0
                       endif
                       if ((tt18 >= 8) .and. &
                           ((tb_18v(i,j) - tb_36v(i,j)) <= 2) .and. &
                           ((tb_23v(i,j) - tb_89v(i,j)) <= 6)) then
                          snow(i,j) = 1
                       endif
                    else
                       ! Precipitation branch
                       snow(i,j) = 0
                       precip(i,j) = 1
                       if (tt18 > 20) then
                          precip(i,j) = 0
                       endif
                       if ((tb_89v(i,j) > 253) .and. (tt18 > 7)) then
                          precip(i,j) = 0
                       endif
                    endif
                 endif
              endif
           end do
        end do
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
      subroutine get_dataset_integer1_3d(file_id, dataset, n, m, slice_index, var2d,ierr)

        ! Defaults
        implicit none

        ! Arguments
        integer(HID_T), intent(in) :: file_id
        character(*), intent(in) :: dataset
        integer, intent(out) :: n
        integer, intent(out) :: m
        integer*4, allocatable, intent(out) :: var2d(:,:)
        integer, intent(in) :: slice_index  ! New parameter
        integer, intent(out) :: ierr

        ! Locals
        integer(HID_T) :: dataset_id, datatype_id
        integer(HID_T) :: dspace_id, mem_space_id
        logical :: link_exists
        integer :: hdferr
        integer(HSIZE_T) :: dims(3), maxdims(3)  ! Changed to 3D
        integer(HSIZE_T) :: start(3), count(3)    ! For hyperslab
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
        
        ! Check size (should be 1 byte for uint8)
        call h5tget_size_f(datatype_id, size, hdferr)
        if (size /= 1) then
            write(LDT_logunit,*)'[ERR] Wrong byte size for ', trim(dataset)
            write(LDT_logunit,*)'[ERR] Expected 1, found ', size
            call h5tclose_f(datatype_id, hdferr)
            call h5dclose_f(dataset_id, hdferr)
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

        ! Check that the rank is 3.
        if (rank .ne. 3) then
           write(LDT_logunit,*) &
                '[ERR] Wrong rank for ', trim(dataset), &
                ', expected 3, found ', rank
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

        ! Check if slice_index is valid
        if (slice_index < 1 .or. slice_index > dims(3)) then
            write(LDT_logunit,*)'[ERR] Invalid slice index for ', trim(dataset)
            call h5sclose_f(dspace_id, hdferr)
            call h5dclose_f(dataset_id, hdferr)
            ierr = 1
            return
        end if

        ! Set up hyperslab selection
        start = [0, 0, slice_index-1]  ! Convert to 0-based indexing
        count = [dims(1), dims(2), 1_HSIZE_T]
        
        call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, start, count, hdferr)
        
        ! Create memory space for 2D slice
        call h5screate_simple_f(2, dims(1:2), mem_space_id, hdferr)


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
        ! Read the dataset slice
        call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, var2d, dims(1:2), hdferr, &
                       mem_space_id, dspace_id)        
        
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
      end subroutine get_dataset_integer1_3d
      
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
          if (allocated(lat89)) deallocate(lat89)
          if (allocated(lon89)) deallocate(lon89)
          if (allocated(tb_time_seconds)) deallocate(tb_time_seconds)
          !if (allocated(scan_qual_flag)) deallocate(scan_qual_flag)
          if (allocated(pixel_qual_flag)) deallocate(pixel_qual_flag)
          !if (allocated(land_ocean_flag)) deallocate(land_ocean_flag)
          m = 0
          n = 0
          ierr = 1
          return
        end subroutine freeall

        subroutine zoom_2d(input, dims_in, output, dims_out) 
           implicit none
           integer, intent(in) :: dims_in(2), dims_out(2)
           real*4, intent(in) :: input(dims_in(1), dims_in(2))
           real*4, intent(out) :: output(dims_out(1), dims_out(2))
           
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
      end subroutine get_amsr_l1r 

END MODULE TOOLSUBS_AMSR
