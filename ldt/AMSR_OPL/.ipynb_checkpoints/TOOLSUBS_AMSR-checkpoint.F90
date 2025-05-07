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
    USE LDT_logMod, only: LDT_logunit, LDT_endrun
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
          use LDT_logMod, only: LDT_logunit, LDT_endrun
        
          ! Arguments
          character(*), intent(in) :: filename
          real*8, allocatable, intent(out) :: tb_time_seconds(:)
          real*4, allocatable, intent(out) :: tb_10v(:,:), tb_10h(:,:), tb_18v(:,:), tb_18h(:,:), tb_23v(:,:), tb_23h(:,:), tb_36v(:,:), &
          tb_36h(:,:), tb_89v(:,:), tb_89h(:,:)
          integer*4, allocatable, intent(out) :: land_water_frac(:,:)
          real*4, allocatable, intent(out) :: lat89(:,:), lon89(:,:)
          real*4, allocatable, intent(out) :: lat(:,:), lon(:,:)
          integer*2, allocatable, intent(out) :: pixel_qual_flag(:,:)
          integer, intent(out) :: m, n, m89, n89 ! m89 & n89 are for the 89GHz band for which lat and lon are provided
          integer :: i, j 
          integer, intent(out) :: ierr
          
          real :: sil, tt18  ! For intermediate snow and precip qual flag calculations
          integer*4, allocatable :: snow(:,:), precip(:,:)
          
#if (defined USE_HDF5)
          ! Locals
          character(100) :: dataset
          integer(HID_T) :: file_id
          integer :: hdferr
          logical :: exists, ishdf5
        
          ierr = 0
          m = 0
          n = 0
          m89 = 0
          n89 = 0
        
          ! Make sure file exists
          inquire(file=trim(filename), exist=exists)
          if (.not. exists) then
             write(LDT_logunit,*)'[ERR] Cannot find file ', trim(filename)
             ierr = 1
             return
          end if
        
          ! Initialize HDF5
          call h5open_f(hdferr)
          if (hdferr /= 0) then
             write(LDT_logunit,*)'[ERR] Cannot initialize HDF5 Fortran interface!'
             call h5close_f(hdferr)
             ierr = 1
             return
          end if
        
          ! Make sure the file is HDF5
          call h5fis_hdf5_f(trim(filename), ishdf5, hdferr)
          if (hdferr /= 0) then
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
          if (hdferr /= 0) then
             write(LDT_logunit,*)'[ERR] Cannot open ', trim(filename)
             call h5close_f(hdferr)
             ierr = 1
             return
          end if
        
          ! ===============================================
          ! STEP 1: Get dimensions from a reference dataset
          ! ===============================================
          
          ! Get dimensions for 89GHz lat/lon data
          call get_dataset_latlon_89(file_id, "Latitude of Observation Point for 89A", &
               n89, m89, lat89, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[ERR] Failed to read 89GHz lat data'
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             ierr = 1
             return
          end if
          
          call get_dataset_latlon_89(file_id, "Longitude of Observation Point for 89A", &
               n89, m89, lon89, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[ERR] Failed to read 89GHz lon data'
             deallocate(lat89)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             ierr = 1
             return
          end if
          
          ! Get dimensions for brightness temperature data (10.7GHz as reference)
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,10.7GHz,H)", &
               n, m, tb_10h, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[ERR] Failed to read 10.7GHz H TB data'
             deallocate(lat89)
             deallocate(lon89)
             call h5fclose_f(file_id, hdferr)
             call h5close_f(hdferr)
             ierr = 1
             return
          end if
          
          ! ===============================================
          ! STEP 2: Allocate remaining arrays with known dimensions
          ! ===============================================
          
          if (.not. allocated(tb_10v)) allocate(tb_10v(n, m))
          if (.not. allocated(tb_18h)) allocate(tb_18h(n, m))
          if (.not. allocated(tb_18v)) allocate(tb_18v(n, m))
          if (.not. allocated(tb_23h)) allocate(tb_23h(n, m))
          if (.not. allocated(tb_23v)) allocate(tb_23v(n, m))
          if (.not. allocated(tb_36h)) allocate(tb_36h(n, m))
          if (.not. allocated(tb_36v)) allocate(tb_36v(n, m))
          if (.not. allocated(tb_89h)) allocate(tb_89h(n, m))
          if (.not. allocated(tb_89v)) allocate(tb_89v(n, m))
          
          if (.not. allocated(lat)) allocate(lat(n, m))
          if (.not. allocated(lon)) allocate(lon(n, m))
          
          if (.not. allocated(snow)) allocate(snow(n, m))
          if (.not. allocated(precip)) allocate(precip(n, m))
          if (.not. allocated(land_water_frac)) allocate(land_water_frac(n, m))
          if (.not. allocated(pixel_qual_flag)) allocate(pixel_qual_flag(n, m))
          
          ! Initialize all arrays
          tb_10v = 0.0
          tb_18h = 0.0
          tb_18v = 0.0
          tb_23h = 0.0
          tb_23v = 0.0
          tb_36h = 0.0
          tb_36v = 0.0
          tb_89h = 0.0
          tb_89v = 0.0
          lat = 0.0
          lon = 0.0
          snow = 0
          precip = 0
          land_water_frac = 0
          pixel_qual_flag = 0
          
          ! ===============================================
          ! STEP 3: Get time data
          ! ===============================================
          
          call get_dataset_scan_time(file_id, "Scan Time", m, tb_time_seconds, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read scan time'
             ! Continue anyway
          end if
          
          ! ===============================================
          ! STEP 4: Read remaining TB fields
          ! ===============================================
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,10.7GHz,V)", &
               n, m, tb_10v, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 10.7GHz V TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,18.7GHz,H)", &
               n, m, tb_18h, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 18.7GHz H TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,18.7GHz,V)", &
               n, m, tb_18v, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 18.7GHz V TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,23.8GHz,H)", &
               n, m, tb_23h, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 23.8GHz H TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,23.8GHz,V)", &
               n, m, tb_23v, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 23.8GHz V TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,36.5GHz,H)", &
               n, m, tb_36h, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 36.5GHz H TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,36.5GHz,V)", &
               n, m, tb_36v, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 36.5GHz V TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,89.0GHz,H)", &
               n, m, tb_89h, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 89.0GHz H TB data'
          end if
          
          call get_dataset_tb_2d(file_id, "Brightness Temperature (res10,89.0GHz,V)", &
               n, m, tb_89v, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read 89.0GHz V TB data'
          end if
          
          ! ===============================================
          ! STEP 5: Read land/water fraction and quality flags
          ! ===============================================
          
          ! Read land/water fraction - use channel 4 (36GHz)
          call get_dataset_land_ocean_flag(file_id, "Land_Ocean Flag 6 to 36", &
               n, m, 4, land_water_frac, ierr)
          if (ierr /= 0) then
             write(LDT_logunit,*)'[WARN] Failed to read land/ocean flag'
             ! Continue anyway
          end if
          
          ! Read pixel quality flags
          !call get_dataset_pixel_quality(file_id, "Pixel Data Quality 6 to 36", &
               !n, m, pixel_qual_flag, ierr)
          !if (ierr /= 0) then
             !write(LDT_logunit,*)'[WARN] Failed to read pixel quality flags'
             ! Continue anyway
          !end if
          
          ! ===============================================
          ! STEP 6: Generate interpolated lat/lon grid for res10 data
          ! ===============================================
          
          ! Interpolate the 89GHz lat/lon grid to the TB grid dimensions
          call zoom_2d(lat89, n89, m89, lat, n, m)
          call zoom_2d(lon89, n89, m89, lon, n, m)
          
          ! ===============================================
          ! STEP 7: Generate snow and precipitation flags
          ! ===============================================
          
          ! Compute snow and precipitation flags based on TB values
          do j = 1, m
             do i = 1, n
                ! Only process data over land
                if (land_water_frac(i,j) >= 50) then
                   ! Ensure all needed values are valid
                   if (tb_18v(i,j) > 0 .and. tb_18h(i,j) > 0 .and. &
                       tb_23v(i,j) > 0 .and. tb_89v(i,j) > 0 .and. &
                       tb_36v(i,j) > 0) then
                      
                      sil = 451.88 - 0.44*tb_18v(i,j) - 1.775*tb_23v(i,j) + &
                            0.00574*tb_23v(i,j)**2 - tb_89v(i,j)
                      tt18 = tb_18v(i,j) - tb_18h(i,j)
                      
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
                endif
             end do
          end do
          
          ! Clean up
          call h5fclose_f(file_id, hdferr)
          call h5close_f(hdferr)
          
          ! Clear error since we reached the end successfully
          ierr = 0
          return
        
        CONTAINS
          ! Internal subroutine for bilinear interpolation of 2D arrays
          subroutine zoom_2d(input, nx_in, ny_in, output, nx_out, ny_out)
            integer, intent(in) :: nx_in, ny_in, nx_out, ny_out
            real*4, intent(in) :: input(nx_in, ny_in)
            real*4, intent(out) :: output(nx_out, ny_out)
            
            real :: x_scale, y_scale, x, y
            integer :: i, j, x1, x2, y1, y2
            real :: dx, dy, c11, c12, c21, c22, f1, f2
            
            ! Compute scaling factors
            x_scale = real(nx_in - 1) / real(nx_out - 1)
            y_scale = real(ny_in - 1) / real(ny_out - 1)
            
            ! Perform bilinear interpolation
            do j = 1, ny_out
               do i = 1, nx_out
                  ! Get input coordinates
                  x = 1.0 + (i-1) * x_scale
                  y = 1.0 + (j-1) * y_scale
                  
                  ! Get surrounding points
                  x1 = max(1, min(int(x), nx_in))
                  x2 = max(1, min(x1 + 1, nx_in))
                  y1 = max(1, min(int(y), ny_in))
                  y2 = max(1, min(y1 + 1, ny_in))
                  
                  ! Interpolation weights
                  dx = max(0.0, min(1.0, x - x1))
                  dy = max(0.0, min(1.0, y - y1))
                  
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

        ! Function for reading 2D brightness temperature fields (uint16 with scale factor 0.01)
        subroutine get_dataset_tb_2d(file_id, dataset, n, m, var2d, ierr)
          use LDT_logMod, only: LDT_logunit
          implicit none
        
          ! Arguments
          integer(HID_T), intent(in) :: file_id
          character(*), intent(in) :: dataset
          integer, intent(out) :: n
          integer, intent(out) :: m
          real*4, allocatable, intent(out) :: var2d(:,:)
          integer, intent(out) :: ierr
        
          ! Locals
          integer(HID_T) :: dataset_id, dataspace_id
          logical :: link_exists
          integer :: hdferr
          integer(HSIZE_T) :: dims(2), maxdims(2)
          integer :: rank
          real*4, allocatable :: temp_data(:,:)  ! Temporary buffer for data
          logical :: already_allocated
        
          ierr = 0
        
          ! Check if dataset exists
          call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
          if (hdferr /= 0 .or. .not. link_exists) then
            write(LDT_logunit,*)'[ERR] Dataset not found: ', trim(dataset)
            ierr = 1
            return
          endif
        
          ! Open the dataset
          call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
          if (hdferr /= 0) then
            write(LDT_logunit,*)'[ERR] Cannot open dataset: ', trim(dataset)
            ierr = 1
            return
          endif
        
          ! Get the dataspace and dimensions
          call h5dget_space_f(dataset_id, dataspace_id, hdferr)
          if (hdferr /= 0) then
            write(LDT_logunit,*)'[ERR] Cannot get dataspace for: ', trim(dataset)
            call h5dclose_f(dataset_id, hdferr)
            ierr = 1
            return
          endif
        
          ! Get the rank
          call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
          if (hdferr /= 0) then
            write(LDT_logunit,*)'[ERR] Cannot get rank for: ', trim(dataset)
            call h5sclose_f(dataspace_id, hdferr)
            call h5dclose_f(dataset_id, hdferr)
            ierr = 1
            return
          endif
        
          ! Check that rank is 2
          if (rank /= 2) then
            write(LDT_logunit,*)'[ERR] Expected 2D dataset, found rank: ', rank
            call h5sclose_f(dataspace_id, hdferr)
            call h5dclose_f(dataset_id, hdferr)
            ierr = 1
            return
          endif
        
          ! Get dimensions
          call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
          if (hdferr < 0) then
            write(LDT_logunit,*)'[ERR] Cannot get dimensions for: ', trim(dataset)
            call h5sclose_f(dataspace_id, hdferr)
            call h5dclose_f(dataset_id, hdferr)
            ierr = 1
            return
          endif
        
          n = int(dims(1))
          m = int(dims(2))
          
          write(LDT_logunit,*)'[DEBUG] Dataset dimensions for ', trim(dataset), ': ', n, 'x', m
          
          ! Check if already allocated with correct dimensions
          already_allocated = allocated(var2d)
          if (already_allocated) then
            if (size(var2d,1) /= n .or. size(var2d,2) /= m) then
              deallocate(var2d)
              already_allocated = .false.
            endif
          endif
          
          ! Allocate if needed
          if (.not. already_allocated) then
            allocate(temp_data(n, m), stat=hdferr)
            if (hdferr /= 0) then
              write(LDT_logunit,*)'[ERR] Memory allocation failed for temp data array'
              call h5sclose_f(dataspace_id, hdferr)
              call h5dclose_f(dataset_id, hdferr)
              ierr = 1
              return
            endif
            
            ! Use a temporary array for reading to avoid potential memory corruption
            temp_data = 0.0
          
            ! Read data into temporary array with native real type
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, temp_data, dims, hdferr)
            if (hdferr /= 0) then
              write(LDT_logunit,*)'[ERR] Cannot read data for: ', trim(dataset)
              deallocate(temp_data)
              call h5sclose_f(dataspace_id, hdferr)
              call h5dclose_f(dataset_id, hdferr)
              ierr = 1
              return
            endif
            
            ! Now allocate the output array and copy data
            allocate(var2d(n, m), stat=hdferr)
            if (hdferr /= 0) then
              write(LDT_logunit,*)'[ERR] Memory allocation failed for output array'
              deallocate(temp_data)
              call h5sclose_f(dataspace_id, hdferr)
              call h5dclose_f(dataset_id, hdferr)
              ierr = 1
              return
            endif
            
            ! Copy data from temporary array and apply scaling if needed
            var2d = temp_data*.01
            
            ! Clean up temporary array
            deallocate(temp_data)
          else
            ! If already allocated with correct dimensions, read directly
            call h5dread_f(dataset_id, H5T_NATIVE_REAL, var2d, dims, hdferr)
            if (hdferr /= 0) then
              write(LDT_logunit,*)'[ERR] Cannot read data for: ', trim(dataset)
              call h5sclose_f(dataspace_id, hdferr)
              call h5dclose_f(dataset_id, hdferr)
              ierr = 1
              return
            endif
          endif
          
          ! Clean up
          call h5sclose_f(dataspace_id, hdferr)
          call h5dclose_f(dataset_id, hdferr)
          
          write(LDT_logunit,*)'[INFO] Successfully read brightness temperature data: ', trim(dataset)
        end subroutine get_dataset_tb_2d
            
            ! Function for reading Land_Ocean Flag (uint8 3D array)
            subroutine get_dataset_land_ocean_flag(file_id, dataset, n, m, layer, var2d, ierr)
              use LDT_logMod, only: LDT_logunit
              implicit none
            
              ! Arguments
              integer(HID_T), intent(in) :: file_id
              character(*), intent(in) :: dataset
              integer, intent(out) :: n
              integer, intent(out) :: m
              integer, intent(in) :: layer        ! Which layer to extract (1-4)
              integer*4, allocatable, intent(out) :: var2d(:,:)
              integer, intent(out) :: ierr
            
              ! Locals
              integer(HID_T) :: dataset_id, dataspace_id, memspace_id
              logical :: link_exists
              integer :: hdferr
              integer(HSIZE_T) :: dims3d(3), maxdims(3), dims(2)
              integer(HSIZE_T) :: start(3), count(3)
              integer :: rank
              integer*1, allocatable :: temp_data(:,:,:)  ! For uint8 data
              integer :: i, j
            
              ierr = 0
            
              ! Check if dataset exists
              call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
              if (hdferr /= 0 .or. .not. link_exists) then
                write(LDT_logunit,*)'[ERR] Dataset not found: ', trim(dataset)
                ierr = 1
                return
              endif
            
              ! Open the dataset
              call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot open dataset: ', trim(dataset)
                ierr = 1
                return
              endif
            
              ! Get the dataspace and dimensions
              call h5dget_space_f(dataset_id, dataspace_id, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot get dataspace for: ', trim(dataset)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
              if (hdferr /= 0 .or. rank /= 3) then
                write(LDT_logunit,*)'[ERR] Expected 3D dataset, found rank: ', rank
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              call h5sget_simple_extent_dims_f(dataspace_id, dims3d, maxdims, hdferr)
              if (hdferr < 0) then
                write(LDT_logunit,*)'[ERR] Cannot get dimensions for: ', trim(dataset)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              write(LDT_logunit,*)'[INFO] Land_Ocean flag dimensions:', &
                  dims3d(1), 'x', dims3d(2), 'x', dims3d(3)
              
              if (layer < 1 .or. layer > dims3d(3)) then
                write(LDT_logunit,*)'[ERR] Invalid layer requested: ', layer, ', max is: ', dims3d(3)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Define hyperslab to extract the requested layer
              start(1) = 0
              start(2) = 0
              start(3) = layer - 1  ! 0-indexed in HDF5
              count(1) = dims3d(1)
              count(2) = dims3d(2)
              count(3) = 1
            
              ! Setup dimensions for the 2D memory space
              n = int(dims3d(1))
              m = int(dims3d(2))
              dims(1) = n
              dims(2) = m
              
              ! Allocate output array
              allocate(var2d(n, m), stat=hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Memory allocation failed for output array'
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Select hyperslab in the file
              call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, start, count, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot select hyperslab'
                deallocate(var2d)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Create memory space for 2D result
              call h5screate_simple_f(2, dims, memspace_id, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot create memory space'
                deallocate(var2d)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Read the selected layer directly into the output array
              call h5dread_f(dataset_id, H5T_NATIVE_INTEGER, var2d, dims, hdferr, &
                            memspace_id, dataspace_id)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot read data for: ', trim(dataset)
                deallocate(var2d)
                call h5sclose_f(memspace_id, hdferr)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Clean up
              call h5sclose_f(memspace_id, hdferr)
              call h5sclose_f(dataspace_id, hdferr)
              call h5dclose_f(dataset_id, hdferr)
              
              write(LDT_logunit,*)'[INFO] Successfully read Land_Ocean flag data, layer: ', layer
            end subroutine get_dataset_land_ocean_flag
            
            ! Function for reading lat/lon data for 89A channel (float32)
            subroutine get_dataset_latlon_89(file_id, dataset, n, m, var2d, ierr)
              use LDT_logMod, only: LDT_logunit
              implicit none
            
              ! Arguments
              integer(HID_T), intent(in) :: file_id
              character(*), intent(in) :: dataset
              integer, intent(out) :: n
              integer, intent(out) :: m
              real*4, allocatable, intent(out) :: var2d(:,:)
              integer, intent(out) :: ierr
            
              ! Locals
              integer(HID_T) :: dataset_id, dataspace_id
              logical :: link_exists
              integer :: hdferr
              integer(HSIZE_T) :: dims(2), maxdims(2)
              integer :: rank
            
              ierr = 0
            
              ! Check if dataset exists
              call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
              if (hdferr /= 0 .or. .not. link_exists) then
                write(LDT_logunit,*)'[ERR] Dataset not found: ', trim(dataset)
                ierr = 1
                return
              endif
            
              ! Open the dataset
              call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot open dataset: ', trim(dataset)
                ierr = 1
                return
              endif
            
              ! Get the dataspace and dimensions
              call h5dget_space_f(dataset_id, dataspace_id, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot get dataspace for: ', trim(dataset)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
              if (hdferr /= 0 .or. rank /= 2) then
                write(LDT_logunit,*)'[ERR] Expected 2D dataset, found rank: ', rank
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
              if (hdferr < 0) then
                write(LDT_logunit,*)'[ERR] Cannot get dimensions for: ', trim(dataset)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              n = int(dims(1))
              m = int(dims(2))
              
              ! Allocate output array
              allocate(var2d(n, m), stat=hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Memory allocation failed for output array'
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Read data directly (already float32, no scaling needed)
              call h5dread_f(dataset_id, H5T_NATIVE_REAL, var2d, dims, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot read data for: ', trim(dataset)
                deallocate(var2d)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Clean up
              call h5sclose_f(dataspace_id, hdferr)
              call h5dclose_f(dataset_id, hdferr)
              
              write(LDT_logunit,*)'[INFO] Successfully read: ', trim(dataset)
            end subroutine get_dataset_latlon_89
            
            ! Function for reading scan time (float64)
            subroutine get_dataset_scan_time(file_id, dataset, n, time_arr, ierr)
                use LDT_logMod, only: LDT_logunit
                implicit none
                
                ! Arguments
                integer(HID_T), intent(in) :: file_id
                character(*), intent(in) :: dataset
                integer, intent(out) :: n
                real*8, allocatable, intent(out) :: time_arr(:)
                integer, intent(out) :: ierr
                
                ! Locals
                integer(HID_T) :: dataset_id, dataspace_id
                logical :: link_exists
                integer :: hdferr
                integer(HSIZE_T) :: dims(1), maxdims(1)
                integer :: rank
                
                ierr = 0
                
                ! Check if dataset exists
                call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
                if (hdferr /= 0 .or. .not. link_exists) then
                    write(LDT_logunit,*)'[ERR] Dataset not found: ', trim(dataset)
                    ierr = 1
                    return
                endif
                
                ! Open the dataset
                call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
                if (hdferr /= 0) then
                    write(LDT_logunit,*)'[ERR] Cannot open dataset: ', trim(dataset)
                    ierr = 1
                    return
                endif
                
                ! Get dimensions
                call h5dget_space_f(dataset_id, dataspace_id, hdferr)
                if (hdferr /= 0) then
                    write(LDT_logunit,*)'[ERR] Cannot get dataspace: ', trim(dataset)
                    call h5dclose_f(dataset_id, hdferr)
                    ierr = 1
                    return
                endif
                
                call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
                if (hdferr /= 0 .or. rank /= 1) then
                    write(LDT_logunit,*)'[ERR] Expected 1D dataset for scan time, found rank: ', rank
                    call h5sclose_f(dataspace_id, hdferr)
                    call h5dclose_f(dataset_id, hdferr)
                    ierr = 1
                    return
                endif
                
                call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
                if (hdferr < 0) then
                    write(LDT_logunit,*)'[ERR] Cannot get dimensions for: ', trim(dataset)
                    call h5sclose_f(dataspace_id, hdferr)
                    call h5dclose_f(dataset_id, hdferr)
                    ierr = 1
                    return
                endif
                
                n = int(dims(1))
                write(LDT_logunit,*)'[INFO] Scan time array length: ', n
                
                ! Allocate output array with correct dimension
                if (allocated(time_arr)) deallocate(time_arr)
                allocate(time_arr(n), stat=hdferr)
                if (hdferr /= 0) then
                    write(LDT_logunit,*)'[ERR] Memory allocation failed for scan time array'
                    call h5sclose_f(dataspace_id, hdferr)
                    call h5dclose_f(dataset_id, hdferr)
                    ierr = 1
                    return
                endif
                
                ! Read data directly (TAI93 values as float64)
                call h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, time_arr, dims, hdferr)
                if (hdferr /= 0) then
                    write(LDT_logunit,*)'[ERR] Cannot read data for: ', trim(dataset)
                    deallocate(time_arr)
                    call h5sclose_f(dataspace_id, hdferr)
                    call h5dclose_f(dataset_id, hdferr)
                    ierr = 1
                    return
                endif
                
                ! Clean up
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                
                write(LDT_logunit,*)'[INFO] Successfully read scan time data'
            end subroutine get_dataset_scan_time
            
            ! Function for reading pixel quality flags (uint8)
            subroutine get_dataset_pixel_quality(file_id, dataset, n, m, var2d, ierr)
              use LDT_logMod, only: LDT_logunit
              implicit none
            
              ! Arguments
              integer(HID_T), intent(in) :: file_id
              character(*), intent(in) :: dataset
              integer, intent(out) :: n
              integer, intent(out) :: m
              integer*2, allocatable, intent(out) :: var2d(:,:)  ! Use int*2 for 0-255 range
              integer, intent(out) :: ierr
            
              ! Locals
              integer(HID_T) :: dataset_id, dataspace_id
              logical :: link_exists
              integer :: hdferr
              integer(HSIZE_T) :: dims(2), maxdims(2)
              integer :: rank
              integer*1, allocatable :: raw_data(:,:)  ! uint8 type
              integer :: i, j
            
              ierr = 0
            
              ! Check if dataset exists
              call h5lexists_f(file_id, trim(dataset), link_exists, hdferr)
              if (hdferr /= 0 .or. .not. link_exists) then
                write(LDT_logunit,*)'[ERR] Dataset not found: ', trim(dataset)
                ierr = 1
                return
              endif
            
              ! Open the dataset
              call h5dopen_f(file_id, trim(dataset), dataset_id, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot open dataset: ', trim(dataset)
                ierr = 1
                return
              endif
            
              ! Get the dataspace and dimensions
              call h5dget_space_f(dataset_id, dataspace_id, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot get dataspace for: ', trim(dataset)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              call h5sget_simple_extent_ndims_f(dataspace_id, rank, hdferr)
              if (hdferr /= 0 .or. rank /= 2) then
                write(LDT_logunit,*)'[ERR] Expected 2D dataset, found rank: ', rank
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, hdferr)
              if (hdferr < 0) then
                write(LDT_logunit,*)'[ERR] Cannot get dimensions for: ', trim(dataset)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
            
              n = int(dims(1))
              m = int(dims(2))
              
              ! Allocate arrays
              allocate(raw_data(n, m), stat=hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Memory allocation failed for raw data'
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              allocate(var2d(n, m), stat=hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Memory allocation failed for output array'
                deallocate(raw_data)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Read the raw uint8 data
              call h5dread_f(dataset_id, H5T_NATIVE_CHARACTER, raw_data, dims, hdferr)
              if (hdferr /= 0) then
                write(LDT_logunit,*)'[ERR] Cannot read data for: ', trim(dataset)
                deallocate(raw_data)
                deallocate(var2d)
                call h5sclose_f(dataspace_id, hdferr)
                call h5dclose_f(dataset_id, hdferr)
                ierr = 1
                return
              endif
              
              ! Convert to int*2 to ensure it can hold 0-255 values
              do j = 1, m
                do i = 1, n
                  var2d(i,j) = int(raw_data(i,j))
                enddo
              enddo
              
              ! Clean up
              deallocate(raw_data)
              call h5sclose_f(dataspace_id, hdferr)
              call h5dclose_f(dataset_id, hdferr)
              
              write(LDT_logunit,*)'[INFO] Successfully read pixel quality flag data'
            end subroutine get_dataset_pixel_quality

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
              if (allocated(land_water_frac)) deallocate(land_water_frac)
              if (allocated(snow)) deallocate(snow)
              if (allocated(precip)) deallocate(precip)
              m = 0
              n = 0
              ierr = 1
              return
            end subroutine freeall
        
#else
          ! Dummy version if LDT was compiled w/o HDF5 support.
          write(LDT_logunit,*) '[ERR] get_amsr_l1r called without HDF5 support!'
          write(LDT_logunit,*) '[ERR] Recompile LDT with HDF5 support and try again!'
          call LDT_endrun()
          ierr = 1
#endif
        
        END SUBROUTINE get_amsr_l1r
END MODULE TOOLSUBS_AMSR        