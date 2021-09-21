program MC_KMeans

   implicit none
   
   integer, parameter :: DP = selected_real_kind(15,307) 

   type image
    integer, allocatable, dimension(:,:) :: pixels
    integer :: rows, cols
    real(kind=DP) :: sx, sy, a, b
   end type image

   type kernel
     integer :: cx, cy
   end type kernel

   type centroid
     type(kernel), allocatable, dimension(:) :: kernels
     integer :: nkernels
   end type centroid
  

   type(image)    :: img
   type(centroid) :: ker

   call read_image(img)
   call calc_centroids(img, ker)

contains
  
  subroutine read_image(img)
    type(image), intent(inout) :: img
    integer :: i, j
    character(len=256) :: buffer
    character(len=:), allocatable :: file_image

    if(command_argument_count() < 2 ) then
      write(*,'(a)') 'Not enough arguments in command line.'
      write(*,'(a)') 'Usage: ./MC_KMeans <input_file> <number_cluster>'
      stop
    end if

    call get_command_argument(1,buffer)
    file_image = trim(adjustl(buffer))

    close(10)
    open(10, file=file_image, status='old')

    read(10,*) img%sx, img%sy
    read(10,*) img%rows, img%cols
    allocate(img%pixels(img%rows,img%cols))

    write(*,*) img%rows, img%cols

    do i = 1, img%rows
       read(10,*) (img%pixels(i,j),j=1,img%cols)
    end do
  end subroutine read_image

  subroutine calc_centroids(img, ker)
    type(image), intent(inout)    :: img
    type(centroid), intent(inout) :: ker
    real(kind=DP) :: dist, dmin

    integer, allocatable, dimension(:,:) :: G
    integer, allocatable, dimension(:) :: sumx, sumy, total
    real(kind=DP), allocatable, dimension(:) :: thres
    real(kind=DP), parameter :: eps = 0.000000001_DP
    integer :: i, j, k, m, x, y, n
    logical :: is = .true.
    real(kind=DP), allocatable, dimension(:) :: areas

    character(len=256) :: buffer
    character(len=:), allocatable :: foo


    call get_command_argument(2,buffer)
    foo = trim(adjustl(buffer))
    read(foo,*) ker%nkernels

    allocate(ker%kernels(ker%nkernels))
    allocate(sumx(ker%nkernels), sumy(ker%nkernels), total(ker%nkernels))
    allocate(thres(ker%nkernels))
    allocate(areas(ker%nkernels))

    allocate(G(img%rows,img%cols))

    do k = 1, ker%nkernels
       ker%kernels(k)%cx = int(rand() * img%cols) + 1
       ker%kernels(k)%cy = int(rand() * img%rows) + 1
    end do 

    G = 0
    thres = 0

    n = 10

    img%a = img%sx/img%cols
    img%b = img%sy/img%rows

  do
       if (.not.is) exit

      do i = 1, img%rows
         do j = 1, img%cols
            dmin = 99999.99_DP   
            do k = 1, ker%nkernels
               dist = sqrt(real((ker%kernels(k)%cx - j)**2 + (ker%kernels(k)%cy - i)**2,kind=DP) )
               if(dist < dmin .and. img%pixels(i,j)==1) then
                 dmin = dist
                 G(i,j) = k
               end if
            end do
         end do
      end do


      sumx  = 0
      sumy  = 0
      total = 0

      n = n+1

      do i = 1, img%rows
         write(n,*) (G(i,j),j=1,img%cols)
      end do

      do i = 1, img%rows
         do j = 1, img%cols
             if (img%pixels(i,j) == 1) then
                 sumx(G(i,j)) = sumx(G(i,j)) + j
                 sumy(G(i,j)) = sumy(G(i,j)) + i
                 total(G(i,j)) = total(G(i,j)) + 1
             end if
         end do
      end do

      do k = 1, ker%nkernels
          sumx(k) = sumx(k) / total(k)
          sumy(k) = sumy(k) / total(k)

          areas(k) = img%a * img%b * total(k)
      end do
   
   thres = 0.0

   is = .false.
   do k = 1, ker%nkernels

     x = ker%kernels(k)%cx
     y = ker%kernels(k)%cy

     thres(k) = (x-sumx(k))**2 +  (y-sumy(k))**2
     thres(k) = real(thres(k),kind=DP)

     thres(k) = sqrt(thres(k))
    
     if (thres(k) > eps) then
        is = .true.
     end if

   end do

    ! update
    do k = 1, ker%nkernels
       ker%kernels(k)%cx  = sumx(k)
       ker%kernels(k)%cy  = sumy(k)
    end do

  end do

  do k = 1, ker%nkernels
     write(*,'(a,2x,i0,f10.2,1x,a,2x,i0)') 'Cluster: ', k, areas(k),'mm^2'
  end do

  end subroutine calc_centroids


end program MC_KMeans
