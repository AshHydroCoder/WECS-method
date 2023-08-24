
program WECS
    implicit none
    integer :: n,i,j
    character(len=500) :: Rainfall_path, Q_path, Area_path
    real, dimension (:,:), allocatable :: R, Q, A, P, log_Q, TR, PR, INV, ID, log_H, H, Q_SIM
    !ID = Intermediate matrix, R = Rainfall matrix, Q = Observed Q, A = Area matrix, P = matrix of log of variables [1 logR logA]
    !TR - Transposed, PR - Product of transpose and original, INV: inverse of the PR, log_H = parameter matrix in log form, H = final parameter matrix, Q_SIM = simulated Q 
    interface
        subroutine Matrix_Trans(In1, In2)
        implicit none
        integer :: nrow, ncol, i, j
        real, dimension (:,:), allocatable :: In1, In2
        end subroutine Matrix_Trans      
    
        subroutine Matrix_Multi(In1,In2,Output)
        implicit none
        integer :: Row1, Col, Col2
        real, dimension (:,:), allocatable :: In1, In2, Output
        end subroutine Matrix_Multi
    
        subroutine MatrixInverse(Input1, Output1)
        implicit none    
        integer :: nrow, ncol
            real, dimension(:,:), allocatable :: Input1, Output1
        end subroutine MatrixInverse    
    end interface
    
    print *, "Delete Quotation Mark when entering file path"
    print *, "Enter File path to Rainfall data"

    open(1, file = "E:\Fortran\Exercies\Assignment\R50.txt")


    open(2,file = "E:\Fortran\Exercies\Assignment\Q50.txt")

    open(3,file = "E:\Fortran\Exercies\Assignment\A.txt")

    print *, "Enter Number of ordinates"
    read *, n
    allocate (R(n,1),Q(n,1),P(n,3),A(n,1), log_Q(n,1), H(3,1), Q_SIM(n,2))
    do i = 1,n
        read(1,*) R(i,1)
    end do
    do i = 1,n
        read(2,*) Q(i,1)
    end do
    do i = 1,n
        read(3,*) A(i,1)
    end do
    do i = 1,n
        log_Q(i,1) = log(Q(i,1))
    end do
!Writing[1  logR    logA] matrix
    do i = 1,n
        P(i,1) = 1
        P(i,2) = log(R(i,1))
        P(i,3) = log(A(i,1))
    end do
    call Matrix_Trans(P,TR)
    call Matrix_Multi(TR,P,PR)
    call Matrix_Multi(TR,log_Q,ID)
    call MatrixInverse(PR,INV)
    call Matrix_Multi(INV,ID,log_H)
    
    H(1,1) = exp(log_H(1,1))
    H(2,1) = log_H(2,1)
    H(3,1) = log_H(3,1)
    print *, "A = ", H(1,1)
    print *, "B = ", H(2,1)
    print *, "G = ", H(3,1)
!Final Results
    do i = 1,n
        Q_SIM(i,1) = Q(i,1)
        Q_SIM(i,2) = H(1,1)*(R(i,1)**H(2,1))*(A(i,1)**H(3,1))
    end do
    

    write(*, '(a, a, a)', advance='no') "Q measured", achar(9), "Q simulated"
    write(*,*)

    
    do i = 1,n
    print *, Q_SIM(i,:)
    end do 
    
    
    end program WECS
    
subroutine Matrix_Trans(In1, In2)
    implicit none
    integer :: nrow, ncol, i, j
    real, dimension (:,:), allocatable :: In1, In2
    nrow = size(In1, dim = 1)
    ncol = size(In1, dim = 2)
    allocate (In2(ncol,nrow))

    do j = 1, ncol
        do i = 1,nrow
            In2(j,i) = In1(i,j)
        end do
    end do
    end subroutine Matrix_Trans

subroutine Matrix_Multi(In1,In2,Output)
    implicit none
    integer :: Row1, Col, Col2, i, j, k
    real, dimension (:,:), allocatable :: In1, In2, Output
    Row1 = size(In1, dim = 1)
    Col = size(In1, dim = 2)
    Col2 = size(In2, dim = 2)
    allocate (Output(Row1,Col2))

    do i = 1,Row1
        do j = 1,Col2
            Output(i,j) = 0.0
            do k = 1,Col
                Output(i,j) = Output(i,j) + In1(i,k)*In2(k,j)
            end do
        end do
    end do
    end subroutine Matrix_Multi
    
    
    
subroutine MatrixInverse(Input1, Output1)
    implicit none
    integer :: nrow, ncol
    real, dimension(:,:), allocatable :: Input1, Output1
    integer :: i, j, k
    real :: X11, Xik
    nrow = size(Input1, dim = 1)
    ncol = nrow
    allocate(Output1(nrow, ncol))
    
    !writing an Identity matrix
    do i = 1, nrow
        do j = 1, ncol
            if (i == j) then
                Output1(i,j) = 1
            else
                Output1(i,j) = 0
            end if
        end do
    end do
    
    do k = 1, nrow
        X11 = Input1(k,k)
        do j = 1, ncol
            Input1(k,j) = Input1(k,j)/X11
            Output1(k,j) = Output1(k,j)/X11
        end do
          
        do i = 1, nrow
            if (i /= k) then
                Xik = Input1(i,k)
                do j = 1, ncol
                  Input1(i,j) = Input1(i,j) - Xik * Input1(k,j)
                  Output1(i,j) = Output1(i,j) - Xik * Output1(k,j)
                end do
            end if
        end do
    end do
    end subroutine MatrixInverse 
