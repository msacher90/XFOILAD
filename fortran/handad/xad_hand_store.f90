!***********************************************************************
!   XFOILAD - storage of the hand-written adjoints (handad/)
!
!   BLSOLV is differentiated in split mode: BLSOLV_FWD solves the Newton
!   system and pushes one record (factors left in VA / VM, blocks VB / VZ,
!   solution, inputs), BLSOLV_BWD pops it.  One record per Newton iteration
!   replayed by the adjoint; the records are kept allocated between calls.
!***********************************************************************
module xad_hand_store
    implicit none
    type blsolv_rec
        integer :: nsys = 0, ivte1 = 0, ivz = 0
        real :: vacc(3), vz(3, 2)
        real, allocatable :: va(:, :, :), vb(:, :, :), vm(:, :, :), d(:, :, :)
        real, allocatable :: va0(:, :, :), vm0(:, :, :), vdel0(:, :, :)
    end type blsolv_rec
    type(blsolv_rec), allocatable, save :: bls_stack(:)
    integer, save :: bls_top = 0
    !> BAKSUB_BWD: columns y = A^-T x_bar and x of every solve with the
    !> factored AIJ, applied at once by LUDCMP_B (A_bar = A_bar - Y X^T)
    real, allocatable, save :: bak_y(:, :), bak_x(:, :)
    integer, save :: bak_m = 0
contains

    !> Empties the stacks (start of an adjoint sweep).
    subroutine xad_hand_reset()
        bls_top = 0
        bak_m = 0
    end subroutine xad_hand_reset

    !> Room for one more BAKSUB column of length n; returns its index.
    function bak_push(n) result(k)
        integer, intent(in) :: n
        integer :: k
        real, allocatable :: t(:, :)
        if (.not. allocated(bak_y)) then
            allocate(bak_y(n, 256), bak_x(n, 256))
        elseif (size(bak_y, 1) < n) then
            deallocate(bak_y, bak_x)
            allocate(bak_y(n, 256), bak_x(n, 256))
            bak_m = 0
        endif
        if (bak_m == size(bak_y, 2)) then
            allocate(t(size(bak_y, 1), 2 * bak_m))
            t(:, 1:bak_m) = bak_y(:, 1:bak_m)
            call move_alloc(t, bak_y)
            allocate(t(size(bak_x, 1), 2 * bak_m))
            t(:, 1:bak_m) = bak_x(:, 1:bak_m)
            call move_alloc(t, bak_x)
        endif
        bak_m = bak_m + 1
        k = bak_m
    end function bak_push

    !> Pushes a record for a system of nsy stations, returns its index.
    function bls_push(nsy) result(k)
        integer, intent(in) :: nsy
        integer :: k
        type(blsolv_rec), allocatable :: tmp(:)
        if (.not. allocated(bls_stack)) allocate(bls_stack(16))
        if (bls_top == size(bls_stack)) then
            allocate(tmp(2 * size(bls_stack)))
            tmp(1:bls_top) = bls_stack(1:bls_top)
            call move_alloc(tmp, bls_stack)
        endif
        bls_top = bls_top + 1
        k = bls_top
        if (bls_stack(k)%nsys /= nsy .or. .not. allocated(bls_stack(k)%vm)) then
            if (allocated(bls_stack(k)%vm)) then
                deallocate(bls_stack(k)%va, bls_stack(k)%vb, bls_stack(k)%vm, bls_stack(k)%d, &
                           bls_stack(k)%va0, bls_stack(k)%vm0, bls_stack(k)%vdel0)
            endif
            allocate(bls_stack(k)%va(3, 2, nsy), bls_stack(k)%vb(3, 2, nsy), bls_stack(k)%vm(3, nsy, nsy), &
                     bls_stack(k)%d(3, 2, nsy), bls_stack(k)%va0(3, 2, nsy), bls_stack(k)%vm0(3, nsy, nsy), &
                     bls_stack(k)%vdel0(3, 2, nsy))
            bls_stack(k)%nsys = nsy
        endif
    end function bls_push

end module xad_hand_store
