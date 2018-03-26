subroutine MPP_TYPE_CREATE_(field, array_of_subsizes, array_of_starts, &
                            newtype_idx)
    MPP_TYPE_, intent(in) :: field(:,:,:)
    integer, intent(in) :: array_of_subsizes(:)
    integer, intent(in) :: array_of_starts(:)
    integer, intent(out) :: newtype_idx

    integer :: newtype      ! MPI datatype ID
    integer :: i            ! Type loop index

    if (.NOT. module_is_initialized) &
        call mpp_error(FATAL, 'MPP_TYPE_CREATE_: You must first call mpp_init.')

    if (current_clock .NE. 0) &
        call SYSTEM_CLOCK(start_tick)

    if (verbose) &
        call mpp_error(NOTE, 'MPP_TYPE_CREATE_: &
                             &using MPI_Type_create_subarray...')

    ! Check if type already exists
    ! NOTE: Something like a hash table (or even a linked list) would be better
    !       here; this is mainly just for convenience, and to resemble existing
    !       methods (eg peset).
    do i = 0, current_datatype_num
        if (datatypes(i)%ndims /= rank(field)) cycle
        if (any(datatypes(i)%sizes /= shape(field))) cycle
        if (any(datatypes(i)%subsizes /= array_of_subsizes)) cycle
        if (any(datatypes(i)%starts /= array_of_starts)) cycle
        if (datatypes(i)%etype /= MPI_TYPE_) cycle

        ! If all parameters match, then the datatype already exists.
        ! We can return its index and exit.
        datatypes(i)%counter = datatypes(i)%counter + 1
        newtype_idx = i
        return
    end do

    ! If the type does not exist, then create a new type
    call MPI_Type_create_subarray( &
        rank(field), &
        shape(field), &
        array_of_subsizes, &
        array_of_starts, &
        MPI_ORDER_FORTRAN, &
        MPI_TYPE_, &
        newtype, &
        error &
    )

    ! Register on the MPI runtime
    call MPI_Type_commit(newtype, error)

    ! Update the datatype list
    ! TODO: Find the lowest empty element in list
    newtype_idx = current_datatype_num + 1
    if (newtype_idx > current_datatype_max)  then
        !call expand_datatypes()    ! TODO: Implement this!
        call mpp_error(FATAL, 'MPP_TYPE_CREATE_: too many types!')
    end if

    ! Append new entry
    allocate(datatypes(newtype_idx)%sizes(rank(field)))
    allocate(datatypes(newtype_idx)%subsizes(rank(field)))
    allocate(datatypes(newtype_idx)%starts(rank(field)))

    ! Populate values
    datatypes(newtype_idx)%counter = 1
    datatypes(newtype_idx)%ndims = rank(field)
    datatypes(newtype_idx)%sizes = shape(field)
    datatypes(newtype_idx)%subsizes = array_of_subsizes
    datatypes(newtype_idx)%starts = array_of_starts
    datatypes(newtype_idx)%etype = MPI_TYPE_
    datatypes(newtype_idx)%id = newtype

    ! Update the size tracker
    current_datatype_num = newtype_idx

    if (current_clock .NE. 0) &
        call increment_current_clock(EVENT_TYPE_CREATE, MPP_TYPE_BYTELEN_)

end subroutine MPP_TYPE_CREATE_
