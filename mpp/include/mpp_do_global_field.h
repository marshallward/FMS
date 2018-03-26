    subroutine MPP_DO_GLOBAL_FIELD_3D_( domain, local, global, tile, ishift, jshift, flags, default_data)
!get a global field from a local field
!local field may be on compute OR data domain
      type(domain2D), intent(in)    :: domain
      integer, intent(in)           :: tile, ishift, jshift
      MPP_TYPE_, intent(in),  contiguous, target :: local(:,:,:)
      MPP_TYPE_, intent(out), contiguous, target :: global(domain%x(tile)%global%begin:,domain%y(tile)%global%begin:,:)
      integer, intent(in), optional :: flags
      MPP_TYPE_, intent(in), optional :: default_data

      integer :: i, j, k, m, n, nd, nwords, lpos, rpos, ioff, joff, from_pe, root_pe, tile_id
      integer :: ke, isc, iec, jsc, jec, is, ie, js, je, nword_me
      integer :: ipos, jpos
      logical :: xonly, yonly, root_only, global_on_this_pe
      MPP_TYPE_ :: clocal ((domain%x(1)%compute%size+ishift)    *(domain%y(1)%compute%size+jshift)    *size(local,3))
      MPP_TYPE_ :: cremote((domain%x(1)%compute%max_size+ishift)*(domain%y(1)%compute%max_size+jshift)*size(local,3))
      integer :: stackuse
      character(len=8) :: text

      ! Alltoallw vectors
      integer, dimension(:), allocatable :: sendcounts(:), recvcounts(:)
      integer, dimension(:), allocatable :: sdispls(:), rdispls(:)
      integer, dimension(:), allocatable :: sendtypes(:), recvtypes(:)
      integer, dimension(3) :: array_of_subsizes, array_of_starts
      integer :: isize, jsize, ihalo, jhalo
      integer :: isg, jsg
      integer :: ierr
      integer, allocatable :: pelist(:)
      integer :: npes
      MPP_TYPE_, dimension(:), pointer :: plocal, pglobal

      ! Cray pointers
      pointer( ptr_local,  clocal  ) 
      pointer( ptr_remote, cremote )

      stackuse = size(clocal(:))+size(cremote(:))
      if( stackuse.GT.mpp_domains_stack_size )then
          write( text, '(i8)' )stackuse
          call mpp_error( FATAL, &
               'MPP_DO_GLOBAL_FIELD user stack overflow: call mpp_domains_set_stack_size('//trim(text)//') from all PEs.' )
      end if
      mpp_domains_stack_hwm = max( mpp_domains_stack_hwm, stackuse )

      ptr_local  = LOC(mpp_domains_stack)
      ptr_remote = LOC(mpp_domains_stack(size(clocal(:))+1))

      if( .NOT.module_is_initialized )call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: must first call mpp_domains_init.' )

      xonly = .FALSE.
      yonly = .FALSE.
      root_only = .FALSE.
      if( PRESENT(flags) ) then
         xonly = BTEST(flags,EAST)
         yonly = BTEST(flags,SOUTH)
         if( .NOT.xonly .AND. .NOT.yonly )call mpp_error( WARNING,  &
                   'MPP_GLOBAL_FIELD: you must have flags=XUPDATE, YUPDATE or XUPDATE+YUPDATE' )
         if(xonly .AND. yonly) then
            xonly = .false.; yonly = .false.
         endif
         root_only = BTEST(flags, ROOT_GLOBAL)
         if( (xonly .or. yonly) .AND. root_only ) then
            call mpp_error( WARNING, 'MPP_GLOBAL_FIELD: flags = XUPDATE+GLOBAL_ROOT_ONLY or ' // &
                 'flags = YUPDATE+GLOBAL_ROOT_ONLY is not supported, will ignore GLOBAL_ROOT_ONLY' )     
            root_only = .FALSE.
         endif
      endif
    
      global_on_this_pe =  .NOT. root_only .OR. domain%pe == domain%tile_root_pe
      ipos = 0; jpos = 0
      if(global_on_this_pe ) then      
         if(size(local,3).NE.size(global,3) ) call mpp_error( FATAL, &
              'MPP_GLOBAL_FIELD: mismatch of third dimension size of global and local')
         if( size(global,1).NE.(domain%x(tile)%global%size+ishift) .OR. size(global,2).NE.(domain%y(tile)%global%size+jshift))then
            if(xonly) then
               if(size(global,1).NE.(domain%x(tile)%global%size+ishift) .OR. &
                   size(global,2).NE.(domain%y(tile)%compute%size+jshift)) &
                  call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: incoming arrays do not match domain for xonly global field.' )
               jpos = -domain%y(tile)%compute%begin + 1
            else if(yonly) then
               if(size(global,1).NE.(domain%x(tile)%compute%size+ishift) .OR. &
                   size(global,2).NE.(domain%y(tile)%global%size+jshift)) &
                  call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: incoming arrays do not match domain for yonly global field.' )
               ipos = -domain%x(tile)%compute%begin + 1
            else
               call mpp_error( FATAL, 'MPP_GLOBAL_FIELD: incoming arrays do not match domain.' )
            endif
         endif
      endif

      if( size(local,1).EQ.(domain%x(tile)%compute%size+ishift) .AND. size(local,2).EQ.(domain%y(tile)%compute%size+jshift) )then
         !local is on compute domain
         ioff = -domain%x(tile)%compute%begin + 1
         joff = -domain%y(tile)%compute%begin + 1
      else if( size(local,1).EQ.(domain%x(tile)%memory%size+ishift) .AND. size(local,2).EQ.(domain%y(tile)%memory%size+jshift) )then
         !local is on data domain
         ioff = -domain%x(tile)%data%begin + 1
         joff = -domain%y(tile)%data%begin + 1
      else
         call mpp_error( FATAL, 'MPP_GLOBAL_FIELD_: incoming field array must match either compute domain or memory domain.' )
      end if

      ke  = size(local,3)
      isc = domain%x(tile)%compute%begin; iec = domain%x(tile)%compute%end+ishift
      jsc = domain%y(tile)%compute%begin; jec = domain%y(tile)%compute%end+jshift

      nword_me = (iec-isc+1)*(jec-jsc+1)*ke

! make contiguous array from compute domain
      m = 0
      if(global_on_this_pe) then
         !z1l: initialize global = 0 to support mask domain
         if(PRESENT(default_data)) then
            global = default_data
         else
#ifdef LOGICAL_VARIABLE
            global = .false.
#else
            global = 0
#endif
         endif

         do k = 1, ke
            do j = jsc, jec
               do i = isc, iec
                  m = m + 1
                  clocal(m) = local(i+ioff,j+joff,k)
                  global(i+ipos,j+jpos,k) = clocal(m) !always fill local domain directly
               end do
            end do
         end do
      else
         do k = 1, ke
            do j = jsc, jec
               do i = isc, iec
                  m = m + 1
                  clocal(m) = local(i+ioff,j+joff,k)
               end do
            end do
         end do
      endif

! if there is more than one tile on this pe, then no decomposition for all tiles on this pe, so we can just return
      if(size(domain%x(:))>1) then
         !--- the following is needed to avoid deadlock.
         if( tile == size(domain%x(:)) ) call mpp_sync_self( )
         return
      end if

      root_pe = mpp_root_pe()

!fill off-domains (note loops begin at an offset of 1)
      if( xonly )then
          nd = size(domain%x(1)%list(:))
          do n = 1,nd-1
             lpos = mod(domain%x(1)%pos+nd-n,nd)
             rpos = mod(domain%x(1)%pos   +n,nd)
             from_pe = domain%x(1)%list(rpos)%pe
             rpos = from_pe - root_pe ! for concurrent run, root_pe may not be 0.
             if (from_pe == NULL_PE) then
                nwords = 0
             else
                nwords = (domain%list(rpos)%x(1)%compute%size+ishift) &
                       * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
             endif
           ! Force use of scalar, integer ptr interface
             call mpp_transmit( put_data=clocal(1), plen=nword_me, to_pe=domain%x(1)%list(lpos)%pe, &
                                get_data=cremote(1), glen=nwords, from_pe=from_pe )
             m = 0
             if (from_pe /= NULL_PE) then
                is = domain%list(rpos)%x(1)%compute%begin; ie = domain%list(rpos)%x(1)%compute%end+ishift
                do k = 1, ke
                   do j = jsc, jec
                      do i = is, ie
                         m = m + 1
                         global(i,j+jpos,k) = cremote(m)
                      end do
                   end do
                end do
             endif
             call mpp_sync_self()  !-ensure MPI_ISEND is done.
          end do
      else if( yonly )then
          nd = size(domain%y(1)%list(:))
          do n = 1,nd-1
             lpos = mod(domain%y(1)%pos+nd-n,nd)
             rpos = mod(domain%y(1)%pos   +n,nd)
             from_pe = domain%y(1)%list(rpos)%pe
             rpos = from_pe - root_pe
             if (from_pe == NULL_PE) then
                nwords = 0
             else
                nwords = (domain%list(rpos)%x(1)%compute%size+ishift) &
                       * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
             endif
           ! Force use of scalar, integer pointer interface
             call mpp_transmit( put_data=clocal(1), plen=nword_me, to_pe=domain%y(1)%list(lpos)%pe, &
                                get_data=cremote(1), glen=nwords, from_pe=from_pe )
             m = 0
             if (from_pe /= NULL_PE) then
                 js = domain%list(rpos)%y(1)%compute%begin; je = domain%list(rpos)%y(1)%compute%end+jshift
                 do k = 1,ke
                    do j = js, je
                       do i = isc, iec
                          m = m + 1
                          global(i+ipos,j,k) = cremote(m)
                       end do
                    end do
                 end do
             endif
             call mpp_sync_self()  !-ensure MPI_ISEND is done.
          end do
      else
         tile_id = domain%tile_id(1)
         nd = size(domain%list(:))

         if(root_only) then
            ! TODO: Trial MPI_alltoallw code
            ! NOTE: Defining types makes the prior clocal/cglobal setup code
            ! unnecessary, so this ought to be done earlier.
            ! If this all works out then we will do all this a bit more
            ! strategically.

            ! Global index
            isg = domain%x(1)%global%begin
            jsg = domain%y(1)%global%begin

            ! Prepare the MPI buffers
            !   Initialize with inert values
            !   TODO: What type to use?
            !         Do I even need to initalise types?
            allocate(sendcounts(0:nd-1))
            allocate(sdispls(0:nd-1))
            allocate(sendtypes(0:nd-1))
            sendcounts(:) = 0
            sdispls(:) = 0
            sendtypes(:) = 0    ! MPI_BYTE

            allocate(recvcounts(0:nd-1))
            allocate(rdispls(0:nd-1))
            allocate(recvtypes(0:nd-1))
            recvcounts(:) = 0
            rdispls(:) = 0
            recvtypes(:) = 0    ! MPI_BYTE

            if (domain%pe .NE. domain%tile_root_pe) then
                ! Rank/PE test
                isize = iec - isc + 1
                jsize = jec - jsc + 1
                ihalo = isc + ioff - 1
                jhalo = jsc + joff - 1

                ! Send one subarray to rank 0
                sendcounts(0) = 1
                sdispls(0) = 0

                array_of_subsizes = [isize, jsize, size(local, 3)]
                array_of_starts = [ihalo, jhalo, 0]

                call mpp_type_create( &
                    local, &
                    array_of_subsizes, &
                    array_of_starts, &
                    sendtypes(0) &
                )
            else
            !if (domain%pe .EQ. domain%tile_root_pe) then
                ! Receive subarrays from all other ranks

                ! Each rank receives one MPI subarray type
                !recvcounts(:) = 1

                !do n = 0, nd - 1
                do n = 1, nd - 1
                    ! Each rank receives one MPI subarray type
                    recvcounts(n) = 1
                    rdispls(n) = 0

                    ! Index to rank
                    rpos = mod(domain%pos + n, nd)
                    is = domain%list(rpos)%x(1)%compute%begin
                    ie = domain%list(rpos)%x(1)%compute%end + ishift
                    js = domain%list(rpos)%y(1)%compute%begin
                    je = domain%list(rpos)%y(1)%compute%end + jshift

                    isize = ie - is + 1
                    jsize = je - js + 1

                    array_of_subsizes = [isize, jsize, ke]
                    array_of_starts = [is - isg, js - jsg, 0]

                    call mpp_type_create( &
                        global, &
                        array_of_subsizes, &
                        array_of_starts, &
                        recvtypes(n) &
                    )
                end do
            end if

            npes = mpp_get_domain_npes(domain)
            allocate(pelist(npes))
            call mpp_get_pelist(domain, pelist)

            plocal(1:size(local)) => local
            pglobal(1:size(global)) => global

            call mpp_alltoall(plocal, sendcounts, sdispls, sendtypes, &
                              pglobal, recvcounts, rdispls, recvtypes, &
                              pelist)

            plocal => null()
            pglobal => null()

            ! Point-to-point method

            !if(domain%pe .NE. domain%tile_root_pe) then
            !   call mpp_send( clocal(1), plen=nword_me, to_pe=domain%tile_root_pe, tag=COMM_TAG_1 )
            !else
            !   do n = 1,nd-1
            !      rpos = mod(domain%pos+n,nd)
            !      if( domain%list(rpos)%tile_id(1) .NE. tile_id ) cycle
            !      nwords = (domain%list(rpos)%x(1)%compute%size+ishift) * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
            !      call mpp_recv(cremote(1), glen=nwords, from_pe=domain%list(rpos)%pe, tag=COMM_TAG_1 )
            !      m = 0
            !      is = domain%list(rpos)%x(1)%compute%begin; ie = domain%list(rpos)%x(1)%compute%end+ishift
            !      js = domain%list(rpos)%y(1)%compute%begin; je = domain%list(rpos)%y(1)%compute%end+jshift

            !      do k = 1,ke
            !         do j = js, je
            !            do i = is, ie
            !               m = m + 1
            !               global(i,j,k) = cremote(m)
            !            end do
            !         end do
            !      end do
            !   end do
            !endif

            !if (domain%pe .EQ. domain%tile_root_pe) then
            !    do n = 1, nd - 1
            !        call mpp_type_free(recvtypes(n))
            !    end do
            !else
            !    call mpp_type_free(sendtypes(0))
            !end if

            ! Cleanup
            deallocate(pelist)
            deallocate(sendcounts, sdispls, sendtypes)
            deallocate(recvcounts, rdispls, recvtypes)
         else
            do n = 1,nd-1
               lpos = mod(domain%pos+nd-n,nd)
               if( domain%list(lpos)%tile_id(1).NE. tile_id ) cycle ! global field only within tile
               call mpp_send( clocal(1), plen=nword_me, to_pe=domain%list(lpos)%pe, tag=COMM_TAG_2 )
            end do
            do n = 1,nd-1
               rpos = mod(domain%pos+n,nd)
               if( domain%list(rpos)%tile_id(1) .NE. tile_id ) cycle ! global field only within tile
               nwords = (domain%list(rpos)%x(1)%compute%size+ishift) * (domain%list(rpos)%y(1)%compute%size+jshift) * ke
               call mpp_recv( cremote(1), glen=nwords, from_pe=domain%list(rpos)%pe, tag=COMM_TAG_2 )
               m = 0
               is = domain%list(rpos)%x(1)%compute%begin; ie = domain%list(rpos)%x(1)%compute%end+ishift
               js = domain%list(rpos)%y(1)%compute%begin; je = domain%list(rpos)%y(1)%compute%end+jshift

               do k = 1,ke
                  do j = js, je
                     do i = is, ie
                        m = m + 1
                        global(i,j,k) = cremote(m)
                     end do
                  end do
               end do
            end do
         endif
      end if

      call mpp_sync_self()
          
      return
    end subroutine MPP_DO_GLOBAL_FIELD_3D_
