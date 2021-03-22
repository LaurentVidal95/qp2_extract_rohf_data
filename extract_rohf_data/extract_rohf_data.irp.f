program extract_rohf_data
  implicit none
  use map_module
  
  BEGIN_DOC
  ! Extract the overlap matrix, the one-e integrals and the number of alpha and beta orbitals.
  ! Also extract the four-index tensor in a condensed writing (by taking symmetries into account).
  END_DOC

  integer                               :: i,j,k,l
  double precision                      :: accu
  double precision                      :: get_ao_two_e_integral
  ! Overlap matrix
  open(1, file = "overlap_matrix.dat")
  do i = 1, ao_num
     write(1,'(100(F16.10,X))')ao_overlap(i,:)
  enddo
  close(1)

  ! Nums_orbitals
  open(2, file = "nums_orbitals.dat")
  write(2, '(3(I4,X))')ao_num,elec_beta_num,elec_alpha_num-elec_beta_num
  close(2)

  ! H_core
  double precision, allocatable :: H_core(:,:)
  allocate(H_core(ao_num,ao_num))
  do j=1,ao_num
     do i = 1,ao_num
        H_core(i,j) = ao_one_e_integrals(i,j)
     enddo
  enddo
  
  open(3, file = "H_core.dat")
  do i = 1, ao_num
     write(3,'(100(F16.10,X))')H_core(i,:)
  enddo
  close(3)

  ! Four index tensor [Non condensed]
  open(4, file = "Four_index_tensor.dat")
  do i=1,ao_num
     print*, "Extracting: ",i
     do j=1,ao_num
        do k=1,ao_num
           do l=1,ao_num
              accu = get_ao_two_e_integral(i,k,j,l,ao_integrals_map)
              write(4,*)i,j,k,l,accu
           enddo
        enddo
     enddo
  enddo
  close(4)
  
  ! ! Four-index tensor [Condensed]
  ! PROVIDE ao_two_e_integrals_in_map
  
  ! integer(map_size_kind)     :: i8
  ! integer                        :: ii(8), jj(8), kk(8), ll(8), k1
  ! integer(cache_map_size_kind)   :: n_elements_max, n_elements
  ! integer(key_kind), allocatable :: keys(:)
  ! double precision, allocatable  :: values(:)
  ! double precision               :: integral
  
  ! ! !$OMP PARALLEL DEFAULT(NONE)                                      &
  ! !    !$OMP PRIVATE(i,j,l,k1,k,integral,ii,jj,kk,ll,i8,keys,values,n_elements_max, &
  ! !    !$OMP  n_elements,coulomb_s_tmp,coulomb_d_tmp,exchange_s_tmp,exchange_d_tmp)&
  ! !    !$OMP SHARED(ao_num,SCF_density_matrix_ao_alpha,SCF_density_matrix_ao_beta,&
  ! !    !$OMP  ao_integrals_map, coulomb_s, coulomb_d, density_mat_d, density_mat_s)
  
  ! call get_cache_map_n_elements_max(ao_integrals_map,n_elements_max)
  ! allocate(keys(n_elements_max), values(n_elements_max))

  ! open(4, file = "Four_index_tensor.dat")
  
  ! ! !$OMP DO SCHEDULE(static,1)
  ! do i8 = 0_8, ao_integrals_map%map_size
  !    n_elements = n_elements_max
  !    call get_cache_map(ao_integrals_map, i8, keys, values, n_elements)
  !    ! keys(n_elements) = i,j,k,l ; values(n_elements) = integral
  !    do k1 = 1,n_elements
  !       ! get the i,j,k,l with all permutations for the key == keys(k1)
  !       call two_e_integrals_index_reverse(kk,ii,ll,jj,keys(k1)) 
  !       i = ii(1); j = jj(1); k = kk(1); l = ll(1);
  !       integral = values(k1)
  !       write(4,*)i,j,k,l,integral
  !    enddo
  ! enddo
  
  ! deallocate(keys,values)

  ! Huckel guess
  call huckel_guess()
  open(5, file = "MOs_huckel_guess.dat")
  do i = 1, ao_num
     write(5,'(100(F16.10,X))')mo_coef(i,:)
  enddo
  close(5)
  

  print*,"############ EXTRACTION DONE  ############"

  mo_label = "Guess"
  call save_mos

  ! !$OMP END PARALLEL
end program extract_rohf_data
