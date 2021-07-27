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
     write(1,'(1000(F16.10,X))')ao_overlap(i,:)
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
     write(3,'(1000(F16.10,X))')H_core(i,:)
  enddo
  close(3)


  logical :: ao_two_e_integral_zero
  double precision :: integral
  
  ! Four index tensor [Non condensed]
  open(4, file = "Four_index_tensor.dat")
  do l=1,ao_num
     do k=1,ao_num
        do j=l,ao_num
           do i=k,ao_num
              if (i>=j) then
                 if( ao_two_e_integral_zero(i,j,k,l) ) cycle
                 integral = get_ao_two_e_integral(i,j,k,l,ao_integrals_map)
                 if (dabs(integral) > ao_integrals_threshold) then
                    write(4,*)i,k,j,l,integral
                 endif
              end if
           enddo
        enddo
     enddo
  enddo
  close(4)

  
  ! Huckel guess
  ! call huckel_guess()
  ! do i = 1, ao_num
  !    write(*,'(1000(F16.10,X))')mo_coef(i,:)
  ! enddo

  

  print*,"############ EXTRACTION DONE  ############"

  ! mo_label = "Guess"
  ! call save_mos

end program extract_rohf_data
