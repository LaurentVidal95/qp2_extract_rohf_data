program extract_rohf_data
  implicit none
  
  BEGIN_DOC
  ! Print the MO coeffs in a file
  END_DOC

  integer                               :: i,j

  open(1, file = "mo_coeffs.dat")
  do i = 1, ao_num
     write(1,'(1000(F16.10,X))')mo_coef(i,:)
  enddo

end program extract_rohf_data
