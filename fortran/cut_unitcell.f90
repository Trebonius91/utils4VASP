!
!    cut_unitcell: cut an arbitrary unit cell from a 
!       given large surface cell
!    Part of VASP4CLINT
!     Julien Steffen, 2023 (julien.steffen@fau.de)
!

program cut_unitcell
implicit none
integer::i,j,k 
real(kind=8)::coord_scale
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
real(kind=8)::a_vec_new(3),b_vec_new(3)
character(len=2)::el_name
character(len=80)::arg_a
integer::natoms,natoms_chosen
integer::ind1,ind2,ind3,ind4
character(len=1),allocatable::active(:,:)
real(kind=8),allocatable::xyz_ref(:,:)
real(kind=8)::pos1(2),pos2(2),pos3(2),pos4(2)
real(kind=8)::vec1(2),vec2(2),origin(3)
real(kind=8)::vec_final(2,2),pos_act(2)
integer::vec_inds(6,2)
integer::readstat
logical::read_span
real(kind=8)::angle,arg
real(kind=8)::pi
real(kind=8)::lambda,mu
real(kind=8)::x_vec(2)
real(kind=8)::unit_mat(2,2),inv_mat(2,2)
real(kind=8)::par_x,par_y,x_act,y_act
real(kind=8)::x_check,y_check,dist
real(kind=8)::xyz_chosen(10000,3)
character(len=1)::active_chosen(10000,3)

pi=3.14159265359d0

write(*,*)
write(*,*) "PROGRAM cut_unitcell: cutting an (almost) arbitrary unit"
write(*,*) " cell from a large given surface."
write(*,*) "A POSCAR file must be given with a large reference unit cell."
write(*,*) " This file must have cartesian coordinates."
write(*,*) "Give the numbers of three atoms that shall span the plane "
write(*,*) " of the cut out unit cell with the following command:"
write(*,*) " -span_atoms=[indices] "
write(*,*) "After execution, a file named POSCAR_cut with the cut unit "
write(*,*)"  cell is written."
!
!     First, read in the given POSCAR of the large surface
!

open(unit=16,file="POSCAR",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file 'POSCAR' could not been found!"
   stop
end if       
read(16,*)
read(16,*) coord_scale
read(16,*) a_vec(:)
read(16,*) b_vec(:)
read(16,*) c_vec(:)
read(16,*) el_name
read(16,*) natoms
read(16,*) 
read(16,*)
allocate(xyz_ref(natoms,3))
allocate(active(natoms,3))
do i=1,natoms
   read(16,*) xyz_ref(i,:),active(i,:)
end do
close(16)


!
!     Second, read indices of chosen atoms from command -span_atoms
!
read_span=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg_a)
   if (trim(arg_a(1:12))  .eq. "-span_atoms=") then
      read(arg_a(13:),*,iostat=readstat) ind1,ind2,ind3
      read_span=.true.
   end if
end do

if (.not. read_span) then
   write(*,*) "Please give the command -span_atoms to define the edged of the unitcell!"
   stop
end if        
if (ind1 .le. 0 .or. ind1 .gt. natoms) then
   write(*,*) "Please give a valid number for the first atom!"
   stop
end if        
if (ind2 .le. 0 .or. ind2 .gt. natoms) then
   write(*,*) "Please give a valid number for the second atom!"
   stop
end if   
if (ind3 .le. 0 .or. ind3 .gt. natoms) then
   write(*,*) "Please give a valid number for the third atom!"
   stop
end if
!
!     Now locate all atoms given in atom_inds.txt     
!

pos1=xyz_ref(ind1,1:2)
pos2=xyz_ref(ind2,1:2)
pos3=xyz_ref(ind3,1:2)


!
!     Calculate the unit cell vectors 
!     v1: atom2-atom1
!     v2: atom3-atom1
!     Origin: atom1
!   

vec1=pos2-pos1
vec2=pos3-pos1

origin(1:2)=pos1
origin(3)=0.d0

write(*,*) "Unit cell vector 1: ", vec1
write(*,*) "Unit cell vector 2: ", vec2
write(*,*) "origin",origin(1:2)
!
!     Determine all atoms that are in the defined unit cell
!     Solve linear system of equations for each atom in the system
!     Subtract position of atom 1 in the unit cell to shift them to the 
!     right position
!
!
!     Calculate inverse matrix of coordinate vectors
!
unit_mat(:,1)=vec1
unit_mat(:,2)=vec2

call matinv2(unit_mat,inv_mat)
natoms_chosen=0
outer: do i=1,natoms
   
   x_vec=matmul(inv_mat,xyz_ref(i,1:2)-origin(1:2))   

!   x_act=xyz_ref(i,1)
!   y_act=xyz_ref(i,2)
!   x_act=x_act-origin(1)
!   y_act=y_act-origin(2)

!   lambda=x_act/vec1(1)-((y_act*vec1(1)-x_act*vec1(2))/(vec1(1)*vec2(2)-&
!               &vec1(2)*vec2(1)))*(vec2(1)/vec1(1))
!   lambda=(vec1(1)*y_act-vec1(2)*x_act)/(vec2(2)*vec1(1)-vec2(1)*vec1(2))
!   mu=(y_act*vec1(1)-x_act*vec1(2))/(vec2(2)*vec1(1)-vec2(1)*vec1(2))

!   if (x_vec(1)+x_vec(2) .lt. 1.0d0 .and. x_vec(1)+x_vec(2) .gt. 0.0d0) then
   do while(x_vec(1) .lt. 1d-8) 
      x_vec(1)=x_vec(1)+1.d0
   end do
   do while(x_vec(2) .lt. 1d-8) 
      x_vec(2)=x_vec(2)+1.d0
   end do
      

   if ((x_vec(1) .ge. 0.d0) .and. (x_vec(1) .lt. 1.0001d0) .and. (x_vec(2) .ge. 0.d0) &
              &  .and. (x_vec(2) .lt. 1.0001d0)) then
      natoms_chosen=natoms_chosen+1
     ! xyz_chosen(natoms_chosen,:)=xyz_ref(i,:)-origin(:)
      do j=1,natoms_chosen-1
         dist=sqrt((xyz_chosen(j,1)-x_vec(1))**2+(xyz_chosen(j,2)-x_vec(2))**2) 
         if (dist .lt. 0.0001d0) then
            natoms_chosen=natoms_chosen-1
            cycle outer
         end if
      end do
      xyz_chosen(natoms_chosen,1)=x_vec(1)
      xyz_chosen(natoms_chosen,2)=x_vec(2)
      xyz_chosen(natoms_chosen,3)=xyz_ref(i,3)/c_vec(3)
      active_chosen(natoms_chosen,:)=active(i,:)
   end if
end do outer
!
!     Write resulting POSCAR with new unit cell
!
open(unit=27,file="POSCAR_cut",status="replace")
write(27,*) "New unit cell, build by cut_unitcell.f90"
write(27,*) coord_scale
a_vec_new(1:2)=vec1
a_vec_new(3)=0
b_vec_new(1:2)=vec2
b_vec_new(3)=0
write(27,*) vec1(1),vec1(2),0.d0
write(27,*) vec2(1),vec2(2),0.d0
write(27,*) 0.d0,0.d0,c_vec(3)
write(27,*) el_name
write(27,*) natoms_chosen
write(27,*) "Selective "
write(27,*) "Cartesian "
do i=1,natoms_chosen
   write(27,*) xyz_chosen(i,1)*a_vec_new(1)+xyz_chosen(i,2)*b_vec_new(1)+&
      &      xyz_chosen(i,3)*c_vec(1),xyz_chosen(i,1)*a_vec_new(2)+xyz_chosen(i,2)*b_vec_new(2)+&
      &      xyz_chosen(i,3)*c_vec(2),xyz_chosen(i,1)*a_vec_new(3)+xyz_chosen(i,2)*b_vec_new(3)+&
      &      xyz_chosen(i,3)*c_vec(3)
end do


close(27)

write(*,*)
write(*,*) "File POSCAR_cut has been written!"
write(*,*) "cut_unitcell exiting normally..."
write(*,*)

end program cut_unitcell

!
!    Taken from: https://fortranwiki.org/fortran/show/Matrix+inversion
!
subroutine matinv2(A,B)
  !! Performs a direct calculation of the inverse of a 2×2 matrix.
  real(kind=8), intent(in) :: A(2,2)   !! Matrix
  real(kind=8), intent(out) :: B(2,2)   !! Inverse matrix
  real(kind=8)             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
  if (abs(detinv) .lt. 1E-8) then
     write(*,*) "Error! The determinant of the coordinate matrix is almost zero!"
     write(*,*) " Please do not choose two parallel unit vectors!"
     stop
  end if        

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * A(2,2)
  B(2,1) = -detinv * A(2,1)
  B(1,2) = -detinv * A(1,2)
  B(2,2) = +detinv * A(1,1)
  return
end subroutine matinv2
