!
!    modify_xdatcar: Perform different tasks on a XDATCAR file,
!      such as translating or multiplying the unit cell, picking
!      certain structures or write an xyz trajectory file
!    Part of VASP4CLINT
!     Julien Steffen, 2023 (julien.steffen@fau.de)
!


program modify_xdatcar
implicit none 
integer::i,j,k,l,m
integer::readstat,openstat,counter,endl
integer::natoms,nelems,xdat_lines,nframes
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
real(kind=8),allocatable::a_vecs(:,:),b_vecs(:,:),c_vecs(:,:)
character(len=2),allocatable::el_names_read(:),el_names(:)
character(len=2),allocatable::at_names(:),at_names2(:)
integer,allocatable::el_nums(:)
real(kind=8)::act_num(3),factor
real(kind=8),allocatable::xyz(:,:,:),xyz2(:,:,:)
real(kind=8)::xlen,ylen,zlen
real(kind=8)::shift_vec(3),act_val,xyz_print(3)
integer::multiply_vec(3),pick_ind,pos_new,multiply_prod
integer::frame_first,frame_last,line_num
integer::read_freq
logical::npt,print_npt
logical::eval_stat(10)
logical::shift_cell,multiply_cell,pick_frame,print_xyz,print_last
character(len=120)::a120,cdum,arg,adum
character(len=220)::a220
character(len=50)::atest

!
!    only print the overview of utils4VASP scripts/programs and stop
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-overview") then
      write(*,*)
      write(*,*) "utils4VASP: Setup and Evaluation of DFT and MLP simulations with VASP"
      write(*,*) "The following scripts and programs are contained:"
      write(*,*) "Setup:"
      write(*,*) " - gen_incar.py: Generate INCAR file templates for several job types"
      write(*,*) " - analyze_poscar.py: Analyze POSCAR, generate KPOINTS and POTCAR"
      write(*,*) " - build_alloy.py: Build regular alloy structures of various shapes"
      write(*,*) " - modify_poscar.py: Modify POSCAR: multiply, transform, shift, insert etc."
      write(*,*) " - cut_unitcell: Generate surface slab unit cells of arbitrary shape"
      write(*,*) " - build_adsorb.py: Place adsorbates on surfaces, set translation, rotation"
      write(*,*) " - split_freq: Divide frequency calculations for large molecules"
      write(*,*) "Evaluation:"
      write(*,*) " - modify_xdatcar: Modify trajectory files: multiply, shift, pick etc."
      write(*,*) " - analyze_bulk: Analyze bulk MD trajectories for RDFs, diffusion etc."
      write(*,*) " - analyze_slab: Analyze slab MD trajectories for RDFs, density etc."
      write(*,*) " - check_geoopt.py: Monitor geometry optimizations with selective dynamics"
      write(*,*) " - manage_cls: Prepare, evaluate core level energy calculations"
      write(*,*) " - eval_bader: Evaluate and visualize Bader charge calculations"
      write(*,*) " - eval_stm: Generate STM images with different settings"
      write(*,*) " - partial_dos: Plot atom and orbital resolved density of states"
      write(*,*) " - manage_neb.py: Setup, monitor and restart NEB calculations"
      write(*,*) "ML-FF:"
      write(*,*) " - mlff_select: Heuristic selection of local reference configurations"
      write(*,*) " - eval_vasp_ml.py: Visualize results of VASP ML-FF on the fly learnings"
      write(*,*) " - vasp2trainset: Generate ML-FF training sets from VASP calculations"
      write(*,*) " - mlp_quality: Determine quality of MLPs for VASP validation set"
      write(*,*) "Management:"
      write(*,*) " - md_long.sh: Automated restart of MD calculations with slurm"
      write(*,*) " - opt_long.sh: Automated restart of geometry optimizations with slurm"
      write(*,*) " - ml_long.sh: Automated restart of VASP ML-FF on the fly learnings"
      write(*,*) " - mace_long.sh: Automated restart of MACE MD trajectories with ASE "
      write(*,*)
      stop
   end if
end do

!
!    Print general information and all possible keywords of the program    
!
write(*,*)
write(*,*) "PROGRAM modify_xdatcar: Modification of XDATCAR trajectory"
write(*,*) " files for arbitrary systems."
write(*,*) "The file XDATCAR must be present."
write(*,*) "The following calculations can be done, called by one of the"
write(*,*) " listed command line arguments:"
write(*,*) " -overview:  print an overview of all scripts and programs in "
write(*,*) "    utils4VASP"
write(*,*) " -shift=a,b,c : Shift each frame of the trajectory by a vector"
write(*,*) "    given in direct coordinates. Example: -shift=0.1,0.0,0.2"
write(*,*) " -multiply=a,b,c : Multiply the unit cell of each frame by some"
write(*,*) "    replications in each of the coordinate directions. Integers"
write(*,*) "    must be given as arguments. Example: -multiply=2,2,1"
write(*,*) " -read_freq=[number]: Only read in and process every nth frame"
write(*,*) "    Example: -read_freq=10 : every 10th frame will be read in and "
write(*,*) "    processed and written out. (DEFAULT: 1)"  
write(*,*) " -pick_frame=[number] : Pick one frame of the XDATCAR file and print"
write(*,*) "    it to a POSCAR file (POSCAR_pick). Example: -pick_frame=283"
write(*,*) " -print_xyz : Print each frame to a xyz trajectory: xdat_mod.xyz"
write(*,*) " -print_last=[number] : Only print the last [number] frames into "
write(*,*) "    the ordered file formate" 
write(*,*) " -print_npt : Print a NVT trajectory in the NPT format, e.g., for "
write(*,*) "    subsequent analysis with TRAVIS."
write(*,*) "One or several of these commands can be executed. The ordering"
write(*,*) " will be the same as the listing of the commands."
write(*,*) "If one job apart pick_frame is ordered, the trajectory will be "
write(*,*) " written to XDATCAR_mod or xdat_mod.xyz"
!
!    PART A: Read in the command line arguments !!!!!!!!!!!!!!!!
!
!    The coordinate shift vector
!
shift_cell=.false.
shift_vec=0.0d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-shift=") then
      read(arg(8:),*,iostat=readstat) shift_vec
      shift_cell=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -shift=..., something went wrong!"
         write(*,*)
      end if        
   end if
end do
!
!    The multiplication of unit cells
!
multiply_cell=.false.
multiply_vec=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-multiply=") then
      read(arg(11:),*,iostat=readstat) multiply_vec
      multiply_cell=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -multiply=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Read in the read in frequency (every n'th frame will be 
!          read in and processed)
!
read_freq=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-read_freq=") then
      read(arg(12:),*,iostat=readstat) read_freq
      pick_frame=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -read_freq=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (read_freq .lt. 1) then
   write(*,*) "The read-in frequency must be at least 1!"
   stop
end if        
!
!     Pick a certain frame from the trajectory
!
pick_frame=.false.
pick_ind=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-pick_frame=") then
      read(arg(13:),*,iostat=readstat) pick_ind
      pick_frame=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -pick_frame=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
!
!     Write the trajectory to a xyz file
!
print_xyz=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-print_xyz") then
      print_xyz=.true.
   end if
end do

!
!     Write the trajectory in the NPT format 
!
print_npt=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-print_npt") then
      print_npt=.true.
   end if
end do

!
!     Only print the last nframes to the trajectory
!
print_last=.false.
frame_last=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-print_last=") then
      read(arg(13:),*,iostat=readstat) frame_last
      print_last=.true.
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -print_last=..., something went wrong!"
         write(*,*)
      end if
   end if
end do


if ((.not. shift_cell) .and. (.not. multiply_cell) .and. (.not. pick_frame) .and. &
           &  (.not. print_xyz) .and. (.not. print_last) .and. (read_freq .eq. 1) .and. &
           &  (.not. print_npt)) then
   write(*,*)
   write(*,*) "Please give at least one of the possible commands!"
   write(*,*)
   stop
end if

!
!    PART B: Read in the XDATCAR file !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    First, determine the number of lines in the XDATCAR file
!
write(*,*)
write(*,*) "Determine number of atoms and frames in XDATCAR..."
call system("wc -l XDATCAR > xdat_length")
open(unit=45,file="xdat_length",status="old")
read(45,*) xdat_lines
close(45)
write(*,*) "completed!"

open(unit=14,file="XDATCAR",status="old",iostat=openstat)
if (openstat .ne. 0) then
   stop "ERROR! The file 'XDATCAR' could not been found!"
end if
read(14,*)
read(14,*) factor   ! the global geometry conversion factor
!
!    Read in the cell dimensions
!
read(14,*) a_vec(1),a_vec(2),a_vec(3)
read(14,*) b_vec(1),b_vec(2),b_vec(3)
read(14,*) c_vec(1),c_vec(2),c_vec(3)
!
!    Read in the elements
!
allocate(el_names_read(10))
el_names_read="XX"
read(14,'(a)') cdum
read(cdum,*,iostat=readstat) el_names_read
nelems=0
do i=1,10
   if (el_names_read(i) .eq. "XX") exit
   nelems=nelems+1
end do
allocate(el_names(nelems),el_nums(nelems))

do i=1,nelems
   el_names(i)=el_names_read(i)
end do
read(14,*) el_nums
!
!    Define the element symbols for all atoms in the system
!
natoms = sum(el_nums)
allocate(at_names(natoms))

!
!    Check if the trajectory is an NVT or NPT trajectory
!
read(14,*)
do i=1,natoms
   read(14,*)
end do
read(14,*) adum
!
!    How does the first line of the second frame begin? If it begins 
!    with Direct or direct, the trajectory is a NVT trajectory
!
if (trim(adum) .eq. "Direct" .or. trim(adum) .eq. "direct") then
   npt=.false.
   write(*,*) "A NVT trajectory (constant unit cell) has been detected."
else
   npt=.true.
   write(*,*) "A NpT trajectory (varying unit cell) has been detected."
end if

counter = 1
do i=1,nelems
   do j=1,el_nums(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do
if (npt) then
   nframes = (xdat_lines - 7)/(natoms+8)
else
   nframes = (xdat_lines - 7)/(natoms+1)
end if
!write(*,*) xdat_lines,natoms,nframes
!
!    If the read frequency is larger than 1, divide the number of 
!       frames by it!
!
if (read_freq .gt. 1) then
   nframes=int(nframes/read_freq)
end if        
!
!    From the print-last command, determine the first frame to be written
!
if (print_last) then
   frame_first=nframes-frame_last+1
else 
   frame_first=1
end if

!
!    Allocate global coordinate arrays
!
allocate(xyz(3,natoms,nframes))
if (npt) then
   allocate(a_vecs(3,nframes))
   allocate(b_vecs(3,nframes))
   allocate(c_vecs(3,nframes))
end if

!
!    Read in the coordinates of the frames and correct for image flags
!

close(14)
open(unit=14,file="XDATCAR",status="old")

eval_stat = .false.
write(*,*)
write(*,*) "Read in the trajectory from XDATCAR..."
line_num=1
do i=1,nframes
!
!    Every 10% of the read in, give a status update
!
   do j=1,10
      if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
         if (.not. eval_stat(j)) then
            write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
            eval_stat(j) = .true.
         end if
      end if
   end do
   read(14,*)
!
!     A NPT trajectory
!
   if (npt) then
      read(14,*)
      read(14,*) a_vecs(:,i)
      read(14,*) b_vecs(:,i)
      read(14,*) c_vecs(:,i) 
      read(14,*)
      read(14,*)
      read(14,*)
!
!     A NVT trajectory
!
   else
      if (i .eq. 1) then 
         read(14,*)
         read(14,*) 
         read(14,*) 
         read(14,*) 
         read(14,*)
         read(14,*)
         read(14,*)
      end if
   end if
   
!
!     A NVT trajectory
!
   do j=1,natoms
      read(14,'(a)') a120
      read(a120,*,iostat=readstat) act_num(:)
!
!   In long trajectories, problems might arise with double-digit negative numbers!
!   (no space between them)   Insert a space again!
!
      if (readstat .ne. 0) then
         a220 = ""
         endl = 1
         do k=1,119
            atest=a120(k:k+1)
            if (index(atest,"-") .ne. 0) then
               write(a220(1:endl+4),'(a,a)') a220(1:endl),"   -"
               endl = endl + 4
            else
               write(a220(1:endl+1),'(a,a)',iostat=readstat) a220(1:endl),atest
               endl = endl + 1
            end if
         end do
    !  write(*,*) "Problem in XDATCAR file!"
         read(a220,*,iostat=readstat) act_num(:)
      end if

      do k=1,3
         do
            if (act_num(k) > 1d0) then
               act_num(k) = act_num(k) - 1d0
            else if (act_num(k) < 0.d0) then
               act_num(k) = act_num(k) + 1d0
            else 
               exit
            end if
         end do
         xyz(k,j,i)=act_num(k)
      end do
   end do
!
!     If read_freq > 1, skip the next read_freq-1 frames
!
   if (read_freq .gt. 1) then
      do j=1,read_freq-1
         if (npt) then
            do k=1,natoms+8
               read(14,*)
            end do
         else
            do k=1,natoms+1
               read(14,*)
            end do
         end if
      end do
   end if
end do
write(*,*) " completed!"
close(14)

write(*,*)
write(*,*) "---------- SETTINGS ---------------------------"
write(*,*) "Number of atoms in the system:",natoms
write(*,*) "Number of frames in the trajectory:",nframes*read_freq
if (read_freq .gt. 1) then
   write(*,'(a,i5,a)') " Every ",read_freq,"th frame will be evaluated."
end if        
write(*,*) "-----------------------------------------------"
write(*,*)
!
!     C: Shift the unit cell along a vector
!
if (shift_cell) then
   eval_stat = .false.
   write(*,'(a,f9.4,a,f9.4,a,f9.4)') " Shift all frames of the XDATCAR along the vector ", &
             & shift_vec(1),",",shift_vec(2),",",shift_vec(3)
   do i=1,nframes
      do j=1,10
         if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do
      do j=1,natoms
         do k=1,3
            act_val=xyz(k,j,i)+shift_vec(k)
            do 
               if (act_val > 1.0d0) then
                  act_val=act_val-1.0d0
               else if (act_val < 0.0d0) then
                  act_val=act_val+1.d0
               else 
                  exit
               end if
            end do
            xyz(k,j,i)=act_val
         end do
      end do
   end do
   write(*,*) "completed!"
   write(*,*)
end if
!
!    D: Multiply the unit cell according to given integer tuple
!
if (multiply_cell) then
   eval_stat = .false.
   write(*,'(a,i4,a,i4,a,i4,a)') " Multipy the unit cell of all frames: ",multiply_vec(1), &
               & " times a, ",multiply_vec(2)," times b, ",multiply_vec(3)," times c"
   multiply_prod=multiply_vec(1)*multiply_vec(2)*multiply_vec(3)
   if (npt) then
      do i=1,nframes
         a_vecs(:,i)=a_vecs(:,i)*multiply_vec(1)
         b_vecs(:,i)=b_vecs(:,i)*multiply_vec(2)
         c_vecs(:,i)=c_vecs(:,i)*multiply_vec(3)
      end do
   else
      a_vec=a_vec*multiply_vec(1)
      b_vec=b_vec*multiply_vec(2)
      c_vec=c_vec*multiply_vec(3)
   end if
   el_nums=el_nums*multiply_prod
   allocate(xyz2(3,natoms*multiply_prod,nframes),at_names2(natoms*multiply_prod))
   do i=1,nframes
      do j=1,10
         if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do
      pos_new=0
      do j=1,natoms
         do k=1,multiply_vec(1)
            do l=1,multiply_vec(2)
               do m=1,multiply_vec(3)
                  pos_new=pos_new+1
                  xyz2(1,pos_new,i)=(xyz(1,j,i)+1.0*real(k))/real(multiply_vec(1))
                  xyz2(2,pos_new,i)=(xyz(2,j,i)+1.0*real(l))/real(multiply_vec(2))
                  xyz2(3,pos_new,i)=(xyz(3,j,i)+1.0*real(m))/real(multiply_vec(3))
                  at_names2(pos_new)=at_names(j)
               end do
            end do
         end do
      end do
   end do
   natoms=natoms*multiply_prod

   write(*,*) "completed!"
   write(*,*)
else 
!
!    If no multiplication shall be done, simply copy contents of xyz to xyz2
!
   allocate(xyz2(3,natoms,nframes))
   allocate(at_names2(natoms))
   xyz2=xyz
   at_names2=at_names
end if
!
!    E: Pick a frame of the trajectory and print it as POSCAR or xyz file
!
if (pick_frame) then
   if (print_xyz) then
      write(*,'(a,i10,a)') " Frame ",pick_ind," will be written to file pos_pick.xyz"
      open(unit=33,file="pos_pick.xyz",status="replace")
      write(33,*) natoms
      write(33,'(a,i10,a)') " Frame ",pick_ind," picked from XDATCAR by modify_xdatcar"
      do j=1,natoms
         if (npt) then
            xyz_print(1)=(xyz2(1,j,pick_ind)*a_vecs(1,pick_ind)+xyz2(2,j,pick_ind)* &
                         & b_vecs(1,pick_ind)+xyz2(3,j,pick_ind)*c_vecs(1,pick_ind))*factor
            xyz_print(2)=(xyz2(1,j,pick_ind)*a_vecs(2,pick_ind)+xyz2(2,j,pick_ind)* &
                         & b_vecs(2,pick_ind)+xyz2(3,j,pick_ind)*c_vecs(2,pick_ind))*factor
            xyz_print(3)=(xyz2(1,j,pick_ind)*a_vecs(3,pick_ind)+xyz2(2,j,pick_ind)* &
                         & b_vecs(3,pick_ind)+xyz2(3,j,pick_ind)*c_vecs(3,pick_ind))*factor
         else
            xyz_print(1)=(xyz2(1,j,pick_ind)*a_vec(1)+xyz2(2,j,pick_ind)*b_vec(1)+ &
                         & xyz2(3,j,pick_ind)*c_vec(1))*factor
            xyz_print(2)=(xyz2(1,j,pick_ind)*a_vec(2)+xyz2(2,j,pick_ind)*b_vec(2)+ &
                         & xyz2(3,j,pick_ind)*c_vec(2))*factor
            xyz_print(3)=(xyz2(1,j,pick_ind)*a_vec(3)+xyz2(2,j,pick_ind)*b_vec(3)+ &
                      & xyz2(3,j,pick_ind)*c_vec(3))*factor
         end if
         write(34,*) at_names2(j),xyz_print(:)
      end do
      close(33)  
      write(*,*) "completed!"
      write(*,*)
   else 
      write(*,'(a,i10,a)') " Frame ",pick_ind," will be written to file POSCAR_pick"
      open(unit=33,file="POSCAR_pick",status="replace")
      write(33,'(a,i10,a)') " Frame ",pick_ind," picked from XDATCAR by modify_xdatcar"
      write(33,*) factor
      if (npt) then
         write(33,'(3f15.6)') a_vecs(:,pick_ind)
         write(33,'(3f15.6)') b_vecs(:,pick_ind)
         write(33,'(3f15.6)') c_vecs(:,pick_ind)
      else
         write(33,'(3f15.6)') a_vec
         write(33,'(3f15.6)') b_vec
         write(33,'(3f15.6)') c_vec
      end if
      do j=1,nelems
         write(33,'(a,a)',advance="no") el_names(j),"  "
      end do
      write(33,*)
      do j=1,nelems
         write(33,'(i6,a)',advance="no") el_nums(j)," "
      end do
      write(33,*)

      write(33,*) "Direct configuration=  ",pick_ind
      do j=1,natoms
         write(33,'(3f15.8)') xyz2(:,j,pick_ind)
      end do
      close(33)
      write(*,*) "completed!"
      write(*,*)
   end if
end if 
!
!    D: Write structure to XYZ trajectory file
!      During this, convert coordinates to cartesian!
!      If not ordered, write a modified XDATCAR file (in direct coordinates) instead
!      Print nothing if the pick_frame option is chosen
!
eval_stat = .false.
if (print_xyz) then
   write(*,*) "Write trajectory in xyz format to file xdat_mod.xyz"
   if (frame_last .ne. 0) then
      write(*,'(a,i10,a)') " Only the last ",frame_last," frames will be written."
   end if
   open(unit=34,file="xdat_mod.xyz",status="replace")
   do i=frame_first,nframes
      do j=1,10
         if (real(i)/real(nframes-frame_first) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do
      write(34,*) natoms
      write(34,*) "Trajectory converted from XDATCAR file via modify_xdatcar"
      do j=1,natoms
         if (npt) then
            xyz_print(1)=(xyz2(1,j,i)*a_vecs(1,i)+xyz2(2,j,i)*b_vecs(1,i)+ &
                         & xyz2(3,j,i)*c_vecs(1,i))*factor
            xyz_print(2)=(xyz2(1,j,i)*a_vecs(2,i)+xyz2(2,j,i)*b_vecs(2,i)+ &
                         & xyz2(3,j,i)*c_vecs(2,i))*factor
            xyz_print(3)=(xyz2(1,j,i)*a_vecs(3,i)+xyz2(2,j,i)*b_vecs(3,i)+ & 
                         & xyz2(3,j,i)*c_vecs(3,i))*factor
         else 
            xyz_print(1)=(xyz2(1,j,i)*a_vec(1)+xyz2(2,j,i)*b_vec(1)+ &
                         & xyz2(3,j,i)*c_vec(1))*factor
            xyz_print(2)=(xyz2(1,j,i)*a_vec(2)+xyz2(2,j,i)*b_vec(2)+ &
                         & xyz2(3,j,i)*c_vec(2))*factor
            xyz_print(3)=(xyz2(1,j,i)*a_vec(3)+xyz2(2,j,i)*b_vec(3)+ &
                         & xyz2(3,j,i)*c_vec(3))*factor
         end if
         write(34,*) at_names2(j),xyz_print(:)
      end do
   end do
   close(34)
   write(*,*) "completed!"
else 
   if (print_npt) then
      if (.not. npt) then
         npt=.true.
         allocate(a_vecs(3,nframes))
         allocate(b_vecs(3,nframes))
         allocate(c_vecs(3,nframes))
         do i=1,nframes
            a_vecs(:,i)=a_vec
            b_vecs(:,i)=b_vec
            c_vecs(:,i)=c_vec
         end do
      end if
      write(*,*) "The -print_npt command is given, the trajectory is written as NpT"
   end if
   if (shift_cell .or. multiply_cell .or. print_last .or. (read_freq .gt. 1) .or. &
                   & print_npt) then
      write(*,*) "Write trajectory in VASP format to file XDATCAR_mod"
      if (frame_last .ne. 0) then
         write(*,'(a,i10,a)') " Only the last ",frame_last," frames will be written."
      end if
      open(unit=34,file="XDATCAR_mod",status="replace")
      do i=frame_first,nframes
         do j=1,10
            if (real(i)/real(nframes-frame_first) .gt. real(j)*0.1d0) then
               if (.not. eval_stat(j)) then
                  write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
                  eval_stat(j) = .true.
               end if
            end if
         end do
         if (npt) then
            write(34,*) "NpT Trajectory written by modify_xdatcar"
            write(34,*) factor
            write(34,'(3f15.6)') a_vecs(:,i)
            write(34,'(3f15.6)') b_vecs(:,i)
            write(34,'(3f15.6)') c_vecs(:,i)
            do j=1,nelems
               write(34,'(a,a)',advance="no") el_names(j),"  "
            end do
            write(34,*)
            do j=1,nelems
               write(34,'(i6,a)',advance="no") el_nums(j)," "
            end do
            write(34,*)
         else if (i .eq. frame_first) then
            write(34,*) "NVT Trajectory written by modify_xdatcar"
            write(34,*) factor
            write(34,'(3f15.6)') a_vec
            write(34,'(3f15.6)') b_vec
            write(34,'(3f15.6)') c_vec
            do j=1,nelems
               write(34,'(a,a)',advance="no") el_names(j),"  "
            end do
            write(34,*)
            do j=1,nelems
               write(34,'(i6,a)',advance="no") el_nums(j)," "
            end do
            write(34,*)
         end if
 
         write(34,*) "Direct configuration=  ",i
         do j=1,natoms
            write(34,'(3f15.8)') xyz2(:,j,i) 
         end do
      end do
      close(34)
      write(*,*) "completed!"
   end if
end if

write(*,*)
write(*,*) "modify_xdatcar has finished all tasks, goodbye!"
write(*,*)

end program modify_xdatcar        
