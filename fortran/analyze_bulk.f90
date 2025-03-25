!
!    analyze_bulk:  reads in a XDATCAR file of a 
!      three-dimensional bulk system and calculates 
!      different dynamical observables like diffusion 
!      coefficients or vibrational density of states
!
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!
module fftw_mod
   use,intrinsic :: iso_c_binding
   include 'fftw3.f03'

   integer::Np
   integer, parameter :: Nmax = 1024
   type(C_PTR) :: plan
   real(kind=8)::factor
end module fftw_mod


program analyze_bulk
use fftw_mod
implicit none
integer::i,j,k,l,m
integer::readstat
logical::calc_msd,calc_vacf,read_dt,npt_traj,print_new
logical::vib_dens
logical::dens_cube
real(kind=8)::a_read(3),b_read(3),c_read(3)
real(kind=8)::a_len,b_len,c_len
real(kind=8),allocatable::a_lens(:),b_lens(:),c_lens(:)
integer::nelems,natoms,nframes,xdat_lines
integer,allocatable::el_nums(:)
integer::frames_skip,nframes_use
integer::gridx,gridy,gridz
integer::stepx,stepy,stepz
real(kind=8)::distance
integer::i1,j1,k1,i_new,j_new,k_new
integer::gridx_act,gridy_act,gridz_act
integer::inc,elnumber
real(kind=8)::a2bohr
character(len=30)::cube_name
real(kind=8),allocatable::dens_3d(:,:,:,:)
real(kind=8),allocatable::xyz(:,:,:)
real(kind=8)::diff_vec(3)
real(kind=8)::time_step
real(kind=8)::msd_act,volume
real(kind=8),allocatable::pos_diff(:)
real(kind=8),allocatable::z_axis(:)
real(kind=8),allocatable::msd_func(:,:),vacf_func(:),vacf_pad(:)
real(kind=8),allocatable::fourier_out(:)
real(kind=8),allocatable::vel_first(:),vel_act(:)
real(kind=8),allocatable::times(:),diff(:)
real(kind=8),allocatable::vector1(:),vector2(:),vector3(:)
real(kind=8),allocatable::vector4(:),vector5(:)
real(kind=8)::vacf_int,norm_factor
character(len=2),allocatable::el_names(:)
real(kind=8),allocatable::rdf_plot(:,:,:)
real(kind=8),allocatable::rdf_sum(:,:,:)
real(kind=8)::pos1(3),pos2(3)
integer::avg_lo,avg_hi
real(kind=8)::avg_diff
real(kind=8)::corrt
character(len=32)::arg
character(len=80)::line,all_els
character(len=2)::el1,el2,el3,el4,el5
logical::eval_stat(10)
real(kind=8)::box_volume
real(kind=8)::dist,pi,rdf_cutoff
real(kind=8)::xlen,ylen,zlen
real(kind=8)::int_act
real(kind=8),allocatable::vacf_plot(:)
integer::frame_first,intervals,frames_part
integer::i_max
integer::vib_inter
real(kind=8)::vib_width
!  RDF calculation
integer::ig,ngr,npart1,npart2
real(kind=8)::nid,r_act,rho,vb
real(kind=8)::slice_step,rdf_binsize
real(kind=8),allocatable::z_dens(:,:),z_dens_tot(:)
real(kind=8),allocatable::time_list(:)
integer::all_tasks,task_act,rdf_bins
logical::calc_rdf,read_time

pi=3.141592653589793238
a2bohr=1.8897259886d0
write(*,*) "PROGRAM analyze_bulk: evaluation of MD trajectories from three-dimensional"
write(*,*) " bulk systems."
write(*,*) "Only a XDATCAR file in VASP format needs to be present."
write(*,*) "Usage: eval_vasp_md -command1 -command2 ..."
write(*,*)
write(*,*) "List of commands:"
write(*,*) "-msd:  the mean square displacement (and self-diffusion coefficient) of up to"
write(*,*) "  three different elements is calculated."
write(*,*) "-vib_dens:  the vibrational density of states (IR spectrum without intensities)"
write(*,*) "        is calculated from the VACF"
write(*,*) "-dt=[time step in fs]: The time step used for MD simulation, in fs."
write(*,*) "-corrt=[time in fs]: Length of the VACF correlation interval to be calculated."
write(*,*) "       (default: 1000 fs)"
write(*,*) "-vib_inter=[number]: The interval in which the IR spectrum is plotted (upper limit)"
write(*,*) "       (default: 6000 cm^-1)"
write(*,*) "-vib_width=[value]: The Gaussian broadening of the VACF IR spectrum. "
write(*,*) "       (default: 0.001, larger value gives less broadening)"
write(*,*) "-readtime: The time of each XDATCAR (in s) will be read in from file 'times.dat'"
write(*,*) "    useful for the evaluation of kinetic Monte Carlo simulations."
write(*,*) "-skip=[steps to skip]: If the first N steps shall be skipped (equilibration.)"
write(*,*) "-npt: If a NpT trajectory (variable volume) shall be analyzed."
write(*,*) "-rdf: Radial distribution functions will be calculated."
write(*,*) "-dens_cube: Write Gaussian cube file with spatially resolved atomic densities."
write(*,*) "-print: A new XDATCAR file (XDATCAR_new) with image flags will be written,"
write(*,*) "        containing only frames after the skipped frames."
write(*,*)
!
!    Use Command line arguments for specification of analysis job
!
calc_msd = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-msd") then
      calc_msd = .true.
      write(*,*) "The mean square displacement (MSD) will be calculated!"
   end if
end do

calc_vacf = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-vib_dens") then
      calc_vacf = .true.
      write(*,*) "The vibrational density of states will be calculated!"
   end if
end do

read_time = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-readtime") then
      read_time = .true.
      write(*,*) "The time of each step will be read in from file 'times.dat'!"
   end if
end do

dens_cube = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-dens_cube") then
      dens_cube = .true.
      write(*,*) "The averaged atomic densities are written to 'densities.cube'!"
   end if
end do

!
!     Number of grid points per dimension for 3D element densities
!
gridx=100
gridy=100
gridz=100




rdf_bins=500
time_step = 0.d0
rdf_cutoff=10d0
rdf_binsize=rdf_cutoff/real(rdf_bins)
if (.not. read_time) then
   read_dt = .false.
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg(1:4))  .eq. "-dt=") then
         read_dt = .true.
         read(arg(5:32),*) time_step
         write(*,*) "The time step shall be:",time_step," fs."
      end if
   end do
   if (.not. read_dt) then
      stop "Please set a time step with the -dt=... flag!"
   end if       
else
   read_dt=.true.
end if

print_new = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-print") then
      print_new = .true.
      write(*,*) "The evaluated part of the trajectory will be written to XDATCAR_new"
   end if
end do


frames_skip = 0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-skip=") then
      read(arg(7:32),*) frames_skip
      write(*,*) "The first ",frames_skip," frames shall be skipped."
   end if
end do
frame_first=1+frames_skip

npt_traj = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-npt") then
      npt_traj = .true.
      write(*,*) "A NpT trajectory will be analyzed (variable volume)."
   end if
end do

calc_rdf=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-rdf") then
      calc_rdf = .true.
      write(*,*) "Radial distribution functions will be calculated."
   end if
end do

corrt = 1000
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-corrt=") then
      read(arg(8:32),*) corrt
      write(*,*) "The VACF correlation time is ",corrt," fs."
   end if
end do

vib_inter = 6000
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-vib_inter=") then
      read(arg(12:32),*) vib_inter
      write(*,*) "The IR plot interval is up to ",vib_inter, " cm^-1."
   end if
end do

vib_width = 0.001
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:11))  .eq. "-vib_width=") then
      read(arg(12:32),*) vib_width
      write(*,*) "The IR Gaussian line broadening prefactor is ",vib_width
   end if
end do






!
!    First, determine the number of lines in the XDATCAR file
!
call system("wc -l XDATCAR > xdat_length")
open(unit=45,file="xdat_length",status="old")
read(45,*) xdat_lines
close(45)


open(unit=15,file="XDATCAR",status="old")
read(15,*)
read(15,*)
read(15,*) a_read(:)
read(15,*) b_read(:)
read(15,*) c_read(:)
read(15,'(a)') line 
all_els=line
!
!    Determine number of elements 
!


read(line,*,iostat=readstat) el1,el2,el3,el4,el5
if (readstat .ne. 0) then
   read(line,*,iostat=readstat) el1,el2,el3,el4
   if (readstat .ne. 0) then
      read(line,*,iostat=readstat) el1,el2,el3
      if (readstat .ne. 0) then
         read(line,*,iostat=readstat) el1,el2
         if (readstat .ne. 0) then
            read(line,*,iostat=readstat) el1
            nelems=1
         else 
            nelems=2
         end if   
      else    
         nelems=3
      end if        
   else      
       nelems=4
   end if         
else 
   nelems=5
end if   

allocate(el_nums(nelems))
allocate(el_names(nelems))
read(15,*) el_nums

natoms=sum(el_nums)

write(*,*) "System setup:"
write(*,*) "Number of atoms:",natoms
write(*,*) "Number of elements:",nelems
write(*,*) "List of elements:"
write(*,*) "  1. ",el1,": ",el_nums(1)," atoms"
el_names(1)=el1
if (nelems .ge. 2) then
   el_names(2)=el2   
   write(*,*) "  2. ",el2,": ",el_nums(2)," atoms"
end if
if (nelems .ge. 3) then
   el_names(3)=el3     
   write(*,*) "  3. ",el3,": ",el_nums(3)," atoms"
end if
if (nelems .ge. 4) then
   el_names(4)=el4
   write(*,*) "  4. ",el4,": ",el_nums(4)," atoms"
end if
if (nelems .ge. 5) then
   el_names(5)=el5
   write(*,*) "  5. ",el5,": ",el_nums(5)," atoms"
end if


!
!    For NVT trajectories: Each frame has only one header line 
!
if (.not. npt_traj) then
   a_len = a_read(1)
   b_len = b_read(2)
   c_len = c_read(3)

   nframes=int((xdat_lines-7)/(natoms+1))

   write(*,*) "Number of frames:",nframes
   

   allocate(xyz(3,natoms,nframes))
   do i=1,nframes
           
      read(15,*)
      do j=1,natoms
         read(15,*,iostat=readstat) xyz(:,j,i)
         if (readstat .ne. 0) then
            stop "Error in reading in XDATCAR? Add the flag -npt?"
         end if        
      end do
   end do
   close(15)

else 
!
!    For Npt trajectories: Read in volume for each frame!
!

   nframes=int((xdat_lines)/(natoms+8))

   write(*,*) "Number of frames:",nframes

   allocate(a_lens(nframes),b_lens(nframes),c_lens(nframes))
   a_lens(1)=a_read(1)
   b_lens(1)=b_read(2)
   c_lens(1)=c_read(3)

   allocate(xyz(3,natoms,nframes))

   do i=1,nframes 
      if (i .gt. 1) then   
         read(15,*)
         read(15,*)
         read(15,*) a_read(:)
         read(15,*) b_read(:)
         read(15,*) c_read(:)
         read(15,*)
         read(15,*)
         read(15,*)
         a_lens(i)=a_read(1)
         b_lens(i)=b_read(2)
         c_lens(i)=c_read(3)
      else 
         read(15,*)
      end if
      do j=1,natoms
         read(15,*) xyz(:,j,i)
      end do
   end do
   close(15)
end if  
!
!    If the times shall be read in from file: read in one 
!        time per frame!
!
if (read_time) then
   open(unit=56,file="times.dat",status="old",iostat=readstat)
   if (readstat .ne. 0) then
      write(*,*) "The file times.dat is not there!"
      stop
   end if
   allocate(time_list(nframes))
   time_list=0.d0
   do i=1,nframes
      read(56,*,iostat=readstat) time_list(i)
      if (readstat .ne. 0) then
         write(*,*) "The file times.dat has a wrong formate on line ",i
         stop
      end if
   end do
end if
 
!
!    Correct for box images 
!
!goto 33
do i=1,nframes-1
   do j=1,natoms
      diff_vec(:) = xyz(:,j,i+1) - xyz(:,j,i)
!
!     Correct the x component
!
      do while (abs(diff_vec(1)) .gt. 0.5d0)
         diff_vec(1)=diff_vec(1)-sign(1.0d0,diff_vec(1))
!         write(*,*) "correct x"
      end do

!
!     Correct the y component
!
      do while (abs(diff_vec(2)) .gt. 0.5d0)
         diff_vec(2)=diff_vec(2)-sign(1.0d0,diff_vec(2))
!         write(*,*) "correct y"
      end do

!
!     Correct the z component
!
      do while (abs(diff_vec(3)) .gt. 0.5d0)
         diff_vec(3)=diff_vec(3)-sign(1.0d0,diff_vec(3))
!         write(*,*) "correct z"
      end do

      xyz(:,j,i+1)=xyz(:,j,i)+diff_vec
   end do
end do
!33 continue

!
!    If desired, print the corrected part of the evaluated trajectory
!

if (print_new) then
   if (.not. npt_traj) then
      open(unit=38,file="XDATCAR_new",status="replace")
      write(38,*) "MD trajectory written by eval_vasp_md"
      write(38,*) 1.0
      write(38,*) a_len,0.0,0.0
      write(38,*) 0.0,b_len,0.0
      write(38,*) 0.0,0.0,c_len
      write(38,*) all_els
      write(38,*) el_nums
      do i=frames_skip+1,nframes
         write(38,*) "Direct configuration= ",i-frames_skip
         do j=1,natoms
            write(38,*) xyz(:,j,i)
         end do
      end do

      close(38)

   else
      open(unit=38,file="XDATCAR_new",status="replace")
      do i=frames_skip+1,nframes
         write(38,*) "MD trajectory written by eval_vasp_md"
         write(38,*) 1.0
         write(38,*) a_lens(i),0.0,0.0
         write(38,*) 0.0,b_lens(i),0.0
         write(38,*) 0.0,0.0,c_lens(i)
         write(38,*) all_els
         write(38,*) el_nums
         write(38,*) "Direct configuration= ",i-frames_skip
         do j=1,natoms
            write(38,*) xyz(:,j,i)
         end do
      end do
      close(38)
   end if        
end if        

!
!    Convert to usual cartesian coordinates 
!
if (.not. npt_traj) then
   do i=1,nframes
      do j=1,natoms
         xyz(1,j,i) = xyz(1,j,i)*a_len
         xyz(2,j,i) = xyz(2,j,i)*b_len
         xyz(3,j,i) = xyz(3,j,i)*c_len     
      end do      
   end do
else 
   do i=1,nframes
      do j=1,natoms
         xyz(1,j,i) = xyz(1,j,i)*a_lens(i)
         xyz(2,j,i) = xyz(2,j,i)*b_lens(i)
         xyz(3,j,i) = xyz(3,j,i)*c_lens(i)
      end do
   end do
end if
!
!    For 3D element densities, determine number of steps around current position
!    One Angstrom around
!
if (.not. npt_traj) then
   stepx=int(gridx/a_len)
   stepy=int(gridy/b_len)
   stepz=int(gridz/c_len)
else 
   stepx=int(gridx/a_lens(nframes))
   stepy=int(gridy/b_lens(nframes))
   stepz=int(gridz/c_lens(nframes))
end if        

!
!    For NPT trajectories, calculate the average volume
!
if (npt_traj) then
   volume=0
   do i=1,nframes
      volume=volume+a_lens(i)*b_lens(i)*c_lens(i)
   end do
   volume=volume/real(nframes)
   write(*,*) "Averaged volume of the system written to vol_avg.dat"
   open(unit=36,file="vol_avg.dat",status="replace")
   write(36,*) volume
   close(36)
else
   volume=a_len*b_len*c_len
end if        
!
!    Calculate the mean square displacement (MSD)
!

allocate(pos_diff(natoms*3))
allocate(z_axis(natoms*3))
allocate(msd_func(nframes-frames_skip,nelems))
allocate(times(nframes-frames_skip))
allocate(diff(nframes-frames_skip))

if (calc_msd) then

   if (nelems .eq. 2) then
      allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3))
   end if        
   if (nelems .eq. 3) then
      allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3),vector3(el_nums(3)*3))
   end if       
   if (nelems .eq. 4) then
      allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3),vector3(el_nums(3)*3),&
               & vector4(el_nums(4)*3))      
   end if 
   if (nelems .eq. 5) then
      allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3),vector3(el_nums(3)*3),&
               & vector4(el_nums(4)*3),vector5(el_nums(5)*3))
   end if
   do i=1,nframes-frames_skip
      if (read_time) then
         times(i) = time_list(i)/1E-15
      else
         times(i)=(i-1)*time_step
      end if
      do j=1,natoms
         do k=1,3
            pos_diff((j-1)*3+k)=xyz(k,j,i+frames_skip)-xyz(k,j,1+frames_skip)
         end do
      end do
      if (nelems .eq. 1) then
         msd_func(i,1)=dot_product(pos_diff,pos_diff)/natoms
      end if   
      if (nelems .eq. 2) then
         vector1 = pos_diff(1:el_nums(1)*3)
         msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
         vector2 = pos_diff(el_nums(1)*3+1:natoms*3)
         msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
      end if        
      if (nelems .eq. 3) then
         vector1 = pos_diff(1:el_nums(1)*3)
         msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
         vector2 = pos_diff(el_nums(1)*3+1:(el_nums(1)+el_nums(2))*3)
         msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
         vector3 = pos_diff((el_nums(1)+el_nums(2))*3+1:natoms*3)
         msd_func(i,3)=dot_product(vector3,vector3)/el_nums(3)
      end if       
      if (nelems .eq. 4) then
         vector1 = pos_diff(1:el_nums(1)*3)
         msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
         vector2 = pos_diff(el_nums(1)*3+1:(el_nums(1)+el_nums(2))*3)
         msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
         vector3 = pos_diff(sum(el_nums(1:2))*3+1:sum(el_nums(1:3))*3)
         msd_func(i,3)=dot_product(vector3,vector3)/el_nums(3)
         vector4 = pos_diff(sum(el_nums(1:3))*3+1:natoms*3)
         msd_func(i,4)=dot_product(vector4,vector4)/el_nums(4)
      end if
      if (nelems .eq. 5) then
         vector1 = pos_diff(1:el_nums(1)*3)
         msd_func(i,1)=dot_product(vector1,vector1)/el_nums(1)
         vector2 = pos_diff(el_nums(1)*3+1:(el_nums(1)+el_nums(2))*3)
         msd_func(i,2)=dot_product(vector2,vector2)/el_nums(2)
         vector3 = pos_diff(sum(el_nums(1:2))*3+1:sum(el_nums(1:3))*3)
         msd_func(i,3)=dot_product(vector3,vector3)/el_nums(3)
         vector4 = pos_diff(sum(el_nums(1:3))*3+1:sum(el_nums(1:4))*3)
         msd_func(i,4)=dot_product(vector4,vector4)/el_nums(4)
         vector5 = pos_diff(sum(el_nums(1:4))*3+1:natoms*3)
         msd_func(i,5)=dot_product(vector5,vector5)/el_nums(5)
      end if
 
   end do
!   write(*,*) natoms,el_nums(1),el_nums(2) 

   open(unit=17,file="msd_plot.dat",status="replace")
   write(17,*) "# time (fs),   MSD (A^2)"
   do i=1,nframes-frames_skip
      if (nelems .eq. 1) then
         write(17,*) times(i),msd_func(i,1)
      end if
      if (nelems .eq. 2) then
         write(17,*) times(i),msd_func(i,1),msd_func(i,2)
      end if
      if (nelems .eq. 3) then
         write(17,*) times(i),msd_func(i,1),msd_func(i,2),msd_func(i,3)
      end if        
      if (nelems .eq. 4) then
         write(17,*) times(i),msd_func(i,1),msd_func(i,2),msd_func(i,3), &
                           & msd_func(i,4)
      end if
      if (nelems .eq. 5) then
         write(17,*) times(i),msd_func(i,1),msd_func(i,2),msd_func(i,3), &
                           & msd_func(i,4),msd_func(i,5)
      end if
      
   end do

   close(17)
!
!     Calculate the diffusion coefficient by averaging the MSD between 
!        0.3 and 0.7 total time 
!     CHANGED: now calculate it based on the last time step!
!
   do k=1,nelems
      do i=1,nframes-frames_skip
         diff(i)=msd_func(i,k)/(6.d0*times(i))
      end do
      diff = diff*(1E-10)**2/(1E-15)
      avg_diff=diff(nframes-frames_skip)

      write(*,*) "calculated_diffusion coefficient, ",trim(el_names(k))," (m^2/s):",avg_diff
      open(unit=19,file="diff_const_MSD_"//trim(el_names(k))//".dat",status="replace")
      write(19,*) avg_diff,"m^2/s"
      close(19)
      
   end do


end if        

!
!    Calculate the vibrational density of states from the velocity autocorrelation
!     function (VACF)
!
!   Number of intervals for which VACF is calculated
!
if (calc_vacf) then

   intervals=int(((nframes-frames_skip)*time_step)/corrt)
   frames_part=int(corrt/time_step)
   write(*,*)
   write(*,*) "Number of intervals for VACF averaging: ",intervals
   write(*,*) "Number of MD frames per interval: ",frames_part
   write(*,*) "Calculate VACF from trajectory..."
   allocate(vel_first(natoms*3))
   allocate(vel_act(natoms*3))
   allocate(vacf_func(frames_part))
   vacf_func=0.d0
   do l=1,intervals
      do i=1,natoms
         do k=1,3
            vel_first((i-1)*3+k)=(xyz(k,i,2+(l-1)*intervals)-xyz(k,i,1+(l-1)*intervals))/time_step
         end do
      end do  
      do i=1,frames_part-1
         do j=1,natoms
            do k=1,3
               vel_act((j-1)*3+k)=(xyz(k,j,i+1+(l-1)*intervals)-xyz(k,j,i+(l-1)*intervals))/time_step
            end do
         end do
         vacf_func(i)=vacf_func(i)+dot_product(vel_first,vel_act)/natoms/3.d0
      end do
   end do

   vacf_func=vacf_func/vacf_func(1)
   open(unit=18,file="vacf_plot.dat",status="replace")
   do i=1,frames_part
      write(18,*) (i-1)*time_step, vacf_func(i)
   end do   
   close(18)

   write(*,*) "Done! VACF written to vacf_plot.dat"
   write(*,*)
   write(*,*) "Perform Fourier transform for vibrational density of states..."
!
!    Now execute the discrete fourier transform to obtain the frequency spectrum
!    Use zero-padding to increase the resolution in the frequency domain
!
   allocate(vacf_pad(1*frames_part))
   allocate(fourier_out(1*frames_part))
   vacf_pad=0.d0
   vacf_pad(1:frames_part)=vacf_func
!   call dfftw_plan_dft_r2c_1d(plan,frames_part*1,vacf_pad,fourier_out,FFTW_ESTIMATE)
!   call dfftw_execute_dft_r2c(plan, vacf_pad,fourier_out)
   call rfft(vacf_pad, frames_part*1)

!
!    Generate a smoothened plot of the VDOS
!    Perform Gaussian broadening of all FFT peaks 
!
   allocate(vacf_plot(vib_inter))
   i_max=int(6.0/33.356*frames_part*time_step)
   do i=1,vib_inter
      int_act=0.d0
      do j=1,i_max
         int_act=int_act+sqrt(vacf_pad(j)*vacf_pad(j))*exp(-vib_width*(real(i)-(j*&
                          & 33.356/(frames_part*time_step)*1000d0))**2)
      end do
      vacf_plot(i)=int_act
   end do


   open(unit=18,file="VDOS.dat",status="replace")
   write(18,*) "# This VDOS has been calculated by analyze_bulk from utils4VASP."
   write(18,*) "# Only frequencies up to ",vib_inter," cm^-1 are plotted."
   write(18,*) "# Wave number (cm^-1)      intensity (a.u.)"
!   do i=1,frames_part*2
!      if ((i*33.356/(frames_part*time_step)*1000d0) .lt. 6000.d0) then
!         write(18,*) i*33.356/(frames_part*time_step)*1000d0,sqrt(vacf_pad(i)*vacf_pad(i))
!      end if
!   end do
   do i=1,vib_inter
      write(18,*) real(i),vacf_plot(i)
   end do
   close(18)
   write(*,*) "Done! VDOS written to tile VDOS.dat"
end if        


if (calc_rdf) then
           write(*,*) "Calculate the RDFs of all element combinations..."
   eval_stat=.false.
   allocate(rdf_plot(rdf_bins,nelems,nelems))
   allocate(rdf_sum(rdf_bins,nelems,nelems))
   rdf_plot=0.d0
   task_act=0
   xlen=a_len
   ylen=b_len
   zlen=c_len
   all_tasks=((nelems**2-nelems)/2+nelems)*(nframes-frame_first)
   do l=1,nelems
      do m=l,nelems
         do i=frame_first,nframes
            if (npt_traj) then
               xlen=a_lens(i)
               ylen=b_lens(i)
               zlen=c_lens(i)
            end if     
            task_act=task_act+1
            !
!    Every 10% of the process, give a status update
!
            do j=1,10
               if (real(task_act)/real(all_tasks) .gt. real(j)*0.1d0) then
                  if (.not. eval_stat(j)) then
                     write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
                     eval_stat(j) = .true.
                  end if
               end if
            end do

            do j=1,el_nums(l)  ! Ni atoms
               do k=1,el_nums(m)  ! Ga atoms
                  if (l .gt. 1) then
                     pos1 = xyz(:,sum(el_nums(1:l-1))+j,i)
                  else
                     pos1 = xyz(:,j,i)
                  end if
                  if (m .gt. 1) then
                     pos2 = xyz(:,sum(el_nums(1:m-1))+k,i)
                  else
                     pos2 = xyz(:,k,i)
                  end if

                  diff_vec=pos1-pos2
!
!     Correct the x component
!
                  do while (abs(diff_vec(1)) .gt. xlen/2.d0)
                     diff_vec(1)=diff_vec(1)-sign(xlen,diff_vec(1))
                  end do

!
!     Correct the y component
!
                  do while (abs(diff_vec(2)) .gt. ylen/2.d0)
                     diff_vec(2)=diff_vec(2)-sign(ylen,diff_vec(2))
                  end do

!
!     Correct the z component
!
                  do while (abs(diff_vec(3)) .gt. zlen/2.d0)
                     diff_vec(3)=diff_vec(3)-sign(zlen,diff_vec(3))
                  end do

                  dist = sqrt((diff_vec(1))**2 + &
                        & (diff_vec(2))**2 + (diff_vec(3))**2)

!
!     Remainder of the calculation taken from Frenkel Smit, page 86
!
                  ig = int(dist/rdf_binsize)
                  if ((ig .le. rdf_bins) .and. (ig .ge. 1)) then
                     rdf_plot(ig,l,m) = rdf_plot(ig,l,m) + 2.d0
                     rdf_sum(ig,l,m)=rdf_sum(ig,l,m)+ 1.d0
                  end if
               end do
            end do
         end do
         do j=1,rdf_bins
            ngr=nframes-frame_first
            npart1=el_nums(l)  ! which of both elements?
            npart2=el_nums(m)
            r_act=rdf_binsize*(real(j)+0.5d0)
            vb=((real(j) + 1.d0)**3-real(j)**3)*rdf_binsize**3
            rho=1.d0/(abs(volume))
            nid=4.d0/3.d0*pi*vb*rho*2d0
            rdf_plot(j,l,m)=rdf_plot(j,l,m)/(ngr*npart1*npart2*nid)
            rdf_plot(j,m,l)=rdf_plot(j,l,m)
            rdf_sum(j,l,m)=rdf_sum(j,l,m)/ngr
            rdf_sum(j,m,l)=rdf_sum(j,l,m)            
         end do
      end do
   end do
   write(*,*) " completed!"
   rdf_plot(1,:,:)=0.d0
   call system("mkdir RDFs")
   call chdir("RDFs")
   do i=1,nelems
      do j=1,nelems
!
!    Write the RDF for the current element combination
!
         open(unit=13,file="RDF_"//trim(el_names(i))//"-"//trim(el_names(j))//".dat", &
                       & status="replace")
         write(13,'(a,a,a,a,a)') "#      Distance (A)        RDF ("&
                       &//trim(el_names(i))//"-"//trim(el_names(j))//") "      
         do k=1,rdf_bins
            write(13,'(2f19.10)') k*rdf_binsize,rdf_plot(k,i,j)
         end do
!
!    Write the integrated RDF (i.e., the number of surrounding atoms) for the 
!    current element combination)
!
         open(unit=13,file="RDF_int_"//trim(el_names(i))//"-"//trim(el_names(j))//".dat", &
                       & status="replace")
         write(13,'(a,a,a,a,a)') "#      Distance (A)       integrated RDF ("&
                       &//trim(el_names(i))//"-"//trim(el_names(j))//") "
         
         do k=1,rdf_bins
            write(13,'(2f19.10)') k*rdf_binsize,sum(rdf_sum(1:k,i,j))/el_nums(i)
         end do

      end do
   end do   
   call chdir("..")

 !  open(unit=14,file="nearest.dat",status="replace")
 !  do i=1,rdf_bins
 !     write(14,*) i*rdf_binsize,neighnum(i)
 !  end do
 !  close(14)

   write(*,*) "RDF plots written to files in folder RDFs."
   write(*,*) "Summed up numbers of surrounding atoms written to "
   write(*,*) "   RDF_int files in folder RDFs."
   write(*,*)
end if

!
!     Calculate the averaged spatially resolved element densities and store
!     them in a Gaussian cube file
!
if (dens_cube) then
   write(*,*) "Calculate the spatially resolved element densities!"
!
!     Determine number grid points per dimension
!
   gridx=100
   gridy=100
   gridz=100
   allocate(dens_3d(nelems,gridx,gridy,gridz))
   dens_3d=0.d0

!
!    Convert back to direct coordinates 
!
   if (.not. npt_traj) then
      do i=1,nframes
         do j=1,natoms
            xyz(1,j,i) = xyz(1,j,i)/a_len
            xyz(2,j,i) = xyz(2,j,i)/b_len
            xyz(3,j,i) = xyz(3,j,i)/c_len
         end do
      end do
   else
      do i=1,nframes
         do j=1,natoms
            xyz(1,j,i) = xyz(1,j,i)/a_lens(i)
            xyz(2,j,i) = xyz(2,j,i)/b_lens(i)
            xyz(3,j,i) = xyz(3,j,i)/c_lens(i)
         end do
      end do
   end if
   
!
!     Remove corrections of box images, project all atoms into central unit cell
!
   do i=1,nframes
      do j=1,natoms
!
!     Correct the x component
!
         do while(xyz(1,j,i) .gt. 1.d0)
            xyz(1,j,i)=xyz(1,j,i)-1.d0 
         end do
         do while(xyz(1,j,i) .lt. 0.d0)
            xyz(1,j,i)=xyz(1,j,i)+1.d0
         end do
!
!     Correct the y component
!
         do while(xyz(2,j,i) .gt. 1.d0)
            xyz(2,j,i)=xyz(2,j,i)-1.d0 
         end do
         do while(xyz(2,j,i) .lt. 0.d0)
            xyz(2,j,i)=xyz(2,j,i)+1.d0
         end do
!
!     Correct the z component
!
         do while(xyz(3,j,i) .gt. 1.d0)
            xyz(3,j,i)=xyz(3,j,i)-1.d0 
         end do
         do while(xyz(3,j,i) .lt. 0.d0)
            xyz(3,j,i)=xyz(3,j,i)+1.d0
         end do
      end do
   end do
!
!     Now assign all atoms in each frame to one of the 3D bins
!
   do i=frame_first,nframes
      inc=0
      do j=1,nelems
         do k=1,el_nums(j)
            inc=inc+1
            gridx_act=ceiling(xyz(1,inc,i)*gridx)
            gridy_act=ceiling(xyz(2,inc,i)*gridy)
            gridz_act=ceiling(xyz(3,inc,i)*gridz)-1
            do i1=-stepx,stepx
               do j1=-stepy,stepy
                  do k1=-stepz,stepz
                     distance=i1**2+j1**2+k1**2
                     i_new=gridx_act+i1
                     if (i_new .gt. gridx) then
                        i_new=i_new-gridx
                     else if (i_new .lt. 1) then
                        i_new=i_new+gridx
                     end if
                     j_new=gridy_act+j1
                     if (j_new .gt. gridy) then
                        j_new=j_new-gridy 
                     else if (j_new .lt. 1) then
                        j_new=j_new+gridy
                     end if
                     k_new=gridz_act+k1
                     if (k_new .gt. gridz) then
                        k_new=k_new-gridz
                     else if (k_new .lt. 1) then
                        k_new=k_new+gridz
                     end if
                     
                     dens_3d(j,i_new,j_new,k_new)= &
                             & dens_3d(j,i_new,j_new,k_new)+exp(-(distance)/2.d0)                     
                  end do 
               end do     
            end do       
         end do
      end do      
   end do
!
!     Write the header of the cube file
!
   do i=1,nelems
      write(cube_name,'(a,a,a)') "density_",trim(el_names(i)),".cube"
      write(*,*) "Write density of ",trim(el_names(i))," to file ",cube_name
      open(unit=45,file=cube_name,status="replace")
      write(45,'(a)') "utils4VASP CUBE FILE"
      write(45,'(a,a)') "Contains spatially resolved density of element:",el_names(i)
      write(45,*) natoms,0d0,0.d0,0d0
      if (.not. npt_traj) then
         write(45,*) gridx,a_len/gridx*a2bohr,0.d0,0.d0
         write(45,*) gridy,0.d0,b_len/gridy*a2bohr,0.d0
         write(45,*) gridz,0.d0,0.d0,c_len/gridz*a2bohr
      else
         write(45,*) gridx,a_lens(1)/gridx*a2bohr,0.d0,0.d0
         write(45,*) gridy,0.d0,b_lens(1)/gridy*a2bohr,0.d0
         write(45,*) gridz,0.d0,0.d0,c_lens(1)/gridz*a2bohr
      end if              
      inc=0
      do j=1,nelems
         do k=1,el_nums(j)
            inc=inc+1
            call elem(el_names(j),elnumber)
            if (.not. npt_traj) then
               write(45,*) elnumber,real(elnumber),xyz(1,inc,nframes)*a_len*a2bohr,xyz(2,inc,nframes)* &
                              & b_len*a2bohr,xyz(3,inc,nframes)*c_len*a2bohr
            else 
               write(45,*) elnumber,real(elnumber),xyz(1,inc,nframes)*a_lens(1)*a2bohr,xyz(2,inc,nframes)* &
                              & b_lens(1)*a2bohr,xyz(3,inc,nframes)*c_lens(1)*a2bohr
            end if             
         end do
      end do     
      write(45,*) 1,48 
!
!    Write the volumetric data in the cube file
!
      do j=1,gridx
         do k=1,gridy
            do l=1,gridz
               write(45,'(f20.10)',advance="no") dens_3d(i,j,k,l)
               if (modulo(l,6) .eq. 5) then
                  write(45,*) " "      
               end if
            end do
            write(45,*) " "
         end do    
      end do
      close(45)   
   end do
end if        


write(*,*) "analyze_bulk exited normally ..."
write(*,*)

end program analyze_bulk 


!
!     subroutine rfft: compute the real fast fourier transform for a
!      given array of data
!      Rewritten 09.10.2023 (updated to Fortran03/C interface)
!    
!     part of EVB
!
subroutine rfft(x,N)
use fftw_mod
implicit none
integer,intent(in) :: N
real(kind=8),intent(inout) :: x(N)
complex(C_DOUBLE_COMPLEX),dimension(N)::ain,aout

ain=x(1:N)

!
!     If this is the first execution, generate the FFT plan (optimized code for 
!     local machine)
!     If the plan already was generated for a different number of beads, destroy
!     it in advance
!
if (N .ne. Np) then
    if (Np .ne. 0) call fftw_destroy_plan(plan)
    plan=fftw_plan_dft_1d(N,ain,aout,FFTW_FORWARD,FFTW_ESTIMATE)
    factor = dsqrt(1.d0/N)
    Np = N
end if
call fftw_execute_dft(plan,ain,aout)

x = factor * real(aout(1:N))

return
end subroutine rfft


!
!     subroutine elem: read in character with element symbol and
!       give out the number
!
!     part of QMDFF
!
subroutine elem(key1, nat)
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CHARACTER(len=*)::KEY1
CHARACTER(len=2)::ELEMNT(107),E

DATA ELEMNT/'h ','he', &
 & 'li','be','b ','c ','n ','o ','f ','ne', &
 & 'na','mg','al','si','p ','s ','cl','ar', &
 & 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu', &
 & 'zn','ga','ge','as','se','br','kr', &
 & 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag', &
 & 'cd','in','sn','sb','te','i ','xe', &
 & 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy', &
 & 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt', &
 & 'au','hg','tl','pb','bi','po','at','rn', &
 & 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx', &
 & 'fm','md','cb','xx','xx','xx','xx','xx'/

nat=0
e='  '
do i=1,len(key1)
   if (key1(i:i).ne.' ') L=i
end do
k=1
DO J=1,L
   if (k.gt.2) exit
   N=ICHAR(key1(J:J))
   if (n.ge.ichar('A') .and. n.le.ichar('Z') ) then
      e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
      k=k+1
   end if
   if (n.ge.ichar('a') .and. n.le.ichar('z') ) then
      e(k:k)=key1(j:j)
      k=k+1
   end if
end do

DO I=1,107
   if (e.eq.elemnt(i)) then
      NAT=I
      RETURN
   END IF
END DO

return
end subroutine elem

