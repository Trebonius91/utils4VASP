!
!     This program reads in a XDATCAR file and calculates 
!     different dynamical observables 
!     (like mean-square-displacements, velocity autocorrelations)
!     from it
!     Finally, diffusion coefficients are calculated as well
!
program eval_md
implicit none
integer::i,j,k,l,m
integer::readstat
logical::calc_msd,calc_vacf,read_dt,npt_traj,print_new
real(kind=8)::a_read(3),b_read(3),c_read(3)
real(kind=8)::a_len,b_len,c_len
real(kind=8),allocatable::a_lens(:),b_lens(:),c_lens(:)
integer::nelems,natoms,nframes,xdat_lines
integer,allocatable::el_nums(:)
integer::frames_skip
real(kind=8),allocatable::xyz(:,:,:)
real(kind=8)::diff_vec(3)
real(kind=8)::time_step
real(kind=8)::msd_act,volume
real(kind=8),allocatable::pos_diff(:)
real(kind=8),allocatable::z_axis(:)
real(kind=8),allocatable::msd_func(:,:),vacf_func(:,:)
real(kind=8),allocatable::vel_first(:),vel_act(:)
real(kind=8),allocatable::times(:),diff(:)
real(kind=8),allocatable::vector1(:),vector2(:),vector3(:)
real(kind=8)::vacf_int
character(len=2),allocatable::el_names(:)
real(kind=8),allocatable::rdf_plot(:,:,:)
real(kind=8)::pos1(3),pos2(3)
integer::avg_lo,avg_hi
real(kind=8)avg_diff
character(len=32)::arg
character(len=80)::line,all_els
character(len=2)::el1,el2,el3
logical::eval_stat(10)
real(kind=8)::box_volume
real(kind=8)::dist,pi,rdf_cutoff
real(kind=8)::xlen,ylen,zlen
integer::frame_first
!  RDF calculation
integer::ig,ngr,npart1,npart2
real(kind=8)::nid,r_act,rho,vb
real(kind=8)::slice_step,rdf_binsize
real(kind=8),allocatable::z_dens(:,:),z_dens_tot(:)
real(kind=8),allocatable::time_list(:)
integer::all_tasks,task_act,rdf_bins
logical::calc_rdf,read_time

pi=3.141592653589793238

write(*,*) "This program calculates dynamical observables from MD trajectories!"
write(*,*) "Only a XDATCAR file in VASP format needs to be present."
write(*,*) "Usage: eval_vasp_md -command1 -command2 ..."
write(*,*)
write(*,*) "List of commands:"
write(*,*) "-msd:  the mean square displacement (and diffusion coefficient) is calculated."
write(*,*) "-vacf:  the velocity autocorrelation function (and diffusion coefficient) "
write(*,*) "        is calculated."
write(*,*) "-dt=[time step in fs]: The time step used for MD simulation, in fs."
write(*,*) "-readtime: The time of each XDATCAR (in s) will be read in from file 'times.dat'"
write(*,*) "    useful for the evaluation of kinetic Monte Carlo simulations."
write(*,*) "-skip=[steps to skip]: If the first N steps shall be skipped (equilibration.)"
write(*,*) "-npt: If a NpT trajectory (variable volume) shall be analyzed."
write(*,*) "-rdf: Radial distribution functions will be calculated."
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
   if (trim(arg)  .eq. "-vacf") then
      calc_vacf = .true.
      write(*,*) "The velocity autocorrelation function (VACF) will be calculated!"
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

allocate(el_nums(nelems))
allocate(el_names(nelems))
read(15,*) el_nums

natoms=sum(el_nums)

write(*,*) nelems
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
if (nelems .eq. 3) then
   el_names(3)=el3     
   write(*,*) "  3. ",el3,": ",el_nums(3)," atoms"
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
   a_lens(1)=a_len
   b_lens(1)=b_len
   c_lens(1)=c_len

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
!   avg_hi= int((nframes-frames_skip)*0.7)
!   avg_lo = int((nframes-frames_skip)*0.3)
!   avg_diff = sum(diff(avg_lo:avg_hi))/(avg_hi-avg_lo)
      avg_diff=diff(nframes-frames_skip)
      if (k .eq. 1) then
         write(*,*) "calculated diffusion coefficient, element 1 (m^2/s):",avg_diff
         open(unit=19,file="diff_const_MSD_el1.dat",status="replace")
         write(19,*) avg_diff,"m^2/s"
         close(19)
      end if 
      if (k .eq. 2) then
         write(*,*) "calculated diffusion coefficient, element 2 (m^2/s):",avg_diff
         open(unit=19,file="diff_const_MSD_el2.dat",status="replace")
         write(19,*) avg_diff,"m^2/s"
         close(19)
      end if
      if (k .eq. 3) then
         write(*,*) "calculated diffusion coefficient, element 3 (m^2/s):",avg_diff
         open(unit=19,file="diff_const_MSD_el3.dat",status="replace")
         write(19,*) avg_diff,"m^2/s"
         close(19)
      end if
   end do


end if        

!
!    Calculate the velocity autocorrelation function (VACF)
!

if (calc_vacf) then
   allocate(vel_first(natoms*3))
   allocate(vel_act(natoms*3))
   allocate(vacf_func(nframes-frames_skip,1))
   do i=1,natoms
      do k=1,3
         vel_first((i-1)*3+k)=(xyz(k,i,2)-xyz(k,i,1))/time_step
      end do
   end do   
   do i=1,nframes-frames_skip-1
      do j=1,natoms
         do k=1,3
            vel_act((j-1)*3+k)=(xyz(k,j,i+1)-xyz(k,j,i))/time_step
         end do
      end do
      vacf_func(i,1)=dot_product(vel_first,vel_act)/natoms/3.d0
   end do

   open(unit=18,file="vacf_plot.dat",status="replace")
   do i=1,nframes-frames_skip-1
      write(18,*) vacf_func(i,1)
   end do   
   close(18)

   open(unit=19,file="vacf_integrate.dat",status="replace")
   vacf_int=0.d0
   do i=1,nframes-frames_skip-1
      vacf_int=vacf_int+vacf_func(i,1)
      write(19,*) vacf_int
   end do
   close(18)

end if        


if (calc_rdf) then
           write(*,*) "Calculate the RDFs of all element combinations..."
   eval_stat=.false.
   allocate(rdf_plot(rdf_bins,nelems,nelems))
   rdf_plot=0.d0
   task_act=0
   xlen=a_len
   ylen=b_len
   zlen=c_len
   all_tasks=((nelems**2-nelems)/2+nelems)*(nframes-frame_first)
   do l=1,nelems
      do m=l,nelems
         do i=1,nframes
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
         end do
      end do
   end do
   write(*,*) " completed!"
   rdf_plot(1,:,:)=0.d0
   open(unit=13,file="rdf.dat",status="replace")
   write(13,'(a)',advance="no") "#      Distance (A)        "
   do i=1,nelems
      do j=1,nelems
         write(13,'(a,a,a,a)',advance="no") trim(el_names(i)),"-",trim(el_names(j)),"       "
      end do
   end do
   write(13,*)
   do i=1,rdf_bins
      write(13,'(f19.10)',advance="no") i*rdf_binsize
      do j=1,nelems
         do k=1,nelems
            write(13,'(f19.10)',advance="no") rdf_plot(i,j,k)
         end do
      end do
      write(13,*)
   end do
   close(13)

 !  open(unit=14,file="nearest.dat",status="replace")
 !  do i=1,rdf_bins
 !     write(14,*) i*rdf_binsize,neighnum(i)
 !  end do
 !  close(14)

   write(*,*) "RDF plot written to 'rdf.dat'."
   write(*,*)
end if


end program eval_vasp_md      
