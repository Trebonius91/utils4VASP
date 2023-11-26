!
!    analyze_scalms: Analyze VASP or LAMMPS simulations of 
!      SCALMS systems from trajectories, preferentially surface slabs
!    Part of VASP4CLINT
!     Julien Steffen, 2023 (julien.steffen@fau.de)
!

program analyze_scalms
implicit none
integer::i,j,k,l,m,r,o,p
integer::xdat_lines,dump_lines
real(kind=8)::factor,rdum1,rdum2
real(kind=8)::xlen,ylen,zlen,zmax
real(kind=8),allocatable::xyz(:,:,:),xyz2(:,:,:) ! all coordinates of the trajectory
character(len=2),allocatable::el_names(:),el_names_read(:) ! symbols of the elements 
integer,allocatable::el_nums(:)  ! numbers of the elements 
integer::nelems ! number of elements in the system
integer::natoms,nframes,frame_first ! number of atoms and frames in the trajectory
real(kind=8)::act_num(3)
character(len=120)::a120
character(len=220)::a220
character(len=1)::atest
character(len=60)::filename,roundname
character(len=32)::arg,cdum,adum
character(len=2)::cls_element
integer::slice_size,frame_round_first,frame_round_last,cls_elem_ind
integer::atom_first,atom_last
integer::readstat,openstat
integer::counter,endl 
integer::cls_rounds
integer::nbins
integer::atom_slices  ! for writeout of different Ni z positions
!  RDF calculation
integer::ig,ngr,npart
real(kind=8)::nid,r_act,rho,vb
real(kind=8)::slice_step
real(kind=8),allocatable::z_dens(:,:)
real(kind=8)::zlo,zhi,zdiff,zstep  ! borders of z-density bins 
character(len=2),allocatable::at_names(:)  ! the element symbols
logical::write_traj,read_dt
logical::calc_rdf,calc_diff
logical::use_reaxff,el_present
logical::surf_tension
logical::skip_xdat
logical::eval_stat(10)
real(kind=8)::rdf_binsize,rdf_range
real(kind=8),allocatable::rdf_plot(:,:,:)
real(kind=8),allocatable::neighnum(:,:,:)
real(kind=8)::pos1(3),pos2(3),diff_vec(3)
real(kind=8)::dist,pi
real(kind=8)::z_shift
real(kind=8)::pres_tensor(3,3),pres_tensor_total(3,3)  ! the pressure tensor for surface tension
real(kind=8)::pres_xx,pres_yy,pres_zz,pres_xy,pres_yz,pres_zx,tension
real(kind=8),allocatable::vector1(:),vector2(:),vector3(:),pos_diff(:)
real(kind=8),allocatable::diff(:),times(:),msd_func(:,:)
real(kind=8)::avg_diff,time_step
real(kind=8)::box_volume
!  Time and accounting
integer::task_act,all_tasks


integer::rdf_bins


write(*,*) "PROGRAM analyze_scalms: Evaluation of VASP ML-FF (or LAMMPS ReaxFF)"
write(*,*) " trajectories for binary and ternary SCALMS."
write(*,*) "The file XDATCAR (or dump.xyz in case of xyz trajectories) must be present!"
write(*,*) "The following command line arguments can/must be given (with - sign!):"
write(*,*) " -reaxff : A ReaxFF xyz trajectory 'dump.xyz' will be analyzed."
write(*,*) "     here, the box dimensions must be given separately in 'box_dims.dat',"
write(*,*) "     with the format 'xlen  ylen  zlen' (one line)"
write(*,*) " -notraj : No file 'trajectory.xyz' shall be written during the analysis."
write(*,*) " -rdf' : The radial distribution functions of the second component around "
write(*,*) "     the first shall be calculated. Then, also the total number of neighbors "
write(*,*) "     up to a certain distance will be calculated."
write(*,*) " -frame_first=[number] : First trajectory frame that shall be evaluated by "
write(*,*) "     the script (e.g., in order to skip equilibration parts) (default: 1)"
write(*,*) " -dens_slices=[number] : In how many parts the z-axis shall be divided for "
write(*,*)"      the evaluation of element densities (default: 501)."
write(*,*) " -z_shift=[value] : The z-coordinates of the frames are shifted by the value,"
write(*,*) "     given in direct coordinates (0 to 1.0)."
write(*,*) " -cls_element=[element]: CLS calculation templates will be generated for"
write(*,*) "     the chosen element."
write(*,*) " -cls_slices=[number]: How many different slices along the z-coordinate where "
write(*,*) "     the atom for which CLS shall be calculated is located (near for far from the "
write(*,*) "     surface of the slab (default: 100)."
write(*,*) " -cls_rounds=[number]: In how many parts the trajectory shall be divided for"
write(*,*) "     CLS template generation in each part (default: 1)."
write(*,*) " -tension : calculates the surface tension averaged over all MD frames."
write(*,*) "     For this, the OUTCAR file needs to be present (MD with ISIF=2)"
write(*,*) " -diffusion : calculates the diffusion coefficient, for each element in the slab"
write(*,*) "     separately, via the mean square displacement (MSD)."
write(*,*) " -timestep=[value] : For diffusion calculations, the time step in fs (consider"
write(*,*) "     a longer step if not every step was written during dynamics!)"

use_reaxff = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-reaxff") then
      use_reaxff = .true.
   end if
end do

!
!    If the trajectory shall not be written out as xyz (only for VASP)
!
if (.not. use_reaxff) then
   write_traj = .true.
end if

do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-notraj") then
      write_traj = .false.
   end if
end do


!
!    Look if Radial Distribution functions shall be calculated 
!
calc_rdf = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-rdf") then
      calc_rdf = .true.
   end if
end do


pi=3.141592653589793238
!
!    For bin packing of abundancies/densities along z-axis
!
nbins = 501
atom_slices=100
!
!    For RDF calculation
!
rdf_bins = 200  ! number of RDF bins 
rdf_range = 8.d0  ! maximum distance for RDF
rdf_binsize = rdf_range/rdf_bins
!    
!    For (optional) diffusion coefficient calculation
! 
time_step=0.d0
!
!    If the evaluation shall start not at the first frame but after 
!     an equilibration phase
!
frame_first = 1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-frame_first=") then
      read(arg(14:),*) frame_first
      write(*,*) "The analysis will start from frame ",frame_first
   end if
end do

z_shift = 0.0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-z_shift=") then
      read(arg(10:),*) z_shift
      write(*,*) "The z-coordinates (direct) will be shifted by ",z_shift
   end if
end do

!
!    The number of slices the z-axis will be divided for element densities
!
nbins=501
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-dens_slices=") then
      read(arg(14:),*) nbins
   end if
end do

!
!    The element for which the CLS shall be calculated
!
cls_element="XX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-cls_element=") then
      read(arg(14:),*) cls_element
   end if
end do

!
!    The number of individual core level shift evaluations
!
cls_rounds = 1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-cls_rounds=") then
      read(arg(13:),*) cls_rounds
   end if
end do

!
!    The number of slices along the z-axis for which CLS values shall be evaluated
!
atom_slices = 100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:12))  .eq. "-cls_slices=") then
      read(arg(13:),*) atom_slices
   end if
end do

!
!    Activates the calculation of surface tension if desired
!
surf_tension=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-tension") then
      surf_tension=.true.
   end if
end do

!
!    Activates the calculation of diffusion coefficients 
!

calc_diff=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-diffusion") then
      calc_diff=.true.
   end if
end do



if (calc_diff) then
   read_dt=.false.
   do i = 1, command_argument_count()
      call get_command_argument(i, arg)
      if (trim(arg(1:10))  .eq. "-timestep=") then
         read_dt = .true.
         read(arg(11:),*) time_step
         write(*,*) "The time step shall be:",time_step," fs."
      end if
   end do
   if (.not. read_dt) then
      stop "Please set a time step with the -timestep=... flag!"
   end if
end if        

!
!    Look if the xyz trajectory file shall be written
!
write_traj = .true.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-notraj") then
      write_traj = .false.
      write(*,*) "The file 'trajectory.xyz' will not be written!"
   end if        
end do   
!
!    Look if Radial Distribution functions shall be calculated 
!
calc_rdf = .false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg)  .eq. "-rdf") then
      calc_rdf = .true.
      write(*,*) "Radial distribution functions will be calculated and "
      write(*,*) "  written to 'rdf.dat'!"
   end if
end do

!
!    If a surface tension calculation is requested and the flag -notraj
!    is turned on (but no other jobs), skip the read in of the XDATCAR 
!    file alltogether!
!
skip_xdat = .false.
if (.not. write_traj .and. surf_tension) then
   write(*,*) "The calculation of surface tension was required as only"
   write(*,*) " job and no trajectory shall be written out." 
   write(*,*) "Therefore, the read in of the XDATCAR file will be skipped!"
   skip_xdat = .true.
end if        


!
!    First, determine the number of lines in the XDATCAR file
!
if (use_reaxff) then
   call system("wc -l dump.xyz > dump_length")
   open(unit=45,file="dump_length",status="old")
   read(45,*) dump_lines
   close(45)
else         
   call system("wc -l XDATCAR > xdat_length")
   open(unit=45,file="xdat_length",status="old")
   read(45,*) xdat_lines
   close(45)
end if
!
!    Open the dump.xyz file (ReaxFF) and read in the frames 
!
if (use_reaxff) then
!
!    Determine cell dimensions from maximum x,y and z values during dynamics
!    read the dimensions from file 'box_dims.dat'

   open(unit=15,file="box_dims.dat",status="old",iostat=openstat)
   if (openstat .ne. 0) then
      stop "ERROR! The file 'box_dims.dat' could not been found!"
   end if
   read(15,*) xlen,ylen,zlen
   close(15)
!
!    Read in the dump.xyz file. If if shall be skipped, jump to the next section
!
   open(unit=14,file="dump.xyz",status="old",iostat=openstat)
   if (openstat .ne. 0) then
      stop "ERROR! The file 'dump.xyz' could not been found!"
   end if        
   read(14,*) natoms
   nframes = (dump_lines)/(natoms+2)
   if (skip_xdat) goto 33
   allocate(xyz(3,natoms,nframes))
!
!    Search first frame for element names 
!
   allocate(at_names(natoms))
   allocate(el_names(10))
   el_names="XX"
   read(14,*)
   do i=1,natoms
      read(14,*) at_names(i),xyz(:,i,1)
   end do   
   
   nelems=0
   do i=1,natoms
      el_present=.false.
      do j=1,i-1
         if (at_names(j) .eq. at_names(i)) el_present=.true.
      end do
      if (el_present) cycle
      nelems=nelems+1
      el_names(nelems)=at_names(i)
   end do

   do i=1,nframes-1
      read(14,*)
      read(14,*)
      do j=1,natoms
         read(14,*) adum,xyz(1,j,i+1),xyz(2,j,i+1),xyz(3,j,i+1)
      end do
   end do   
!
!    Determine cell dimensions from maximum x,y and z values during dynamics
!    read the dimensions from file 'box_dims.dat'

   open(unit=15,file="box_dims.dat",status="old",iostat=openstat)
   if (openstat .ne. 0) then
      stop "ERROR! The file 'box_dims.dat' could not been found!"
   end if        
   read(15,*) xlen,ylen,zlen
   close(15)
!
!    Box volume for surface tension calculation
!
   box_volume=xlen*ylen*zlen
!
!     We assume that the bulk is located in the lower half of the simulation
!     cell. If atoms go through the lower x-y surface z-values near 1,
!     move them to values close below zero for better appearance
!
!   zmax=zlen*0.8

!   do i=1,nframes-1
!      do j=1,natoms
!         if (xyz(3,j,i) .gt. zmax) then
!            xyz(3,j,i)=xyz(3,j,i)-zlen
!         end if     
!      end do
!   end do

!
!     Apply the z-shift if defined 
!
   do i=1,nframes-1
      do j=1,natoms
         xyz(3,j,i)=xyz(3,j,i)+z_shift
      end do
   end do


!
!    Open the XDATCAR file and read in the frames 
!
        
else         
   open(unit=14,file="XDATCAR",status="old",iostat=openstat)
   if (openstat .ne. 0) then 
      stop "ERROR! The file 'XDATCAR' could not been found!"
   end if        
   read(14,*)
   read(14,*) factor   ! the global geometry conversion factor 
!
!    Read in the cell dimensions: Cubic unit cell is always assumed!
!
   read(14,*) xlen,rdum1,rdum2
   read(14,*) rdum1,ylen,rdum2
   read(14,*) rdum1,rdum2,zlen
!
!    Box volume for surface tension calculation
!
   box_volume=xlen*ylen*zlen
!
!    Read in the elements (so far, only two different allowed)
!
   allocate(el_names_read(10))
   el_names_read="XX"
   read(14,'(a)') cdum
   read(cdum,*,iostat=readstat) el_names_read
   nelems=0
   do i=1,10
      if (el_names_read(i) .eq. "XX") exit
         nelems=nelems+1
      if (el_names_read(i) .eq. cls_element) then
         cls_elem_ind=nelems
      end if        
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

   counter = 1
   do i=1,nelems
      do j=1,el_nums(i)
         at_names(counter) = el_names(i)
         counter = counter +1
      end do
   end do

   nframes = (xdat_lines - 7)/(natoms+1)
!
!    If the XDATCAR file shall not be read in, skip the rest
!
   if (skip_xdat) then
      close(14)
      goto 33
   end if

   allocate(xyz(3,natoms,nframes))
!
!    Read in the coordinates of the frames and correct for image flags
!
   eval_stat = .false.
   write(*,*)
   write(*,*) "Read in the trajectory from XDATCAR..."
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
                  write(a220(1:endl+1),'(a,a)') a220(1:endl),atest
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
               else 
                  exit
               end if   
            end do
            do 
               if (act_num(k) < 0d0) then
!
!    Special case: move atoms near the lower border to negative values
!
!
                  if (act_num(k) >= -0.2d0) then
                     exit 
                  end if        
                  act_num(k) = act_num(k) + 1d0   
               else      
                  exit
               end if   
            end do
!
!     We assume that the bulk is located in the lower half of the simulation
!     cell. If atoms go through the lower x-y surface z-values near 1, 
!     move them to values close below zero for better appearance
!

            if (act_num(k) .gt. 0.9d0) then
               act_num(k) = act_num(k)-1.d0
            end if        

            if (k .eq. 1) then
               xyz(k,j,i) = act_num(k)*xlen
            end if
            if (k .eq. 2) then   
               xyz(k,j,i) = act_num(k)*ylen
            end if   
            if (k .eq. 3) then
               xyz(k,j,i) = act_num(k)*zlen
            end if   
         end do
      end do
   end do
   write(*,*) " completed!"
   close(14)
end if
allocate(z_dens(nbins,nelems))
write(*,*)
write(*,*) "---------- SETTINGS ---------------------------"
write(*,*) "Number of atoms in the system:",natoms
write(*,*) "Number of frames in the trajectory:",nframes
if (calc_rdf) then
   write(*,*) "Radial distribution functions will be calculated."
else
   write(*,*) "No radial distribution functions will be calculated."
end if
if (write_traj) then
   write(*,*) "The trajectory will be written to 'trajectory.xyz'."
else
   write(*,*) "No xyz trajectory will be written."
end if

write(*,*) "The first ",frame_first," frames will be skipped!"
write(*,*) "Number of slices along z-axis for element densities:",nbins
write(*,'(a,a,a,i3,a)') " CLS will be calculated for element: ",cls_element," (index: ",cls_elem_ind,")"
write(*,*) "Number of slices along z-axis for CLS calculations:",atom_slices
write(*,*) "Number of trajectory parts for CLS calculations:",cls_rounds
write(*,*) "-------------------------------------------------"
write(*,*)


!
!    Determine lowest and highest z-values in the coordinates 
!
zlo = minval(xyz(3,:,:))
zhi = maxval(xyz(3,:,:))
!write(*,*) "z-borders were defined manually!"
!zlo = -2.63d0
!zhi = 39.31d0
write(*,*) "Lowest z-value during dynamics: ",zlo
write(*,*) "Highest z-value during dynamics: ",zhi
write(*,*)

zdiff = zhi-zlo
zstep = zdiff/(nbins-1)
!
!    Write the trajectory in xyz format to file 
!
if (write_traj) then
   write(*,*)
   write(*,*) "Write the trajectory of the system to 'trajectory.xyz'..."     
   open(unit=15,file="trajectory.xyz",status="replace")
   eval_stat=.false.
   do i=1,nframes 
      do j=1,10
         if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do

      write(15,*) natoms 
      write(15,*) "Frame No.",i
      do j=1,natoms 
         write(15,*) at_names(j), xyz(:,j,i)
      end do
   end do
   close(15)
   write(*,*) " completed!"
   write(*,*)
end if
!
!    Calculate the RDFs of all elements if desired 
!
if (calc_rdf) then
   write(*,*) "Calculate the RDFs of all element combinations..."
   eval_stat=.false.
   allocate(rdf_plot(rdf_bins,nelems,nelems))
   allocate(neighnum(rdf_bins,nelems,nelems))
   rdf_plot=0.d0
   neighnum=0.d0
   task_act=0
   all_tasks=((nelems**2-nelems)/2+nelems)*(nframes-frame_first)
   do l=1,nelems
      do m=l,nelems
         do i=frame_first,nframes
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
                     pos1 = xyz(:,sum(el_nums(1:l-1))+j,i+frame_first)
                  else 
                     pos1 = xyz(:,j,i)
                  end if   
                  if (m .gt. 1) then
                     pos2 = xyz(:,sum(el_nums(1:m-1))+k,i+frame_first)
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
            npart=el_nums(l)  ! which of both elements?
            r_act=rdf_binsize*(real(j)+0.5d0)
            vb=((real(j) + 1.d0)**3-real(j)**3)*rdf_binsize**3
            rho=npart/(abs(box_volume))
            nid=4.d0/3.d0*pi*vb*rho 
            rdf_plot(j,l,m)=rdf_plot(j,l,m)/(ngr*npart*nid)
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
!
!    Evaluate the distribution of atoms along the z axis
!
write(*,*) "Calculate element density distributions along z-axis..."
z_dens = 0.d0
if (use_reaxff) then
   do i=frame_first,nframes
      do l=1,natoms
         do j=1,nbins-1
            if ((xyz(3,l,i) .ge. zlo+(j-1)*zstep) .and. (xyz(3,l,i) &
                 & .le. zlo+j*zstep)) then
               do k=1,nelems
                  if (at_names(l) .eq. el_names(k)) then
                     z_dens(j,k) = z_dens(j,k) + 1.d0
                  end if        
               end do       
            end if 
         end do   
      end do
   end do
else         
   do i=frame_first,nframes
      counter = 1
      do j=1,nelems
         do k=1,el_nums(j)
            do l=1,nbins-1
               if ((xyz(3,counter,i) .ge. zlo+(l-1)*zstep) .and. (xyz(3,counter,i) & 
                            & .le. zlo+l*zstep)) then
                  z_dens(l,j) = z_dens(l,j) + 1.d0

               end if             
            end do
            counter = counter + 1
         end do
      end do
   end do
end if
write(*,*) " completed!"
!
!    Write the density profile to file 
!
open(unit=16,file="dens_elems.dat",status="replace")
write(16,'(a)',advance="no") " # z-coordinate    "
do i=1,nelems
   write(16,'(a,a)',advance="no") el_names(i),"      "
end do  
write(16,*)
do i=1,nbins-1
   write(16,*) zlo+(i-0.5d0)*zstep,z_dens(i,:)/(nframes*xlen*ylen*zstep)
end do
close(16)

!open(unit=17,file="active_positions.dat",status="replace")
!
!    track positions of individual Ni atoms in the bulk in z-direction
!
!do i=frame_first,nframes
!   do j=el_nums(1)+1,sum(el_nums)
!      write(17,'(f15.7)',advance="no") xyz(3,j,i)      
!   end do
!   write(17,*)
!end do

!close(17)
slice_step = zdiff/(atom_slices)

!
!     Calculates the averaged surface tension of the system if desired
!     Opens the OUTCAR file, reads the components of the pressure tensor,
!     averages them and finally calculates the pressure along the z axis
!
33 continue
if (surf_tension) then
   write(*,*)
   write(*,*) "Calculate surface tension ..." 
   task_act=0
   all_tasks=(nframes-frame_first)
   eval_stat=.false.
   if (use_reaxff) then



   else        
      open(unit=18,file="OUTCAR",status="old",iostat=openstat)
      if (openstat .ne. 0) then
         stop "ERROR! The file 'OUTCAR' could not been found!"
      end if
      tension=0.d0
      pres_tensor_total=0.d0
      do 

         pres_tensor=0.d0

         read(18,'(a)',iostat=readstat) a120 
         if (readstat .ne. 0) exit
         if (index(a120,"stress matrix after NEB project") .ne. 0) then
!
!    Every 10% of the process, give a status update
!
            task_act=task_act+1
            do j=1,10
               if (real(task_act)/real(all_tasks) .gt. real(j)*0.1d0) then
                  if (.not. eval_stat(j)) then
                     write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
                     eval_stat(j) = .true.
                  end if
               end if
            end do
                 
            read(18,*) pres_xx,pres_xy,pres_zx
            read(18,*) pres_xy,pres_yy,pres_yz
            read(18,*) pres_zx,pres_yz,pres_zz
            pres_tensor(1,1)=pres_xx
            pres_tensor(2,2)=pres_yy
            pres_tensor(3,3)=pres_zz
            pres_tensor(1,2)=pres_xy
            pres_tensor(2,3)=pres_yz
            pres_tensor(1,3)=pres_zx
            pres_tensor(3,1)=pres_tensor(1,3)
            pres_tensor(3,2)=pres_tensor(2,3)
            pres_tensor(2,1)=pres_tensor(1,2)
       !     pres_tensor=-pres_tensor
         end if 
         if (task_act .gt. frame_first) then
            pres_tensor_total=pres_tensor_total+pres_tensor
         end if    
      end do
     
!
!     Now calculate surface tension from formula
!
      pres_tensor_total=pres_tensor_total/(nframes-frame_first)
      tension=(pres_tensor_total(3,3)-0.5d0*(pres_tensor_total(1,1)+&
                  & pres_tensor_total(2,2)))*1d8*1602.1766d0/box_volume
      tension=0.5d0*zlen*1d-10*tension
      write(*,*) " completed!"
      write(*,*) "The surface tension is (N/m) or (J/m^2)",tension
      write(*,*) 
   end if

!
!     If the skip-xdat option is activated, jump directly to the end 
!
   if (skip_xdat) goto 34
end if


!
!     Calculates the diffusion coeffients for each element (self-diffusion)
!
if (calc_diff) then
   write(*,*) "Calculate the diffusion coefficients of all elements..."
   allocate(xyz2(3,natoms,nframes))
   allocate(pos_diff(natoms*3))
   allocate(msd_func(nframes-frame_first,nelems))
   allocate(times(nframes-frame_first))
   allocate(diff(nframes-frame_first))

!
!    Correct for box images: reintroduce the images for diffusion coefficients! 
!
   write(*,*) "Part1: Correct for box images ..."
   xyz2=xyz
   do i=1,nframes-1
      do j=1,natoms
!     Correct the x component
!
         do while ((xyz2(1,j,i+1) - xyz2(1,j,i)) .gt. xlen/2d0)
            xyz2(1,j,i+1)=xyz2(1,j,i+1)-xlen
         end do
         do while ((xyz2(1,j,i+1) - xyz2(1,j,i)) .lt. -xlen/2d0)
            xyz2(1,j,i+1)=xyz2(1,j,i+1)+xlen
         end do
!
!     Correct the y component
!
         do while ((xyz2(2,j,i+1) - xyz2(2,j,i)) .gt. ylen/2d0)
            xyz2(2,j,i+1)=xyz2(2,j,i+1)-ylen
         end do
         do while ((xyz2(2,j,i+1) - xyz2(2,j,i)) .lt. -ylen/2d0)
            xyz2(2,j,i+1)=xyz2(2,j,i+1)+ylen
         end do
!
!     Correct the z component
!
         do while ((xyz2(3,j,i+1) - xyz2(3,j,i)) .gt. zlen/2d0)
            xyz2(3,j,i+1)=xyz2(3,j,i+1)-zlen
         end do
         do while ((xyz2(3,j,i+1) - xyz2(3,j,i)) .lt. -zlen/2d0)
            xyz2(3,j,i+1)=xyz2(3,j,i+1)+zlen
         end do
      end do
   end do
   write(*,*) " ... done"

   write(*,*) "Part 2: Calculate mean square displacements..."
   if (nelems .eq. 1) then
      allocate(vector1(el_nums(1)*3))
   else if (nelems .eq. 2) then
      allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3))
   else if (nelems .eq. 3) then
      allocate(vector1(el_nums(1)*3),vector2(el_nums(2)*3),vector3(el_nums(3)*3))
   else 
      stop "Currently only up to three elements possible for diffusion coefficients!"
   end if        
   do i=1,nframes-frame_first
      times(i)=(i-1)*time_step
      do j=1,natoms
         do k=1,3
            pos_diff((j-1)*3+k)=xyz2(k,j,i+frame_first)-xyz2(k,j,1+frame_first)
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
   write(*,*) " ... done"

   open(unit=17,file="msd_plot.dat",status="replace")
   write(17,*) "# time (fs),   MSD (A^2)"
   do i=1,nframes-frame_first
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
!    Write the trajectory in xyz format (centered to unit cell) to file
!
   write(*,*) "Write trajectory centered to unit cell to 'traj_center.xyz' ..."
   open(unit=17,file="traj_center.xyz",status="replace")
   eval_stat=.false.
   do i=1,nframes
      do j=1,10
         if (real(i)/real(nframes) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do

      write(17,*) natoms
      write(17,*) "Frame No.",i
      do j=1,natoms
         write(17,*) at_names(j), xyz2(:,j,i)
      end do
   end do
   close(17)
   write(*,*) " completed!"
!
!     Calculate the diffusion coefficient by calculating the MSD 
!     based on the last time step!
!
   do k=1,nelems
      do i=1,nframes-frame_first
         diff(i)=msd_func(i,k)/(6.d0*times(i))
      end do
      diff = diff*(1E-10)**2/(1E-15)
      avg_diff=diff(nframes-frame_first)
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
   write(*,*) "Diffusion coefficient calculation finished!"
   write(*,*) "Coefficients written to diff_const_MSD_el*.dat"
end if        

!
!     Additionally, write POSCAR files for Core Level energy calculations 
!     of the chosen active atom species at different positions if wanted
!
if (cls_element .ne. "XX") then 
   write(*,*) "Write POSCARs for example calculatons..."
   atom_first=sum(el_nums(1:cls_elem_ind-1))+1
   atom_last=sum(el_nums(1:cls_elem_ind))
   call system("mkdir core_levels")
   slice_size=int((nframes-frame_first)/cls_rounds)
   rounds: do r=1,cls_rounds
      if (r .le. 9) then
          write(roundname,'(a,i1)') "round",r
      else if (r .le. 99) then
          write(roundname,'(a,i2)') "round",r
      else        
          write(roundname,'(a,i3)') "round",r
      end if    
      call system("mkdir core_levels/" //roundname) 

      open(unit=18,file="core_levels/" // trim(roundname) // "/active_examples.xyz",status="replace")
      open(unit=19,file="core_levels/" // trim(roundname) // "/active_examples.XDATCAR",status="replace")
!
!     Write header of XDATCAR file
!
      write(19,*) "Example positions of active atoms with different z-coordinates"
      write(19,*) factor
      write(19,*) xlen,0.0,0.0
      write(19,*) 0.0,ylen,0.0
      write(19,*) 0.0,0.0,zlen
      do i=1,nelems
         write(19,'(a,a)',advance="no")"  ", el_names(i)
      end do
      write(19,*)
      write(19,*) el_nums



      frame_round_first=frame_first+(r-1)*slice_size
      frame_round_last=frame_first+r*slice_size   
      write(*,*) "CLS round ",r,": From frame ", frame_round_first," to ",frame_round_last
      slices: do i=1,atom_slices
         if (i .le. 9) then
            write(filename,'(a,a,a,i1)') "core_levels/",trim(roundname),"/POSCAR_slice",i
         else if (i .le. 99) then
            write(filename,'(a,a,a,i2)') "core_levels/",trim(roundname),"/POSCAR_slice",i
         else 
            write(filename,'(a,a,a,i3)') "core_levels/",trim(roundname),"/POSCAR_slice",i
         end if

         frames: do j=frame_round_first,frame_round_last
            do k=atom_first,atom_last
               if ((xyz(3,k,j) .ge. zlo+(i-1)*slice_step) .and. (xyz(3,k,j) &
                           & .le. zlo+i*slice_step)) then
                  write(18,*) natoms
                  write(18,'(a,i4,a,i4,a,f13.6,a,f13.6)') "Layer ",i," of ", & 
                           & atom_slices,": z=",zlo+(i-1)*slice_step," to ", &
                           & zlo+i*slice_step
                  write(19,'(a,i4,a,i4,a,f13.6,a,f13.6)') "Direct Layer ",i," of ", &
                           & atom_slices,": z=",zlo+(i-1)*slice_step," to ", &
                           & zlo+i*slice_step
!
!     Write header of POSCAR file for Core Level shift
!
                  open(unit=20,file=filename)
                  write(20,*) "Active atom at z = ",xyz(3,k,j),", input for Core Level calculation."
                  write(20,*) factor
                  write(20,*) xlen,0.0,0.0
                  write(20,*) 0.0,ylen,0.0
                  write(20,*) 0.0,0.0,zlen
                  do l=1,nelems
                     write(20,'(a,a)',advance="no")"  ", el_names(l)
                  end do
                  write(20,*) "   ",el_names(cls_elem_ind)
                  do l=1,nelems-1
                     write(20,'(i10)',advance="no") el_nums(l)
                  end do
                  write(20,'(i10)',advance="no") el_nums(nelems)-1
                  write(20,*) 1
                  write(20,*) "Direct"
                  do l=1,natoms
                     write(18,*) at_names(l), xyz(:,l,j)
                     write(19,*) xyz(1,l,j)/xlen,xyz(2,l,j)/ylen,xyz(3,l,j)/zlen
                     if (l .ne. k) then
                        write(20,*) xyz(1,l,j)/xlen,xyz(2,l,j)/ylen,xyz(3,l,j)/zlen
                     end if
                  end do 
                  write(20,*) xyz(1,k,j)/xlen,xyz(2,k,j)/ylen,xyz(3,k,j)/zlen
                  close(20)   
                  cycle slices
               end if            
            end do
         end do frames
         write(*,*) "Warning: No Active atom found in slice No.",i 
      end do slices
      close(18)
      close(19)
   end do rounds
   write(*,*) " completed!"
   write(*,*) 
end if
write(*,*)
!if (write_traj) then
!   write(*,*) "File 'trajectory.xyz' with xyz trajectory written."
!end if
write(*,*) "File 'dens_elems' with element densities in z-direction written."
if (cls_element .ne. "XX") then
   write(*,*) "POSCAR files for Core level shifts written to folder 'core_levels/round...'."
   write(*,*) "Files 'active_examples.xyz/.XDATCAR' with active atom positions written "
   write(*,*) "  to core_levels/round...."
end if
34 continue
write(*,*)
write(*,*) "analyze_scalms ended normally..."
write(*,*)
end program analyze_scalms


subroutine replace_text (s,text,rep,outs)
CHARACTER(*)        :: s,text,rep
CHARACTER(LEN(s)+100) :: outs     ! provide outs with extra 100 char len
INTEGER             :: i, nt, nr

outs = s ; nt = LEN_TRIM(text) ; nr = LEN_TRIM(rep)
DO
   i = INDEX(outs,text(:nt)) ; IF (i == 0) EXIT
   outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
END DO

return 
end subroutine replace_text
