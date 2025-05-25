!
!    eval_stm: Plots time-dependent radial distribution functions from
!      MD trajectories (XDATCAR files). Time-dependent processes like 
!      intermetallic phase formations for other phase transitions can
!      be observed in real time.
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!


program td_rdf
implicit none 
integer::i,j,k,m
integer::xdat_lines,natoms,nframes,nelems,average
integer::readstat,nframes2,eval_ind
integer::nbins
integer::shift,endl,ig,npart
integer,allocatable::el_nums(:),el_inds(:)
character(len=32)::arg
character(len=100)::line
character(len=120)::a120
character(len=220)::a220
character(len=1)::atest
character(len=2)::el1,el2,el3
real(kind=8)::rdf_binsize,r_cut
real(kind=8)::dist,pi
real(kind=8),allocatable::rdf_plot(:,:)
real(kind=8)::diff_vec(3),pos1(3),pos2(3),act_num(3)
real(kind=8)::a_vec(3),b_vec(3),c_vec(3)
real(kind=8)::a_vec_avg(3),b_vec_avg(3),c_vec_avg(3)
real(kind=8),allocatable::a_vecs(:,:),b_vecs(:,:),c_vecs(:,:)
real(kind=8),allocatable::xyz(:,:,:),a_read(:),b_read(:),c_read(:)
real(kind=8)::volume
real(kind=8)::vb,rho,nid,ngr,r_act
logical::npt_traj

pi=3.141592653589793238

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
write(*,*) "PROGRAM td_rdf: calculation of time-dependent radial-"
write(*,*) " distributio functions from XDATCAR trajectories."
write(*,*) " For this, the given trajectory is subdivided into a "
write(*,*) " number of parts, for each those, a seprate RDF plot is"
write(*,*) " generated. Currently, this is done for one element "
write(*,*) " (index given by user)"
write(*,*) " The resulting time-dependent RDF is written to a 3D plot file."
write(*,*) " Only a XDATCAR file is needed for the analysis!"
write(*,*)
write(*,*) "The following command line arguments can be given::"
write(*,*) " -overview:  print an overview of all scripts and programs in utils4VASP"
write(*,*) "  -elem=[number]: Index of analyzed element (order in POSCAR)"
write(*,*) "  -avg=[number]: No. of MD frames to average for analysis (def. 10)"
write(*,*) "  -rcut=[value]: Cutoff distance for RDF calculation (Angs.) (def. 10.0)"
write(*,*) "  -nbins=[number]: Number of bins along distance (def. 500)"
write(*,*) "  -npt: If a NpT trajectory with variable volume was calculated."
write(*,*)

!
!    Use Command line arguments for specification of analysis job
!

npt_traj=.false.
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:4))  .eq. "-npt") then
      write(*,*) "A NpT trajectory will be evaluated!"
      npt_traj=.true.
   end if
end do


eval_ind=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-elem=") then
      read(arg(7:32),*) eval_ind
      write(*,*) "Index of element that will be evaluated:",eval_ind
   end if
end do
if (eval_ind .lt. 1) then
   write(*,*) "Please give a valid element index with -elem=[number]!"
   stop
end if

average=10
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:5))  .eq. "-avg=") then
      read(arg(6:32),*) average
   end if
end do
write(*,*) "MD steps to be averaged for histogram:",average

r_cut=10.0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-rcut=") then
      read(arg(7:32),*) r_cut
   end if
end do
write(*,*) "Cutoff for the RDF calculation:",r_cut

nbins=500
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-nbins=") then
      read(arg(8:32),*) nbins
   end if
end do
write(*,*) "Number of bins along distance for RDF:",nbins

rdf_binsize=real(r_cut/nbins)
!
!    First, determine the number of lines in the XDATCAR file
!
call system("wc -l XDATCAR > xdat_length")
open(unit=45,file="xdat_length",status="old")
read(45,*) xdat_lines
close(45)

allocate(a_read(3),b_read(3),c_read(3))
open(unit=15,file="XDATCAR",status="old")
read(15,*)
read(15,*)
read(15,*) a_read(:)
read(15,*) b_read(:)
read(15,*) c_read(:)
read(15,'(a)') line
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
read(15,*) el_nums

natoms=sum(el_nums)
allocate(el_inds(natoms))
el_inds(1:el_nums(1))=1
if (nelems .ge. 2) then
   el_inds(el_nums(1):el_nums(1)+el_nums(2))=2
end if
if (nelems .ge. 3) then
   el_inds(el_nums(1)+el_nums(2):)=3
end if

write(*,*)
write(*,*) "System setup:"
write(*,*) "Number of atoms:",natoms
write(*,*) "Number of elements:",nelems
write(*,*) "List of elements:"
write(*,*) "  1. ",el1,": ",el_nums(1)," atoms"
if (nelems .ge. 2) then
   write(*,*) "  2. ",el2,": ",el_nums(2)," atoms"
end if
if (nelems .eq. 3) then
   write(*,*) "  3. ",el3,": ",el_nums(3)," atoms"
end if

!
!    For NVT trajectories: Each frame has only one header line 
!
write(*,*)
write(*,*) "Read in trajectory from XDATCAR ..."
if (.not. npt_traj) then
   a_vec(:) = a_read(:)
   b_vec(:) = b_read(:)
   c_vec(:) = c_read(:)

   nframes=int((xdat_lines-7)/(natoms+1))

   write(*,*) "Number of frames:",nframes

   allocate(xyz(3,natoms,nframes))
   do i=1,nframes
      read(15,*)
      do j=1,natoms
         read(15,'(a)') a120
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

         xyz(:,j,i)=act_num(:)
!
!     Remove image flags!
!
         do k=1,3
            do while(xyz(k,j,i) .gt. 1.d0)
               xyz(k,j,i) = xyz(k,j,i) -1.d0
            end do
            do while(xyz(k,j,i) .lt. -1.d0)
               xyz(k,j,i) = xyz(k,j,i) +1.d0
            end do
         end do

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

   allocate(a_vecs(3,nframes),b_vecs(3,nframes),c_vecs(3,nframes))
   a_vecs(:,1)=a_read(:)
   b_vecs(:,1)=b_read(:)
   c_vecs(:,1)=c_read(:)

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
         a_vecs(:,i)=a_read(:)
         b_vecs(:,i)=b_read(:)
         c_vecs(:,i)=c_read(:)
      else
         read(15,*)
      end if
      do j=1,natoms
         read(15,'(a)') a120
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

         xyz(:,j,i)=act_num(:)
!
!     Remove image flags!
!
         do k=1,3
            do while(xyz(k,j,i) .gt. 1.d0)
               xyz(k,j,i) = xyz(k,j,i) -1.d0
            end do
            do while(xyz(k,j,i) .lt. -1.d0)
               xyz(k,j,i) = xyz(k,j,i) +1.d0
            end do
         end do
      end do
   end do
   close(15)
end if
!
!    Convert to usual cartesian coordinates 
!
!if (.not. npt_traj) then
!   do i=1,nframes
!      do j=1,natoms
!         xyz(1,j,i) = xyz(1,j,i)*a_len
!         xyz(2,j,i) = xyz(2,j,i)*b_len
!         xyz(3,j,i) = xyz(3,j,i)*c_len
!      end do
!   end do
!else
!   do i=1,nframes
!      do j=1,natoms
!         xyz(1,j,i) = xyz(1,j,i)*a_lens(i)
!         xyz(2,j,i) = xyz(2,j,i)*b_lens(i)
!         xyz(3,j,i) = xyz(3,j,i)*c_lens(i)
!      end do
!   end do
!end if
write(*,*) "... done!"

!
!    Number of trajectory parts
!
nframes2=int(nframes/average)

!
!    Calculate the radial distribution function for each trajectory part
!
write(*,*)
write(*,*) "Calculate RDF function ..."

if (eval_ind .gt. 0) then
   shift=sum(el_nums(1:eval_ind-1))
else 
   shift=0
end if
allocate(rdf_plot(nbins,nframes2))
rdf_plot=0.d0

write(*,*) "Number of parts/RDF plots: ",average
!
!    For NVT trajectories
!
if (.not. npt_traj) then
   do m=1,nframes2
      write(*,*) "Parts evaluated: ",m," of ",nframes2
      do i=(m-1)*average+1,m*average
         do j=1,el_nums(eval_ind)-1  ! Ni atoms
            do k=j+1,el_nums(eval_ind)  ! Ga atoms
               pos1 = xyz(:,j+shift,i)
               pos2 = xyz(:,k+shift,i)


               diff_vec=pos1-pos2
!
!     Correct the x component
!
               do while (abs(diff_vec(1)) .gt. 0.5d0)
                  diff_vec(1)=diff_vec(1)-sign(1.0d0,diff_vec(1))
               end do

!
!     Correct the y component
!
               do while (abs(diff_vec(2)) .gt. 0.5d0)
                  diff_vec(2)=diff_vec(2)-sign(1.0d0,diff_vec(2))
               end do

!
!     Correct the z component
!
               do while (abs(diff_vec(3)) .gt. 0.5d0)
                  diff_vec(3)=diff_vec(3)-sign(1.0d0,diff_vec(3))
               end do

               diff_vec(:)=diff_vec(1)*a_vec(:)+diff_vec(2)*b_vec(:)+diff_vec(3)*c_vec(:)

               dist = sqrt((diff_vec(1))**2 + &
                     & (diff_vec(2))**2 + (diff_vec(3))**2)

!  Remainder of the calculation taken from Frenkel Smit, page 86
               ig=int(dist/rdf_binsize)
               if (ig .le. nbins) then
                   rdf_plot(ig,m)=rdf_plot(ig,m)+2.d0
               end if        
            end do
         end do
      end do
      do j=1,nbins
         ngr=real(average)
         npart=el_nums(eval_ind)
         r_act=rdf_binsize*(real(j)+0.5d0)
         vb=((real(j)+1.0)**3-real(j)**3)*rdf_binsize**3
         call trip_prod(a_vec(:),b_vec(:),c_vec(:),volume)
         rho=npart/(abs(volume))

         nid=4.d0/3.d0*pi*vb*rho
         rdf_plot(j,m)=rdf_plot(j,m)/(ngr*npart*nid)
      end do
   end do

!
!    For NpT trajectories
!
else 
   do m=1,nframes2
  !    write(*,*) "Parts evaluated: ",m," of ",nframes2
  !    write(*,*) "hhh",(m-1)*average+1,m*average
      a_vec_avg=0.d0
      b_vec_avg=0.d0
      c_vec_avg=0.d0
      do i=(m-1)*average+1,m*average
!
!     Average the sizes of the unit cell
!
         a_vec_avg(:)=a_vec_avg(:)+a_vecs(:,i)
         b_vec_avg(:)=b_vec_avg(:)+b_vecs(:,i)
         c_vec_avg(:)=c_vec_avg(:)+c_vecs(:,i)

         do j=1,el_nums(eval_ind)-1  ! Ni atoms
            do k=j+1,el_nums(eval_ind)  ! Ga atoms
               pos1 = xyz(:,j+shift,i)
               pos2 = xyz(:,k+shift,i)
               

               diff_vec=pos1-pos2
!
!     Correct the x component
!
               do while (abs(diff_vec(1)) .gt. 0.5d0)
                  diff_vec(1)=diff_vec(1)-sign(1.0d0,diff_vec(1))
               end do

!
!     Correct the y component
!
               do while (abs(diff_vec(2)) .gt. 0.5d0)
                  diff_vec(2)=diff_vec(2)-sign(1.0d0,diff_vec(2))
               end do

!
!     Correct the z component
!
               do while (abs(diff_vec(3)) .gt. 0.5d0)
                  diff_vec(3)=diff_vec(3)-sign(1.0d0,diff_vec(3))
               end do

               diff_vec(:)=diff_vec(1)*a_vecs(:,i)+diff_vec(2)*b_vecs(:,i)+&
                             & diff_vec(3)*c_vecs(:,i)
               dist = sqrt((diff_vec(1))**2 + &
                     & (diff_vec(2))**2 + (diff_vec(3))**2)
!  Remainder of the calculation taken from Frenkel Smit, page 86
               ig=int(dist/rdf_binsize)
               if ((ig .le. nbins) .and. (ig .ge. 1)) then
                   rdf_plot(ig,m)=rdf_plot(ig,m)+2.d0
               end if
            end do
         end do
      end do
      a_vec_avg(:)=a_vec_avg(:)/real(average)
      b_vec_avg(:)=b_vec_avg(:)/real(average)
      c_vec_avg(:)=c_vec_avg(:)/real(average)
      do j=1,nbins
         ngr=real(average)
         npart=el_nums(eval_ind)
         r_act=rdf_binsize*(real(j)+0.5d0)
         vb=((real(j)+1.0)**3-real(j)**3)*rdf_binsize**3
         call trip_prod(a_vec_avg(:),b_vec_avg(:),c_vec_avg(:),volume)
         rho=npart/(abs(volume))
         nid=4.d0/3.d0*pi*vb*rho
  !       write(*,*) "plot",rdf_plot(j,m),ngr,npart,nid,vb,rho    
         rdf_plot(j,m)=rdf_plot(j,m)/(ngr*npart*nid)
      end do
   end do
end if  
write(*,*) "... done!"

!
!    Print time-dependent RDF plot to file 
!
write(*,*)
write(*,*) "Print RDF to file ..."

open(unit=13,file="td_rdf.dat",status="replace")
write(13,*) "# time (fs)    dist(Angs)     RDF "
do i=1,nframes2
   do j=1,nbins
      write(13,*) i,(real(j)+0.5)*rdf_binsize,rdf_plot(j,i)
   end do
   write(13,*)
end do
write(*,*) "... done!"
write(*,*)
write(*,*) "File 'td_rdf.dat' with RDF frames has been written."
write(*,*) "td_rdf exited normally."
write(*,*)


end program td_rdf


!
!    Calculation of triple product of vectors (needed for volumes for 
!       RDF renormalizations
!
subroutine trip_prod(vec1,vec2,vec3,triple)
implicit none 
real(kind=8)::vec1(3),vec2(3),vec3(3)
real(kind=8)::triple

triple=((vec1(2)*vec2(3)-vec1(3)*vec2(2))*vec3(1))- &
      &  ((vec1(1)*vec2(3)-vec1(3)*vec2(1))*vec3(2))+ &
      &  ((vec1(1)*vec2(2)-vec1(2)*vec2(1))*vec3(3))

return
end subroutine trip_prod        
