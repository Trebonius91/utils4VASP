!
!    partial_dos: evaluates DOSCAR files and prints out the DOS 
!      plots for one or several elements, atom indices or orbitals
!    Part of VASP4CLINT
!     Julien Steffen, 2023 (julien.steffen@fau.de)
!
program partial_dos
implicit none 
integer::i,j,k
integer::readstat,counter,ispin
integer::nelems,natoms,natoms_test,nelems_choice
character(len=120)::cdum,adum
character(len=2),allocatable::el_list(:),el_list_read(:)
character(len=5),allocatable::orb_list_read(:),orb_list(:)
character(len=100)::arg
character(len=2),allocatable::el_names_read(:),el_names(:),at_names(:)
integer,allocatable::el_nums(:)
real(kind=8)::enmax,enmin,efermi,fdum
real(kind=8),allocatable::e_dos(:),dos_tot(:),dos_tot_up(:),dos_tot_down(:)
real(kind=8),allocatable::dos_partial(:,:,:)
real(kind=8),allocatable::dos_part_out(:),dos_part_out_up(:),dos_part_out_down(:)
integer::npoints,ninds_choice,norb_choice
integer,allocatable::ind_list(:),ind_list_read(:)
logical::dos_elements,dos_indices,act_read,el_all,orb_all,use_orbital
logical,allocatable::orb_used(:)

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
write(*,*) "PROGRAM partial_dos: Evaluation of VASP density of states (DOS)"
write(*,*) " calculations. Element- or atom-resolved DOS can be printed out,"
write(*,*) " further, one or several certain orbitals can be chosen."
write(*,*) "The POSCAR and DOSCAR files of a VASP calculation need to be "
write(*,*) " present, the calculation must be done with the LORBIT=11 command."
write(*,*) "The following command line arguments can/must be given:"
write(*,*) " -overview:  print an overview of all scripts and programs in utils4VASP"
write(*,*) " -element=[list of elements] : which element shall be analyzed."
write(*,*) "    one or several can be given. State 'all' if all atoms shall be "
write(*,*) "    incorporated into the DOS evaluation. Examples: -element=Pt or "
write(*,*) "    -element=Pt,Ga,In or -element=all" 
write(*,*) " -index=[one or more indices] : If certain atoms shall be analyzed."
write(*,*) "    One or more indices can be given. Examples: -index=103 or "
write(*,*) "    -index=13,103,201"
write(*,*) " -orbital=[name(s)] : Which orbitals shall be analyzed for the given "
write(*,*) "    elements or atom indices. Either 'all' for all orbitals or "
write(*,*) "    one or several descriptors need to be given, possible are:"
write(*,*) "    s,px,py,pz,dxy,dyz,dz2,dxz,dx2y2. As shortcuts can further be "
write(*,*) "    used: p (all p-orbitals), d (all d-orbitals).  "
write(*,*) "    Examples: -orbital=all or -orbital=s or -orbital=p,dxy,dxz"
write(*,*) "If both -element=all and -orbital=all are given, only the global DOS"
write(*,*) " will be written, since then the whole wavefunction is chosen."
write(*,*) "The Fermi energy will always be subtracted to build the x-axis. If"
write(*,*) " you want to plot with without the correction, use the Fermi energy"
write(*,*) " printed in the headers of the files."
write(*,*)
!
!    Read in and process the command line arguments
!

!
!    The elements whose DOS shall be calculated
!
dos_elements=.false.
el_all=.false.
allocate(el_list_read(10))
el_list_read="XX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-element=") then
      read(arg(10:),*,iostat=readstat) el_list_read
      dos_elements=.true.
   end if
end do

if (dos_elements) then
   nelems_choice=0
   do i=1,10
      if (el_list_read(i) .eq. "XX") exit
!
!    If all elements are chosen (all)
!
      if (el_list_read(i) .eq. "al") then
         el_all = .true.
         exit
      end if
      nelems_choice=nelems_choice+1
   end do

   allocate(el_list(nelems_choice))
   do i=1,nelems_choice
      el_list(i)=el_list_read(i)
   end do
   if (el_all) then
      write(*,*) "All atoms in the system were chosen for the DOS evalulation."
   else 
      write(*,*) "The following elements were chosen for the DOS evaluation:"
      do i=1,nelems_choice
         write(*,*) " - ",el_list(i)
      end do
   end if
end if


!
!     The atom indices whose DOS shall be calculated
!

dos_indices=.false.
allocate(ind_list_read(1000))
ind_list_read=0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:7))  .eq. "-index=") then
      read(arg(8:),*,iostat=readstat) ind_list_read
      if (dos_elements) then
         write(*,*) "Both -element=... and -index=... flags were given."
         write(*,*) " Please use only one of them at a time!"
         stop
      end if
      dos_indices=.true.
   end if
end do

if (dos_indices) then
   ninds_choice=0

   do i=1,1000
      if (ind_list_read(i) .eq. 0) exit
      ninds_choice=ninds_choice+1
   end do
 
   allocate(ind_list(ninds_choice))
   do i=1,ninds_choice
      ind_list(i)=ind_list_read(i)
   end do   

   write(*,*) "The following atoms were chosen for the DOS evaluation:"
   do i=1,ninds_choice
      write(*,*) " - No. ",ind_list(i)
   end do
 
end if

if ((.not. dos_indices) .and. (.not. dos_elements)) then
   write(*,*) "Please give either a -element=... or -index=... selection!"
   stop
end if
   
!
!    The orbital angular momentum specifiers 
!
orb_all=.false.
allocate(orb_list_read(10))
orb_list_read="XXXXX"
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:9))  .eq. "-orbital=") then
      read(arg(10:),*,iostat=readstat) orb_list_read
      use_orbital=.true.
   end if
end do
if (.not. use_orbital) then
   write(*,*) "Please give a -orbital=... selection!"
   stop
end if

norb_choice=0
do i=1,10
   if (orb_list_read(i) .eq. "XXXXX") exit
!
!    If all elements are chosen (all)
!
   if (orb_list_read(i) .eq. "all") then
      orb_all = .true.
      exit
   end if
   norb_choice=norb_choice+1
end do

allocate(orb_list(norb_choice))
do i=1,norb_choice
   orb_list(i)=orb_list_read(i)
end do

allocate(orb_used(9)) 
orb_used = .false.
if (orb_all) orb_used = .true.
do i=1,norb_choice
   if (orb_list(i) .eq. "s") orb_used(1) = .true.
   if (orb_list(i) .eq. "p") orb_used(2:4) = .true.
   if (orb_list(i) .eq. "d") orb_used(5:9) = .true.
   if (orb_list(i) .eq. "px") orb_used(2) = .true.
   if (orb_list(i) .eq. "py") orb_used(3) = .true.
   if (orb_list(i) .eq. "pz") orb_used(4) = .true.
   if (orb_list(i) .eq. "dxy") orb_used(5) = .true.
   if (orb_list(i) .eq. "dyz") orb_used(6) = .true.
   if (orb_list(i) .eq. "dz2") orb_used(7) = .true.
   if (orb_list(i) .eq. "dxz") orb_used(8) = .true.
   if (orb_list(i) .eq. "dx2y2") orb_used(9) = .true. 
end do
write(*,*) "The following orbitals were chosen for the DOS evaluation:"
write(*,*) "  (T: used, F: not used)"
write(*,*) " - s     : ",orb_used(1)
write(*,*) " - px    : ",orb_used(2)
write(*,*) " - py    : ",orb_used(3)
write(*,*) " - pz    : ",orb_used(4)
write(*,*) " - dxy   : ",orb_used(5)
write(*,*) " - dyz   : ",orb_used(6)
write(*,*) " - dz2   : ",orb_used(7)
write(*,*) " - dxz   : ",orb_used(8)
write(*,*) " - dx2y2 : ",orb_used(9)
write(*,*)


!
!    First, open the POSCAR file for elements and atom indices
!

open(unit=28,file="POSCAR",status="old",iostat=readstat)
if (readstat .ne. 0) then
   stop "The POSCAR file is not there!"
end if

read(28,*)
read(28,*)
read(28,*)
read(28,*)
read(28,*)
!
!    Determine the element symbols and their numbers
!
allocate(el_names_read(10))
el_names_read="XX"
read(28,'(a)') cdum

read(cdum,*,iostat=readstat) el_names_read
nelems=0

do i=1,10
   if (el_names_read(i) .eq. "XX") exit
   nelems=nelems+1
end do
!
!    Define permanent element symbol and number arrays and the total
!    number of atoms in the system
!
allocate(el_names(nelems),el_nums(nelems))
do i=1,nelems
   el_names(i)=el_names_read(i)
end do
read(28,*) el_nums

natoms = sum(el_nums)

!
!    Define individual elements for each atom
!

allocate(at_names(natoms))
counter = 1
do i=1,nelems
   do j=1,el_nums(i)
      at_names(counter) = el_names(i)
      counter = counter +1
   end do
end do

close(28)

!
!    Second, read in the DOSCAR completely
!

open(unit=29,file="DOSCAR",status="old",iostat=readstat)
if (readstat .ne. 0) then
   stop "The DOSCAR file is not there!"
end if
read(29,*) natoms_test
if (natoms_test .eq. natoms) then
   write(*,*) "Check compatibility of POSCAR and DOSCAR ..."
   write(*,*) " ...  success!"
else 
   write(*,*) "Check compatibility of POSCAR and DOSCAR ..."
   write(*,*) "POSCAR has ", natoms, " atoms, DOSCAR has ",natoms_test," atoms!"
   stop "Check POSCAR and DOSCAR!" 
end if
read(29,*) 
read(29,*) 
read(29,*) 
read(29,*) 
read(29,*) enmax,enmin,npoints,efermi
write(*,*) "General information:"
write(*,'(a,f16.8,a)') "  - Minimum DOS energy: ",enmin," eV"
write(*,'(a,f16.8,a)') "  - Maximum DOS energy: ",enmax," eV"
write(*,'(a,f16.8,a)') "  - Fermi energy: ",efermi, " eV"
write(*,*)
allocate(e_dos(npoints))

!
!     Determine if the calculation was spin-polarized or not and read in the 
!     total DOS (always written at the beginning of DOSCAR)
!
write(*,*) "Read in the total DOS ..."
ispin=1
do i=1,npoints
   if (i .eq. 1) then
      read(29,'(a)') adum
      read(adum,*,iostat=readstat) fdum,fdum,fdum,fdum,fdum
      if (readstat .eq. 0) then
         write(*,*) "You have done a spin-polarized calculation (ISPIN=2)!"
         write(*,*) "Two columns will be written for the DOS: up and down."
         ispin=2
         allocate(dos_tot_up(npoints),dos_tot_down(npoints))
         read(adum,*) e_dos(1),dos_tot_up(1),fdum,dos_tot_down(1) 
      else 
         write(*,*) "You have done no spin-polarized calculation (ISPIN=1)!"
         write(*,*) "The DOS will be written into one column."
         ispin=1
         allocate(dos_tot(npoints))
         read(adum,*) e_dos(1),dos_tot(1)
      end if
      cycle
   end if
   if (ispin .eq. 2) then
      read(29,*) e_dos(i),dos_tot_up(i),fdum,dos_tot_down(i)
   else if (ispin .eq. 1) then
      read(29,*) e_dos(i),dos_tot(i)
   end if
end do
write(*,*) " ...  success!"
!
!     Write the total DOS to a file 
!
open(unit=30,file="dos_total.dat",status="replace")
write(30,*) "# DOS calculated by partial_dos "
write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"

if (ispin .eq. 1) then
   write(30,*) "#   energy(eV)             DOS"
else if (ispin .eq. 2) then 
   write(30,*) "#   energy(eV)             DOS(up)            DOS(down)"
end if
do i=1,npoints
   if (ispin .eq. 1) then
      write(30,*) e_dos(i)-efermi,dos_tot(i)
   else if (ispin .eq. 2) then
      write(30,*) e_dos(i)-efermi,dos_tot_up(i),dos_tot_down(i)
   end if
end do
close(30)
write(*,*) "Total DOS of the system written to file 'dos_total.dat'."

!
!    If spin polarization is present, write an averaged DOS 
!     (e.g., for effective calculation of DOS overlaps)
!
if (ispin .eq. 2) then
   open(unit=30,file="dos_tot_sum.dat",status="replace")
   write(30,*) "# DOS calculated by partial_dos (alpha and beta summed up) "
   write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"
   write(30,*) "#   energy(eV)             DOS"
   do i=1,npoints
      write(30,*) e_dos(i)-efermi,dos_tot_up(i)+dos_tot_down(i)
   end do 
   close(30)
   write(*,*) "Summed total DOS (alpha+beta) written to file 'dos_tot_sum.dat'."
end if
        
if (el_all .and. orb_all) then
   write(*,*)
   write(*,*) "All atoms and all orbitals are chosen, the partial DOS will be skipped."
   goto 22
end if 
!
!     Now read the individual DOS of all atoms and their respective orbitals
!
write(*,*)
write(*,*) "Read in the partial DOS ..."
if (ispin .eq. 1) then
   allocate(dos_partial(natoms,npoints,9))
else if (ispin .eq. 2) then
   allocate(dos_partial(natoms,npoints,18))
end if
do i=1,natoms
   read(29,*,iostat=readstat) fdum
   if (readstat .ne. 0) then
      stop "Missing data in DOSCAR? Did you add the LORBIT=11 keyword?"
   end if
   do j=1,npoints
      read(29,*,iostat=readstat) fdum,dos_partial(i,j,:)
      dos_partial(i,j,:)=abs(dos_partial(i,j,:))
      if (readstat .ne. 0) then
         stop "Missing data in DOSCAR? Did you add the LORBIT=11 keyword?"
      end if
   end do
end do
close(29)
write(*,*) " ...  success!"

!
!     Now select the DOS data according to element,index and orbital
!
write(*,*)
write(*,*) "Sum up partial DOS according to your element/index/orbital choice."
if (ispin .eq. 1) then
   allocate(dos_part_out(npoints))
   dos_part_out=0.d0
else if (ispin .eq. 2) then
   allocate(dos_part_out_up(npoints),dos_part_out_down(npoints))
   dos_part_out_up=0.d0
   dos_part_out_down=0.d0
end if

do i=1,natoms
   act_read=.false.
   if (dos_indices) then
      do j=1,ninds_choice
         if (i .eq. ind_list(j)) then
            act_read=.true.
         end if
      end do
   else if (dos_elements) then
      if (el_all) act_read=.true.
      do j=1,nelems_choice
         if (at_names(i) .eq. el_list(j)) then
            act_read=.true.
         end if
      end do
   end if
   if (act_read) then
      do j=1,npoints
         if (ispin .eq. 1) then
            if (orb_all) then
               dos_part_out(j)=dos_part_out(j)+sum(dos_partial(i,j,:))
            else 
               do k=1,9
                  if (orb_used(k)) then 
                     dos_part_out(j)=dos_part_out(j)+dos_partial(i,j,k)
                  end if
               end do
            end if
         else if (ispin .eq. 2) then
            if (orb_all) then
               do k=1,9     
                  dos_part_out_up(j)=dos_part_out_up(j)+dos_partial(i,j,k*2-1)
                  dos_part_out_down(j)=dos_part_out_down(j)+dos_partial(i,j,k*2)
               end do 
            else
!
!     The upspin component
!
               do k=1,9
                  if (orb_used(k)) then 
                     dos_part_out_up(j)=dos_part_out_up(j)+dos_partial(i,j,k*2-1)
                  end if
               end do
!
!     The downspin component
!
               do k=1,9
                  if (orb_used(k)) then
                     dos_part_out_down(j)=dos_part_out_down(j)+dos_partial(i,j,k*2)
                  end if
               end do
            end if
         end if
      end do
   end if
end do

write(*,*) " ...  success!"
if (ispin .eq. 1) then
   if (sum(dos_part_out) .lt. 0.0001d0) then
      write(*,*)
      write(*,*) "The partial DOS is zero! Maybe no valid element or index chosen?"
      write(*,*)
   end if
else if (ispin .eq. 2) then
   if (sum(dos_part_out_up) .lt. 0.0001d0) then
      write(*,*)
      write(*,*) "The partial DOS is zero! Maybe no valid element or index chosen?"
      write(*,*)
   end if
end if
!
!     Write the partial DOS to a file 
!
open(unit=30,file="dos_partial.dat",status="replace")
write(30,*) "# DOS calculated by partial_dos "
write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"

if (ispin .eq. 1) then
   write(30,*) "#   energy(eV)             DOS  "
else if (ispin .eq. 2) then
   write(30,*) "#   energy(eV)             DOS(up)            DOS(down) "
end if
do i=1,npoints
   if (ispin .eq. 1) then
      write(30,*) e_dos(i)-efermi,dos_part_out(i)
   else if (ispin .eq. 2) then
      write(30,*) e_dos(i)-efermi,dos_part_out_up(i),dos_part_out_down(i)
   end if
end do
close(30)
write(*,*) "Partial DOS of the system written to file 'dos_partial.dat'."

!
!    If spin polarization is present, write an averaged DOS 
!     (e.g., for effective calculation of DOS overlaps)
!
if (ispin .eq. 2) then
   open(unit=30,file="dos_part_sum.dat",status="replace")
   write(30,*) "# DOS calculated by partial_dos (alpha and beta summed up) "
   write(30,*) "# Fermi energy has been subtracted (E_fermi = ",efermi," eV)"
   write(30,*) "#   energy(eV)             DOS"
   do i=1,npoints
      write(30,*) e_dos(i)-efermi,dos_part_out_up(i)+dos_part_out_down(i)
   end do
   close(30)
   write(*,*) "Summed partial DOS (alpha+beta) written to file 'dos_part_sum.dat'."
end if 

22 continue
write(*,*)
write(*,*) "partial_dos terminated normally."
write(*,*)

end program partial_dos

