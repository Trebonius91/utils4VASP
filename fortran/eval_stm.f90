!
!    eval_stm: evaluates the results of a STM calculation 
!      (partial charge density, to be evaluated by the Tersoff-Hamann
!      approach). A PARCHG file is read in and depending on the 
!      user keywords, a constant current or constant height STM
!      picture is evaluated. Both a bitmap file and a plottable 
!      data file with a gnuplot file to plot it are given as 
!      results.
!    Part of VASP4CLINT
!     Julien Steffen, 2024 (julien.steffen@fau.de)
!

program eval_stm
implicit none 
integer::readstat
integer::i,j,k,l,m,n
integer::k_act,l_act,m_act
integer::natoms,nelems
integer::grida,gridb,gridc
integer::nx,ny,nz
integer::stm_gridx,stm_gridy
integer::near_x,near_y,near_z
integer::repeat_x,repeat_y
integer::include_x,include_y,include_z
real(kind=8)::scale
real(kind=8)::a_fac,cutoff,g_width
real(kind=8)::x_len,y_len,z_len
real(kind=8)::stm_height,isos_dens
real(kind=8)::x_plot_min,x_plot_max
real(kind=8)::y_plot_min,y_plot_max
real(kind=8)::dist_act,xtoy
real(kind=8)::weight_sum,weight
real(kind=8)::prev_z,prev_z_init
real(kind=8)::def_z_shift,z_step
real(kind=8)::stm_act
real(kind=8)::coord_mat(3,3)
real(kind=8)::act_pos(3)
real(kind=8)::vec_act(3)
real(kind=8)::diff_vec(3)
real(kind=8),allocatable::xyz(:,:)
real(kind=8),allocatable::cdens(:,:,:)
real(kind=8),allocatable::c_coord(:,:,:,:)
real(kind=8),allocatable::stm_dat(:,:)
real(kind=8),allocatable::stm_pos(:,:,:)
logical::eval_stat(10)
character(len=120)::cdum,arg
character(len=30)::stm_mode
character(len=2),allocatable::el_names_read(:),el_names(:)
character(len=2),allocatable::at_names(:)
integer,allocatable::el_nums(:)

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
write(*,*) "PROGRAM eval_stm: Plots a STM image from a VASP"
write(*,*) " PARCHG file based on the Tersoff-Hamann approach."
write(*,*) "Constant height and constant current images are possible."
write(*,*) "The only needed input is a PARCHG file for a certain"
write(*,*) " energy range below the Fermi level, corresponding to an"
write(*,*) " experimental STM tunneling voltage."
write(*,*) "The following command line arguments can be given:"
write(*,*) "  -overview:  print an overview of all scripts and programs "
write(*,*) "    in utils4VASP"
write(*,*) "  -mode=[height or current] : Determines if the constant"
write(*,*) "    height or constant current STM mode shall be used."
write(*,*) "  -pos=[value] : For constant height STMs: the position of"
write(*,*) "    imaginary tip along the z-axis of the unit cell."
write(*,*) "    The value needs to be between 0 and the height of the"
write(*,*) "    unit cell! (in Angstroms)"
write(*,*) "  -dens=[value] : For constant current STMs : The charge "
write(*,*) "    density at which the tip shall be located at each grid"
write(*,*) "    point along x and y, resulting in a z-value (height)."
write(*,*) "    Reasonable values are between 0.01 and 0.1"
write(*,*) "  -repeat_x=[number] : For the generated picture: the number"
write(*,*) "    of unit cell repetitions along x-axis to generate a larger  "
write(*,*) "    picture which shows the periodicity. (DEFAULT: 1)"
write(*,*) "  -repeat_y=[number] : For the generated picture: the number"
write(*,*) "    of unit cell repetitions along y-axis to generate a larger  "
write(*,*) "    picture which shows the periodicity. (DEFAULT: 1)"
write(*,*) "  -grid_x=[number] : The number of grid points along x-axis "
write(*,*) "    in the unit cell at which the STM picture shall be "
write(*,*) "    calculated. (DEFAULT: 100)"
write(*,*) "  -grid_y=[number] : The number of grid points along y-axis "
write(*,*) "    in the unit cell at which the STM picture shall be "
write(*,*) "    calculated. (DEFAULT: 100)"
write(*,*) "  -gauss_width=[value] : The full width of half maximum of the"
write(*,*) "    Gaussian that is used to smooth the image by including the "
write(*,*) "    charge densities around an evaluation point by weighting them"
write(*,*) "    with the Gaussian function value. The larger the value, the"
write(*,*) "    smoother/smeared the picture will be. (DEFAULT: 0.3 Angs.)"
write(*,*) 
!
!    For the constant current mode: the moving up of the tip for each 
!    new row position (in Angstroms)
!
def_z_shift=2.d0
!
!    Read command line arguments
!
!    If the constant height or the constant current mode 
!      shall be used
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-mode=") then
      read(arg(7:),*,iostat=readstat) stm_mode
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -mode=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "height") then
   write(*,*) "The constant height mode will be evaluated."
else if (trim(stm_mode) .eq. "current") then
   write(*,*) "The constant current mode will be evaluated."
else 
   write(*,*) "Please give either the constant height mode (-mode=height) "
   write(*,*) " or the constant current mode (-mode=current)!"
   stop
end if        
!
!    If constant height mode: read in the height (z/c coordinate)
!
stm_height=-1.d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:5))  .eq. "-pos=") then
      read(arg(6:),*,iostat=readstat) stm_height
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -pos=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "height") then
   if (stm_height .lt. 0.d0) then
      write(*,*) "Please give the height (z-coordinate) at which the STM"
      write(*,*) " tip shall scan the surface!"
      stop
   else
      write(*,'(a,f13.6,a)') " The STM tip will scan the surface at z=", &
                      & stm_height," Angstroms" 
   end if        
end if        
!
!    If constant current mode: read in the isos density 
!
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:6))  .eq. "-dens=") then
      read(arg(7:),*,iostat=readstat) isos_dens
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -dens=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
if (trim(stm_mode) .eq. "current") then
   if (isos_dens .lt. 0.d0) then
      write(*,*) "Please give a positive density (e/Ang^3) to which the "
      write(*,*) " shall be moved down in each grid point!"
      stop
   else
      write(*,'(a,f10.6,a)') " The local density to which the STM tip will be dropped is ", &
                      & isos_dens," e/Ang.^3"
   end if
end if

!
!     The number of STM image repetitions along the x axis
!
!     Default (only actual unit cell)
!
repeat_x=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-repeat_x=") then
      read(arg(11:),*,iostat=readstat) repeat_x
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -repeat_x=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i3)') " Number of unit cell repetitions along x: ",repeat_x

!
!     The number of STM image repetitions along the y axis
!
!     Default (only actual unit cell)
!
repeat_y=1
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:10))  .eq. "-repeat_y=") then
      read(arg(11:),*,iostat=readstat) repeat_y
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -repeat_y=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i3)') " Number of unit cell repetitions along y: ",repeat_y
!
!     The number of grid points for STM printout along the x axis 
!
!     Default
!
stm_gridx=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-grid_x=") then
      read(arg(9:),*,iostat=readstat) stm_gridx
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -grid_x=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i6)') " Number of plotted STM grid points (unit cell) along x: ",stm_gridx

!
!     The number of grid points for STM printout along the y axis
!
!     Default
!
stm_gridy=100
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:8))  .eq. "-grid_y=") then
      read(arg(9:),*,iostat=readstat) stm_gridy
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -grid_y=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,i6)') " Number of plotted STM grid points (unit cell) along y: ",stm_gridy

!
!     The width at half maximum (FWHM) of the Gaussian used for 
!       smoothing the picture
!
!     Default
!
g_width=0.3d0
do i = 1, command_argument_count()
   call get_command_argument(i, arg)
   if (trim(arg(1:13))  .eq. "-gauss_width=") then
      read(arg(14:),*,iostat=readstat) g_width
      if (readstat .ne. 0) then
         write(*,*)
         stop "Check the command -gauss_width=..., something went wrong!"
         write(*,*)
      end if
   end if
end do
write(*,'(a,f12.7)') " Width of the smooting Gaussian at half height (Angs.): ",g_width
!
!    Determine the properties and the function of the smoothing Gaussian
!
a_fac=dlog(2.d0)/((g_width/2d0)**2)
cutoff=sqrt(dlog(100d0)/a_fac)

write(*,*) 
write(*,*) "Read in the PARCHG file ..."
!
!    Open file PARCHG and read in relevant data
!
open(unit=15,file="PARCHG",status="old",iostat=readstat)
if (readstat .ne. 0) then
   write(*,*) "The file PARCHG is not there!"
   stop
end if        
!
!    Read the unit cell definition
!
read(15,*)
read(15,*) scale
read(15,*) coord_mat(1,:)
read(15,*) coord_mat(2,:)
read(15,*) coord_mat(3,:)
!
!    The lengths of coordinate axes
!
x_len=sqrt(dot_product(coord_mat(1,:),coord_mat(1,:)))
y_len=sqrt(dot_product(coord_mat(2,:),coord_mat(2,:)))
z_len=sqrt(dot_product(coord_mat(3,:),coord_mat(3,:)))
!
!    Read in the elements
!
allocate(el_names_read(10))
el_names_read="XX"
read(15,'(a)') cdum
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
read(15,*) el_nums

!
!    Define the element symbols for all atoms in the system
!
natoms = sum(el_nums)
allocate(at_names(natoms))
k=0
do i=1,nelems
   do j=1,el_nums(i)
      k=k+1
      at_names(k)=el_names(i)
   end do
end do
!
!    Read in the atomic coordinates
!
allocate(xyz(3,natoms))
read(15,*)
do i=1,natoms
   read(15,*) xyz(:,i)
end do
!
!    Read in the number of grid points
!
read(15,*)
read(15,*) grida,gridb,gridc
allocate(cdens(grida,gridb,gridc))
!
!    Read in the partial charge density 
!    Use compactified read command for array
!
read(15,*) (((cdens(nx,ny,nz),nx=1,grida),ny=1,gridb),nz=1,gridc)
close(15)
write(*,*) "  done!"
write(*,'(a)') " Length of coordinate axes (Angstroms):"
write(*,'(a,f12.6,a,f12.6,a,f12.6,a)') "   x:",x_len,", y:",y_len,", z:",z_len
write(*,'(a,i6,a,i6,a,i6,a,i9,a)') " Number of partial density grid points:"
write(*,'(a,i6,a,i6,a,i6,a,i9,a)') "   x:",grida, &
              & ", y:",gridb,", z:",gridc, "   (total:",grida*gridb*gridc,")"
write(*,*)
!
!    Calculate x,y,z cartesian coordinates for all grid points, for the 
!    subsequent determination of distances to the current sampling point
!
allocate(c_coord(3,grida,gridb,gridc))
do i=1,gridc
   do j=1,gridb
      do k=1,grida
         vec_act(1)=real((k-1))/real(grida)
         vec_act(2)=real((j-1))/real(gridb)
         vec_act(3)=real((i-1))/real(gridc)
         c_coord(:,k,j,i)=vec_act(:)!matmul(vec_act,coord_mat)
      end do
   end do
end do    
!
!    Allocate data array for 2D STM picture
!
allocate(stm_dat(stm_gridy,stm_gridx))
allocate(stm_pos(2,stm_gridy,stm_gridx))
stm_dat=0.d0
!
!    Determine number of grid points to be included into the Gaussian
!    smoothin region
!    Along z only half the length! (since we do not probe into depth)
!
include_x=nint(real(cutoff)/real(x_len)*real(grida))
include_y=nint(real(cutoff)/real(y_len)*real(gridb))
include_z=nint(real(cutoff)/real(z_len)*real(gridc)/2)
write(*,*) "Generate the STM picture at all plot grid points ...  "
!
!    A: CONSTANT HEIGHT MODE
!
if (trim(stm_mode) .eq. "height") then
!
!    Loop over all picture gridpoints in the x-y (or a-b) plane and calculate 
!    the Gaussian-weighted interpolation of nearby grid points
!
   eval_stat = .false.
   do i=1,stm_gridx
!
!    Every 10% of the read in, give a status update
!
      do j=1,10
         if (real(i)/real(stm_gridx) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do

      do j=1,stm_gridy
         act_pos(1)=real((i-1))/real(stm_gridx)
         act_pos(2)=real((j-1))/real(stm_gridy)
         act_pos(3)=stm_height/z_len
!
!    Determine nearest point in density grid
!
         near_x=nint(act_pos(1)*grida)
         near_y=nint(act_pos(2)*gridb)
         near_z=nint(act_pos(3)*gridc)
       !  act_pos(:)=matmul(act_pos,coord_mat)
         
!
!    Now loop over all other grid points within the proposed cutoff of the 
!    Gaussian smearing function
!
         weight_sum=0
         do k=near_x-include_x,near_x+include_x
            if (k .lt. 1) then
               k_act=grida+k
            else if (k .gt. grida) then
               k_act=k-grida
            else     
               k_act=k
            end if     
            do l=near_y-include_y,near_y+include_y
               if (l .lt. 1) then
                  l_act=gridb+l
               else if (l .gt. gridb) then
                  l_act=l-gridb     
               else 
                  l_act=l
               end if     
               do m=near_z-include_z,near_z+include_z
                  if (m .lt. 1) then
                     m_act=gridc+m
                  else if (m .gt. gridc) then
                     m_act=m-gridc
                  else
                     m_act=m
                  end if  
!
!    Now sum over all grid points within the chosen cutoff for consideration
!    and add the components together
!
                !  write(*,*) k_act,l_act,m_act
                  diff_vec(:)=act_pos(:)-c_coord(:,k_act,l_act,m_act)       
                  
                  do n=1,3
                     if (diff_vec(n) .gt. 0.5d0) then
                        diff_vec(n) = diff_vec(n)-0.5d0     
                     end if  
                     if (diff_vec(n) .lt. 0.d0) then
                        diff_vec(n) = diff_vec(n)+0.5d0 
                     end if        
                  end do

                  diff_vec=matmul(diff_vec,coord_mat)
                  dist_act=sqrt(dot_product(diff_vec,diff_vec))
                  weight=exp(-g_width*dist_act*dist_act)
                  weight_sum=weight_sum+weight
                  stm_dat(j,i)=stm_dat(j,i)+weight*cdens(k_act,l_act,m_act)
               end do
            end do
         end do      
         stm_dat(j,i)=stm_dat(j,i)/weight_sum
         act_pos=matmul(act_pos,coord_mat)
         stm_pos(:,j,i)=act_pos(1:2)

      end do
   end do    
end if
!
!    B: CONSTANT CURRENT MODE
!
if (trim(stm_mode) .eq. "current") then
!
!    Loop over all picture gridpoints in the x-y (or a-b) plane and calculate 
!    the Gaussian-weighted interpolation of nearby grid points
!    For constant current, start for the first point of each row at 0.8 times the 
!    maximum z-value of the cell and scan below until a value equal or higher the 
!    desired current (local charge density) is reached, similar to a STM tip that
!    is moved down to a surface
!    Within each row, from one point to the next, only scan the next +/- 1 Angstrom
!    relative to the previous point, assuming that the charge density won't change too
!    much between two gridpoints
!
!    The default movement along z for scanning the correct height, a fraction of the 
!    z-grid density within the given PARCHG file.
!
   eval_stat = .false.
   z_step=1.d0/(gridc*3d0)
   do i=1,stm_gridx
!
!    Every 10% of the read in, give a status update
!
      do j=1,10
         if (real(i)/real(stm_gridx) .gt. real(j)*0.1d0) then
            if (.not. eval_stat(j)) then
               write(*,'(a,i4,a)')  "  ... ",j*10,"% done "
               eval_stat(j) = .true.
            end if
         end if
      end do
 
      if (i .eq. 1) then
         prev_z=0.8d0-def_z_shift/z_len
      end if
!
!    Now loop trough the remaining points of the current row and start with the init z value
!    of this row
!
!      prev_z=prev_z_init
      do j=1,stm_gridy
         act_pos(1)=real((i-1))/real(stm_gridx)
         act_pos(2)=real((j-1))/real(stm_gridy)
         act_pos(3)=prev_z+def_z_shift/z_len
!
!    Scan the starting position for a new grid point
!
         do
            if (act_pos(3) .ge. 1.d0) then
               write(*,*) "The tip would be located outside the cell!"
               write(*,*) "Increase the -dens parameter if possible!"
               stop
            end if        
!
!    Determine nearest point in density grid
!
            near_x=nint(act_pos(1)*grida)
            near_y=nint(act_pos(2)*gridb)
            near_z=nint(act_pos(3)*gridc)
!
!    Now loop over all other grid points within the proposed cutoff of the 
!    Gaussian smearing function
!
            weight_sum=0.d0
            stm_act=0.d0
            do k=near_x-include_x,near_x+include_x
               if (k .lt. 1) then
                  k_act=grida+k
               else if (k .gt. grida) then
                  k_act=k-grida
               else
                  k_act=k
               end if
               do l=near_y-include_y,near_y+include_y
                  if (l .lt. 1) then
                     l_act=gridb+l
                  else if (l .gt. gridb) then
                     l_act=l-gridb
                  else
                     l_act=l
                  end if
                  do m=near_z-include_z,near_z+include_z
                     if (m .lt. 1) then
                        m_act=gridc+m
                     else if (m .gt. gridc) then
                        m_act=m-gridc
                     else
                        m_act=m
                     end if
!
!    Now sum over all grid points within the chosen cutoff for consideration
!    and add the components together
!
                     diff_vec(:)=act_pos(:)-c_coord(:,k_act,l_act,m_act)
   
                     do n=1,3
                        do while (abs(diff_vec(n)) .gt. 0.5d0)
                           diff_vec(n) = diff_vec(n) -sign(1.0d0,diff_vec(n))
                        end do
                     end do
                     diff_vec=matmul(diff_vec,coord_mat)
                     dist_act=sqrt(dot_product(diff_vec,diff_vec))
                     weight=exp(-g_width*dist_act*dist_act)
                     weight_sum=weight_sum+weight
                     stm_act=stm_act+cdens(k_act,l_act,m_act)
                  end do
               end do
            end do
!
!    Leave the scanning loop if either the corrent density has been found or 
!    the tip reached the ground of the unit cell
!
            stm_act=stm_act/weight_sum
            if ((stm_act .ge. isos_dens) .or. (act_pos(3) .lt. 0.01d0)) then
               act_pos=matmul(act_pos,coord_mat)
               stm_pos(:,j,i)=act_pos(1:2)
               stm_dat(j,i)=act_pos(3)
               prev_z=act_pos(3)/z_len
               exit
            end if
            act_pos(3)=act_pos(3)-z_step
         end do
      end do   
   end do
end if
        
write(*,*) " done!"
!
!    Write plot file with STM data
!
write(*,*) "Write file with 2D plot data for STM image (stm_plot.dat) ..."
open(unit=37,file="stm_plot.dat",status="replace")
write(37,*) "# This data has been generated by the program eval_stm."
write(37,*) "# x-coord     y-coord     intenstiy  "
x_plot_min=100.d0
x_plot_max=-100.d0
y_plot_min=100.d0
y_plot_max=-100.d0
do k=1,repeat_x
   do i=1,stm_gridx
      do l=1,repeat_y
         do j=1,stm_gridy
            vec_act(1)=real(k-1)
            vec_act(2)=real(l-1)
            vec_act(3)=1.d0
            vec_act=matmul(vec_act,coord_mat)
            vec_act(1)=vec_act(1)+stm_pos(1,j,i)
            vec_act(2)=vec_act(2)+stm_pos(2,j,i)
            if (vec_act(1) .lt. x_plot_min) then
               x_plot_min=vec_act(1)
            else if (vec_act(1) .gt. x_plot_max) then
               x_plot_max=vec_act(1)
            end if
            if (vec_act(2) .lt. y_plot_min) then
               y_plot_min=vec_act(2)
            else if (vec_act(2) .gt. y_plot_max) then
               y_plot_max=vec_act(2)
            end if
            write(37,*) vec_act(1), vec_act(2), stm_dat(j,i)
         end do
      end do 
      write(37,*)  
   end do   
end do   
close(37)
write(*,*) " done!"
!
!    Write gnuplot file for the plot
!    Automatically align the plotted png picture to the proportion of 
!     x and y axes lengths
!
write(*,*) "Write gnuplot file for 2D plot of STM image (stm_plot.gnu) ..."
xtoy=(x_plot_max-x_plot_min)/(y_plot_max-y_plot_min)
open(unit=38,file="stm_plot.gnu",status="replace")
write(38,*) "# This file generates a png picture from 'stm_plot.dat'."
write(38,*) "set terminal png size 2400,",nint(2000/xtoy)," lw 3.5 font 'Helvetica,46'"
write(38,*) "set output 'stm_plot.png'"
write(38,*) "set encoding iso_8859_1"
write(38,*) "set lmargin screen 0.10"
write(38,*) "set rmargin screen 0.80"
write(38,*) "set tmargin screen 0.13"
write(38,*) "set bmargin screen 0.95"
write(38,*) "set xrange [",x_plot_min,":",x_plot_max,"]"
write(38,*) "set yrange [",y_plot_min,":",y_plot_max,"]"
write(38,*) "set xlabel 'x-coordinate ({\305})"
write(38,*) "set ylabel 'y-coordinate ({\305})"
if (trim(stm_mode) .eq. "height") then
   write(38,*) "set cblabel 'local density (e/{\305}^3)' offset -.8,0"
end if
if (trim(stm_mode) .eq. "current") then
   write(38,*) "set cblabel 'height ({\305})' offset -.8,0"
end if       
write(38,*) "set cblabel offset 0.8,0" 
write(38,*) "set palette gray"
write(38,*) "set pm3d map interpolate 0,0"
write(38,*) "splot 'stm_plot.dat' with pm3d notitle"
close(38)

write(*,*) " done!"
write(*,*) "Execute gnuplot to obtain image (stm_plot.png) ..."
call system("gnuplot stm_plot.gnu")
write(*,*) " done!"
write(*,*)
!
!    Write second gnuplot file for plot without axis labels or tics
!
write(*,*) "Write gnuplot file for 2D plot wihtout axis labels (stm_plot_blank.gnu) ..."
xtoy=(x_plot_max-x_plot_min)/(y_plot_max-y_plot_min)
open(unit=38,file="stm_plot_blank.gnu",status="replace")
write(38,*) "# This file generates a png picture from 'stm_plot.dat'."
write(38,*) "set terminal png size 2400,",nint(2400/xtoy)," lw 3.5 font 'Helvetica,46'"
write(38,*) "set output 'stm_plot_blank.png'"
write(38,*) "set encoding iso_8859_1"
write(38,*) "set lmargin screen 0.0"
write(38,*) "set rmargin screen 1.0"
write(38,*) "set tmargin screen 0.0"
write(38,*) "set bmargin screen 1.0"
write(38,*) "set xrange [",x_plot_min,":",x_plot_max,"]"
write(38,*) "set yrange [",y_plot_min,":",y_plot_max,"]"
write(38,*) "unset xtics"
write(38,*) "unset ytics"
write(38,*) "unset cbtics"
write(38,*) "unset colorbox"

write(38,*) "set palette gray"
write(38,*) "set pm3d map interpolate 0,0"
write(38,*) "splot 'stm_plot.dat' with pm3d notitle"
close(38)

write(*,*) " done!"
write(*,*) "Execute gnuplot to obtain image without labels (stm_plot_blank.png) ..."
call system("gnuplot stm_plot_blank.gnu")
write(*,*) " done!"
write(*,*)
write(*,*) "eval_stm has finished all tasks, goodbye!"
write(*,*)

end program eval_stm        
