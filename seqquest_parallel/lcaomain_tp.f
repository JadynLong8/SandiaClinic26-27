       PROGRAM SEQQUEST
csubc:start
csubc To turn Quest into subroutine
csubc  (1) activate this section of code
csubc  (2) activate "Lsubbed" flag in utl*.f file
csub      IMPLICIT NONE
csub      REAL*8  engytotl
csub      write(*,'(a)') 'SEQQUEST wrapper'
csub      call SUBQUEST( engytotl )
csub      write(*,'(a,f20.8)') 'SEQQUEST energy=',engytotl
csub      STOP
csub      END
csub      SUBROUTINE SUBQUEST( engytotl )
csubc:end
c
c***********************************************************************
c
c PROGRAM: SEQQUEST
c
c PURPOSE: Gaussian-based "lcao", total energy&forces, pseudopotential,
c          density functional, electronic structure code
c          for both clusters and periodic systems, neutral or charged
c
c WRITTEN: Peter A. Schultz
c          Sandia National Laboratories
c          paschul@sandia.gov
c
c          Based on a code originally written by Peter J. Feibelman
c
c Documentation:  http://www.cs.sandia.gov/~paschul/Quest/
c           and:  http://dft.sandia.gov/Quest/
c
c***********************************************************************
c
c Copyright 2007 Sandia Corporation.
c Under the terms of Contract DE-AC04-94Al85000, there is a non-exclusive
c license for use of this work by or on behalf of the U.S. Government.
c Export of this program may require a license from the U.S. Government.
c
c NOTICE:
c
c For five (5) years from 09/20/2007, the United States Government is granted
c for itself and others acting on its behalf a paid-up, nonexclusive,
c irrevocable worldwide license in this data to reproduce, prepare derivative
c works, and perform publicly and display publicly, by or on behalf of the
c Government. There is provision for the possible extension of the term of
c this license. Subsequent to that period or any extension granted, the
c United States Government is granted for itself and others acting on its
c behalf a paid-up, nonexclusive, irrevocable worldwide license in this data
c to reproduce, prepare derivative works, distribute copies to the public,
c perform publicly and display publicly, and to permit others to do so.
c The specific term of the license can be identified by inquiry made to
c Sandia Corportation or DOE.
c
c NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF
c ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES ANY
c WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE
c ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
c OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY
c OWNED RIGHTS.
c
c Any licensee of This software has the obligation and responsibility to abide
c by the applicable export control laws, regulations, and general prohibitions
c relating to the export of technical data. Failure to obtain an export control
c license or other authority from the Government may result in criminal
c liability under U.S. laws.
c
c***********************************************************************
c
c Thanks to (in alphabetical order):
c  Arthur H. Edwards (AHE)     - post-processing, mpi/clusters, band structure
c  Dean McCoullough (DMM)      - simple task-parallelism
c  Ann E. Mattsson (AEM)       - spin-gga stress, density fcnal development(BLYP)
c  Jonathan E. Moussa (JEM)    - HSE semi-local functionals
c  Richard P. Muller (RPM)     - cell optimization, molecular dynamics
c  Andrew C. Pineda  (ACP)     - reduced memory k-parallelism, parallel-parallel solver
c  David B. Raczkowski (DBR)   - real eigensolvers; eigfcn density; level-3 blas
c  Jamil Tahir-Kheli (JTK)     - external electric fields
c  Aidan P. Thompson (APT)     - NEB transition state finder
c  Renee M. Van Ginhoven (RMV) - expand multi-image (NEB) methods, SLIC
c
c***********************************************************************
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      CHARACTER  version*5,         lastup*11
      DATA       version /'2.68e'/, lastup/'29Mar20-PAS'/
c
c-----------------------------------------------------------------------
c
c REVISION HISTORY: (see lcao.history)
c 
c v2.68/Jun-2018: Expanded(bonds/slinky)/unified slics & idle proc solver bug (Aug17); git
c   (a) 29Dec19-PAS - compiler warnings (nonscfk dim(vpsrad,makelies); broydata(blsav),
c                     FLSPECS->FLNSTUFF commons w/flstuff
c   (b)  2Jan20-PAS - reduced mem prhoij (to fix nstate=norb memory bug)
c   (c)  7Jan20-PAS - more minor bug-fixes
c   (d) 29Jan20-PAS - oopsy in slic-updated GUPDATE for NEB
c   (e) 29Mar20-PAS - to do partial (npop<nstate) read of ivecfl, add READB1
c v2.67/Aug-2015: label-selected atoms; graceful neb stall; bvl init; bugfixes
c   (a) 24Dec15-PAS - print out gamma point eigenspectrum in bands (fix bug)
c   (b) 25Apr16-PAS - LNREADTYP iatm>99 bugfix
c v2.66/24Apr15-PAS: reduced disk I/O; geom file input; ges(k) off
c v2.65/Aug13-PAS: reset defaults; hse semilocal(JEM); vdw corrections (ulg,D2);
c              enhanced cell constraints and optimization
c    (a) 24Sep13-PAS: NEB+LMCC fix: propagate lmcc reference file to all image masters
c                     compute enthalpy for pressure calculations.
c    (b)  5Feb14-ACP&PAS: fix dscal arg in cellopt (fixed-volume); ACP: fix bands symbols
c    (c)  9May14-PAS: fix vdwdata input bugs (2): ijvdw args, and errant goto
c    (d) 15Dec14-PAS: fix bugs to band structure (nblocksz->nkblocksz), eV scale twice
c v2.64/Jul13: band structure(AHE), k-parallel-parallel(ACP), cleanup(ACP&PAS)
c v2.63/May12-ACP&PAS: k-parallel scf and k-point distributed memory; dble-lyr
c    (a) 16Jan13-PAS: bugfix parallel solver memory dist for some small 4-proc cases
c    (b) 26Feb13-PAS: bugfix wkmem parallel stripe v. square(ijk)
c v2.62/Feb12-PAS: stripe-parallel 2/3-center routines; vgfix/cell-volume constraints;
c                  basic double-layer boundary conditions
c v2.61/Mar08-PAS: spin-optimization; grid parallelism; MD/NVT thermostats
c    (a) 17Apr08-PAS: bugfix: nlat1c void with floaters at origin
c    (b)  2May08-PAS: bugfix: parallel populations code argument
c    (c) 17May08-PAS: bugfix: memory allocation in (parallel) spin code
c    (d) 24May08-PAS: tweak mem allocs to allow larger (parallel) spin calcs
c    (e)  2Jun08-PAS: bugfix: fix to 261(a)
c    (f) 28Aug08-PAS: bugfix: fix to 261(a) (refix an old bug reinserted in a)
c    (g)  4Sep08-PAS: bugfix: GGA and left-handed vectors (bug reinserted 2.55a)
c    (h) 23Sep08-PAS: bugfix: memory allocator in wkmem
c    (i) 30Sep08-PAS: bugfix: pop uninit variable, EF pop error
c    (j) 26Jan09-PAS: bugfix: spsteps input bug; rhosave-as-subroutine
c    (k) 29May10-PAS: bugfix: dosdump, lcao.post file repairs for prop1e
c v2.60/Jun07-PAS: wrap more mpi; in-core broyden history; merge serial-tp,
c    SLIC minimization (RMV);
c    (a)  2Oct07-PAS: Sandia Copyright notice.
c    (b) 16Nov07-PAS: nebsetup/mprecv bug fixed
c    (c) 10Dec07-PAS: portability issues
c    (d) 20Dec07-PAS: dangling BLACS_GRIDINIT now exited
c    (e) 28Jan08-PAS: parallel- griddft uninit file unit rewind
c    (f) 12Feb08-PAS: indexing error in nearby correction
c v2.59/Jan06-PAS: task-parallelism; changes to maximize consistency
c    with serial version; wall timer; defect/closed-shell sampling;
c    update NEB (RMV); Quest as subroutine capability; mat-files (idhamfl);
c    master-only matrix files.
c    (a) 17Jul07-ACP/PAS: tp/pgrid0xc - initialized variable
c    (b)  9Apr07-AHE/PAS: Art's DOSDUMP/post-proc update
c
c-----------------------------------------------------------------------
c
c To do:
c   - poscar input
c   - complete HSE implementation
c   - rework energy analysis (kpt), memory usage
c   - compact/distribute 3N^2 matrix bottlenecks
c   - change grid field default to distrbuted boxed rather than unboxed
c   - "tight" memory packing
c   - clean up density matrix usage
c   - new geometry specification/optimization (Z-matrix?)
c   - qdb training set output
c
c   - vdw configuration into separate file
c   - global graded output controls (master-obly, and silent)
c   - new nearby correction weighting scheme
c   - atom blocking, and install order-N code in production code
c   - clean up cell optimization (rscale,cartesian atomic coords,rsym)
c   - clean out do termination lines with operations
c   - unrestricted spin-DFT for net zero spin systems (e.g.,enable AFM)
c   - lmcc for initial iteration of charged system
c   - xc 1-ctr corrections
c   - develop new radial grids
c   - split command and input i/o
c   - new optimized fft routine
c   - body-centered xc grids (make symmetry work off-grid)
c   - complete boxing grids: fast stuff
c   - rework v1ctrij (split v1ctrij?=vl1cpot,vl1cmat)
c
c-----------------------------------------------------------------------
c
c     INPUT ROUTINES FOR THIS LCAO ELECTRONIC STRUCTURE PROGRAM
c
c   Input of command options is done in the routine RUNOPT
c      - the default is to do (scf) iterations and not do forces
c
c   Input of structure/setup is done in the routine SETDATA
c      - data that defines the structure/symmetry/grids of the system
c
c   Input of atomic pseudopotentials/bases done in ATMDAT1
c
c   Input of run (scf) options is done in the routine RUNDATA
c      - modify parameters affecting the scf calcalation
c
c   Input of geometry relaxation options is done in GEOMDAT
c      - constraints/criteria for the relaxation of the structure
c
c   Input of cell optimization options is done in CELLDAT
c      - modify parameters controlling optimization of cell shape
c
c   Input for control of nudged elastic band (NEB) done in NEBDAT
c      - setup for the NEB transition state finder
c
c   Input for control of molecular dynamics (MD) done in QMDDAT
c      - setup for molecular dynamics - MD
c
c-----------------------------------------------------------------------
c
c           NOTES FOR THIS LCAO ELECTRONIC STRUCTURE PROGRAM:
c
c   Currently the atom coordinates given in the input are shifted in
c   the direction along each non-periodic direction by half of the
c   unit cell primitive vector in that direction.  For example, in a
c   slab calculation, each atom coordinate is shifted by half the
c   length of the third primitive vector (assumed to be the z-axis).
c   Thus if an atom is specified to lie at z=0 in the input data,
c   the program automatically shifts to a value of z=0.5*rprim(3,3).
c   For a molecular calculation, the point (0,0,0) will be moved to
c   the center of the unit cell.  This origin can be shifted in input.
c   For internal use, the code automatically projects the atoms outside
c   the primary unit cell (the parallelpiped defined by the given
c   primitive lattice vectors) into the cell along periodic directions.
c   Atoms outside the primary cell in a non-periodic direction will
c   trigger an error.
c
c-----------------------------------------------------------------------
c
c             *************** DIMENSIONING ***************
c
      INCLUDE    'lcao.params'
c
c-----------------------------------------------------------------------
c
c Allocate space for big work array
c MACHINE-DEPENDENCE:
c  Methods for allocating memory for arena wk():
c   (1) blank common - SGI and most others, including pgi
c   (2) explicit declaration - if blank common fails (some linux)
c   (3) dynamic allocation (thx, C. Bauschlicher) - if both fail
CM1   COMMON  //  wk(maxwkd)
CM2   DIMENSION   wk(maxwkd)
      ALLOCATABLE :: wk(:) ! and uncomment "CM3" below to allocate
c
c-----------------------------------------------------------------------
c
c  Control flags:
      LOGICAL  dotests
      LOGICAL  dosetup,doiters,doforce,dorelax,do_cell,do_neb,do_md
      LOGICAL  doflag
      LOGICAL  DO_VDW
      LOGICAL  Lijkmat
      LOGICAL  do_eigvecs, do_eigpops
      LOGICAL  do_kppsolve, do_psolv
      LOGICAL  do_post
      LOGICAL  do_bands, do_optic, do_nonscf,Lnonscf
      LOGICAL  do_geom, opt_spin
      LOGICAL  dosplit,doblas3, dolocxc
      LOGICAL  dospars,doeigen,doprjct,dorspsy,dodefct
c
      LOGICAL  do_gga, do_spin
      LOGICAL  do_field, do_pvcalc
      LOGICAL  dokdefct
      LOGICAL  dolvatm, dokgrid, redusym
      LOGICAL  anycore, lmcc_on
      LOGICAL  any1ctr,anyxc1c
      LOGICAL  engydone
ctp:
CTP      DATA     lvlengy / 1 /
      DATA     lvlengy / 4 /
c
      CHARACTER  fcnal*6
      CHARACTER  filenm*128
c cutoffs:
      DIMENSION  convs(8)
c
c-----------------------------------------------------------------------
c
c  Task-parallel data structures (should be fused into serial mem):
      PARAMETER  ( maxnodes =  512 )
      INTEGER    icluster
      DOUBLE PRECISION  gapp
      DIMENSION  icluster(2*maxnodes)
      DIMENSION  gapp(maxnodes)
c  parallel-parallel solver data structures (PEIGSOLVKP/PSCHROEDKP)
      INTEGER    masterlist
      DIMENSION  masterlist(maxnodes)
c
c  External electric field:
      DIMENSION  efield(3)
c  Arrays for box stuff:
      DIMENSION  ndbox(2,2,3), inboxs(nrecd)
      DIMENSION  nearatm(3,natmd)
c   ... these are only used for scf fast density code (non-standard 2D):
      DIMENSION  mngr3(norbd),mxgr3(norbd),ngr(norbd)
c
c  Atom definition data arrays:
      PARAMETER  (natmnm=12)
      CHARACTER*(natmnm)  typnm(ntypd), atmnm(natmd)
      DIMENSION  itypa(natmd), numshl(ntypd),lshel(nshld,ntypd),
     $  nala(nshld,ntypd),ala(nald,nshld,ntypd),cala(nald,nshld,ntypd)
      DIMENSION  nocc(ntypd),locc(nshld,ntypd),
     $ nalsh(nshld,ntypd),alsh(nald,nshld,ntypd),calsh(nald,nshld,ntypd)
      DIMENSION  alamin(ntypd), norba(ntypd),ioffa(natmd)
      DIMENSION  znuc(ntypd),occsh(nshld,ntypd),occij(nshld,nshld,ntypd)
      DIMENSION  atmass(ntypd), atengy(ntypd)
cv vdw arrays
      DIMENSION  c6vdw(ntypd), r0vdw(ntypd)
      DIMENSION  c6ijvdw(ntypd*ntypd), r0ijvdw(ntypd*ntypd)
      DIMENSION  lmxnlp1(ntypd),almnnl(ntypd)
      DIMENSION  nrad(ntypd),radmsh(nrd,ntypd),radwt(nrd,ntypd),
     $  nrps(ntypd),vpsrad(nrd,4,ntypd)
      DIMENSION  nrcor(ntypd),corden(nrd,ntypd),cordrv(nrd,ntypd)
c  ... these only used for scf fast density code (non-standard)
      DIMENSION  lmx1ctr(ntypd),almnv1c(ntypd),
     $  alxcmin(ntypd),nrxcfit(ntypd),nalfxc(ntypd),alfxc(nalfxcd,ntypd)
c
c  Coordinate arrays:
c     nchrgd = number of independent LMCC charges
      PARAMETER  (nchrgd=2)
      DIMENSION  rprim(3,4),rscale(3)
      DIMENSION  dorig(3),orig(3), rsym(3)
      DIMENSION  origws(3),origws0(3)
      DIMENSION  rchrg(3*nchrgd),rchrg0(3*nchrgd)
      DIMENSION  ratm(3,natmd),ratm0(3,natmd), atveloc(3,natmd)
      DIMENSION  ratma(3,natmd),ratmb(3,natmd)
      DIMENSION  hh(3,3),gprim(3,3)
      DIMENSION  rlat(3,nlatd)
      DIMENSION  veck(3,nkd),veckg(3,nkd),wtk(nkd),nku(3)
c
c  Arrays needed for real space symmetrization of charge density:
      DIMENSION  isymtyp(nsymd),syminfo(6*nsymd)
      DIMENSION  rmatsym(12*nsymd),isymgrd(nsymd)
      DIMENSION  dsym(5,5,3,nsymd),naofsym(natmd*nsymd)
c
c  Array to track various grid/analytic reference atom integrals:
      DIMENSION  xc0val(12)
c
c  Integrals of local density, grid v. analytic
      DIMENSION  rhoatom(3*natmd), defatom(2*natmd,2)
      DIMENSION  efermi(2),egap(2),elecno(2)
      DIMENSION  nstbulk(2),nedefct(2)
c
c The following arrays now taken out of big work space(v2.31):
c      DIMENSION  s1atom(noad,noad,nkd,natmd)
c      DOUBLE COMPLEX    orbint(norbd,nkd),orbofk(natmd,norbd,nkd)
c  Bloch phase factors:
c      DIMENSION  coskr(nkd,nlatd),sinkr(nkd,nlatd)
c  eigenvalue problem arrays:
c      DIMENSION  eigval(nstd,nkd),eigpop(nstd,nkd)
c      DIMENSION  angsav(mxang)
c                 40095 =   30375   +  1620   +   8100
      PARAMETER  (mxang = (125*3*81 + 4*5*9*9 + 4*25*9*9) )
c
c Radial arrays ...
c
c  ... potentials of various flavors
      DIMENSION  vatrad(nrd,ntypd),vesrad(nrd,ntypd),
     $           vxcrad(nrd,ntypd),excrad(nrd,ntypd)
c
c  ... scratch arrays for radial integrals (eg. bessel functions)
c mxwkrbx must be at least 14(?) times biggest possible box.
c See defbox for biggest box (216), and grid0xc for need (14x).
c I give it 16x to provide some flexibility in future grid0xc use
c wkatv space is also in grid0xc, and needs 11x natm as of 2.59
      EQUIVALENCE  (wkrbx,wkrad,wkatv)
      PARAMETER  (maxvbox=16,mxbsiz=216, mxwkrbx=maxvbox*mxbsiz)
      PARAMETER  (maxvrad=27,            mxwkrad=maxvrad*nrd)
      PARAMETER  (maxvala=11, maxatv=maxvala*natmd)
      DIMENSION  wkrbx(mxwkrbx), wkrad(nrd,maxvrad), wkatv(maxatv)
      DIMENSION  wspl(nrd,2)
c
c  Small scratch space (wksmX):
      EQUIVALENCE  (wksm, wksmo,wksma,wksml, wksmg,wksmr,wksmc,wksmk)
      DIMENSION  wksm(46*naldsq)
      DIMENSION  wksmo(norbd,9), wksma(natmd,24), wksml(nlatd,2)
      DIMENSION  wksmg(9,25,3), wksmr(nrd,2), wksmc(noad,noad)
      DIMENSION  wksmk(2,nkd)
c  fft arrays: wft(4*n1r+15 + 4*n2r+15 + 4*n3r+15 + 2*max(n1r,n2r,n3r))
      DIMENSION  wft(14*maxr1d+45)
c  Work array for angular stuff needed in non-local integrals:
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
c
c Spin optimization arrays (2*{spinold, new, ges, slope})
c
      DIMENSION  spindata(4,2)
c
c Arrays for non-SCF (band structure) calculations:
c
cxxx264: dimensions for non-scf should be moved to lcao.params/wkmem
      PARAMETER  ( nbranchd=12, nknscfd=1200 )
      DIMENSION  veckn(3,nknscfd), wtkn(nknscfd),vecknsm(3,nkd)
      DIMENSION  wtknsm(nkd),dk(nknscfd)
      DIMENSION  nbr(nbranchd)
cxxx264: eigenvals is O(N^2) memory, now taken out of wk() before nonscfk
C      DIMENSION  eigenvals(norbd,nkd,2)
c Arrays/variables for band structure calculations
      CHARACTER*3  bvl
      CHARACTER*2  bsymb(nbranchd)
c
c Force and stress arrays:
c
      DIMENSION  frctot(3,natmd), fdefct(3)
      DIMENSION  strtot(3,3), strlast(3,3), ucstepfac(3,3)
      DIMENSION  watmg(natmd)
c
      DIMENSION  frc1(3,natmd),frc2(3,natmd),frc3(3,natmd)
      DIMENSION  str1(3,3),str2(3,3),str3(3,3)
      DIMENSION  strxtrnl(3,3)
c
c Relaxation arrays:
c
      PARAMETER  (mxgeom=20)
      DIMENSION  gengy(mxgeom), ucengy(mxgeom)
c
c NEB arrays:
c
      DIMENSION  imgstat(2*nimgd)
      DIMENSION  engyneb(nimgd)
      DIMENSION  igrouplist_neb(nimgd+1)
      LOGICAL    Lantikink_neb,Loptall_neb,Lclimb_neb
      LOGICAL    Lnoninteract
c
c Various force/relaxation constraints
      DIMENSION  iforce(natmd)
      DIMENSION  iatframe(3)
c SLIC: projected coordinate variables and arrays
      DIMENSION  islics(4,natmd), vslics(3,natmd), islicats(natmd)
c
c Logical flag to produce lmcc/dble-lyr debug diagnostics
      LOGICAL    Ldevtest
      DATA       Ldevtest / .false. /
c
c1c#####################################################################
c
c Non-standard dynamic 1-center expansion of scf fast density and pots
c
c  ... coefficients for one-center expansion of scf fast density:
      DIMENSION  dcountr(10,natmd),rh1coef(10,noad,natmd)
c  ... scratch arrays (v1ctrij/blkvofn) to construct one-center terms:
      DIMENSION  orb(nlatd),grorb(3,nlatd),grgrorb(3,3,nlatd),
     $  wf1tmp1(3,nlatd),wf1tmp2(3,3,nlatd),
     $  wf1tmp3(3,nlatd),wf1tmp4(nlatd),wf1tmp5(nlatd)
c  ... bookkeeping for v1ctrij
      DIMENSION  ieqi(neq3d),ieqj(neq3d),ieqn(neq3d),ipermut(neq3d**2),
     $  vradsr(7,10,neq3d),v1csum(neq3d,49),velem(neq3d),
     $  pradsr(7,10,neq3d),p1csum(neq3d,49),pelem(neq3d),
     $  xradsr(7,10,neq3d),x1csum(neq3d,49),xelem(neq3d)
c  ... arrays for xc-fit fit code (not active, development only):
      PARAMETER  (nangd=12)
      DIMENSION  angpts(3*nangd),angwts(nangd),angylm(nangd,9)
      DIMENSION  nfityp(ntypd),radwf(nangd,nrd,noad,nfitd)
      DIMENSION  orbj(noad),grorbj(3,noad),grgrorj(3,3,noad),
     $           gofac(nangd),ggofac(nangd)
      DIMENSION  alsq(nalfxcd,nalfxcd,3,ntypd),
     $  rholoc(nrd,nangd),vxcloc(nrd,nangd),excloc(nrd,nangd),
     $  vlm(nrd),elm(nrd),rfac(nrd,nalfxcd),
     $  ipvt(nalfxcd,3,ntypd),vg(nalfxcd,9),eg(nalfxcd,9),
     $  cv(10,nalfxcd,natmd),ce(10,nalfxcd,natmd)
c1c#####################################################################
c
c Set up various cutoffs for code:
c
c  Cuts: gaus->exp argument; grid->wfn grid range; gwfn->wfn range;
c       slow-> extend wfn range for slow; ii->1c pot:wfn range;
      COMMON  /CUTS/ cutgaus,cutgrid,cutgwfn,cutslow,cutii,alfast
c
c  Zero density cutoff for xc potential calculation:
      DATA  xcrhocut / 1.d-21 /
c  Atomic xc fast function cutoff defining values:
      DATA  xcalcut,xcfac / 0.45d0,0.025d0 /
c        xcalcut > 0 => gaussian alpha =  xcalcut
c        xcalcut < 0 => gaussian alpha = -xcalcut*alamin(iatm)
c        xcfac = extending function (see below)
c        To recover original pjf heuristic, use (-1.,0.05)
c        Recommend "optimal" set =(.45,.025)
c          - 15Jun94-PAS
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Module control data
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c NB: none of this data is to be used directly.  It is only to be used
c     through the appropriate accessor functions within the "modules"
c
c  Timing module data:
      COMMON  /TIMSTUFF/  timetotl,timeuser,timesys,timelim,
     $ tw,tw0,twoff, itimer
C      DATA  itimer / 0 /
c
c  File i/o system (used in FLxxxx routines)
      COMMON  /FLUSTUFF/ IRDU,IWRU,IERRU,IFLU1,IFLU2,iflstat(100)
      COMMON  /FLNSTUFF/ ndirnm,njobnm,dirname,jobname
      CHARACTER*128  jobname,dirname
c
c  Big I/O wrappers
      COMMON  /BIGIO/ maxreclen
c
c  Parallel/message-passing data
      COMMON  /MPSTUFF/ nodetot(4),nodeid(4),node0id(4),icommid(4)
      COMMON  /KPSTUFF/ kpartyp
c
c  Define dft flavor types (used/set in DFxxxx routines) ...
      COMMON  /DFTTYPS/ i_fcnal, i_gga,i_spin,
     $ i_capz,i_pw91,i_pbe,i_blyp,i_am05,i_hse,i_hse2
c   ... and set defaults to undefined values
C      DATA  i_fcnal / -1 /, i_gga / -1 /, i_spin / -1 /
c
c  Geometry relaxation types (used in Gxxxxx routines)
      COMMON  /GEOMTYPS/ ig_type, ig_steep,ig_broy,ig_md0,
     $ ig_dmd,ig_dmdat,ig_asd,ig_lbfgs,ig_cg, ig_qmd
c   ... set default to undefined value
C      DATA  ig_type / -1 /
c
c  Geometry constraint/local-coord (SLICS) types
      COMMON  /SLICTYPS/ igc_fixed, igc_free,
     $ gc_slic, igc_slic1,igc_slic2,
     $ igc_frame,igc_bond,igc_vgp,igc_vgd,igc_vga,igc_vgl,igc_slink
      COMMON  /SLICPARS/ sl_pars(9), isl_fank
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  CONSTANTS INITIALIZATION
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  To have consistent values of some constants, put in common:
      COMMON  /CONST/ pi,rtpi,cutexp
c
c  Some physical constants worth having around:
      DATA  toev / 13.605805d0 /, toang / 0.5291771d0 /
c
c  Convenient numerical constants:
      DATA  zero,half,one,two,three / 0.d0,0.5d0,1.d0,2.d0,3.d0 /
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * * * * * * * *                       * * * * * * * * * * * *
c * * * * * * * * * * * *   EXECUTABLE CODE   * * * * * * * * * * * * *
c  * * * * * * * * * * *                       * * * * * * * * * * * *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
CM3: 3rd method for allocating arena workspace: dynamic allocation
      ALLOCATE( wk(maxwkd) )
c
c Pre-initialize some module data, to "undefined" values:
      itimer = 0
      i_fcnal = -1
      i_gga   = -1
      i_spin  = -1
      ig_type = -1
      call VDWINIT( 0 )
c
c Initialize code:
      call STARTX
c
c Get the standard default output unit:
      call FLGETIWR( IWR )
c announce this is slics version
      write(IWR,*) 'version 2.68: slics foundation version'
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c            Announce arrival, and identify the culprits:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write(IWR,9000)  version,lastup
 9000 format(1x, 72('=')
     $  /4x,'Program: SEQQUEST - LCAO-PS-DFT electronic structure code'
     $  /4x,'Version:  ',a,'     Last revised: ',a
     $  /4x,'Copyright 2007 Sandia Corporation.',
C     $  /4x,'Under the terms of Contract DE-AC04-94Al85000, there is a'
C     $      ' non-exclusive
C     $  /4x,'license for use of this work by or on behalf of the U.S.'
C     $      ' Government.'
C     $  /4x,'Export of this program may require a license from the',
C     $      ' U.S. Government.'
     $ //4x,'Written by Peter A. Schultz (Sandia National Laboratories)'
     $  /4x,'Based on a code originally written by Peter J. Feibelman'
     $ //4x,'Documentation:  http://www.cs.sandia.gov/~paschul/Quest/'
     $  /4x,'           and  http://dft.sandia.gov/Quest/'
     $  /4x,'Questions, comments to:  paschul@sandia.gov'
     $  /1x, 72('=') )
c
      write(IWR,*   ) ' Maximum dimension parameters for this run are:'
 9011 format(3x, a,i10 )
 9012 format(3x, a,i10, 5x, a,i10 )
 9010 format(1x, 72('-') )
      write(IWR,*   ) ' '
      write(IWR,9011) 'maxwkd (r8 workspace)   =', maxwkd
      write(IWR,9011) 'maxr1d (grid dimension) =', maxr1d
      write(IWR,*   ) ' '
      write(IWR,9011) 'norbd (basis fcns) =', norbd
      write(IWR,9012) 'natmd (atoms)      =', natmd
     $              , 'ntypd (atom types) =', ntypd
      write(IWR,9012) 'nshld (shell/atm)  =', nshld
     $              , 'noad  (bf/atm)     =', noad
      write(IWR,9012) 'nald  (gauss/shell)=', nald
     $              , 'nrd   (rad pts/atm)=', nrd
      write(IWR,9012) 'nkd   (kpts/ibz)   =', nkd
     $              , 'nsymd (sym/group)  =', nsymd
      write(IWR,9011) 'nimgd (NEB images) =', nimgd
      write(IWR,*   ) ' '
      write(IWR,9011) 'nlatd (lattice images)=',nlatd
      write(IWR,9011) 'neq3d (do not ask) =',neq3d
c
      call IOGETRECL( mxreclen)
      write(IWR,9011) 'Default record length (r8) for i/o=',mxreclen
      write(IWR,9010)
c
c Maximimum workspace in wkrad/wkrbx
      maxwkr = MAX( mxwkrbx, mxwkrad )
c Need to pass dimension, but do not like passing parameter, hence ...
      nang = nangd
c
c  Set up constants to be carried in common: 18Dec92-PAS
c     pi = 3.14159265358979d0
      pi = 4.d0 * ATAN( one )
      rtpi = SQRT( pi )
      call EXPUND( cutexp )
c
cACP variables that must be initialized to zero to avoid compiler problems.
      call MKZERO( 9,str1 )
      call MKZERO( 9,str2 )
      call MKZERO( 9,str3 )
c
      call TIMER('Before input     ')
      write(IWR,*) ' '
c
c  Format for diagnostics:
 9020 format(1x,a,3f20.10)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Read input data
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Open input data file 'lcao.in' (with reserved std input unit IRD):
      call FLGETIRD( IRD )
      call FLOPENU( IRD, 'lcao.in', 'OLD', 'FOR' )
      IDAT = IRD
c      call FLOPENA( IDAT, 'in' )
c
c Open status/restart (ASCII) file:
      call FLOPENA( istatfl, 'stat' )
c
c Set default to spinless lda (Ceperley/Alder + Perdew/Zunger)
      call DFTSET( 'lda   ' )
c
c Get job options ...
c
c  ... set job options defaults:
      call DAT0OPT( lvlout, mkrhopt, ion_opt, kparopt,
     $ dotests, dosetup,doiters,doforce,dorelax,do_cell, do_pvcalc,
     $ opt_spin, do_md, do_neb, do_post,
     $ do_bands,bvl, do_optic, do_nonscf,
     $ dosplit,doblas3, redusym, dolocxc,
     $ dospars,doeigen,doprjct,dorspsy,dodefct,
     $ do_kppsolve,do_psolv )
c
c  ... and get job options from input:
      call RUNOPT( IDAT,IWR, lvlout, mkrhopt, kparopt,
     $ dotests, dosetup,doiters,doforce,dorelax,do_cell, do_pvcalc,
     $ opt_spin, do_md, do_neb, do_post,
     $ do_bands, do_optic, do_nonscf,
     $ dosplit,doblas3, redusym, dolocxc,
     $ dospars,doeigen,doprjct,dorspsy,dodefct,
     $ do_kppsolve,do_psolv )
c
cpas: this is temporary kludge to longer-term dealing with geom engine
      do_geom = ( dorelax .or. do_cell .or. do_neb .or. do_md )
      madeges = 0
c
c Get input data for setup phase ...
c
      call SETDATA( IDAT,IWR, do_gga,do_spin,
     $ c6vdw,r0vdw, c6ijvdw,r0ijvdw,
     $ do_field, efield, dielec,
     $ dokdefct, nedefct,
     $ ndim, n1r,n2r,n3r, nearopt,near2c,
     $ iatmfmt, natmnm,typnm,atmnm,
     $ natm,natmd,ntyp,ntypd,nshld,nald, noad,
     $ itypa, numshl,lshel,nala,ala,cala,  occsh, norba, znuc,
     $ atmass,atengy,
     $ dolvatm, ratm, rprim,dorig,rscale,
     $ nk,nkd, veckg,wtk,wtkscal,nku,nku0,dokgrid,
     $ nsymd,nsymi,ndimsy, rsym,isymtyp,syminfo,
     $ spinpol, elchrg,rchrg,ion_opt,
     $ any1ctr, lmx1ctr,almnv1c,lmxnlp1,almnnl,
     $ nrd,nrad,nrps,radmsh,radwt,vpsrad,
     $ anycore, nrcor,corden,cordrv,
     $ anyxc1c,nfityp,alxcmin,nalfxc,nrxcfit,nalfxcd,alfxc,alsq,ipvt,
     $ wk(1) )
c
c  ... do some arithmetic needed to process run phase input:
      call INSETUP( IWR, dokdefct,
     $ natm,ntyp, norb,norbd,  nstbulk,nedefct,
     $ elchrg, nocc0k, do_spin,opt_spin,spinpol, nspin,elecno,spinfac,
     $ itypa, ioffa, norba, znuc )
c
c Get data controlling run/scf phase and beyond ...
c
c ... set up coordinate origins:
      call SETORIG( dorig,orig,rsym, origws,rchrg,
     $ ndim, dolvatm,rprim,rscale )
c
      if( Ldevtest )then
        write(IWR,1341) 'DBG: postseto, rchrg=',rchrg
 1341   format(1x,a,6f10.4)
      endif
c
c  ... set some defaults for run phase:
      call DAT0RUN( doforce,do_post, dokdefct,icloseocc,
     $ itstart,itstop,meth_bl,nhiste,igesopt,  nocc0k,norb,nstate,
     $ convgr,convs,convsl,convii, conv1c, alfast,
     $ opt_spin, meth_sp,conv_sp,nstep_sp,blend_sp, istart_sp,
     $ edegen,etemp, scfblnd,scfbl2,scfconv, idosolv, scut )
c
c  ... set some defaults for geometry relaxation:
      call DAT0GEOM( do_neb,
     $ natm,nafrc,iforce,ifdefct,
     $ nhistg,gblend,gconv,tstep,igstart,igstop,
     $ nslics,    iatframe )
c
c  ... set defaults for molecular dynamics run
      if( do_md ) call DAT0QMD
c
      if( do_cell )then
c      ... and defaults for cell-opt phase:
        call DAT0CELL(
     $   do_cell, ndim, iuctype, icellscheme,strsbroy,
     $   iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $   strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $   strlast )
      endif
c
c  ... and get run-phase options from input file:
      call RUNDATA( IDAT,IWR,  doflag,
     $ itstart,itstop,meth_bl,nhiste,igesopt,
     $ convgr,convs,convsl,convii, conv1c, alfast,
     $ icloseocc,
     $ edegen,etemp, scfblnd,scfbl2,scfconv, nstate,norb,nocc0k,
     $ opt_spin, meth_sp,conv_sp,nstep_sp,blend_sp,
     $ dorelax, natm,nafrc,iforce,ifdefct,
     $          atmnm,natmnm,
     $          nslics,islics,islicats,vslics,
     $          iatframe,
     $          nhistg,gblend,gconv,tstep,igstart,igstop,
     $ do_cell, do_pvcalc, ndim, iuctype, icellscheme,strsbroy,
     $          iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $          strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $          strlast,
     $ idosolv, scut,
     $ nbranchd,nbranch,bvl,bsymb )
c
c     Recap all slic constraints
      call MTBUFF
      call SLICSUM( IWR, natm,
     $ nslics,islics,islicats,vslics,    iatframe )
      call MTBUFF
c
c  ... get NEB information, if we have an NEB calculation:
c
c     Non-NEB calculation gets named image #0
      image = 0
      if( do_neb )then
cpas:xxx   should run a status/sync step here?
c
        nimg_max = nimgd
        call NEBSETUP( IDAT,IWR,IWRNEB,  do_neb, image, igstart,
     $   nimg_neb, nimg_max, ea_neb,eb_neb, spring_neb,
     $   Lantikink_neb,Loptall_neb,Lclimb_neb,Lnoninteract,
     $   ndim, dolvatm, rprim,orig,rscale,
     $   natm, imgstat, engyneb, ratm,ratma,ratmb, wksma,
     $   maxwkd, wk(1), igrouplist_neb, ngroup_neb )
c
cpas::xxx
        igstep = igstart
        igeom = 1
        igstat = -1
        call GSTAMP( IWR, igstep,igeom,igstat,image,
     $   iatmfmt, ntyp,natm, natmnm, itypa, atmnm,typnm, ratm )
c
      endif
c
c  A good time to print out the local processor allocation
      call MPNODES( nodes )
      if( nodes .gt. 1 ) call MPSTAMP( IWR )
c
c I will prevent this calculation going forward if there are
c too many nodes.  This restriction may ultimately be lifted,
c but at this point, this is not a constraint because the
c code does not scale beyond at most 16 task-parallel nodes.
      if( nodes.gt.maxnodes )then
        write(IWR,'(a,i8)') 'Max # core dimension, maxnodes=',maxnodes
        call STOPXERR( 'too many tp nodes' )
      endif
c
c  ... and update with status in restart file:
      call RDSTAT( IWR, istatfl,
     $ dosetup, doiters, doforce, do_geom, do_cell, opt_spin,
     $ iucstart,igstart,itstart, engytotl,excengy,esnsengy,
     $ istart_sp, elecno, spindata,
     $ rprim, natm, ratm, frctot,strtot )
cpas: force restart as first step cell:
      iucstart = 1
c
cqmd:start - where to put MD initialization?
      if( do_md )then
c       Initialize the q-md module
        call QMDINIT( IWR,natm,itypa,atmass,ratm,igstart )
c       Determine total number of gsteps in MD
        call QMDSTEPSGET( igstop )
      endif
cqmd:end
c
c Tell broyden/elec that this is initial transit
      itbroye = 0
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Check for incompatible options, and obvious bad things
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c   non-scf (bands or optical) should be alone, no force or relax
      if( do_nonscf )then
        if( doforce .or. do_geom )then
          write(IWR,*) 'bands/optical incompatible with'//
     &                 ' force/relax/cell/neb/md'
          call STOPXERR( 'bands/optical + any relax not compatible' )
        endif
      endif
c
      if( do_md )then
        if( do_cell .or. do_neb )then
          write(IWR,*) 'Dynamics incompatible with relax/cell/neb'
          call STOPXERR( 'MD+relax/cell/neb not compatible' )
        endif
      endif
c
      if( do_field )then
        do  idim=1,ndim
          edotlat = efield(1)*rprim(1,idim) + efield(2)*rprim(2,idim)
     $            + efield(3)*rprim(3,idim)
          if( ABS( edotlat ) .gt. 1.d-6 )then
            write(IWR,*) 'E-field cannot be in periodic direction'
            call STOPXERR( 'E-field in periodic direction' )
          endif
        enddo
      endif
c
      if( do_cell )then
        if( .not. dolvatm )then
c         Reconciling atomic coordinates can be a problem.
C     $   call STOPXERR( 'Cell+lat - lattice coords for do cell' )
          do  idim=1,ndim
            if( ABS( rscale(idim) - one ) .gt. 1.d-8 )
     $       call STOPXERR( 'scale+cell - incompatible' )
          enddo
        endif
        if( do_neb ) call STOPXERR( 'cell+NEB - cannot do both' )
        if( do_gga .and. do_spin )
C     $   call STOPXERR( 'spgga+cell - spin+GGA stress not coded' )
     $   write(IWR,*) '>>>>> spin+gga stress done by A.E. Mattsson'
        if( elchrg .ne. zero )
     $   call STOPXERR( 'cell+charge - cannot do both' )
      endif
c
c Non-standard dynamic 1-center code is strictly spinless lda
      if( any1ctr .and. (do_gga .or. do_spin) )then
c       Would need to gga/spin-ify v1ctrij and blkovfn.
        call STOPXERR( '1ctr+gga/cannot do 1ctr AND (gga.or.spin)' )
      endif
c
c Check that flat density sphere countercharge only with molecule:
      if( IABS(ion_opt).eq.1.and.elchrg.ne.zero.and.ndim.ne.0 )then
c       Flat spherical countercharge can only work for molecule
        call STOPXERR( 'ion_opt=1/not allowed' )
      endif
c
c Check 1D grid dimensioning:
      mxnr = MAX( n1r,n2r,n3r )
      if( mxnr .gt. maxr1d )then
c       fft arrays will fail
        write(IWR,*) '>>>>> grid dim=',mxnr,' exceeds maxr1d =',maxr1d
        call STOPXERR( 'maxr1d  /grid 1d dim too large' )
      endif
c
c Leave records of (settable) states of calculation
      call IOGETRECL( mxreclen )
      write(IWR,*) 'Maximum record length (r8) for i/o=',mxreclen
      call DFTGETTYP( fcnal )
      write(IWR,*) 'Flavor of dft that will be used for this run=',fcnal
      if( DO_VDW() )then
c       Leave a record of what is being used for vdw
        call VDWGETNAME( fcnal )
        write(IWR,*) 'VDW C6 corrections are ON, method= ',fcnal
c       Dump out the full ff (pair) parammeters, as they are never quoted:
        if( lvlout.gt.0 )
     $  call VDWSTAMP( IWR, ntyp, c6vdw,r0vdw, c6ijvdw,r0ijvdw )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                Assign io-units and open working files
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call MPNODE0( master )
      call MPNODE( iproc )
c
c >>>>> Files generated during setup and used in run phase:
c
c  iset0fl = data from setup for use in iteration phase
      call FLOPENB( iset0fl, 'set0' )
c
      if( iproc .eq. master )then
c
c  igrd0fl = grid sums over spherical atoms
        call FLOPENB( igrd0fl, 'grd0' )
c  idrhfl0 = derivatives over reference density (sum of spherical atoms)
        if( do_gga) call FLOPENB( idrhfl0, 'drh0' )
c  iovlpfl = overlap matrix
        call FLOPENB( iovlpfl, 's'    )
c  iham0fl = iteration independent part of Hamiltonian
        call FLOPENB( iham0fl, 'h0'   )
c  ivxc0fl = iteration-independent reference xc potential of Ham
        call FLOPENB( ivxc0fl, 'h0vxc'  )
c  iexc0fl = iteration-independent reference xc energy-density matrix Ham
        call FLOPENB( iexc0fl, 'h0exc'  )
c  ih0tfl  = iteration-independent kinetic energy part of Ham
        call FLOPENB( ih0tfl , 'h0t'  )
c  ih0nlfl  = iteration-independent non-local pseudopotential part of Ham
        call FLOPENB( ih0nlfl , 'h0nl'  )
c
      endif
c
c >>>>> Run phase files:
c
      if( iproc .eq. master )then
c
c  igridfl = iteration grid data
        call FLOPENB( igridfl, 'grd'  )
c  idrhofl = density derivatives
        if( do_gga ) call FLOPENB( idrhofl, 'drho' )
c  idrhsfl = additional derivatives needed for spin
        if( do_gga .and. do_spin ) call FLOPENB( idrhsfl, 'drhs' )
c  igrhofl = density gradients
        if( do_gga ) call FLOPENB( igrhofl, 'grho' )
c
c  idhamfl = running delta-Ham file (h0+dh=Htot)
        call FLOPENB( idhamfl, 'dh'   )
c  ivecfl  = eigenvector file
        call FLOPENB( ivecfl , 'vec'  )
c  idmatfl,iematfl = density matrix, energy-weight dmat
        call FLOPENB( idmatfl, 'dmat' )
        call FLOPENB( iematfl, 'emat'  )
        if( do_spin )then
c         Spin-resolved matrices
          call FLOPENB( idmatsfl, 'dmats' )
          call FLOPENB( iematsfl, 'emats'  )
        endif
c
      endif
c1c#####################################################################
c Files for non-standard scheme for "fast" 1-center density and pots
      if( any1ctr )then
c  i1ctrfl = (NON-STD) 1-ctr fast fcns on grid
        call FLOPENB( i1ctrfl, '1ctr' )
c  ilocfl  = (NON-STD) fast atom radial density reference
        call FLOPENB( ilocfl , 'loc'  )
c  i1cgrfl = (NON-STD) scf fast xc-potentials on grid
        call FLOPENB( i1cgrfl, '1cgr' )
      endif
c1c#####################################################################
c
c  igrdref = reference system delta-density for defect-LMCC:
c  Check for presence of grid file from reference system (for LMCC):
      igrdref = 0
      if( ( ( ndim.eq.3 .and. elchrg.ne.zero ) .and. ion_opt.ne.0 )
     $    .and. ( iproc .eq. master ) )then
c       Attempt to get LMCC reference rho/pot if charged bulk defect
        filenm = 'lmcc.grdref'
        call FLEXIST( iflexist, filenm , 'UNF' )
        if( iflexist .le. 0 )then
c         Try old convetional name ...
          filenm = 'lcao.grdref'
          call FLEXIST( iflexist, filenm , 'UNF' )
        endif
        if( iflexist.gt. 0 )then
          write(IWR,*) 'Found LMCC reference density file=',filenm
          call FLOPENF( igrdref, filenm, 'OLD', 'UNF' )
        else
          write(IWR,*) '***** WARNING: no LMCC reference density found'
        endif
      endif
c
c Open history file for post-processing (master-only):
      ihistfl = -1
      if( iproc .eq. master ) call FLOPENA( ihistfl, 'hist' )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Construct default integral cutoff parameters:
      cutgaus = cutexp
      cutgrid = -LOG( convgr )
c
c  Check/setup data for xc1c local fit development code:
      if( anyxc1c )then
        call DAT0ANG( nang, angwts,angpts )
      else
        ilocfl = 0
      endif
c
c Set number of total points in grid field:
      nptr = n1r*n2r*n3r
c Set boundary potential
      ibndpot = ndim
c Partition array wft() for fft workspace:
c  need (4*n1r+15);(4*n2r+15);(4*n3r+15);(2*max(n1r,n2r,n3r))
      ift1 = 1
      ift2 = ift1 + 4*n1r + 15
      ift3 = ift2 + 4*n2r + 15
      ifta = ift3 + 4*n3r + 15
      call CFFTI( n1r, wft(ift1) )
      call CFFTI( n2r, wft(ift2) )
      call CFFTI( n3r, wft(ift3) )
c
c Break up unit cell of grid points into boxes:
c
      call DEFBOX( IWR, ndim, n1r,n2r,n3r, nbox,mxrbox, ndbox )
c
      if( 2*nbox.gt.nrecd )then
        write(IWR,*) '>>>>> nbox =',nbox,', but only have nrecd=',nrecd
        call STOPXERR( 'nrecd-up/not enough box records (nrecd)' )
      endif
c
c     Number of points in biggest box
      mxrb = mxrbox
c
c Set up basis-shell relationship:
c
      call BFSETUP( IWR, natm,ntyp,nshld,nald,
     $ norb,norbd,
     $ numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ znuc,occsh,occij,alamin )
c
      if( dosplit )then
c       Apply split-basis transformation:
        call BFSPLIT( ntyp,nshld,nald, wk(1),occsh,occij,
     $   numshl,lshel,nala,ala,cala )
cpas:   Transformed basis should be written to the output file
        write(IWR,*) '***** Note: basis transformed to split basis form'
      endif
c
c Build full symmetry group specified in input ...
c
      call SYMBUILD( IWR, nsym, nsymi,nsymd,  ndimsy,
     $ isymtyp,syminfo,  rprim, rmatsym )
c
c      ... and preset size of reduced symmetry groups to out-of-bounds
      nsyma = -1
      nsymg = -1
c
      call QHIST0( ihistfl, 'INIT',
     $ ndim,natm,ntyp,itypa, natmnm,typnm,atmnm )
c
      if( lvlout.gt.2 ) call TIMER('After data       ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    Cell step returns here
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  itconv = have potential from previous converged geom (0=n,>0=y)?
      itconv = 0
c
      iucstep = 1
      iuciter = iucstart
      iucstart = 1
  100 continue
c
      if( do_cell ) call QHIST( ihistfl, 'UCST', iucstep, rdata, wksm )
      call QHIST( ihistfl, 'CELL', ndata, rdata, rprim )
c
c  Calculate: gprim= recip. lattice vectors, hh=mesh intervals vectors
c
      call MESH( n1r,n2r,n3r,rprim, weight,hh, gprim )
c
c  Calculate bulk polarization energy in response to charge:
c
      engypol = zero
      if( ndim.eq.3 .and. ABS( elchrg ).ne.zero )then
        call BULKPOL( IWR, ndim,elchrg,dielec,rprim, engypol )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    Geometry step returns here
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  igeom = total number of geometry steps performed in this run
c  igstep = step number in current broyden blend geometry sequence
c
      igeom = 0
      igstep = igstart
      igstart = 1
  200 continue
      igeom = igeom + 1
c
c  Generate internal coordinates from input atomic positions ...
      call RINTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $ ratm,ratm0 )
c   ... and project into primary unit cell if needed
      call UCATCHK( IWR, ndim,natm, dolvatm, itypa,znuc,
     $ rprim,orig,rscale,  ratm0,ratm )
c
      if( do_neb )  call QHIST( ihistfl, 'IMAG', image , rdata, wksm )
      if( do_geom ) call QHIST( ihistfl, 'GSTE', igstep, rdata, wksm )
      call QHIST( ihistfl, 'COOR', natm , rdata, ratm )
c
c  Generate internal coords of LMCC stuff:
      if( ion_opt.ne.0 )then
c       Generate internal coords of wigner-seitz origin ...
        call RINTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $   origws,origws0 )
c       ... and project into primary unit cell if needed
        call UCATCHK( IWR, ndim,  1 , dolvatm, itypa,znuc,
     $   rprim,orig,rscale,  origws0,origws )
c
c       Convert and check neutralizing charge coordinates ...
c
      if( Ldevtest )then
        write(IWR,1341) 'DBG: pre-Rint, rchrg=',rchrg
        write(IWR,1341) 'DBG: pre-Rint, rchrg0=',rchrg0
      endif
        call RINTRNL( ndim, nchrgd, dolvatm, rprim,orig,rscale,
     $   rchrg,rchrg0 )
c        ... and project into primary unit cell if needed
c
      if( Ldevtest )then
        write(IWR,1341) 'DBG: postRint, rchrg=',rchrg
        write(IWR,1341) 'DBG: postRint, rchrg0=',rchrg0
      endif
        call UCATCHK( IWR, ndim, 1     , dolvatm, itypa,znuc,
     $   rprim,orig,rscale,  rchrg0,rchrg )
c       This is a less-than-optimal means to do convert,
c       but must be this way until fullness of b.c. stuff is defined
CX        do  i=4,3*nchrgd
        do  i=1,3*nchrgd
          rchrg(i) = rchrg0(i)
        enddo
      endif
c
      if( Ldevtest )then
        write(IWR,1341) 'DBG: postsetg, rchrg=',rchrg
        write(IWR,1341) 'DBG: postsetg, rchrg0=',rchrg0
      endif
c
c  Find/check atom images under symmetry operations
c  Reduce symmetry from input default if "redusym" is set to true.
c
      nsyma = nsym
      call SYMAGES( IWR, redusym, rsym,rprim, natm,ratm,itypa,
     $ natmnm,atmnm,
     $ nsym ,nsyma, naofsym, rmatsym )
c
c  Set number of symmetries to that used by the system:
      if( nsyma .ne. nsym )then
        write(IWR,*)  'Size of atomic symmetry group (no id)=',nsyma
        nsym = nsyma
      endif
c
c  Check WS origin and countercharge positions for symmetry:
      if( ion_opt.ne.0 )then
        do  jr=1,3
          origws0(jr) = origws(jr)
        enddo
        call DCOPY( 3*nchrgd, rchrg,1, rchrg0,1 )
        call SYMVEC( 1,origws, rsym, nsym, rmatsym )
        call SYMVEC( nchrgd,rchrg , rsym, nsym, rmatsym )
        delo = zero
        delr = zero
        do  jr=1,3
          delo = delo + (origws(jr)-origws0(jr))**2
          delr = delr + (rchrg(jr)-rchrg0(jr))**2
        enddo
        if( delo .gt. 1.d-8 )then
          write(IWR,*) '***** WARNING: WS location symmetry looks wrong'
        endif
        if( elchrg.ne.zero .and. delr .gt. 1.d-8 )then
          write(IWR,*) '***** WARNING: charge location symtry looks bad'
c          call STOPXERR( 'lmcc+sym/charge or WS position symmetry off' )
        endif
      endif
c
      if( Ldevtest )then
        write(IWR,1341) 'DBG: postseti, rchrg=',rchrg
        write(IWR,1341) 'DBG: postseti, rchrg0=',rchrg0
      endif
c
c  Next we check the mesh symmetry by symmetrizing dummy grid field.
c  Reduce symmetry from atomic symmetry if redusym is set to true.
c  Set input number of symmetries to number of atomic symmetries,
c  unless grid symmetry already reduced and ordered:
c
      nsymg0 = nsyma
      if( nsymg.ge.0 ) nsymg0 = nsymg
c
      i01p = 1
      call MKZERO( nptr, wk(i01p) )
c
      call SYMRHO( IWR, redusym,ndimsy, nsymg0,nsymg,
     $ rsym, rmatsym, isymgrd, rprim, n1r,n2r,n3r, hh,
     $ wk(i01p), wk(i01p+nptr) )
c -->  rho-i    iwrk-2s
c
      if( nsymg .ne. nsymg0 )then
c       The symmetry of the grid is reduced from that of the atoms
c       and need to reorder some arrays the first time through:
        call SYMORDR( natm,nsyma,nsymg, isymgrd,rmatsym,naofsym )
      endif
c
c     And now the symmetry should *never* reduce again:
      redusym = .false.
c
c Generate Bloch vector grid and fold into IBZ:
c
      if( igeom.eq.1 )then
        ikdefct = 0
        if( dokgrid )then
c         Only if we have a grid to generate, and only if we have not
c         already generated a kgrid in previous geom for same cell.
c          ... generate raw Bloch vector grid:
          call KGRID( ndim, rprim,
     $     nku0,nku(1),nku(2),nku(3),nkfull,wk(1) )
c          ... fold into IBZ using symmetry:
          call SYMIBZ( ndim, nkfull,wk(1),wk(1+3*nkfull),
     $     nkd,nk,veck,wtk,wtkscal,
     $     rprim, nsym, rmatsym )
c         Convert to lattice (fractional) units:
          call KLATT( ndim, rprim, nk,veck, veckg )
c         For defect k-sampling, make sure that
c         gamma-pt is in the list:
          if( dokdefct )then
            call KGAMMA( IWR, ikgamma, nk,nkd, wtk, veck,veckg )
            ikdefct = ikgamma
          endif
        else
c         Set total number of kgrid to actual number used:
          nkfull = nk
        endif
c       Write the k-vectors to be used to the output file:
        call KWRITE( IWR, ndim, rprim, nk,nkfull,wtk,wtkscal,veckg )
      endif
c     Convert k-points to Cartesian units:
      call KCART( ndim, rprim, nk,veckg, veck )
c
c  This sets whether I have real or complex basis orbitals
      ncplx = 2
      if( ndim.eq.0 .or. ( nk.eq.1 .and.
     $  (veck(1,1)**2+veck(2,1)**2+veck(3,1)**2).eq.zero ) ) ncplx = 1
c
c
      if( igesopt.lt.0 )then
c       ges not set by user, choose default ges option
        if( ncplx .eq. 1 )then
c         Gamma-pt works, turn ges on
          igesopt = IABS( igesopt )
          write(IWR,*) 'SCF guessing turned on, option=',igesopt
        else
c         Complex fails if atom crosses cell boundary, turn guess off
          igesopt = 0
          write(IWR,*) 'SCF guessing turned off for complex'
        endif
      endif
c
c  This decides whether to promote this to k-parallel, and sets up the module
      call KPSETUP( nk, kparopt )
c
      call KPFLAG_GET( kparopt )
      if( kparopt.ne.2 )then
c       Not k-parallel, disable parallel-parallel solver
        do_kppsolve = .false.
      else 
c      See if KPSETUP gave us a processor distribution that benefits
c      from using the parallel-parallel solver.
        call MPNODES( nodes )
        call MPNODES_K( nodes_k )
c       If all processors are in the same k-point communicator, then
c       there is no point to using the k-parallel-parallel solver.
        if( nodes_k.eq.nodes ) do_kppsolve = .false.
cxxx264  Unless there is at least 4 procs in every parallel solve, turn off kpp
cxxx264  (if proc allocations diff from kp-allocs, can get different kpp criterion)
        if( (nodes/nk) .lt. 4 ) do_kppsolve = .false.
c       Turn off kpp, if problems are too small (adapting acp criterion)
C        if( norb .lt. 200*(nodes/nk+1) - acp criterion, seems mysterious
        if( norb .lt. 800 ) do_kppsolve = .false.
      endif
c
c Compute, filter, and order lattice vectors in real space:
c
      call MKCUTS( convs(1), convsl, convii )
      call LATVEC( ndim, nlatd,nlat,nlat1c,nlat2c,
     $ natm,ntyp,nshld,nald, itypa, numshl,nala,ala,
     $ ratm,rprim, rlat,wksml(1,1),wksml(1,2) )
c
      if( ndim.eq.0 ) nlat = 1
c
c  Calculate range of grid orbitals, determine which boxes they
c  have non-negligible amplitude in, and number of orbs/box.
c
      maxrec = nrecd/nbox
      call BOXRNG( IWR,lvlout, maxrec,          mngr3,mxgr3,ngr,
     $ ndim, n1r,n2r,n3r, natm,ntyp, nshld,nald,
     $ itypa,numshl,lshel,nala,norba,ala, alamin, ratm,rprim,
     $ hh, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1), boxrad,
     $ wk(i01p) )
c-->   rmesh-4bs
c
c  Compute "npopdm", the number of occupied orbitals at which a
c  density-matrix scheme for computing the grid density becomes
c  cheaper than an eigenfunction scheme (strictly operations ct).
c
      call EF2DM( IWR,lvlout, doblas3, nocc0k,npopdm, nbox,inboxs )
c
c  Set the scheme for computing grid density (DM=dmat or EF=eigenfcn)
      if( mkrhopt .eq. 0 )then
c       Default to DM ...
        irhoopt = 1
c        ... and switch to EF scheme if few enough occupied orbitals:
        if( nocc0k .lt. npopdm ) irhoopt = 2
      elseif( mkrhopt .eq. 1 .or. mkrhopt .eq. 2 )then
c       User has forced one of the schemes:
        irhoopt = mkrhopt
      else
        call STOPXERR( 'mkrhopt must be 0,1, or 2 (gridrho scheme)' )
      endif
      if( lvlout.gt.2 )
     $ write(IWR,*) 'Grid density scheme (1=DM;2=EF)=',irhoopt
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c             Carve up temp space in big work array
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Determine arrays sizes for this problem:
      call KPMINE( nk, nkloc, nk0 )
      mat = norb*norb
      nmat = nk*mat
      nmatkp = nkloc*mat
      nmats = nmat*nspin
      nmatskp = nmatkp*nspin
c
c  The routine WKMEM assigns all space in the big work array wk()
c  This routine should be consulted before attempting to use any
c  space in wk().
c
      i00 = 1
c
c     These flags must be set up for scf-phase of the calculation here
c      (note coexisting do_nonscf flag, deciding if nonscf afterwards)
      do_eigvecs = .true.
      Lnonscf = .false.
c
      call WKMEM( IWR, lvlout, maxwkd,
     $ do_gga,do_spin,doblas3, Lijkmat,
     $ idosolv, do_kppsolve,do_psolv, do_eigvecs, Lnonscf,
     $ mat,nmat,nmatkp, nptr, nhiste,
     $ ncplx,nk,nlat, norb,nstate, noad,natm,ntyp,itypa,norba,
     $ mxang, ndbox,nbox,inboxs,mxrb,
     $ near2c,irhoopt, idwf,  memwf,memwfw, iwf,iwfw,
     $ icoskr,isinkr, is1at,
     $ ioint,ioofk, ieigval,ieigpop,
     $ idbr, imblnde, memblnde,
     $ i00, i01s, i01,i02,i03,i04,
     $ i05,i06,i07,i08,i09,i10,i11,i12 )
c
c  Initialize checks for memoery overrun into blend history data:
      call MEMCHKIBR( -1, 0, 'Initialize', wk(imblnde), savbr )
c
c Compute Bloch phase factors:
c
      call MKEIKL( nk,nlat, veck,rlat, wk(icoskr),wk(isinkr) )
c
c Locate nearest mesh point to each atom in list
c
      call NEARBY( nearatm, natm, itypa, norba,
     $ ratm, rprim,hh,n1r,n2r,n3r,nptr, mngr3,mxgr3,ngr,
     $ wk(i01s), wk(i01s+3*nptr) )
c -->  r-3s      rsq-s
c
c Write out geometry/basis parameters to disk file:
c
      if( do_post )then
        nsv1 = nshld*ntyp
        nsv2 = nald*nshld*ntyp
        call DATSAVE(
     $   norb, natm,ntyp,nshld,nald,noad,
C     $   itypa, numshl,lshel,nala,ala,cala, znuc,occsh,
     $   itypa, nocc,locc,nalsh,alsh,calsh, znuc,occsh,
     $   ratm, rprim, nlat,rlat, hh, n1r,n2r,n3r,
     $   orig, rscale,
     $   wk(i01s),wk(i01s+nsv1),wk(i01s+2*nsv1),wk(i01s+2*nsv1+nsv2) )
c -->    lshl2-s  nal2-s        al2-s           cal2-s
      endif
c
      if( lvlout.gt.2 ) call TIMER('Aft-preliminaries')
c
      if( dotests )then
        call STOPXOK( 'tests-ok/system passes input and memory tests' )
        goto 999
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * *  Start setup routines/big setup phase if-block  * * * * *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( dosetup )then
c
c  Clear out reference atom xc integrals
      call MKZERO( 12, xc0val )
c
c  Construct, for use in local potential matrix routines, *local* atomic
c  potentials made up of the nuclear coulomb attraction minus Zval/r,
c  plus cutoff atom xc potential, on same radial mesh used for the
c  non-local pseudopotential.  Also compute Exc(r), xc energy density.
c
      call VLOCPOT( IWR, do_gga, xcalcut,xcfac,xcrhocut,
     $ exc0rad,exc0slo,exccrad, strxc0,strxcc,  rhoatom,
     $ wksma(1,1),wksma(1,2),wksma(1,3),wksma(1,4),wksma(1,5),
     $ natm,ntyp,nshld,nald,
     $ itypa, nocc,locc,nalsh,alsh,calsh, znuc,occsh,lmxnlp1,
     $ nrad,nrd, radmsh,radwt, nrcor, corden,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad,
     $ maxwkr, wkrad, wspl(1,1), wspl(1,2) )
c
c Store some reference atom (semi-)analytic xc terms:
      xc0val(2) = exc0rad
      xc0val(6) = exccrad
      xc0val(8) = exc0slo
      xc0val(10) = strxc0
      xc0val(11) = strxcc
c
      if( lvlout.gt.2 ) call TIMER('After vlocpot    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Compute iteration-independent (semi-)analytic matrix elements
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Tell "mat0set" how much mem is available in wk():
      maxwk = maxwkd - i01s + 1
c
      if( Lijkmat )then
c       We use the square (i,j,k) matrices
      call MAT0SET1( IWR, lvlout,
     $ iovlpfl,iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $ atm0engy, mat,nmat, convs,convsl,convii,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad, wspl,
     $ ratm, rlat, wk(icoskr),wk(isinkr), wk(is1at),
     $ wksm, gntcomb, maxwk,mxang,
     $ wk(i01s) )
c -->  wkmat-s(see wkmem)
      else
c       We use the stripe-parallel matrix
      call MAT0SET( IWR, lvlout,
     $ iovlpfl,iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $ atm0engy, mat,nmat, convs,convsl,convii,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad, wspl,
     $ ratm, rlat, wk(icoskr),wk(isinkr), wk(is1at),
     $ wksm, gntcomb, maxwk,mxang,
     $ wk(i01s) )
c -->  wkmat-s(see wkmem)
      endif
c
      if( lvlout.gt.1 ) call TIMER('After mat0set    ')
c
      if( DO_VDW() )then
c       Compute energies (and forces/stresses) due to vdW corrections
        engyvdw = zero
c
        call VDWFF( engyvdw, frc1,str1, c6vdw,r0vdw,c6ijvdw,r0ijvdw,
     $   natm,ntyp,nlat,nlat2c, itypa,znuc,alamin,
     $   ratm,rlat )
c
        write(IWR,'(2x,a,f20.10)')  'vdW ENERGY(Ry)=', engyvdw
        atm0engy = atm0engy + engyvdw
        if( lvlout.gt.2 ) call TIMER( 'After vdW term' )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         Do grid part of iteration-independent calculations
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Let routine know how much memory is available in wk():
      maxwk = maxwkd - i01s + 1
c
      call GRID0SET( IWR, lvlout, igrd0fl,idrhfl0,igrhofl,
     $ do_gga,do_spin, anycore,  xc0val,rhoatom,
     $ xcalcut,xcfac,xcrhocut,  maxwk,maxwkr,mxrb,
     $ natm,ntyp,nshld,nald, nlat, nlat1c, nptr, n1r,n2r,n3r,
     $ itypa, numshl,lshel,nala,ala,cala,
     $        nocc,locc,nalsh,alsh,calsh,  znuc, occsh,
     $ nrd,nrcor,radmsh,corden,
     $ ratm,rlat, nearatm, hh, weight,
     $ ndbox,nbox,nrec,inrec,inboxs,boxrad,
     $ wksma,wspl,wkrbx,
     $ wk(i01s) )
c -->  wk-5s(lda)/14s(gga)
c1c#####################################################################
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c         Non-standard "fast" dynamic 1-center scheme prep
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( any1ctr )then
c
c       Grid stuff for dynamic "fast" density scheme (non-standard)
c
        REWIND( unit=i1ctrfl )
        mxwk2 = maxwkd - (i01s+nptr-1)
        call SETGR1C( i1ctrfl,
     $   natm,ntyp,nshld,nald, n1r,n2r,n3r,nptr,nlat1c,mxwk2,
     $   itypa, numshl,lshel,nala,ala,cala, ratm,hh,rlat,
     $   lmx1ctr,almnv1c,  mngr3,mxgr3,ngr,
     $   wk(i01s),wk(i01s+10*nptr) )
c -->    wk-10s   wk-25s
c
        if( lvlout.gt.2 ) call TIMER('After setgr1c    ')
c
        if( anyxc1c )then
c
c         Dynamic 1-center XC radial atom setup (development code).
c
          REWIND( unit=ilocfl )
          call SETXC1C( ilocfl,
     $     natm,ntyp,nshld,nald, noad,nrd,nfitd, nlat1c, nang,
     $     itypa, numshl,lshel,nala,ala,cala,nocc,locc,nalsh,alsh,calsh,
     $     occsh,nfityp, nrad,radmsh, angpts,angylm, radwf,  ratm,rlat,
     $     wk(i01s),        wksmr )
c -->      wk(9*nrd*nang)-s wksm(2*nrd)-s
c
          if( lvlout.gt.2 ) call TIMER('After setxc1c    ')
c
        else
          ilocfl = 0
        endif
c
      endif
c1c#####################################################################
c
c Compute analytic parts of 2-ctr nearby construction defects:
c
      if( near2c.eq.2 )then
c
        call MKZERO( 2*norb*nk, wk(ioint) )
c
        call FTORB( wk(ioint),  norb,nk, natm,ntyp,nshld,nald,
     $   itypa, numshl,lshel,nala,ala,cala, norba, veck )
c
        call WF1COFK( wk(ioofk),  norb,nk, natm,ntyp,nshld,nald, nlat1c,
     $   itypa,numshl,lshel,nala,ala,cala,norba, ratm,rlat,
     $   wk(icoskr),wk(isinkr) )
c
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Iteration-independent Hamiltonian complete
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Save setup data for iteration phase
c
      iprgm = 242
      jprgm = 1
      nxcd = 12
      call WRSETD( iset0fl, iprgm,jprgm, nxcd,
     $ norb,nk, natm,ntyp,nrd, noad, nlat,nlat1c,nlat2c,nlat3c,
     $ atm0engy, xc0val,
     $ vatrad,vesrad,vxcrad,excrad, gntcomb,
     $ rhoatom,nearatm, near2c, wk(is1at), wk(ioint),wk(ioofk),
     $ nang,nfitd, angylm,nfityp,radwf )
c
c Record status as having completed setup phase:
      call WRSTAT( IWR, istatfl, 'SETUP ',
     $ do_geom, do_cell, opt_spin,
     $ iucstep,igstep,iterscf, engytotl,excengy,esnsengy,
     $ istart_sp, elecno, spindata,
     $ rprim, natm, ratm )
c
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * * * * * *  End big setup phase if-block * * * * * * * * * *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call TIMER('AFTER SETUP PHASE')
c
      write(IWR,*)
c
c *********************************************************************
c **********               end of setup phase               ***********
c *********************************************************************
c **********              start run/scf phase               ***********
c *********************************************************************
c
      if( .not. dosetup )then
c
c       Retrieve setup data for iteration phase:
c
        nxcd = 12
        call RDSETD( iset0fl, iprgm,jprgm, nxcd,
     $   norb,nk, natm,ntyp,nrd, noad, nlat,nlat1c,nlat2c,nlat3c,
     $   atm0engy, xc0val,
     $   vatrad,vesrad,vxcrad,excrad, gntcomb,
     $   rhoatom,nearatm, near2c, wk(is1at), wk(ioint),wk(ioofk),
     $   nang,nfitd, angylm,nfityp,radwf )
c
c       And on subsequent geoms, we will want to do setup calculation:
        dosetup = .true.
      endif
c
c Generate and keep grid orbitals if we have space:
c
      if( idwf.gt.0 )then
        call GRIDWF( ncplx,idwf,
     $   nk, natm,ntyp, nshld,nald, nlat1c,
     $   itypa,numshl,lshel,nala,ala,cala, alamin,ratm,rlat,
     $   wk(icoskr),wk(isinkr),
     $   hh, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1), boxrad,
     $   wk(iwf), wk(iwf),
c -->    wvfcns   wf
     $   wkrbx )
c -->    wkbox(10B)-s
c
        if( lvlout.gt.2 ) call TIMER('After gridwf-i   ')
c       Initialize memory check of grid orbital data
        call MEMCHKIWF( idwf, 0, 'gridwf', wk(iwf), wfiwf )
      endif
c
      if( .not. doiters )then
c       Restarting after scf, reset doiters for later geoms
        doiters = .true.
        itstart = 0
        itconv = 1
        goto 2000
      endif
c
c Subtract grid overlaps from analytic to get 1-ctr "defect" overlaps:
c
      call GRDOVLP( ncplx, idwf, near2c, wk(is1at),wk(ioint),
     $ norb,nk, natm,ntyp, nshld,nald, noad, nlat1c,
     $ itypa,numshl,lshel,nala,ala,cala, alamin,norba,
     $ ratm,rlat,wk(icoskr),wk(isinkr),  veck, wksmk,
     $ weight, hh, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1), boxrad,
     $ wk(iwf),wk(iwf), wkrbx(1),    wkrbx(1+10*mxrb),
c -->  wvfcns  wf       wkbox(10B)-s r(3B)-s
     $ wk(i01),      wk(i01+2*nk*mxrb) )
c -->  eikr(2B*nk)-s otmp(MAX(noad^2*nk*natm,2*norb*nk))-s
c
      if( lvlout.gt.2 ) call TIMER('After grdovlp    ')
      call MEMCHKIWF( idwf, 1, 'grdovlp   ', wk(iwf), wfiwf )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Iterate to self-consistency or itstop, whichever comes first.
c    Branch here according to how starting.
c      1> itstart=0: iterate from start of setup spherical atoms
c      2> itstart>0: continue blend from previous calculation
c      3> itstart<0: restart blend from previous calculation, treating
c                    input potential as first guess potential matrix
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      itdone = 0
      if( dorelax .and. igesopt.eq.1 ) itstart = -itconv
c
      if( itstart.eq.0 .or. itstart.gt.0 )then
c       Initial iteration, or restart as continued scf sequence ...
        iterscf = itstart
        itoff = itstart
      elseif( itstart.lt.0 )then
c       Restart with previous potential as first iter ...
        iterscf = 1
        itoff = 1
      else
        call STOPXERR( 'code err - impossible itstart value' )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     On iteration 0, need to complete Hamiltonian for overlapping
c     spherical reference grid matrix elements of long xc-potential
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Open the (image-specific) guess file ...
      call GESOPEN( image, iges0fl, igesfl, igesopt, nk )
c
c Spin optimization returns here
      istep_sp = istart_sp
 1100 continue
      istep_sp = istep_sp + 1
c
      if( iterscf .eq. 0 ) goto 1500
c
c  Restart scf cycle with potential from a previous run, on "idhamfl".
c  If itstart less than zero, restart blend procedure (as if this were
c  the first guess potential matrix).
c
C      if( ncplx.eq.1 .and. nk.eq.1 )then
      if( ncplx.eq.1 .and. nk.eq.1 .and. nprocs.eq.1 )then
cxxx:   This read of gamma-pt Ham is likely redundant, too, esp. w/spin
c       Retrieve old dHam from idhamfl if gamma pt (else take later)
        call MPNODE( iproc )
        call MPNODE0( master )
        call MPCOMM( icomm )
        if( iproc .eq. master )then
          REWIND( unit=idhamfl )
          call H0READ( idhamfl, mat,nk, wk(i01) )
C          if( do_spin ) call H0READ( idhamfl, mat,nk, wk(i01+nmat) )
        endif
Ccxxx:   The following bcasts of full gamma-pt Hams WAS surely redundent
C        call MPBCAST8( master, nmat, wk(i01), icomm )
C        if( do_spin ) call MPBCAST8( master, nmat, wk(i01+nmat), icomm )
      endif
c
      if( lvlout.gt.2 ) call TIMER('Before iterations')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Iteration cycles begin here
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 1200 continue
      write(IWR,9100)  itdone + itoff
 9100 format( 1x, 25('*'),'  iteration',i4,' begins here  ',25('*') / )
c
c Check if enough time left for another iteration
c
      if( itdone.eq.1 )then
        call GETTLFT( startit )
      elseif( itdone.eq.2 )then
        call GETTLFT( afterit )
        timeit = startit - afterit
      else
        call GETTLFT( timerem )
        if( timerem .lt. 0.95*timeit )then
          write(IWR,9110) timeit,timerem
 9110     format(1x,'Time/iteration=',f16.3/1x,'Time remaining=',f16.3)
          call STOPXOK( 'time out/not enough for another iteration' )
          goto 999
        endif
      endif
c
c  Given input delta-Hamiltonian matrix in "01", "XXXsolv" solves
c  the Schroedinger equation and sums to produce a density matrix
c
      if( iproc .eq. master )then
        REWIND( unit=idhamfl )
      endif
c
      call MEMCHKIWF( idwf, 1, 'b4 solvers', wk(iwf), wfiwf )
      call MEMCHKIBR( idbr, 1, 'b4 solvers', wk(imblnde), savbr )
c
      if( doeigen )then
c
        if( ncplx.eq.1 )then
c
c         Invoke eigensolver for real symmetric matrix
c
c         Get the space (AND dimensioning of the distributed matrices)
ctp:start
CTP          matdist = mat
          call PSOLVSIZES( -1,norb,norb,nlocrow,nloccol,lwork )
          matdist = nlocrow*nloccol
cxxx      If proc is not in pdsyg solver contect, nloc might by a problem?
ctp:end
c         Determine top of open memory here, avoiding wf/broyden space
c         Entire free space up to this used for WORK space in solver routine
          memceil  = maxwkd
          if( idwf .gt. 0 ) memceil = MIN( memceil, (iwf-1) )
          if( idbr .gt. 0 ) memceil = MIN( memceil, (imblnde-1) )
cxxx      This memory alloc for persolv is UGLY - needs to be cleaned up
cxxx      3 matdist data plus (WORKspace for solver | 4th matdist temp}
cxxx      Also, followed by tempio(mat), myst follow 4th matdist
          msolv = (i01-1) + 3*matdist
          lwork8 = memceil - msolv
          i01t = i01 + 4*matdist
c
          if( lvlout.gt.3 ) write(IWR,*) 'DSYG work mem=',lwork8
          if( lwork8 .lt. matdist )then
            write(IWR,*) '***** ERROR: lwork for dsyg too small'
            write(IWR,*) '>>>>> This is a coding error. lwork=',lwork8
            call STOPXERR( 'DSYG memory-prep code error' )
          endif
ctp:
CTP          call ERLSOLV(
          call PERLSOLV( icluster, gapp, nlocrow, nloccol,
     $     idosolv,scut, nspin,
     $     idhamfl,iham0fl,iovlpfl,
     $     ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $     icloseocc, doblas3,
     $     etemp,edegen, elecno, efermi,egap, esumsp,
     $     norb,nstate, wtk,  wk(ieigval),wk(ieigpop),wksmo(1,8),npop,
     $     wksmo(1,1),wksmo(1,7),lwork8,
c -->      iwork-is6o ifail-iso
     $     wk(i01),   wk(i01t) )
c -->      vmat-(4s)  tempio-(s)
c
          if( lvlout.gt.2 ) call TIMER('After erlsolv    ')
          call MEMCHKIWF( idwf, 1, ' erlsolv  ', wk(iwf), wfiwf )
          call MEMCHKIBR( idbr, 1, ' erlsolv  ', wk(imblnde), savbr )
c
        elseif( ncplx.eq.2 )then
c
c         Invoke general eigensolver for complex Hermitian matrices
c
          do_eigvecs = .true.
          do_eigpops = .true.
ctp:
CTP          call EIGSOLV(
          call PEIGSOLV( masterlist,do_kppsolve, do_psolv,
c263          call PEIGSOLV0( do_psolv,
     $     icluster, gapp,
     $     do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $     idhamfl,iham0fl,iovlpfl,
     $     ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $     icloseocc, dokdefct, doblas3, ikdefct,nstbulk,
     $     etemp,edegen, elecno, efermi,egap, esumsp,
     $     norb,nstate, nk, veck,wtk,
     $     wk(ieigval),wk(ieigpop), npop, wksmo,
     $     wk(i01) )
c -->      wmat-6s+
c
          if( lvlout.gt.2 ) call TIMER('After eigsolv    ')
          call MEMCHKIWF( idwf, 1, ' eigsolv  ', wk(iwf), wfiwf )
          call MEMCHKIBR( idbr, 1, ' eigsolv  ', wk(imblnde), savbr )
c
        endif
c
c       Put out summary of eigenspectrum
        lprnt = 0
        if( lvlout.gt.1 ) lprnt = lprnt + 1
        if( lvlout.gt.4 ) lprnt = lprnt + 1
        if( iterscf.eq.1 .or. iterscf.eq.itstart .or.
     $      iterscf.eq.itstop )then
          if( istep_sp.lt.2 ) lprnt = lprnt + 1
        endif
c
        call EVPRNT( IWR, lprnt, nspin, efermi,egap, npop,nocc0k,
     $   nstate,nk, wtk,wk(ieigval),wk(ieigpop) )
c
      else
c
c       Invoke alternate solver
c
        call STOPXERR( 'solver--unknown option' )
c
      endif
c
      if( lvlout.gt.2 )then
        write(IWR,*) ' '
        write(IWR,9020) 'Single-particle energy sum, esumsp= ', esumsp
        write(IWR,*) 'Number of populated states=',npop
      endif
c
      if( iterscf.eq.1 .and. itstart.eq.0 .and. lvlout.gt.2 )then
        eharris = eharris + esumsp + atm0engy
        write(IWR,9020) 'Harris Fcnal TOTAL ENERGY = ', eharris
      endif
c
c Evaluate density on regular grid:
c
      if( iproc .eq. master )then
        REWIND( unit=ivecfl )
        REWIND( unit=idmatfl )
        if( do_spin ) REWIND( unit=idmatsfl )
      endif
      idmfile = idmatfl
      i05i06 = i05
      jeigpop = ieigpop
      do  ispin=1,nspin
        call MKZERO( nptr, wk(i05i06) )
c
        if( ispin.eq.1 .or. elecno(ispin).ne.zero )then
          call GRIDRHO( ivecfl,idmfile,
     $     irhoopt,ncplx,idwf,doblas3,
     $     nstate,npop, norb,nk,
     $     natm,ntyp, nshld,nald, nlat1c,
     $     itypa,numshl,lshel,nala,ala,cala, alamin, ratm,rlat,
     $     wtk, wk(icoskr),wk(isinkr), wk(jeigpop),wksmo,
     $     weight,hh,ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1),boxrad,
     $     wk(iwf),wk(iwf), wkrbx(1),    wkrbx(1+mxrb), wk(iwfw),
c -->      wvfcns  wf       rhobox(1B)-s wkbox(10Br8)-s wkbi(wkmem)-s
     $     wk(i01), wk(i02),  wk(i05i06),  nptr,wk(i07) )
c -->      dmat-2s  dmatrd-2k rho-io            rho_tp-s/tp
c         Note: wk(i07) overlaps iwfw, but scheduled to not conflict
        endif
c
c        ... shift to down-spin electrons, if we have them:
        idmfile = idmatsfl
        i05i06 = i06
        jeigpop = ieigpop + nk*nstate
      enddo
c
      call MEMCHKIWF( idwf, 1, 'gridrho   ', wk(iwf), wfiwf )
      call MEMCHKIBR( idbr, 1, 'gridrho   ', wk(imblnde), savbr )
c
c  ... retrieve rho-zero
      if( iproc .eq. master ) REWIND( unit=igrd0fl )
      call MPBCREADBIG( igrd0fl, nptr, wk(i10) )
c
      i05i06 = i05
      i08i09 = i08
      do  ispin=1,nspin
c        ... move grid density out of boxes
        call FROMBOX( n1r,n2r,n3r, wk(i08i09), ndbox, wk(i05i06) )
c        ... subtract rho0 to get del-rho
        call DAXPY( nptr, -spinfac, wk(i10),1, wk(i08i09),1 )
c
c        ... symmetrize delta-density on real space mesh (in-place):
        call SYMRHO( IWR, redusym,ndimsy, nsymg ,nsymg,
     $   rsym, rmatsym, isymgrd, rprim, n1r,n2r,n3r, hh,
     $   wk(i08i09), wk(i01) )
c -->    rho-io      iwrk-2s
c
        if( lvlout.gt.1 )then
c         Find out how many electrons moved around:
          call GRIDCHG( delchg,dum,weight,nptr, wk(i08i09) )
          if( nspin.eq.1 )then
            write(IWR,*) ' '
            write(IWR,9020) 'delta movement in electrons =', delchg
          elseif( ispin.eq.1 )then
            write(IWR,*) ' '
            write(IWR,9020) 'delta movement in up-elecs =', delchg
          else
            write(IWR,9020) 'delta movement in dn-elecs =', delchg
          endif
        endif
c
c         ... and shift to dn-spin electrons, if we have them
        i05i06 = i06
        i08i09 = i09
      enddo
c
      if( lvlout.gt.1 .and. nspin.gt.1 )then
c       Get count of total electron movement(only up/dn separate above)
        call DCOPY( nptr, wk(i08),1, wk(i07),1 )
        call DAXPY( nptr, one, wk(i09),1, wk(i07),1 )
        call GRIDCHG( delchg,dum,weight,nptr, wk(i07) )
        write(IWR,9020) 'delta movement in electrons =', delchg
      endif
c
      if( lvlout.gt.2 ) call TIMER('After gridrho    ')
c
      if( do_gga )then
c       Construct gga gradients from reference+scf derivatives
c
        call GRADRHO( do_spin, idrhfl0,idrhofl,idrhsfl,igrhofl,
     $   n1r,n2r,n3r, nptr, gprim,
     $   wft(ift1),wft(ift2),wft(ift3),wft(ifta),
     $   wk(i08),    wk(i09),    wk(i01),wk(i11) )
c -->    rho[up]-io  rho[dn]-io  wkr-7s  wkrd-1s
c
        if( lvlout.gt.2 ) call TIMER('After gradrho    ')
        call MEMCHKIWF( idwf, 1, 'gradrho   ', wk(iwf), wfiwf )
        call MEMCHKIBR( idbr, 1, 'gradrho   ', wk(imblnde), savbr )
      endif
c
c Compute grid charge defects for atoms using nearby construction:
c
      if( iproc .eq. master )then
        REWIND( unit=idmatfl )
        if( do_spin ) REWIND( unit=idmatsfl )
      endif
      idmfile = idmatfl
      do  ispin=1,nspin
        call MKZERO( natm*nspin, defatom(1,ispin) )
        call ATDEFCT( IWR, idmfile,
     $   ncplx,znuc,rhoatom,defatom(1,ispin),
     $   norb,nk, natm,ntyp, noad, itypa, norba,
     $   near2c, wk(is1at), wk(ioint),wk(ioofk), wksma,
     $   wk(i01),      wk(i01+mat) )
c -->    dmat(mat)-s   dmatloc-s(nkloc*mat)
c       Switch to dn-spin dmat
        idmfile = idmatsfl
      enddo
c
      if( lvlout.gt.2 ) call TIMER('After atdefct    ')
c
c MEM: 08=delrho/09=delrhdn/10=rho0
      call DCOPY( nptr, wk(i08),1, wk(i12),1 )
c MEM: 08=delrho/09=delrhdn/10=rho0/12=delrho
c
c Initialize total Hamiltonian on grid:
      call MKZERO( nmatkp, wk(i01) )
c
c1c#####################################################################
c Following variables used in energy even if no 1-ctr performed
      ves1ctr = zero
      ves1c0  = zero
      exc1c = zero
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   Non-standard treatment of "fast" local non-reference (delta)-rho
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( any1ctr )then
c
c       Store rho0 in 05 and delrho in 06 for safe-keeping
        call DCOPY( nptr, wk(i10),1, wk(i05),1 )
        call DCOPY( nptr, wk(i08),1, wk(i06),1 )
c       "v1ctrij/blkvofn" constructs 1-ctr pots that reproduce the
c         rapid near-nucleus spatial variations of the delta-density
c       "v1ctrij" evaluates materix elements of fast es and xc pot'ls
c       "rho1ctr" removes from grid countercharge matching fast rho
c       "xconmsh" removes from grid fast part of xc-potentials
c
c       Set up sym info for symmetrization of local density matrix
c       in put-and-take treatment of rapidly varying local density:
c
        nsyms = nsym + 1
        call ROTMAT( nsymi,nsyms, isymtyp,syminfo, dsym )
c
c       Initialize vxc1c hamiltonian:
        call MKZERO( nmatkp, wk(i04) )
c
c       Compute fast local pots, and all matrix elements over fast pots
c
        call V1CTRIJ( idmatfl,ilocfl, conv1c,
     $   norb,nk, natm,ntyp,nshld,nald, noad, nlat,
     $   wk(i01),wk(i04),wk(i02),  dcountr,rh1coef,
c -->    v-o     xp-o    dmatrd-2o
     $   exc1c, cv,ce,
     $   itypa, numshl,lshel,nala,ala,cala, alamin, norba,ioffa, occij,
     $   lmx1ctr,almnv1c, nrad,nrd, radmsh,radwt,
     $   ratm,rlat, veck, wk(icoskr),wk(isinkr),
     $   neq3d,ieqn,ieqi,ieqj,ipermut,
     $   vradsr,v1csum,velem, pradsr,p1csum,pelem,xradsr,x1csum,xelem,
     $   wkrad(1,1),wkrad(1,8), wkrad(1,15),wkrad(1,16),wkrad(1,17),
     $                          wkrad(1,18),wkrad(1,19),wkrad(1,20),
     $   wk(i07), wk(i07+nrd*10*natm), wk(i07+2*nrd*10*natm),
c -->    ves1c    pot1c                xen1c
     $   nsyms,dsym,naofsym,
     $   alxcmin,nalfxc,nalfxcd,alfxc,alsq,ipvt, nrxcfit,vxcrad,excrad,
     $   nang, angylm,angpts,angwts, nfitd,nfityp,radwf,
     $   rholoc,vxcloc,excloc, vlm,elm,rfac, vg,eg,
     $   orb,grorb,grgrorb, wf1tmp1,wf1tmp2,wf1tmp3,wf1tmp4,wf1tmp5,
     $   orbj,grorbj,grgrorj, gofac,ggofac,
     $   wksmr, wk(i07+3*nrd*10*natm) )
c -->           dmatwk-s(nlat*noad**2)
c
        if( (i07 + 3*nrd*10*natm + nlat*noad**2) .gt. i11 )then
          call STOPXERR( 'v1ctrij scratch stomps stored memory' )
        endif
c
c       MEM: 01,02=v1cij,vxc1cij/03,04=dmatk/05=rho0/06=delrho
        call DAXPY( nmatkp, -one, wk(i01),1, wk(i04),1 )
        call DSCAL( nmatkp, -one, wk(i04),1 )
c       MEM: 01,04=v1cij,es1cij/02,03=dmatk/05=rho0/06=delrho
c
c       Take the needed traces here:
        call DOENG1C( ves1ctr,ves1c0,
     $   nk,norb,natm,ntyp,nshld, itypa,numshl,lshel, occij,wtk,
     $   wk(i04), wk(i02) )
c -->    es1cij-i dmat-2i
c
        if( lvlout.gt.3 )then
          write(IWR,*) ' '
          write(IWR,9020) '1-center es/es0 energy =',ves1ctr,ves1c0
          write(IWR,9020) '1-center xc-energy correction = ', exc1c
        endif
c
        if( lvlout.gt.2 ) call TIMER('After v1ctrij    ')
c
c       MEM: 01=v1cij/02,03=dmatk/05=rho0/06=delrho
c
c       Copy the del-rho(r) just calculated for use in rho1ctr:
        call DCOPY( nptr, wk(i06),1, wk(i12),1 )
c
c       Subtract fast local 1-center density from the density on grid
c       Remainder density should be slow, and better treated via fft
c
        REWIND( unit=i1ctrfl )
        call RHO1CTR( IWR, i1ctrfl,i1ctrfl,i1ctrfl,
     $   norb,nk, natm,ntyp,nshld,nald, noad,
     $   itypa, numshl,lshel,nala,ala,cala,  norba,ioffa, occij,
     $   lmx1ctr,almnv1c, nsyms,dsym,naofsym,
     $   hh,veck, n1r,n2r,n3r,nptr, weight, mngr3,mxgr3,ngr,
     $   dcountr,rh1coef, wksmc,
     $   wk(i12), wk(i02),   wk(i07), wk(i04) )
c -->    rho-io   dmatofk-2i r-3s     cnrho-s
c
        call GRIDCHG( slochg,dum,weight,nptr, wk(i12) )
        if( lvlout.gt.1 )
     $   write(IWR,9020) 'slow movement in electrons = ', slochg
c
        if( lvlout.gt.2 ) call TIMER('After rho1ctr    ')
c
c       MEM: 01=v1cij/05=rho0/06=delrho/12=slorho
c
c       Blank initial "fast" xc-potentials:
        call MKZERO( nptr, wk(i03) )
        call MKZERO( nptr, wk(i04) )
c
        if( anyxc1c )then
c
c         Remove from grid pot fast local 1-ctr xc-potentials
c         Strictly development code!!!
c
          call XCONMSH( natm,ntyp,nshld, nlat1c, n1r,n2r,n3r,nptr,
     $     itypa, numshl,lshel, ratm,hh,rlat, mngr3,ngr,
     $     nfityp,almnv1c,alxcmin,lmx1ctr,
     $     alfxc,nalfxc,cv,ce,nalfxcd,
     $     wk(i03),  wk(i04),  wk(i11), wk(i07) )
c -->      vxcmsh-io excmsh-io expfac-s r-3s
c
c         Save 1ctr grid fields to 1cgrfl ...
          REWIND( unit=i1cgrfl )
c          ... first on icgrfl: -vxc1c
          call WRITBIG( i1cgrfl, nptr, wk(i03) )
c          ... second on icgrfl: -exc1c
          call WRITBIG( i1cgrfl, nptr, wk(i04) )
c
          if( lvlout.gt.2 ) call TIMER('After xconmsh    ')
c         End "fast" xc 1-center
        endif
c
c       Move rho0 and delta-density back to assigned spaces:
        call DCOPY( nptr, wk(i05),1, wk(i10),1 )
        call DCOPY( nptr, wk(i06),1, wk(i08),1 )
c
c       End "fast" 1-center (non-standard)
      endif
c1c#####################################################################
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   Manipulation of electrostatics for supercell charge and moments
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Add dn to up-rho to get full delta-rho/sp-scaling
      if( do_spin ) call DAXPY( nptr, one, wk(i09),1, wk(i12),1 )
c
      call INTOBOX( n1r,n2r,n3r, wk(i12), ndbox, wk(i04) )
c
c MEM: 01=v1cij/08=delrho[up]/09=delrhdn/10=rho0/12=slorho
c MEM: 04=slorho:B
c
c Extract moments/potentials for supercell es with LMCC:
c
      if( Ldevtest )then
        write(IWR,1341) 'DBG: pre-lmcc, rchrg=',rchrg
        write(IWR,1341) 'DBG: pre-lmcc, rchrg0=',rchrg0
      endif
      call LMCC( IWR, igrdref,  lvlout, lmcc_on, ion_opt,
     $ esion,vacplus, eslmref, elchrg, ndim,nsym,
     $ dolvatm, rprim,rscale, orig,rsym,origws,rchrg, rmatsym,
     $ n1r,n2r,n3r,nptr, hh,weight, ndbox,nbox,
     $ wk(i04), wk(i05), wk(i06), wk(i07), wk(i11) )
c -->  wk1-i    wk2-o    wk3-o    wk4-3s   wk5-s
c      slorho   rho_lmcc v_lmcc
c
      if( Ldevtest )then
        write(IWR,1341) 'DBG: postlmcc, rchrg=',rchrg
        write(IWR,1341) 'DBG: postlmcc, rchrg0=',rchrg0
C        write(IWR,9020) 'DBG: postlmcc, esion =',esion
        write(IWR,*) 'DBG: lmcc_on=', lmcc_on
      endif
c
      if( lmcc_on )then
c
c       Transfer LMCC density, and potential, out of boxes:
        call FROMBOX( n1r,n2r,n3r, wk(i04), ndbox, wk(i05) )
        call FROMBOX( n1r,n2r,n3r, wk(i11), ndbox, wk(i06) )
c
c       Take LMCC density out of the slow rho for fft:
        call DAXPY( nptr, -one, wk(i04),1, wk(i12),1 )
c
        if( Ldevtest )then
          qslogrd = 0.d0
          do  i=1,nptr
            qslogrd = qslogrd + wk(i12-1+i)
          enddo
          qslogrd = weight*qslogrd
          write(IWR,9020) 'DBG: post-lmcc SLAB CHARGE=',qslogrd
          if( lvlout .gt. 1 )then
            call MOMENT( IWR, elchrg, str1  ,str2   , ndim,
c -->                                 dipole,rdipole
     $       origws, rprim, n1r,n2r,n3r, hh, weight, ndbox,nbox,
     $       wk(i12),    wk(i07) )
c
            call SYMVEC( 1, dipole, vec0, nsym, rmatsym )
            call SYMVEC( 1,rdipole, rsym, nsym, rmatsym )
c
            write(IWR,9020) 'DBG: post-lmcc SLAB DIPOLE=',dipole
            write(IWR,9020) 'DBG: post-lmcc SLAB RDIPOLE=',rdipole
          endif
c         End lmcc diagnostics
        endif
c
        write(IWR,*) 'Using LMCC treatment of supercell electrostatics'
        if( lvlout.gt.2 ) call TIMER('After lmcc       ')
c
      endif
c MEM: 01=v1cij/04=rholmcc/08=drho[up]/09=drhodn/10=rho0/
c MEM: 11=veslmcc/12=slorho
c
c Compute magnitudes of reciprocal lattice vectors
c
      call GVECMAG( n1r,n2r,n3r,nptr,gprim, wk(i07) )
c -->                                       gvecmag-o
c
c Compute slow-Ves on grid from slowly-varying grid density
c
      call VESSLO( IWR, ibndpot, n1r,n2r,n3r,nptr, weight,
     $ wft(ift1),wft(ift2),wft(ift3),wft(ifta),
     $ wk(i12),  wk(i02),  wk(i07) )
c -->  rhoslo-io espot-2os gvecmag-i
c
c     Move Ves into its appropriate space:
      call DCOPY( nptr, wk(i02),1, wk(i03),1 )
c     Need to get total delrho (=delrhup+delrhdn) somewhere ...
c      ... copy delrh(up) into a space
      call DCOPY( nptr, wk(i08),1, wk(i07),1 )
c      ... and add delrh(dn) into it if spin
      if( do_spin ) call DAXPY( nptr, one, wk(i09),1, wk(i07),1 )
c
c Compute contributions to electrostatic energy from grid terms
c
      if( lmcc_on )then
c
c       Compute energy terms for LMCC electrostatics:
c
        call LMCCINT( IWR, igrdref,  lvlout, lmcc_on,
     $   esnsgr, esion,eslmref,  ndim, nptr, weight,
     $   wk(i12), wk(i03),  wk(i04),  wk(i11),  wk(i07) )
c -->    rhofft-i vesfft-i  rholmcc-i veslmcc-i delrho-io
c <--    rhoslo-o vesslo-o
c
        if( lvlout.gt.1 )
     $   write(IWR,9020) 'LMCC esns total grid integral=',esnsgr
c
      else
c
c       Compute standard electrostatic grid energy terms:
c
        esnslo = half*weight*DDOT( nptr, wk(i12),1, wk(i03),1 )
        if( lvlout.gt.3 .and. any1ctr )
     $   write(IWR,9020) 'esns grid int of slo-rho:slo-ves =',esnslo
c       For esns, need int[del-rho:slo-ves], NOT int[slo-rho:slo-ves]
        esnsgr = half*weight*DDOT( nptr, wk(i07),1, wk(i03),1 )
c
        if( lvlout.gt.3 )
     $   write(IWR,9020) 'esns grid int of del-rho:slo-ves =',esnsgr
c
      endif
c
c  For Vxc trace calculation, need to take full Ves trace ...
c      ... rho0:Ves
      ves0xc = weight*DDOT( nptr, wk(i10),1, wk(i03),1 )
c      ... combined with delrho:Ves
      vesdxc = weight*DDOT( nptr, wk(i07),1, wk(i03),1 )
      vesgrxc = ves0xc + vesdxc
      if( lvlout.gt.3 )
     $write(IWR,9020) 'esns grid integral for Vxc trace =     ',vesgrxc
c
      if( lvlout.gt.2 ) call TIMER('After vesslo     ')
c
c MEM: 01=v1cij/03=sloves/07=delrho/08=delrhup/09=delrhdn/10=rho0/
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c          Iteration 0 enters code at this point:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 1500 continue
c
      if( iterscf .eq. 0 )then

        write(IWR,*) ' '
        write(IWR,9100)  iterscf
c
c       Initialize arrays for 0-th iteration ...
c
        vesgrxc = zero
        if( iproc .eq. master ) REWIND( unit=igrd0fl )
c        Retrieve rho[atom] from igrd0fl ...
        call MPBCREADBIG( igrd0fl, nptr, wk(i10) )
c        ... initialize current total delta-density to zero:
        call MKZERO( nptr, wk(i07) )
c        ... initialize delta (up-spin) density to zero:
        call MKZERO( nptr, wk(i08) )
c        ... and dn-spin density:
        if( do_spin ) call MKZERO( nptr, wk(i09) )
c        Initialize (blank) delta-Hamiltonian
        call MKZERO( nmatkp, wk(i01) )
c        Blank slo-Ves(r):
        call MKZERO( nptr, wk(i03) )
c
        if( lvlout.gt.1 ) call TIMER('Before iterations')
c
      endif
c
c MEM: 01=vij1c/03=sloves/07=delrho/08=delrhup/09=delrhdn/10=rho0
c
c Obtain electric field potential, and take relevant integrals:
c
      engydfld = zero
      engy0fld = zero
      if( do_field )then
c
c       Generate potential due to external electric field ...
        call MKZERO( 3, wksm )
        call MKZERO( nptr, wk(i04) )
        call FIELDPOT( efield, wksm, n1r,n2r,n3r,hh, wk(i04) )
c
CCCCCC  This was test code ---- to fix a problem? --- I forget
CCCc       Set average field from E-field to zero
CCC        call GRIDAVE( -1 ,n1r,n2r,n3r,1,wk(i04),wksm,wft(ifta),
CCC     $      'E-field' )
CCC        write(IWR,9020) 'Average potential from E-field=',wksm(1)
CCC        wksm(1) = -wksm(1)
CCC        call DAXPY( nptr, one, wksm(1),0, wk(i04),1 )
c
c        ... put it into the grid electrostatic potential
        call DAXPY( nptr, one, wk(i04),1, wk(i03),1 )
c        ... take the energy integral with the delta-density:
        engydfld = weight*DDOT( nptr, wk(i07),1, wk(i04),1 )
c        ... take the energy integral with the reference rho0:
        engy0fld = weight*DDOT( nptr, wk(i10),1, wk(i04),1 )
c        ... add into the grid integral used to create Vxc engy:
        vesgrxc = vesgrxc + engydfld + engy0fld
c
        write(IWR,9020)'Field energy, INT:V[E]*drho (Ry)=',engydfld
        write(IWR,9020)'Field energy, INT:V[E]*drho (eV)=',engydfld*toev
        write(IWR,9020)'Field energy, INT:V[E]*rho0 (Ry)=',engy0fld
        write(IWR,9020)'new esns grid int for Tr(Vxc) =',vesgrxc
      endif
c
c Compute xc-potential(s) over the grid
c
      call GRIDDFT( IWR, igrd0fl,igrhofl,
     $ i1cgrfl, anyxc1c,
     $ anycore, do_gga,do_spin, iterscf,lvlout,
     $ nptr, weight,xcrhocut, exctot,excslo,corslow,eharris,
     $ wk(i10),wk(i02),wk(i08), wk(i04),wk(i06),wk(i07),
c -->  rho0-i  rhoc-s  drhoa-io dvxca-o dexca-o exc0-s
     $                 wk(i09), wk(i05),wk(i11),
c -->                  drhob-io dvxcb-o dexcb-o
     $ wk(i01), wk(i12) )
c -->  wksp1-s  wksp2-s
c
c     Adding slow-Ves into slow-Vxc, to get full slow potential on grid
      call DAXPY( nptr, one, wk(i03),1, wk(i04),1 )
      if( do_spin ) call DAXPY( nptr, one, wk(i03),1, wk(i05),1 )
c
c MEM: 01=vij/03=sloves/04=vsloa/06=dexca/07=exc0/08=drhoa
c MEM: (spin) 09=drhob/05=vslob/11=dexcb
c
      if( lvlout.gt.2 ) call TIMER('After griddft    ')
c
      if( iterscf .gt. 0 .and. iproc .eq. master )then
c        Write to igridfl quantities to be used later in the code
c         (for force/stress calculation)
        REWIND( unit=igridfl )
c        ... full-delta-rho (non-spin); delta-rho-up/dn(spin) to igridfl
        call WRITBIG( igridfl, nptr, wk(i08) )
        if( do_spin ) call WRITBIG( igridfl, nptr, wk(i09) )
c        ... slow-delta-electrostatic potential to igridfl
        call WRITBIG( igridfl, nptr, wk(i03) )
c        ... slow-delta-pot(xc+es) to igridfl, total(nonspin);up/dn(spin)
        call WRITBIG( igridfl, nptr, wk(i04) )
        if( do_spin ) call WRITBIG( igridfl, nptr, wk(i05) )
c        End of igridfl
        if( lvlout.gt.2 ) call TIMER('After igridfl    ')
      endif
c
c If doing spin, initialize both spin Hamiltonians
      if( do_spin ) call MKZERO( 2*nmatkp, wk(i01) )
c
      if( nearopt.gt.0 )then
c
c       Compute values of grid potentials "nearby" atoms:
        mxwkwt = nptr
        i01i02 = i01
        i04i05 = i04
        i06i11 = i06
        do  ispin=1,nspin
          call POTNEAR( nearopt, natm,
     $     ratm,hh,rprim, mxwkwt ,n1r,n2r,n3r,nptr, nearatm,
     $     wksmo(1,ispin), wksmo(1,3), wksmo(1,3+ispin),
c -->      vslnear-o       vesnear-o   excnear-o
     $     wk(i04i05), wk(i03), wk(i06i11), wk(i07) )
c -->      vslo-io     vesslo-i excslo-i    wtnr-s
c
c         Compute correction to matrix elements from nearby grid defect:
          call VSLOFIX( near2c, wk(is1at), wk(ioint),wk(ioofk),
     $     norb,nk, natm,ntyp,nshld, noad,  itypa, numshl,lshel,
     $     wksmo(1,ispin), ncplx,
c -->      vslnear-i
     $     wk(i01i02) )
c -->      v-io
c
c       ... switch to dn-spin electrons, if we have them
          i01i02 = i01 + nmatkp
          i04i05 = i05
          i06i11 = i11
        enddo
c
        if( lvlout.gt.2 ) call TIMER('After vslofix    ')
      else
c       Blank nearby arrays, so that doengy routine gets zeroes:
        call MKZERO( natm, wksmo(1,1) )
        call MKZERO( natm, wksmo(1,2) )
        call MKZERO( natm, wksmo(1,3) )
        call MKZERO( natm, wksmo(1,4) )
        call MKZERO( natm, wksmo(1,5) )
      endif
c
c Compute matrix elements of (slow-varying) grid potential:
c
      i01i02 = i01
      i04i05 = i04
c     Position Ham scratch directly behind stored Hams:
      i03s = i01 + nmatskp
      checkb4 = wk(i05)
      do  ispin=1,nspin
        call INTOBOX( n1r,n2r,n3r, wk(i04i05), ndbox, wk(i06) )
c
        call VSLOMAT( ncplx,idwf,doblas3,
     $   norb,nk, natm,ntyp, nshld,nald, nlat1c,
     $   itypa,numshl,lshel,nala,ala,cala, alamin, ratm,rlat,
     $   wk(icoskr),wk(isinkr),
     $   weight, hh, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1),boxrad,
     $   wk(iwf),wk(iwf), wkrbx(1),     wkrbx(1+mxrb), wk(iwfw),
c -->    wvfcns  wf       vslobox(1B)-s wkbox(10Br8)-s wkbi(lots)-s
     $   wk(i01i02), wk(i03s), wk(i06) )
c -->    vijk-io     vkij-s    vslo-i
c
c        ... switch to dn-spin electrons, if we have them
        i01i02 = i01 + nmatkp
        i04i05 = i05
      enddo
      checkaft = wk(i05)
c     At this point, grid potential in i04 is assumed to be lost
c     Check that our scratch space in i03s did not run over (into i05):
      if( checkb4 .ne. checkaft )then
c       Oops, we have a problem, i03s ran over when it should not have
c       Put this check in because this is a repeat offender
        write(IWR,*) '***** ERROR: code memory error, scratch overrun'
        write(IWR,*) '>>>>> this is a coding error, notify author'
        write(IWR,*) '>>>>> (mvslomat is incorrect in wkmem)'
        call STOPXERR( 'mvslomat scratch memory coding error' )
      endif
      call MEMCHKIWF( idwf, 1, 'vslomat   ', wk(iwf), wfiwf )
      call MEMCHKIBR( idbr, 1, 'vslomat   ', wk(imblnde), savbr )
c
      if( lvlout.gt.2 ) call TIMER('After vslomat    ')
      write(IWR,*) ' '
c
      engydone = .false.
      if( iterscf .gt. 0 )then
        if( lvlout .ge. lvlengy )then
c
c         Compute total energy of system and output energy analysis
c
          engydone = .true.
          call DOENGY( engytotl, IWR,
     $     idmatfl,idmatsfl, iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $     nspin, nearopt,near2c,
     $     atm0engy,esnsgr,vesgrxc,excengy,esnsengy, ves1ctr,ves1c0,
     $     do_field,engydfld,engy0fld, engypol,
     $     xc0val,exctot,excslo,exc1c, corslow,
     $     norb, nk, natm,ntyp,nshld,
     $     itypa,numshl,lshel,znuc,occij,wtk,
     $     rhoatom, defatom, defatom(1,2),
     $     wksmo(1,1),wksmo(1,3),wksmo(1,4), norbd,
c -->      vslnear-i  vesnear-i  excnear-i
     $     wk(i01), wk(i01+nmatskp), wk(i01+nmatskp+mat) )
c -->      v-io     dmat-ks             h0mat-ks
c
          if( lvlout.gt.2 ) call TIMER('After energies   ')
          call MEMCHKIWF( idwf, 1, 'doengy-1  ', wk(iwf), wfiwf )
          call MEMCHKIBR( idbr, 1, 'doengy-1  ', wk(imblnde), savbr )
        else
          write(IWR,9020)'TOTAL single particle eigenvalue sum= ',esumsp
        endif
      endif
c
      write(IWR,*) ' '
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c      Hamiltonian complete - do blend/check for self consistency
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  At this point the output delta-Hamiltonian matrix is in "01".
c  If iterscf=0, use as input for iteration 1, otherwise, blend with
c  earlier potential matrices to obtain input for next iteration.
c
      call MPNODES( nprocs )
      call MPNODE0( master )
      call MPNODE( iproc )
      call MPCOMM( icomm )
c  Set up space in arrays used for blending.  I have my current
c  dHams in i01 (and i02 if spin), and need to allocated the space
c  for the blend on the combined spin Hams.
      nmatsbl = nmats
      nprocbl = nprocs
      call KPFLAG_GET( kparopt )
      if( nprocs.gt.1 .and. kparopt.eq.2 )then
        nmatsbl = nmatskp
        call MPNODES_K( nprocbl )
      endif
      nmatsloc = ( nmatsbl + nprocbl - 1 ) / nprocbl
      i01bl = i01
      i02bl = i01bl + nmatsbl
      i03bl = i02bl + nmatsloc
      i04bl = i03bl + nmatsloc
c
      iscfstat = 0
      if( iterscf .eq. 0 )then
c
c       FIRST ITER, USE INITIAL HAMILTONIAN WITH GUESS
c
        call GESSER( iges0fl,igesfl, igesopt,
     $   do_geom, madeges, mat, nk, nspin,
     $   wk(i01bl),  wk(i02bl) )
c -->    vij(sp2)-io wk1(sp2)-s
c
        call MEMCHKIWF( idwf, 1, 'gesser    ', wk(iwf), wfiwf )
        call MEMCHKIBR( idbr, 1, 'gesser    ', wk(imblnde), savbr )
c
        if( lvlout.gt.1 ) call TIMER('After 0 iteration')
c
      else
c
c       * * * * * * * * * * * * *
c       CHECK FOR SCF CONVERGENCE
c       * * * * * * * * * * * * *
c
c       Have output Ham, get input Ham and take difference for convergene check
        call MPNODE0_K( master_k )
        call MPNODE_k( iproc_k )
        if( iproc.eq.master .or. iproc_k.eq.master_k )then
c         This assumes data is k-parallel AND that master is a k-master too.
          if( iproc.eq.master ) REWIND( unit=idhamfl )
          call H0READKP( idhamfl, mat,nk, wk(i02bl),
     $                                    wk(i02bl+nmatskp) )
          if( do_spin )
     $      call H0READKP( idhamfl, mat,nk, wk(i02bl+nmatkp),
     $                                      wk(i02bl+nmatskp) )
c         Construct v(out)-v(in) difference
          call DAXPY( nmatskp, -one, wk(i02bl),1, wk(i01bl),1 )
          nmatscf = nmatskp
        else
          nmatscf = 0
        endif
c
        call RUNSTEER( scfconv, gconv,  cellconv,
     $                 itstop,  igstop, maxucstep )
c
        call SCFCHECK( IWR, 'scfchange', iscfstat,
     $   iterscf,itstop,scfconv, chgmax,chgrms,  nmatscf,
     $   wk(i01bl) )
c
        call MEMCHKIWF( idwf, 1, 'scfcheck  ', wk(iwf), wfiwf )
        call MEMCHKIBR( idbr, 1, 'scfcheck  ', wk(imblnde), savbr )
c
        if( lvlout.gt.2 ) call TIMER('After scfcheck   ')
c
c       * * * * * * * * * * * * *
c       EVALUATE ENERGY FINAL ITER
c       * * * * * * * * * * * * *
c
        if( iscfstat.gt.0 .and. .not. engydone )then
c         Compute total energy of system and output energy analysis
c
c         Reconstruct output Ham after scfcheck
          if( iproc.eq.master .or. iproc_k.eq.master_k )then
c           k-master k-parallel
            call DAXPY( nmatskp, one, wk(i02bl),1, wk(i01bl),1 )
          endif
c
c         ***** Evaluate energy *****
c
          call DOENGY( engytotl, IWR,
     $     idmatfl,idmatsfl, iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $     nspin, nearopt,near2c,
     $     atm0engy,esnsgr,vesgrxc,excengy,esnsengy, ves1ctr,ves1c0,
     $     do_field,engydfld,engy0fld, engypol,
     $     xc0val,exctot,excslo,exc1c, corslow,
     $     norb, nk, natm,ntyp,nshld,
     $     itypa,numshl,lshel,znuc,occij,wtk,
     $     rhoatom, defatom, defatom(1,2),
     $     wksmo(1,1),wksmo(1,3),wksmo(1,4), norbd,
c -->      vslnear-i  vesnear-i  excnear-i
     $     wk(i01), wk(i01+nmatskp), wk(i01+nmatskp+mat) )
c -->      v-io     dmat-ks             h0mat-ks
c
c         Reconstruct the v(out)-v(in) difference again ...
          if( iproc.eq.master .or. iproc_k.eq.master_k )then
            if( iproc.eq.master) REWIND( unit=idhamfl )
            call H0READKP( idhamfl, mat,nk, wk(i02bl),
     $                                      wk(i02bl+nmatskp) )
            if( do_spin )
     $        call H0READKP( idhamfl, mat,nk, wk(i02bl+nmatkp),
     $                                        wk(i02bl+nmatskp) )
            call DAXPY( nmatskp, -one, wk(i02bl),1, wk(i01bl),1 )
          endif
c
          if( lvlout.gt.2 ) call TIMER('After energies   ')
        endif
c
c       * * * * * * * * * * * * *
c       DO BLEND [Modified broyden: D.D. Johnson, PRB 38, 12807 (1988)]
c       * * * * * * * * * * * * *
c
        call MEMCHKIWF( idwf, 1, 'b4 ebroy  ', wk(iwf), wfiwf )
        call MEMCHKIBR( idbr, 1, 'b4 ebroy  ', wk(imblnde), savbr )
c
        if( iscfstat .lt. 2 )then
c         Not converged: Blend to get input to next scf iteration
c
c         Distribute vectors over processors:
          if( nprocs .gt. 1 )then
            call KPFLAG_GET( kparopt )
            if( kparopt.eq.2 )then
c             Distribute k-parallel matrices from k-masters
              call MPNODES_K( nproc_k )
              call MPNODE0_K( master_k )
              call MPNODE_K( iproc_k )
              call MPCOMM_K( icomm_k )
              nmatsbl = nmatskp
              nprocbl = nproc_k
              masterbl = master_k
              iprocbl = iproc_k
              icommbl = icomm_k
            else
c             Distribute full matrices from master
              nmatsbl = nmats
              nprocbl = nprocs
              masterbl = master
              iprocbl = iproc
              icommbl = icomm
            endif
            call MPSPLITR8( nprocbl, masterbl, iprocbl,
     $           nmatsbl,nmatsloc, wk(i01bl), icommbl )
            call MPSPLITR8( nprocbl, masterbl, iprocbl,
     $           nmatsbl,nmatsloc, wk(i02bl), icommbl )
          else
c           Single processor, local allotment is all ...
            nmatsloc = nmats
          endif
cxxx  Test code for installing new ASD blend
          if( meth_bl.eq.1 )then
            write(IWR,*)  'electronic meth_bl(BROYDEN)=',meth_bl
            call EBROYDDJ( IWR, 'elec', iterscf,itbroye,nhiste,scfblnd,
     $       chgmax,nmats,nmatsloc, memblnde,
     $       wk(i01bl), wk(i02bl), wk(i03bl), wk(i04bl), wk(imblnde) )
c -->        vv1-i(dv)  vv2-io(v)  vv3-s      vv4-s      wbroy-io
c
          elseif( meth_bl.eq.2 )then
            write(IWR,*)  'electronic meth_bl(ASD)=',meth_bl
c           Ugly, hard code asd blend meth=3 into this routine:
            call EBLENDASD( IWR, 'elec', iterscf, memblnde,
     $       nmatsloc,  scfblnd,
c     $       ratmg,vatmg,frcg, watmg )
     $       wk(i01bl),wk(i02bl),wk(i03bl),wk(i04bl), wk(imblnde) )
c -->        ri-i      fi-i/ri-o vi-s      wt-s
          else
            write(IWR,*)  'electron meth_bl(*****ERROR*****)=',meth_bl
            call STOPXERR( 'Invalid value for Meth_bl' )
          endif
c
          if( nprocbl .gt. 1 )then
c           Merge resulting guess vector onto appropriate master(s) again ...
            call MPMERGER8( nprocbl, masterbl, iprocbl,
     $           nmatsbl, 1   ,nmatsloc, wk(i02bl), icommbl )
          endif
          call MEMCHKIWF( idwf, 1, 'ebroyddj  ', wk(iwf), wfiwf )
          call MEMCHKIBR( idbr, 0, 'ebroyddj  ', wk(imblnde), savbr )
c
        endif
c
c       Move updated Hamiltonian into right place:
        call DCOPY( nmatskp, wk(i02bl),1, wk(i01bl),1 )
c
        if( lvlout.gt.1 ) call TIMER('After blending   ')
c
        write(IWR,*) ' '
c
c       Converged SCF, go on ...
        if( iscfstat .eq. 2 ) goto 1800
c
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    Not SCF converged
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call MEMCHKIWF( idwf, 1, 'notSCF    ', wk(iwf), wfiwf )
      call MEMCHKIBR( idbr, 1, 'notSCF    ', wk(imblnde), savbr )
c Write next input Hamiltonian to file:
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      if( iproc.eq.master .or. iproc_k.eq.master_k )then
        if( iproc.eq.master)  REWIND( unit=idhamfl )
        call H0WRITKP( idhamfl, mat,nk, wk(i01bl), wk(i02bl) )
        if( do_spin )
     $  call H0WRITKP( idhamfl, mat,nk, wk(i01bl+nmatkp), wk(i02bl) )
      endif
c
c Record status that this iteration complete:
      call WRSTAT( IWR, istatfl, 'ITERATION',
     $ do_geom, do_cell, opt_spin,
     $ iucstep,igstep,iterscf+1, engytotl,excengy,esnsengy,
     $ istep_sp, elecno, spindata,
     $ rprim, natm, ratm )
c
      if( iscfstat .eq. 0 )then
c       Not converged, go on to next iteration:
c
        iterscf = iterscf + 1
        itdone = itdone + 1
        call MTBUFF
        goto 1200
c
      elseif( iscfstat .eq. 1 )then
c       Not converged, and we ran out of iterations
c
        write(IWR,9020)'FINAL *NOT* CONVERGED ENERGY (Ry)=',engytotl
        if( engypol .ne. zero )
     $  write(IWR,9020)'FINAL *NOT* CONVERGED ENERGY w/polarization=',
     $   engytotl - ABS( engypol )
c
        call DOEBIND( IWR, engytotl, engyref, engybind,
     $   toev, natm, itypa, atengy )
c
        write(IWR,9180)  iterscf, iterscf+1
 9180   format(1x,'>>>>> itstop=',i4,', self-consistency NOT achieved'
     $        /1x,'      to continue on, set first iteration=',i4 )
c
        if( do_neb )then
c         Must tell neb_master, to be able to close less badly
          write(IWR,*) 'Warning NEB master that this scf failed'
          call NEBSTATX( IWRNEB, 'scf', image )
        endif
c
        call TIMER('AFTER ITERATIONS ')
        call STOPXOK( 'NO CONVG' )
        goto 999
c
      else
        call STOPXERR( 'code err - iscfstat must be 0 or 1 here' )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c      SCF Converged, do some final output and go on to next phase
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 1800 continue
c
c     Converged, commit the good news to the output file
      itconv = iterscf
      write(IWR,9181)  iterscf
 9181 format(1x,'>>>>> self-consistency achieved, iteration',i4/)
c
c  Spin polarization convergence check, and optimization cycle
c
      if( opt_spin .or. ( do_spin .and. lvlout.gt.2 ) )then
c
c       Write out the fermi level here, at least
        lprnt = 1
c       With higher output level, put out occupied eigenstates
        if( lvlout.gt.3 ) lprnt = 2
        call EVPRNT( IWR, lprnt, nspin, efermi,egap, npop,nocc0k,
     $   nstate,nk, wtk,wk(ieigval),wk(ieigpop) )
        write(IWR,*) ' '
c
c       Compute joint fermi level, elec(sp) with this eigenspectrum
        call SPINOPT( nspin, efspin, spinold,spinnew,
     $   etemp, elecno,  nstate,nk, wtk, wk(ieigval),
     $   wk(i02bl) )
c -->    wklvl(6*n*k*sp)-s
c
c       If doing an optimization, decide if converged, next step
        if( opt_spin )then
          call QHIST( ihistfl, 'ESPI', ndata, engytotl, elecno )
c       
c         Is spin converged, out of steps, done?
c     
          call SPUPDATE( istat_sp, istep_sp, spindata,
     $     elecno, spinold,spinnew,spinges,
     $     nstep_sp,meth_sp,conv_sp,blend_sp )
c
          write(IWR,9182) istep_sp,spinold,spinnew,spinges
 9182     format(1x,'Spin step',i3,'; spinpol: inp,out,ges=',3f12.6)
c
          if( istat_sp .eq. 0 .or. istat_sp .eq. 1 )then
c           If spin-pol not converged, commit current step to status:
            call WRSTAT( IWR, istatfl, 'SPIT',
     $       do_geom, do_cell, opt_spin,
     $       icellstep,igstep+1,iterscf, engytotl,excengy,esnsengy,
     $       istep_sp, elecno, spindata,
     $       rprim, natm, ratm )
          endif
c
          if( istat_sp .eq. 0 )then
c           Not converged spin polarization, go try again ...
            write(IWR,9183)  spinold, engytotl
 9183       format(
     $       /1x,'SPINPOL=',f16.6,', CURRENT SPINSCF ENERGY=',f20.10)
            call TIMER('After spin-scf   ')
            write(IWR,*)
c           Start with the last Ham (prevent GESSER w/iterscf.ne.0):
            itstart = -1
            iterscf = 1
            itoff = 0
            goto 1100
c
          elseif( istat_sp .eq. 1 )then
c           Not converged spin, and out of steps
            write(IWR,9183)  spinold, engytotl
            write(IWR,9186)  istep_sp, spinnew
 9186       format(1x,'Out-of-steps,  step=',i4,', spin NOT converged',
     $       /1x,'; Next Spin Polarization=',f16.6)
            call TIMER('AFTER NOSPINCONV ')
            call STOPXOK( 'NO SPIN CONV' )
            goto 999
c    
          elseif( istat_sp .eq. 2 )then
c           Converged spin polarization
            write(IWR,9184)  spinold,engytotl
 9184       format(/1x,'SPINOPT Converged',
     $       /1x,'SPINPOL=',f16.6,', FINAL SPINSCF ENERGY=',f20.10)
c
            call QHIST( ihistfl, 'SP-F', ndata, spinold, elecno )
            call WRSTAT( IWR, istatfl, 'SP-F',
     $       do_geom, do_cell, opt_spin,
     $       icellstep,igstep+1,iterscf, engytotl,excengy,esnsengy,
     $       istep_sp, elecno, spindata,
     $       rprim, natm, ratm )
          endif
c
        endif
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c               SCF+spinopt is converged at this point
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Write out eigenspectrum
      if( doeigen .and. iproc .eq. master )then
c       Summarize final eigenspectrum
        lprnt = 1
        if( lvlout.gt.1 ) lprnt = lprnt + 1
        if( lvlout.gt.4 ) lprnt = lprnt + 1
        if( ndim.eq.2 .and. ion_opt.eq.2 .and. vacplus .ne. zero )then
c         Slab calculation, neutral, with dipole ...
          write(IWR,9020) 'Final (+z) VACUUM level(Ry)=', vacplus
          write(IWR,9020) 'Final (-z) VACUUM level(Ry)=',-vacplus
        endif
c
        if( nspin.eq.1 )then
          write(IWR,9191) 'Final GAP(eV)=',egap(1)*toev
        else
          write(IWR,9192) 'Final GAP(eV)=',egap(1)*toev,egap(2)*toev
        endif
 9191   format(/1x,a,f8.4)
 9192   format(/1x,a,f8.4,' (up),', f8.4,' (dn)')
c
        call EVPRNT( IWR, lprnt, nspin, efermi,egap, npop,nocc0k,
     $   nstate,nk, wtk,wk(ieigval),wk(ieigpop) )
        write(IWR,*) ' '
c
        if( do_post )then
c         Output data needed for density of states calculation
          call FLOPENA( ipostfl, 'post' )
          REWIND( ipostfl )
          call DOSDUMP( ipostfl,
     $     nspin, norb,nstate,nk,ncplx, natm,ntyp,nshld,nald,
     $     itypa,numshl,lshel,
     $     ala,cala,nala,
     $     etemp, efermi, wk(ieigval),wk(ieigpop), wtk,veck, typnm )
          call FLCLOSE( ipostfl )
        endif
c
      endif
c
      call GESSET( iges0fl,igesfl, igesopt,
     $ do_geom, madeges, mat,nk,nspin,
     $ wk(i01bl),  wk(i02bl) )
c -->  vij(sp2)-io wk1(sp2)-s
c
c Now finished with guess file, close it
      call GESCLOSE( iges0fl, igesfl )
c
c Save scf-Ham result for scf-guess in next geometry
      call MPNODE( iproc )
      call MPNODE0( master )
      call MPCOMM( icomm )
c
c If we asked for do_post, transfer grid density to output file
      if( do_post .and. iproc .eq. master )then
c       Save converged grid delta density to data file:
        REWIND( unit=igridfl )
        call READBIG( igridfl, nptr, wk(i08) )
        if( do_spin )then
          call READBIG( igridfl, nptr, wk(i09) )
          call DAXPY( nptr, one, wk(i09),1, wk(i08),1 )
        endif
        call RHOSAVE( nptr,wk(i08) )
      endif
c
c Communicate electrostatic potential out ...
      if( do_post .and. iproc .eq. master )then
c       This read (del-Ves) assumes the delta-rho's read immediately above ...
        call READBIG( igridfl, nptr, wk(i09) )
        if( ndim .eq. 2 )then
c         write the density (do not read the density, do output to IWR) ...
          ipotfil = -1
          ioutfil = IWR
          call ZPOTOUT( ipotfil,ioutfil,
     $     iatmfmt, ntyp,natm, natmnm, itypa, atmnm,typnm, ratm,
     $     nspin, hh, n1r,n2r,n3r,nptr,
     $     wft(ifta), wk(i08), wk(i09) )
c -->      wkz(n3r)   delrho-i delves-i
        endif
      endif
c
      write(IWR,*)
      if( do_neb )then
        write(IWR,9022)
     $   'IMAGE',image,' CONVERGED ENERGY (Ry)=',engytotl
 9022   format(1x,a,i3,a,3f20.10)
      else
        write(IWR,9020) 'FINAL CONVERGED ENERGY (Rydberg) =',engytotl
cqmd: where to get energy for MD module?
        if ( do_md ) call QMDEPSET( engytotl )
        if( engypol .ne. zero )
     $  write(IWR,9020) 'FINAL CONVERGED ENERGY w/polarization=',
     $   engytotl - ABS( engypol )
      endif
c
      call DOEBIND( IWR, engytotl, engyref, engybind,
     $ toev, natm, itypa, atengy )
c
c Record status as converged, with energies:
      call QHIST( ihistfl, 'ESCF', ndata, engytotl, wksm )
      call WRSTAT( IWR, istatfl, 'SCF-ENERGIES',
     $ do_geom, do_cell, opt_spin,
     $ iucstep,igstep,iterscf, engytotl,excengy,esnsengy,
     $ istart_sp, elecno, spindata,
     $ rprim, natm, ratm )
c
      call TIMER('AFTER ITERATIONS ')
c
c *********************************************************************
c **********             end of iteration phase             ***********
c *********************************************************************
c **********          start force and stress phase          ***********
c *********************************************************************
c
 2000 continue
c Use a bcast to do a sync/SYNC here, to let k-masters catch up
      call KPSYNC()
c
      if( do_nonscf )then
c       Branch to non-scf (band structure) phase of code
c       For now, do_nonscf implies no force/relax/cell
        goto 6000
      endif
c
      if( .not. doforce )then
        if( .not. do_geom )then
          call QHIST( ihistfl, 'END ' , ndata, rdata, wksm )
          call STOPXOK( 'no force calc requested' )
          goto 999
        endif
        goto 3000
      endif
c
c  Check for output of piece-wise force contributions:
      IWRF = 0
      if( lvlout.gt.3 ) IWRF = IWR
c
c  Initialize force/stress to zero:
      call MKZERO( 3*natm, frctot )
      call MKZERO( 9, strtot )
c
c  Recover analytic atomic xc stress integrals.
      strxc0 = xc0val(10)
      strxcc = xc0val(11)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Do contributions to force/stress from grid-based pieces
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  Partial evaluation of e.s. contribution to stress is accomplished in
c  routine "dvstres." This requires slow charge density on "igridfl."
c
      if( iproc .eq. master ) REWIND( unit=igridfl )
c Read delta charge density from "igridfl"
      call MPBCREADBIG( igridfl, nptr, wk(i08) )
      call DCOPY( nptr, wk(i08),1, wk(i07),1 )
      if( do_spin )then
c        ... and combine with delrho(dn) from "igridfl"
        call MPBCREADBIG( igridfl, nptr, wk(i09) )
        call DAXPY( nptr, one, wk(i09),1, wk(i07),1 )
      endif
c
c MEM: 07=delrho 08=drho/up 09=drhodn
c
      if( ndim.gt.0 )then
c
        call DVSTRES( IWR, ndim, str1,
     $   n1r,n2r,n3r,nptr,weight, gprim,
     $   wft(ift1),wft(ift2),wft(ift3),wft(ifta),
     $   wk(i07),  wk(i03), wk(i05) )
c -->    rhoslo-io espot-2s rhofft-2s
c
c  Symmetrize this stress calculation:
        call SYMSTR( ndim, str1, nsym, rmatsym )
c  Add contributions to stress due to change in size of the unit cell:
        strduc = excengy + esnsengy
        do  id=1,ndim
          str2(id,id) = strduc
        enddo
        do  id=1,ndim
          strtot(id,id) = strtot(id,id) + str1(id,id) + str2(id,id)
        enddo
c
        if( lvlout.gt.2 )then
          write(IWR,9020)'dvstres stress terms=',(str1(id,id),id=1,ndim)
          write(IWR,9020) 'other vol stress=',(str2(id,id),id=1,ndim)
        endif
c
        if( lvlout.gt.2 ) call TIMER('After dvstres    ')
      endif
c
c MEM: 07=delrho 08=drho/up 09=drhodn
c
      if( iproc .eq. master ) REWIND( unit=igrd0fl )
c
c If periodic calculation and GGA, add GGA-specific term to stress:
c
      if( do_gga .and. ndim.gt.0 )then
c
c       Set up densities ...
c        ... get rho0 from igrd0fl ...
        call MPBCREADBIG( igrd0fl, nptr, wk(i03) )
        if( anycore )then
c        ... get rho_core from igrd0fl:
          call MPBCREADBIG( igrd0fl, nptr, wk(i04) )
c          ... and combine with ref density in 03:
          call DAXPY( nptr, one, wk(i04),1, wk(i03),1 )
        endif
c
c       Construct xc-density (full, or spin-polarized)
        if( do_spin )then
c         Add half reference+core density into up and dn delta-density
          call DAXPY( nptr, half, wk(i03),1, wk(i08),1 )
          call DAXPY( nptr, half, wk(i03),1, wk(i09),1 )
        else
c         Put reference+core density into delta-density
          call DAXPY( nptr, one , wk(i03),1, wk(i08),1 )
        endif
c
c       Compute the GGA-specific term in stress:
        call GGASTRESS( str1, do_gga, do_spin,
     $   idrhfl0,idrhofl,idrhsfl,igrhofl,
     $   ndim, nptr,weight,xcrhocut,
     $   wk(i08),  wk(i09),  wk(i02),   wk(i03),     wk(i04),
c -->    rhoa(r)-i rhob(r)-i g_dxc(r)-s g2_dxc2(r)-s gg_d2dn(r)-s
     $   wk(i05),  wk(i10) )
c -->    drhoa-3s  drhob-3s
c
c       Symmetrize stress (should not be needed, but completeness):
        call SYMSTR( ndim, str1, nsym, rmatsym )
        do  jd=1,ndim
          do  id=1,ndim
            strtot(id,jd) = strtot(id,jd) + str1(id,jd)
          enddo
        enddo
c
        if( lvlout.gt.2 )then
          write(IWR,*)
          write(IWR,*) 'GGA derivative correction to stress:'
          do  jd=1,ndim
            write(IWR,9020)  '    ',(str1(id,jd),id=1,ndim)
          enddo
        endif
c
        call TIMER('After ggastress  ')
c
      else
c
c       Set location within reference grid file:
c      ... skip rho0 from igrd0fl:
        call MPBCREADSKP( igrd0fl, nptr )
c        ... skip rho_core from igrd0fl:
        if( anycore ) call MPBCREADSKP( igrd0fl, nptr )
c
      endif
c
c MEM: everything destroyed.
c
c        ... retrieve Vxc0 from igrd0fl ...
      call MPBCREADBIG( igrd0fl, nptr, wk(i03) )
c
c  Read slow Ves part of potential, also saved on "igridfl"
      call MPBCREADBIG( igridfl, nptr, wk(i02) )
c  Read slow Vslo full potential on mesh, saved on "igridfl"
      call MPBCREADBIG( igridfl, nptr, wk(i04) )
      if( do_spin) call MPBCREADBIG( igridfl, nptr, wk(i05) )
c
c MEM: 02=Vesslo/03=Vxc0/04=Vslo-up/05=Vslo-dn
c
      if( anycore )then
c       Construct total Vxc potential ...
c        ... start with delta-potential (Ves+delVxc)
        call DCOPY( nptr, wk(i04),1, wk(i06),1 )
        if( do_spin )then
c          ... add in dn-spin delta-potential to up-spin potential
          call DAXPY( nptr, one, wk(i05),1, wk(i06),1 )
c          ... split in half since we have both up/dn core together
          call DSCAL( nptr, half, wk(i06),1 )
        endif
c        ... remove es potential to get del-Vxc ...
        call DAXPY( nptr,-one, wk(i02),1, wk(i06),1 )
c        ... and add Vxc0 to get full Vxc...
        call DAXPY( nptr, one, wk(i03),1, wk(i06),1 )
c        ... and put into boxes
        call INTOBOX( n1r,n2r,n3r, wk(i06), ndbox, wk(i01) )
c
c       Compute force on core charge from xc-potential
c
        call CORFORC( IWR, ndim, frc1, str1,
     $   natm,ntyp,nrd, nlat1c,  itypa,numshl, ratm,rlat,
     $   hh,weight, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1),boxrad,
     $   nrcor, radmsh,corden, wspl(1,1),wspl(1,2),
     $   wk(i01),  wk(i06), wkrbx )
c -->    vxcgrd-io r-3s     wkb(2B)-s
c
        call SYMFRC( IWRF, 'corforc', str1, frc1, frc3, fdefct,
     $   ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
        call DAXPY( 3*natm, one, frc1,1, frctot,1 )
        do  jd=1,ndim
          do  id=1,ndim
            strtot(id,jd) = strtot(id,jd) + str1(id,jd)
          enddo
        enddo
c
        if( lvlout.gt.2 ) call TIMER('After corforc    ')
      endif
c
c Reconstruct grid fields needed for vslofrc (eslofrc?):
c
      if( .not. dolocxc )then
c       Cutoff reference atom Vxc will *not* be taken care of
c       in the 3ctr routines, and, hence, must be done in vslofrc
c       Add Vxc0 into del-Vxc to get total Vxc (up and dn)
        call DAXPY( nptr, one, wk(i03),1, wk(i04),1 )
        if( do_spin )
     $  call DAXPY( nptr, one, wk(i03),1, wk(i05),1 )
      endif
c
c  Move es grid potential into boxes:
      call INTOBOX( n1r,n2r,n3r, wk(i02), ndbox, wk(i01) )
c
c Compute phantom force from spherical xc atoms, and grid es force:
c
      call ESLOFRC( IWR, do_gga, ndim,
     $ frc1,str1, frc2,str2, frc3,str3,
     $ xcalcut,xcfac,xcrhocut,  natm,ntyp,nshld,nald, nlat1c,
     $ itypa, numshl,lshel,nala,ala,cala, nocc,locc,nalsh,alsh,calsh,
     $ znuc,occsh,  nrd,nrcor,radmsh,corden, wspl(1,1),wspl(1,2),
     $ ratm,rlat,
     $ hh, weight, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1), boxrad,
     $ wkrbx,     wksma,
c -->  wkb(10B)-s frctp
     $ wk(i01), wk(i06) )
c -->  ves-io   r-3s
c
      write(IWR,*)
      write(IWR,9020) 'atom xc0frc stress term=', strxc0/three
c
      call SYMFRC( IWRF, 'xc0frc', str2, frc2, wksma, fdefct,
     $ ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
      call DAXPY( 3*natm, one, frc2,1, frctot,1 )
c
      if( anycore )then
        write(IWR,9020) 'atom xc0core stress term=', strxcc/three
c
        call SYMFRC( IWRF, 'xc0core', str3, frc3, wksma, fdefct,
     $   ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
        call DAXPY( 3*natm, one, frc3,1, frctot,1 )
      endif
c
      call SYMFRC( IWRF, 'eslofrc', str1, frc1, wksma, fdefct,
     $ ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
      call DAXPY( 3*natm, one, frc1,1, frctot,1 )
      do  jd=1,ndim
        do  id=1,ndim
          strtot(id,jd) = strtot(id,jd) + str1(id,jd)
     $                    + str2(id,jd) + str3(id,jd)
        enddo
        strtot(jd,jd) = strtot(jd,jd) + (strxc0+strxcc)/three
      enddo
c
      if( lvlout.gt.2 ) call TIMER('After eslofrc    ')
c
c  Set up scratch memory for "vslofrc", making sure do not
c  stomp on memory I care about.
c
      memdph = (12*(ncplx*nk*mxrb)+1) / 2
      if( idwf .eq. 0 )then
c       On-the-fly grid orbitals, try *after* iwf, as iwfw and iwf
c       are packed back-to-back (see "wkmem")
        idph = iwf + memwf
        lastdph = idph + memdph - 1
        if( lastdph .gt. maxwkd )then
          write(IWR,*) '>>>>> maxwkd: have, need=',maxwkd,lastdph
          call STOPXERR( 'vslofmem/out of memory for grid forces' )
        endif
      elseif( idwf .eq. 1 )then
c       Stored grid orbitals, tuck between iwfw and iwf, or after iwf:
        idph = iwfw + memwfw
        lastdph = idph + memdph - 1
        if( lastdph .gt. (iwf-1) )then
c         After iwf ...
          idph = iwf + memwf
          lastdph = idph + memdph - 1
          if( lastdph .gt. maxwkd )then
            write(IWR,*) '>>>>> maxwkd: have, need=',maxwkd,lastdph
            call STOPXERR( 'vslofmem/out of memory for grid forces' )
          endif
        endif
        call MEMCHKIWF( idwf, 1, 'b4 forces ', wk(iwf), wfiwf )
      else
        call STOPXERR( 'idwf-frc/B4 vslofrc, idwf is not (0 or 1 )' )
      endif
c
      if( iproc .eq. master )then
        REWIND( unit=idmatfl )
        if( do_spin ) REWIND( unit=idmatsfl )
      endif
      idmfile = idmatfl
      i04i05 = i04
c     Position Ham scratch directly behind stored Hams:
      i03s = i01 + nmatskp
      do 2300 ispin=1,nspin
        if( ispin.eq.2 .and. elecno(ispin) .eq. zero ) goto 2300
c       Move grid potential into boxes ...
        call INTOBOX( n1r,n2r,n3r, wk(i04i05), ndbox, wk(i06) )
c
c       Compute forces from grid potentials
c
        call VSLOFRC( idmfile,  ndim,ncplx,idwf,doblas3,
     $   nafrc,iforce, frc1,str1,
     $   norb,nk, natm,ntyp, nshld,nald, nlat1c,
     $   itypa,numshl,lshel,nala,ala,cala, alamin,
     $   ratm, rlat, wtk,wk(icoskr),wk(isinkr),
     $   weight, hh, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1),boxrad,
     $   wkrbx(1),   wkrbx(1+3*mxrb), wkrbx(1+4*mxrb), wksma,
c -->    rmesh(3B)-s vbox(1B)-s       wkbox(7B)-s      frctp(3a)
     $   wk(iwf),wk(iwf),
c -->    wvfcns  wf
     $   wk(iwfw),           wk(idph),
c -->    wkbi(nc*nk*nob*B)r4 dphi(12*ncplx*nk*B)r4
     $   wk(i03s)  , wk(i01),     wk(i06) )
c -->    dmatrd-ks   dmat-2s(r4)  pot-i
c
        call SYMFRC( IWRF, 'vslofrc', str1, frc1, frc3, fdefct,
     $   ndim, natm, nsyma, rprim, rmatsym, naofsym )
c
        call DAXPY( 3*natm, one, frc1,1, frctot,1 )
        do  jd=1,ndim
          do  id=1,ndim
            strtot(id,jd) = strtot(id,jd) + str1(id,jd)
          enddo
        enddo
c        .... switch to down-spin, if we have it
        idmfile = idmatsfl
        i04i05 = i05
 2300 continue
c
      if( lvlout.gt.2 ) call TIMER('After vslofrc    ')
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c        Do contributions to force from (semi-)analytic pieces
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( DO_VDW() )then
c       Compute forces/stresses from vdW-FF correction
        call VDWFF( engyvdw, frc1,str1, c6vdw,r0vdw,c6ijvdw,r0ijvdw,
     $   natm,ntyp,nlat,nlat2c, itypa,znuc,alamin,
     $   ratm,rlat )
c
        IWRvdw = IWR
        call SYMFRC( IWRvdw, 'vdW-FF ', str1, frc1, frc3, fdefct,
     $   ndim, natm, nsyma, rprim, rmatsym, naofsym )
        call DAXPY( 3*natm, one, frc1,1, frctot,1 )
        do  jd=1,ndim
          do  id=1,ndim
            strtot(id,jd) = strtot(id,jd) + str1(id,jd)
          enddo
        enddo
        call TIMER( 'After vdW forces' )
      endif
      
c Tell "mat0frc" how much memory it has to work with in dmat()
      maxwk = maxwkd - i01 + 1
c
c Tell the following whether to worry about spin-dmat:
      nspinf = 1
      if( do_spin .and. elecno(2).ne.zero ) nspinf = 2
c
      call MAT0FRC( IWR,IWRF, idmatfl,idmatsfl,iematfl,iematsfl,
     $ lvlout, nspinf, ndim, dolocxc,
     $ mat,nmat, convs,convsl,convii,
     $ nafrc,iforce, frctot,strtot,fdefct, frc1,frc2,frc3, str1,str2,
     $ norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $ itypa, numshl,lshel,nala,ala,cala,
     $ nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $ alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $ nrd, nrad,nrps,radmsh,radwt,
     $ vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad,
     $ ratm, rlat, wk(icoskr),wk(isinkr), nsyma,rprim,rmatsym,naofsym,
     $ wksm, naldsq, gntcomb, maxwk,mxang,
c -->  wksm(46*naldsq)-s
     $ wk(i01) )
c -->  dmat-s
c      spin: [nmatpr+{nk:nmatstr+2mat,1k:nmatstr+mat,gamma:2mat}]
c    nospin: [nmatpr+nmatstr+mat]
c
c Output final total force of this geometry:
c
      call SYMFRC( IWR, 'total', strtot, frctot,frc3, fdefct,
     $ ndim, natm, 0, rprim, rmatsym, naofsym )
c
      engypv = zero
      if( ndim.eq.3 )then
c       Put out stress tensor in condensed form:
        write(IWR,2990) strtot(1,1), strtot(2,2), strtot(3,3),
     $                  ( strtot(1,2)+strtot(2,1) ) / two,
     $                  ( strtot(2,3)+strtot(3,2) ) / two,
     $                  ( strtot(3,1)+strtot(1,3) ) / two
 2990   format(1x,'STR(Ry)=',6f10.6)
        engypv = (strtot(1,1)+strtot(2,2)+strtot(3,3)) / three
      endif
c
c Save status as force done:
      call QHIST( ihistfl, 'FORC' , natm , rdata, frctot )
      call WRSTAT( IWR, istatfl, 'FORCE',
     $ do_geom, do_cell, opt_spin,
     $ iucstep,igstep,iterscf, engytotl,excengy,esnsengy,
     $ istart_sp, elecno, spindata,
     $ strtot, natm, frctot )
c
      if( ndim.eq.3 .and. do_pvcalc )then
c       We want an enthalpy, E+PV
        enthalpy = engytotl - engypv
        write(IWR,9020)'TOTAL Pressure X Volume (Ry)=    ',-engypv
        write(IWR,9020)'TOTAL ENTHALPY (E+PV) (Ry)=      ',enthalpy
      endif
c
      call TIMER('AFTER FORCE RTNS ')
c
c *********************************************************************
c **********         end of force calculation phase         ***********
c *********************************************************************
c **********       start of geometry relaxation phase       ***********
c *********************************************************************
c
 3000 continue
      if( .not. do_geom )then
        call QHIST( ihistfl, 'END ' , ndata, rdata, wksm )
        call STOPXOK( 'FORCE DONE' )
        goto 999
      endif
c
      write(IWR,*) ' '
c Reset control flags for follow-on geoms:
cpas: redundancy in control flag reset needs to be resolved
      dosetup = .true.
      doiters = .true.
      doforce = .true.
      itconv = 1
      itstart = 0
      istart_sp = 0
c
c If we are doing scf-guessing, change scf blend for follow-on geoms:
      if( igesopt .eq. 2 ) scfblnd = scfbl2
c
      if( do_cell .and. .not. (dorelax .or. do_neb ) )then
c       Must transform coordinates appropriately, and skip
c       atomic relaxation ...
cpas:   Must resolve restart/status stuff, too ...
c
c      ... first transform atomic coordinates into external frame:
        call REXTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $   ratm,ratm0 )
c
c   ... and also must transform charge/ws origin vectors:
        if( ion_opt.ne.0 )then
          do  j=1,3
            origws0(j) = origws(j)
C            rchrg0(j) = rchrg(j)
          enddo
          call DCOPY( 3*nchrgd, rchrg,1, rchrg0,1 )
          call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $     origws,origws0 )
          call REXTRNL( ndim,  nchrgd, dolvatm, rprim,orig,rscale,
     $     rchrg,rchrg0 )
        endif
c
        goto 5000
      endif
c
c Keep track of recent total energies (not used right now):
c
      if( igeom.le.mxgeom )then
        gengy(igeom) = engytotl
      else
        do  ig=1,mxgeom-1
          gengy(ig) = gengy(ig+1)
        enddo
        gengy(mxgeom) = engytotl
      endif
c
c Install force defect into specified atom, if viable:
c
      call FDEFECT( IWR, ifdefct, natm, frctot, fdefct,
     $ nsyma,rmatsym,naofsym )
c
      if( do_neb )then
c
c       ***** NEB multiple-image relaxation *****
c
        call NEBRUN( IWR,IWRNEB, inebstat, ierrg,
     $   do_neb, image, engytotl,
     $   igstep,igstop,nhistg, gconv,gblend, tstep,
     $   nimg_neb, nimg_max, ea_neb,eb_neb, spring_neb,
     $   Lantikink_neb,Loptall_neb,Lclimb_neb,Lnoninteract,
     $   itypa,atmass, nafrc, iforce,
     $   nslics,islics,islicats,vslics,    iatframe,
     $   ndim, dolvatm, rprim,orig,rscale, nsym,nsyma, rmatsym, naofsym,
     $   natm, imgstat, engyneb, frctot, ratm,ratma,ratmb, wksma,
     $   maxwkd, wk(1),
c -->            tmp()-s
     $   igrouplist_neb, ngroup_neb )
c
        igstat = inebstat
c
        if( inebstat.gt.0 )then
          call STOPXOK( 'END NEB' )
          goto 999
        elseif( inebstat.lt.0 )then
          call STOPXERR( 'ERROR IN NEB' )
          goto 999
        endif
c
c       End multiple-image branch of geometry update
      else
c
c       ***** Single geometry relaxation *****
c
        call RUNSTEER( scfconv, gconv,  cellconv,
     $                 itstop,  igstop, maxucstep )
c
        if( iproc.eq.master )then
c         Only master computes update,, exeryone else finds out later
c
c         Impose any frame constraint:
c
          call FIXFRAME( natm, iatframe,
     $     ratm0, frctot )
c
c       Set up geom/force array to be optimized:
c
          call GSELECT( natm, nafrc, iforce, nrlxdim,
     $     nslics,islics,islicats,vslics,
     $     wksma(1,19),wksma(1,22),
c -->      ratmP(3a)-o frcP(3a)-s
     $     itypa, atmass,  natmg,wksma(1,16),
c -->                            nxatm(a)-o
     $     ratm0,frctot, wksma(1,1), wksma(1,4), wksma(1,7) )
c -->            frc-i   ratmg(3a)-o frcg(3a)-o  amassg(3a)-o
c
          if( do_md )then
c            Molecular dynamics, test for termination
             call QMDCONVCHK( igstat )
          else
c           Geometry minimization, test for geometry convergence.
c
            call CONVCHK( IWR, 'force', igstat, igstep,igstop, gconv,
     $       chgmax,chgrms,  3,nafrc,  wksma(1,4) )
c -->                                  frcg(3a)-io
            call QHIST( ihistfl, 'FMAX' , ndata, chgmax, wksm )
            call QHIST( ihistfl, 'FRMS' , ndata, chgrms, wksm )
          endif
c
          if( igstat.lt.2 )then
c           We are not yet converged ...
c
c           Blend to get new geometry
c
            call GRELAX( IWR, gblend, tstep, engytotl,
     $       igstep,nhistg, nrlxdim, chgmax,
     $       natmg,watmg,wksma(1,16),
c -->                    nxatm(a)-i
     $       wksma(1,1),  wksma(1,4), atveloc,     wksma(1,7),
c -->        ratmg(3a)-io frcg(3a)-i  vatmg(3a)-io amassg(3a)-i
     $       wksma(1,10),wksma(1,13), maxwkd,wk(1) )
c -->        wk3(3a)-s   wk4(3a)-s
c
c           Take new geometry, check symmetry, and update
c           the master coordinate vector
c
            call GUPDATE( ierrg, ratm0,
     $       natm, iforce, ndim,nsym,nsyma, rprim, rmatsym, naofsym,
     $       nslics,islics,islicats,vslics,    iatframe,
     $       wksma(1,19),
c -->        ratmpP(3a)-i
     $       wksma(1,1),  wksma(1,4), wksma(1,7),  wksma(1,10) )
c -->        ratmg(3a)-is ratmx(3a)-s dratm0(3a)-s dratmsy(3a)-s
c
            if( ierrg .ne. 0 ) igstat = -1
          endif
c
c         Master has computed convergence, and new coordinates ...
        endif
c
        if( nodes.gt.1 ) then
c          ... tell everyone the news:
          call MPBCASTI( master, 1, igstat, icomm )
          if( igstat .lt. 0 ) call STOPXERR( 'geom update error' )
          lenmsg = 3*natm
          call  MPBCAST8( master, lenmsg, ratm0, icomm )
        endif
        if( igstat .lt. 0 ) call STOPXERR( ' geom update error' )
c
c       Now we take our new/converged coordinates, and put back
c       into the user external frame ...
c
        call REXTRNL( ndim,natm, dolvatm, rprim,orig,rscale,
     $   ratm,ratm0 )
c
c       End single-image branch of geometry update.
      endif
c
c   ... and also must transform charge/ws origin vectors:
      if( ion_opt.ne.0 )then
        do  j=1,3
          origws0(j) = origws(j)
C          rchrg0(j) = rchrg(j)
        enddo
        call DCOPY( 3*nchrgd, rchrg,1, rchrg0,1 )
        call REXTRNL( ndim,  1 , dolvatm, rprim,orig,rscale,
     $   origws,origws0 )
        call REXTRNL( ndim, nchrgd, dolvatm, rprim,orig,rscale,
     $   rchrg,rchrg0 )
      endif
c
c  If doing QMD: output the geometry and other info now:
c
      call MPNODE0( master )
      call MPNODE( iproc )
      if( do_md .and. iproc .eq. master )then
c       Only dump this if you are master processor (contention for lcao.vxyz)
        call QMDPRINT( ratm )
      endif
c
c  Write out new geometry to user listing file:
c
      call GSTAMP( IWR, igstep,igeom,igstat,image,
     $ iatmfmt, ntyp,natm, natmnm, itypa, atmnm,typnm, ratm )
c
c Save status/new geometry:
      call WRSTAT( IWR, istatfl, 'GEOMETRY',
     $ do_geom, do_cell, opt_spin,
     $ iucstep,igstep,iterscf, engytotl,excengy,esnsengy,
     $ istart_sp, elecno, spindata,
     $ rprim, natm, ratm )
c
c Are we converged, out of steps, done, or do we continue?
c
      if( igstat.eq.0 )then
c       Not converged geometry, go do scf/force with next geometry
        call TIMER('After geometry   ')
        if( igstep.le.igstop ) goto 200
c
      elseif( igstat.eq.1 )then
c       Not converged geometry, and out of steps
        write(IWR,9491)  engytotl
 9491   format(/1x,'FINAL *NOT* RELAXED ENERGY =',f20.10 / )
c
      elseif( igstat.eq.2 )then
c       Converged geometry
        write(IWR,9492)  engytotl
 9492   format(/1x,'FINAL RELAXED ENERGY =',f20.10 / )
c
      endif
c
      if( .not. do_cell )then
c
c       See if further instructions for relaxation:
c
        dorelax = .false.
        call GEOMDAT( IWR,IDAT, dorelax,
     $   natm,nafrc,iforce,ifdefct,
     $   nhistg,gblend,gconv,tstep,igstart,igstop,
     $   nslics,islics,islicats,vslics,    iatframe )
c
        if( dorelax )then
          igstep = 1
          igstat = 0
          goto 200
        endif
c
        call TIMER('AFTER GEOMETRY   ')
c
        call QHIST( ihistfl, 'END ' , ndata, rdata, wksm )
        call STOPXOK( 'RELAX DONE' )
        goto 999
      endif
c
c *********************************************************************
c **********           end atom relaxation phase            ***********
c *********************************************************************
c **********         start of cell relaxation phase         ***********
c *********************************************************************
c
 5000 continue
      if( .not. do_cell ) goto 6000
c
c Keep track of recent total cell energies (not used right now):
c
      if( iucstep.le.mxgeom )then
        ucengy(iucstep) = engytotl
      else
        do  iuc=1,mxgeom-1
          ucengy(iuc) = ucengy(iuc+1)
        enddo
        ucengy(mxgeom) = engytotl
      endif
c
      call QHIST( ihistfl, 'ECEL' , ndata, engytotl, wksm )
c
c Do the cell optimization:
c
      call RUNSTEER( scfconv, gconv,  cellconv,
     $               itstop,  igstop, maxucstep )
c
      call CELLOPT( iucstat, ihistfl,
     $ do_cell, ndim, iuctype, icellscheme,strsbroy,
     $ iucstep,iuciter,maxucstep,nhistuc, cellconv,strxtrnl,
     $ strnmax, ucstepfac,cellstepdecr,cellstepincr, ucblend,
     $ strlast, engytotl,
     $ rprim, orig,rscale, strtot, nsym,rmatsym,
     $ nk,veck, dolvatm, natm,ratm,ratm0 )
c
c Record new cell vectors to output file:
c
      call UCSTAMP( IWR, iucstat, iucstep,iuciter, rprim )
c
c Decide the next step:
c
      if( iucstat.eq.0 )then
c       Not converged, go to next cell step:
c        ... save status/new geometry:
        icellstep = iucstep
        if( iuciter.gt.0 ) icellstep = iuciter
        icellstep = icellstep + 1
cpas:   Right now, just pick up cell vectors, with no history
        icellstep = 1
        call WRSTAT( IWR, istatfl, 'CELL',
     $   do_geom, do_cell, opt_spin,
     $   icellstep,igstep+1,iterscf, engytotl,excengy,esnsengy,
     $   istart_sp, elecno, spindata,
     $   rprim, natm, ratm )
c
        call TIMER('After cellopt    ')
        write(IWR,*)
        goto 100
c
      elseif( iucstat.eq.1 )then
c       Not converged cell, and out of steps
        write(IWR,9591)  engytotl
 9591   format(/1x,'*NOT* RELAXED CELL ENERGY =',f20.10 / )
c
      elseif( iucstat.eq.2 )then
c       Converged cell
        write(IWR,9592)  engytotl
 9592   format(/1x,'FINAL RELAXED CELL ENERGY =',f20.10 / )
        if( ndim.eq.3 .and. do_pvcalc )then
c         We want an enthalpy, E+PV
          write(IWR,9020)'TOTAL Pressure X Volume (Ry)=      ',-engypv
          write(IWR,9020)'TOTAL RELAXED ENTHALPY (E+PV) (Ry)=',enthalpy
        endif
c
      endif
c
      call TIMER('AFTER CELLOPT    ')
c
c *********************************************************************
c **********           end cell optimization phase          ***********
c *********************************************************************
c **********           start band structure phase           ***********
c *********************************************************************
c
 6000 continue
      if( .not. do_nonscf ) goto 7000
c
      call TIMER( 'Starting non-scf code' )
c
      write(IWR,*)
     $ '>>>>> AT NON-SCF: do_nonscf, do_optic, do_bands =',
     $                    do_nonscf, do_optic, do_bands
c
c     Non-SCF portion of code (Art Edwards' band structure code)
c
c     We first need the new set of k-points. We will branch to two
c     different routines. For now, we will only opt for do_bands=.true.
c     and put the gamma point into the calculation.
c
      if( do_optic )then
        call STOPXERR( 'Optical properties not (yet) implemented' )
      endif
c
      if( do_bands )then
        call TIMER( 'Starting bandkpts' )
        if( bvl .eq. 'nil' )then
c         Forgot to give a band input!
          write(IWR,*) 'Band structure, but no bands input'
          call STOPXERR( 'do bands, but no input bands' )
        endif
c
c       Output data needed for band structure
        scfefermi = efermi(1)
        efev = efermi(1)*toev
        call FLOPENA( IBNDFL, 'bands' )
        REWIND( IBNDFL )
c
c       Generate new sampling k-vectors for bands.
        call BANDKPTS( veckn,veck,wtkn,nknscf,nknscfd,nkd,bvl,
     $       nbranch, nbranchd,bsymb,rprim,gprim,dk,nbr,IWR )
c
        write(IWR,*) 'After bandkpts: total band k-pts nknscf =',nknscf
        write(IWR,*) 'Here are the new k-vectors (??? UNITS ???)'
        do  ik=1,nknscf
          write(IWR,*) (veckn(ic,ik), ic=1,3), wtkn(ik)
        enddo
        nstate = MIN( nocc0k+6, norb )
        call TIMER( 'After bandkpts' )
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c           Outer loop for non-self-consistent calculations
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The strategy is to place small blocks of veckn into a smaller array
c     that can be sent transparently into the appropriate routines.
c
c     We know processing up to nk k-points in parallel works to this point.
c     For small problems, we can probably do blocks up to nkd.
c     For large problems, nk probably strains limits of what can be done.
c     For now, we use the smaller of nk and nknscf until we
c     come up with logic to pick the largest safe block size.
c     (or have it as input?)
c
c  Now we set up the old-style while loop over these blocks:
c   nkblocksz = the number of k-points to be processed in one shot.
      nkblocksz = MIN( nk, nknscf )
c   we certainly can do (at least) 1k/proc:
      nkblocksz = MAX( nkblocksz, nprocs )
c    ... but can only do what is available ...
      if( nkblocksz .gt. nk ) nkblocksz = nk
      write(IWR,*) 'nk, nkblocksz, nknscf (initial)=',
     $              nk, nkblocksz, nknscf
c   iblockoffset = offset of current block of k-points in list of k-points from bandkpts
      iblockoffset = 0
c   nleft = the number of k-points left to do. It is decremented in the loop below.
      nleft = nknscf
c
c  Set up storage of band eigenvalues
      nbands = nstate
      ieigval = 1
      neigval = nbands*nkblocksz
      iwk1nonscf = ieigval + neigval
c
c  Return here for old-style while loop, same as F90 do ... while (nleft.gt.0)
c
 6100 continue
        write(IWR,*) 'nonscf k-blocks loop: nk, nkblocksz, nknscf=',
     $                                      nk, nkblocksz, nknscf
c
c       Load a block of k-points to process.
        do  jvec=1,nkblocksz
          do  i=1,3
            vecknsm(i,jvec) = veckn(i,iblockoffset+jvec)
          enddo
          wtknsm(jvec) = wtkn(iblockoffset+jvec)
        enddo
c
c       Set up space for band eigeenvalues (nkblocksz might have changed)
        neigval = norb*nkblocksz
        call MKZERO( neigval, wk(ieigval) )
c
c       Solve for these k-points non-self-consistently.
c
        call NONSCFK( ala,alamin,almnnl,alsh,atm0engy,boxrad, 
     $       cala,calsh,convs,convii,convsl,do_gga,
     $       do_spin,doblas3,dokdefct,do_bands,do_optic,scfefermi,
     $       edegen,efermi,egap,elecno,etemp,
     $       excrad,gntcomb,hh,icloseocc,icluster,
     $       idhamfl,idmatfl,idmatsfl,iematfl,iematsfl,
     $       idosolv,
     $       igrd0fl,igridfl,iham0fl, ivxc0fl, iexc0fl, ih0tfl,
     $       ih0nlfl, ikdefct,inboxs,
     $       inrec,iovlpfl,irhoopt,ivecfl,itypa,
     $       IWR,lmxnlp1,locc,lshel,lvlout,
     $       Lijkmat,
     $       mat,mxang,maxvrad,maxwkd,mxrb,mxwkrbx,n1r,n2r,n3r,nala,
     $       nald,naldsq,nalsh,natm,natmd,nbox,
     $       ncplx,ndbox,nearatm,nearopt,
     $       near2c,nhiste, nkblocksz,
     $       nlat,nlat1c,nlat2c,nlat3c,nlatd,
     $       noad,nocc,norb,norba,norbd,nptr,nrad,nrd,nrec,nrecd,nrps,
     $       nshld,nspin,nbands, nstbulk,ntyp,ntypd,numshl, occ,
     $       occij,occsh, radmsh,radwt, ratm,rlat,
     $       rprim,scut,vatrad,vecknsm,
     $       vesrad,vpsrad,vxcrad, weight,
     $       wk(ieigval), wk, iwk1nonscf, wkrad,
     $       wkrbx,wksm,wksmo,wspl,
     $       wtknsm,znuc, masterlist,
     $       do_kppsolve,do_psolv )
c
       call TIMER( 'after nonscfk' )
c
c       Write out eigenvalues and Fermi level of this block of k-points
c       for band structure plot if requested
        if( do_bands )then
          iev0 = ieigval - 1
          call DSCAL( neigval, toev, wk(ieigval),1 )
          do  ik1=1,nkblocksz
            write(IBNDFL,9610) dk(ik1+iblockoffset),efev,
C     $           (eigenvals(iv,ik1,1)*toev,iv=1, nbands)
     $           ( wk(iev0+iev), iev=1,nbands )
            iev0 = iev0 + nbands
          enddo
        endif
c       30*nbranchd - huh?
 9610   format( 360( 1x, f11.5 ) )
c
c       Figure out how many more k-points to do, and then go do them ...
        nleft = nleft - nkblocksz
        if( nleft .gt. 0 )then
c         Still have work, increment offset by old block size
          iblockoffset = iblockoffset + nkblocksz
c         Update block size
          nkblocksz = MIN( nleft, nkblocksz )
          goto 6100
        endif
c
c  Non-scf loop done, wrap up BAND file
      if( do_bands )then
C        CLOSE( IBNDFL,status='keep' )
        call FLCLOSE( IBNDFL )
        call TIMER( 'Band structure done' )
      endif
c
      call TIMER( 'Non-scf code done' )
c
c *********************************************************************
c **********            end band structure phase            ***********
c *********************************************************************
c
c *********************************************************************
c **********               End of main program              ***********
c *********************************************************************
c
 7000 continue
c
      call QHIST( ihistfl, 'END ' , ndata, rdata, wksm )
      call STOPXOK( 'QUEST DONE' )
c
  999 continue
c
      STOP
      END
********************************************************************************
c NONSCF section
c
c 10Jan13-PAS: Important notes on this NONSCFK routine
c  NONSCFK contains serial/mpi-dependent calls (to eigesolver routines)
c  In addition, this duplicates much of the main code logic in building H
c  Hence, the NONSCFK routine is appended to main routine.
c
********************************************************************************
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> NONSCFK
c
c
      subroutine NONSCFK( ala,alamin,almnnl,alsh,atm0engy,boxrad, 
     $ cala,calsh,convs,convii,convsl,do_gga,
     $ do_spin,doblas3,dokdefct,do_bands,do_optic,scfefermi,
     $ edegen,efermi,egap,elecno,etemp,
     $ excrad,gntcomb,hh,icloseocc,icluster,
     $ idhamfl,idmatfl,idmatsfl,iematfl,iematsfl,
     $ idosolv,
     $ igrd0fl,igridfl,iham0fl, ivxc0fl, iexc0fl, ih0tfl,
     $ ih0nlfl, ikdefct,inboxs,
     $ inrec,iovlpfl,irhoopt,ivecfl,itypa,
     $ IWR,lmxnlp1,locc,lshel,lvlout,
     $ Lijkmat,
     $ mat,mxang,maxvrad,maxwkd,mxrb,mxwkrbx,n1r,n2r,n3r,nala,
     $ nald,naldsq,nalsh,natm,natmd,nbox,
     $ ncplx,ndbox,nearatm,nearopt,
     $ near2c,nhiste,nk,
     $ nlat,nlat1c,nlat2c,nlat3c,nlatd,
     $ noad,nocc,norb,norba,norbd,nptr,nrad,nrd,nrec,nrecd,nrps,
     $ nshld,nspin,nstate, nstbulk,ntyp,ntypd,numshl, occ,
     $ occij,occsh, radmsh,radwt, ratm,rlat,
     $ rprim,scut,vatrad,vecknsm,
     $ vesrad,vpsrad,vxcrad, weight,
     $ eigenvals, wk, iwk1nonscf, wkrad,
     $ wkrbx,wksm,wksmo,wspl,
     $ wtknsm,znuc, masterlist, 
     $ do_kppsolve,do_psolv )
c---------------------------------------------------------------
c Purpose: control routine for non-scf (post-scf) calculation of ham(k)
c
c Written: Arthur H. Edwards, for 2.64 insertion
c     This is the control routine for non-SCF calculations. We will require repartition
c     of the main workspace, a recalculation using MAT0SET, ...
c     Remember to dimension everything that is brought down.
c
c Revision history:
c  03May12-ACP: Merged into 2.64 (17Jan13-PAS: renamed/reintegrated)
c  02May12-ACP: Cleaned up: Now use alternate entry point to solvers to skip calculation
c               of eigenvectors when not needed. (Previous revision computed eigenvectors
c               even when they were not needed. (Only do_optic needs them and that is not implemented yet).
c  26Apr12-ACP: Bug fixes: original version only worked if passed all the k-points.
c               It failed in POP where the set of k-points needed to have weights
c               that summed to one. Also had a compiler dependent bug relating to
c               rewinding of a file (run-time error). This version also redistributes
c               processors depending on the number of k-points being processed. To
c               take full advantage of this a new routine is needed for the calling
c               calling program to compute the maximum number of k-points that will
c               run in available memory.
c  20Dec10-ACP: Fully merged into 2.62dev. Uses k-parallel versions
c               of routines developed for 2.62dev. 
c  2009-10-AHE: Written for 2.61b.
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
c
c NOTES: here veck is sent down from vecknsm
c
c  Task-parallel data structures (should be fused into serial mem):
c     Dimensioning. Be careful that the dimensions of past arrays
c     match those in lcaomain.f
c
      LOGICAL   do_gga, do_spin, doblas3,dokdefct
      LOGICAL   Lijkmat
      LOGICAL   do_kppsolve, do_psolv
      LOGICAL   do_bands, do_optic, Lnonscf
      LOGICAL   do_eigvecs, do_eigpops
      DIMENSION ala(nald,nshld,ntypd),alamin(ntypd),almnnl(ntypd)
      DIMENSION alsh(nald,nshld,ntypd),cala(nald,nshld,ntypd)
      DIMENSION calsh(nald,nshld,ntypd),convs(8), efermi(2),egap(2)
      DIMENSION elecno(2),excrad(nrd,ntypd),hh(3,3), itypa(natmd)
c      DIMENSION iwork(norb,*)
      DIMENSION lmxnlp1(ntypd),locc(nshld,ntypd),lshel(nshld,ntypd)
      DIMENSION nala(nshld,ntypd),nalsh(nshld,ntypd),nocc(ntypd)
      DIMENSION norba(ntypd),nrad(ntypd),nrps(ntypd),nstbulk(2)
      DIMENSION numshl(ntypd),occsh(nshld,ntypd)
      DIMENSION occij(nshld,nshld,ntypd),radmsh(nrd,ntypd)
      DIMENSION radwt(nrd,ntypd),ratm(3,natmd),rlat(3,nlatd),rprim(3,4)
      DIMENSION vatrad(nrd,ntypd),vesrad(nrd,ntypd),vxcrad(nrd,ntypd)
      DIMENSION vpsrad(nrd,4,ntypd)
      DIMENSION wk(maxwkd),wkrad(nrd,maxvrad),wkrbx(mxwkrbx)
      DIMENSION wksm(46*naldsq)
cAHE
c    Remove dimensioning for h0scf and h0nscf
      DIMENSION h0scf(mat),h0nscf(mat)
c
      DIMENSION wksmo(norbd,9),wspl(nrd,2),wtknsm(nk),znuc(ntypd)
c  Work array for angular stuff needed in non-local integrals:
      DOUBLE COMPLEX    gntcomb(9,5,5,2,9)
      PARAMETER  ( maxnodes =  512 )
      INTEGER    icluster
      DOUBLE PRECISION  gapp
      DIMENSION  icluster(2*maxnodes)
      DIMENSION  gapp(maxnodes)
c  parallel-parallel solver data structures (PEIGSOLVKP/PSCHROEDKP)
      INTEGER    masterlist
      DIMENSION  masterlist(maxnodes)
c  Arrays for box stuff:
      DIMENSION  ndbox(2,2,3), inboxs(nrecd)
      DIMENSION  nearatm(3,natmd)
c
c     Arrays for non-SCF calculations
c
      DIMENSION vecknsm(3,nk),eigenvals(nstate,nk,2)
c
c*****************************************************************
c >>>> Executable code
c*****************************************************************
c
c     do_bands does not require eigenvectors
c     do_optic does require eigenvectors
c     Neither require the solver to compute populations.
      do_eigvecs = do_optic
      do_eigpops = .false.
      if( nk .eq. 1 )then
c       A single point: not band but probably gamma, likely want occupations, too ...
        do_eigpops = .true.
      endif
      Lnonscf = .true.
c
      if( lvlout.gt.2 )then
        write(IWR,*)'NONSCFK: processing new block of non-scf k-points'
        write(IWR,*)'NONSCFK do_bands, do_optic=',
     $                       do_bands, do_optic
        write(IWR,*)'NONSCFK do_eigvecs, do_eigpops=',
     $                       do_eigvecs, do_eigpops
        call flush(IWR)
      endif
c     Free the k-parallel communicator as we don't necessarily have the same number
c     of k-points as in the main code or the previous loop.
CTP:  Comment out MPFREECOMK for serial code.
      call MPFREECOMMK()
      call MPNODE0( master )
      call MPNODE( iproc )
c     Set up new k-communicators.
      kparopt = 0
      call KPSETUP( nk, kparopt )
      call KPFLAG_GET( kparopt )
      if( kparopt.ne.2 )then
c       Not k-parallel, disable parallel-parallel solver
        do_kppsolve = .false.
      else 
c       See if KPSETUP gave us a processor distribution that benefits
c       from using the parallel-parallel solver.
        call MPNODES( nodes )
        call MPNODES_K( nodes_k )
c       If all processors are in the same k-point communicator, then
c       there is no advantage to using the parallel-parallel solver.
        if( nodes_k .eq. nodes ) do_kppsolve = .false.
cxxx264  Unless there is at least 4 procs in every parallel solve, turn off kpp
cxxx264  (if proc allocations diff from kp-allocs, can get different kpp criterion)
        if( (nodes/nk) .lt. 4 ) do_kppsolve = .false.
c       Turn off kpp, if problems are too small (adapting acp criterion)
C        if( norb .lt. 200*(nodes/nk+1) - acp criterion, seems mysterious
        if( norb .lt. 800 ) do_kppsolve = .false.
      endif
c     Determine arrays sizes for this problem:
      call KPMINE( nk, nkloc, nk0 )
c      mat = norb*norb
      nmat = nk*mat
      nmatloc = nkloc*mat
      nmatskloc = nmatloc*nspin
      i00 = iwk1nonscf
c
c extra args : Lijkmat, nmatloc, nkloc, nhiste, ntyp, itypa, norba, i13
c
      call WKMEM( IWR, lvlout, maxwkd,
     $ do_gga,do_spin,doblas3, Lijkmat,
     $ idosolv, do_kppsolve, do_psolv, do_eigvecs, Lnonscf,
     $ mat,nmat,nmatloc, nptr, nhiste,
     $ ncplx,nk,nlat, norb,nstate, noad,natm,ntyp,itypa,norba, 
     $ mxang, ndbox,nbox,inboxs,mxrb,
     $ near2c,irhoopt, idwf,  memwf,memwfw, iwf,iwfw,
     $ icoskr,isinkr, is1at,
     $ ioint,ioofk, ieigval,ieigpop,
     $ idbr, imblnde, memblnde,
     $ i00, i01s, i01,i02,i03,i04,
     $ i05,i06,i07,i08,i09,i10,i11,i12 )
c
c Compute Bloch phase factors:
c
      call MKEIKL( nk,nlat, vecknsm,rlat, wk(icoskr),wk(isinkr) )
c
c Locate nearatm mesh point to each atom in list
c
c      call NEARBY( nearatm, natm, itypa, norba,
c     $ ratm, rprim,hh,n1r,n2r,n3r,nptr, mngr3,mxgr3,ngr,
c     $ wk(i01s), wk(i01s+3*nptr) )
c
c I'm skipping routines that appear to be k-independent
c
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Compute iteration-independent (semi-)analytic matrix elements
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c Tell "mat0set" how much mem is available in wk():
      maxwk = maxwkd - i01s + 1
c
c extra args: ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,ncplx
c
      if( Lijkmat )then
c       We use the square (i,j,k) matrices
        call MAT0SET1( IWR, lvlout,
     $   iovlpfl,iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $   atm0engy, mat,nmat, convs,convsl,convii,
     $   norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $   itypa, numshl,lshel,nala,ala,cala,
     $   nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $   alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $   nrd, nrad,nrps,radmsh,radwt,
     $   vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad, wspl,
     $   ratm, rlat, wk(icoskr),wk(isinkr), wk(is1at),
     $   wksm, gntcomb, maxwk,mxang,
     $   wk(i01s) )
c -->  wkmat-s(see wkmem)
      else
c       We use the stripe-parallel matrix
        call MAT0SET( IWR, lvlout,
     $   iovlpfl,iham0fl,ivxc0fl,iexc0fl,ih0tfl,ih0nlfl,
     $   atm0engy, mat,nmat, convs,convsl,convii,
     $   norb,nk,ncplx, natm,ntyp,nshld,nald, noad, nlat,nlat2c,nlat3c,
     $   itypa, numshl,lshel,nala,ala,cala,
     $   nocc,locc,nalsh,alsh,calsh, occsh,occij,
     $   alamin, norba, occ,znuc, lmxnlp1,almnnl,
     $   nrd, nrad,nrps,radmsh,radwt,
     $   vpsrad, vatrad,vesrad,vxcrad,excrad, wkrad, wspl,
     $   ratm, rlat, wk(icoskr),wk(isinkr), wk(is1at),
     $   wksm, gntcomb, maxwk,mxang,
     $   wk(i01s) )
c -->  wkmat-s(see wkmem)
      endif
c size: 3*nmatpr+matstr*ncplx+mat
c
      if( iproc.eq.master ) REWIND( iham0fl )
c
c Put bloch functions on grid
c
c Generate and keep grid orbitals if we have space:
c
      if( idwf.gt.0 )then
         call GRIDWF( ncplx,idwf,
     $        nk, natm,ntyp, nshld,nald, nlat1c,
     $        itypa,numshl,lshel,nala,ala,cala, alamin,ratm,rlat,
     $        wk(icoskr),wk(isinkr),
     $        hh, ndbox,nbox,nrec,inrec,inboxs,
     $        inboxs(nbox+1), boxrad,
     $        wk(iwf), wk(iwf),
c-------->    wvfcns   wf
     $        wkrbx)
c-------->    wkbox(10B)-s
      endif
c
c Subtract grid overlaps from analytic to get 1-ctr "defect" overlaps:
c
      call GRDOVLP( ncplx, idwf, near2c, wk(is1at),wk(ioint),
     $ norb,nk, natm,ntyp, nshld,nald, noad, nlat1c,
     $ itypa,numshl,lshel,nala,ala,cala, alamin,norba,
     $ ratm,rlat,wk(icoskr),wk(isinkr),  vecknsm, wksmk,
     $ weight, hh, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1), boxrad,
     $ wk(iwf),wk(iwf), wkrbx(1),    wkrbx(1+10*mxrb),
c -->  wvfcns  wf       wkbox(10B)-s r(3B)-s
     $ wk(i01),      wk(i01+2*nk*mxrb) )
c -->  eikr(2B*nk)-s otmp(MAX(noad^2*nk*natm,2*norb*nk))-s
c
      if( iproc .eq. master ) REWIND( unit=igrd0fl )
c        Retrieve rho[atom] from igrd0fl ...
      call MPBCREADBIG( igrd0fl, nptr, wk(i10) )
c        Retrieve total delta-density (or up-spin) from igridfl
      if( iproc .eq. master ) REWIND( unit=igridfl )
      call MPBCREADBIG( igridfl, nptr, wk(i08) )
      call DCOPY( nptr, wk(i08),1, wk(i07),1 )
c     ... and dn-spin density:
      if( do_spin ) call MPBCREADBIG( igridfl, nptr, wk(i09) )
      call MPBCREADBIG(igridfl,nptr,wk(i03))
      call MPBCREADBIG(igridfl,nptr,wk(i04))
      if( do_spin ) call MPBCREADBIG( igridfl, nptr,wk(i05) )
c     We need vexslo
c     PAS: No, we don't, don't use energy here, we can simply pass zeros
c       However, also must spoof minority spin excslo in i11, if have spin
C      REWIND(iexcslo)
C      call READBIG(iexcslo,nptr,wk(i06))
      call MKZERO( nptr, wk(i06) )
      if( do_spin ) call MKZERO( nptr, wk(i11) )
c        Initialize delta-Hamiltonian
      call MKZERO( nmatloc, wk(i01) )
c
c     We are, at this point, including only parts of the spin-dependent
c     code. When we want to do spin-bands, we will have to return.
c     Also, we are not including external fields
c
c If doing spin, initialize both spin Hamiltonians
      if( do_spin ) call MKZERO( 2*nmatloc, wk(i01) )
c
      if( nearopt.gt.0 )then
c
c       Compute values of grid potentials "nearby" atoms:
        mxwkwt = nptr
        i01i02 = i01
        i04i05 = i04
        i06i11 = i06
        do  ispin=1,nspin
          call POTNEAR( nearopt, natm,
     $     ratm,hh,rprim, mxwkwt ,n1r,n2r,n3r,nptr, nearatm,
     $     wksmo(1,ispin), wksmo(1,3), wksmo(1,3+ispin),
c -->      vslnear-o       vesnear-o   excnear-o
     $     wk(i04i05), wk(i03), wk(i06i11), wk(i07) )
c -->      vslo-io     vesslo-i excslo-i    wtnr-s
c
c         Compute correction to matrix elements from nearby grid defect:
c
          call VSLOFIX( near2c, wk(is1at), wk(ioint),wk(ioofk),
     $     norb,nk, natm,ntyp,nshld, noad,  itypa, numshl,lshel,
     $     wksmo(1,ispin), ncplx,
c -->      vslnear-i
     $     wk(i01i02) )
c -->      v-io
c
c       ... switch to dn-spin electrons, if we have them
          i01i02 = i01 + nmatloc
          i04i05 = i05
          i06i11 = i11
        enddo
      else
c       Blank nearby arrays, so that doengy routine gets zeroes:
        call MKZERO( natm, wksmo(1,1) )
        call MKZERO( natm, wksmo(1,2) )
        call MKZERO( natm, wksmo(1,3) )
        call MKZERO( natm, wksmo(1,4) )
        call MKZERO( natm, wksmo(1,5) )
      endif
c
c Compute matrix elements of (slow-varying) grid potential:
c
      i01i02 = i01
      i04i05 = i04
c     Position Ham scratch directly behind stored Hams:
      i03s = i01 + nspin*nmatloc
cACP      if(nspin.eq.1) then
cACP        i03s=i01+nmatloc
cACP      else
cACP        i03s=i02+nmatloc
cACP      endif
      do  ispin=1,nspin
        call INTOBOX( n1r,n2r,n3r, wk(i04i05), ndbox, wk(i06) )
c
        call VSLOMAT( ncplx,idwf,doblas3,
     $   norb,nk, natm,ntyp, nshld,nald, nlat1c,
     $   itypa,numshl,lshel,nala,ala,cala, alamin, ratm,rlat,
     $   wk(icoskr),wk(isinkr),
     $   weight, hh, ndbox,nbox,nrec,inrec,inboxs,inboxs(nbox+1),boxrad,
     $   wk(iwf),wk(iwf), wkrbx(1),     wkrbx(1+mxrb), wk(iwfw),
c -->    wvfcns  wf       vslobox(1B)-s wkbox(10Br8)-s wkbi(lots)-s
     $   wk(i01i02), wk(i03s), wk(i06) )
c -->    vijk-io     vkij-s    vslo-i
c        nmatloc     2*nmatloc(r*4) nptr
c
c        ... switch to dn-spin electrons, if we have them
        i01i02 = i01 + nmatloc
c        i01i02 = i01 + nmat
        i04i05 = i05
      enddo
      i01bl = i01      
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( master )
      if( iproc .eq. master )then
        REWIND( unit=idhamfl )
      endif
      call H0WRITKP( idhamfl, mat,nk, wk(i01bl),
     &     wk(i01bl+nmatloc+nmatloc) )
      if( do_spin ) then
        call H0WRITKP( idhamfl, mat,nk, wk(i01bl+nmatloc),
     &       wk(i01bl+nmatloc+nmatloc) )
      endif
      if( iproc.eq.master )then
        REWIND( unit=idhamfl )
      endif
      call MTBUFF
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Invoke generalized eigensolvers for complex Hermitian matrices
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
ctp:
CTP        call EIGSOLV(
        call PEIGSOLV( masterlist,do_kppsolve, do_psolv,
     $   icluster, gapp,
     $   do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $   idhamfl,iham0fl,iovlpfl,
     $   ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $   icloseocc, dokdefct, doblas3, ikdefct,nstbulk,
     $   etemp,edegen, elecno, efermi,egap, esumsp,
     $   norb,nstate, nk, vecknsm,wtknsm,
     $   wk(ieigval),wk(ieigpop), npop, wksmo,
     $   wk(i01) )
c -->    wmat-6s+
c
      if( nk .eq. 1 )then
c       If single point, then not bands, probing single pt (gamma, probably)
        lprnt = 2
        call EVPRNT( IWR, lprnt, nspin, efermi,egap, npop,nocc0k,
     $   nstate,nk, wtk,wk(ieigval),wk(ieigpop) )
      endif
c
c     Put eigenvalues into eigenvals for convenience
c
      i1 = ieigval
      do  ispin=1,nspin
        do  ik=1,nk
          do  ist=1,nstate
            eigenvals(ist,ik,ispin) = wk(i1)
            i1 = i1 + 1
          enddo
        enddo
      enddo
c
c Stop here if no eigenvectors were computed above.
      if( .not. do_eigvecs ) RETURN
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  Calculations needing eigvenvectors, like (future) optical calculations
c  would follow from this point.  *** Not yet implemented ***
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  scfefermi, passed to this routine, contains the fermi level from the 
c  SCF calculation. You can use this to compute populations for the new 
c  k-vectors using the appropriate distribution function.
c  
 3555 format(/2x,'natm = ',i3,2x,'ntyp = ',i3,2x,'ityp = ',i3,2x)
 3556 format(3(3x,f9.5))
 3557 format(/2x,'mxwkwt=',2x,i6,2x,'n1r,n2r,n3r',3(3x,i4))
 3558 format(3(3x,i4))
 3559 format(5(3x,f9.5))
c
c    That's all, Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000 characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
