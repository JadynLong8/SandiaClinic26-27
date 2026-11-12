c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STARTX
c
c
      subroutine STARTX
c---------------------------------------------------------------
c Purpose: initialize system stuff the code needs to run
c          ***** MACHINE DEPENDENT *****
c
c Written: Peter A. Schultz,  6-December-2001, for v2.51
c
c Revision history:
c 30Apr12-PAS/2.63: More mptools to be wrapped
c 19Nov08-PAS/2.62: retool mptools for k-parallel communicators
c 15Jun05-PAS/2.58: tools to turn code into subroutine
c 28Oct02-PAS/2.54c: polish up path stuff; more Cplant-proofing
c  4Oct02-PAS/2.54: MPI stuff for NEB images enabled
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      COMMON  /UTLCONFIG/ isubcall
      INTEGER  isubcall
c
      LOGICAL    Lsubbed
      LOGICAL    Lcplant, Ljobname, build_path
      CHARACTER  chtmp*128, jobnm*128, prepath*10
      CHARACTER  pathtyp*6, chnode*6
      CHARACTER  blank
      DATA  blank / ' ' /
c
c >>>>> EXECUTABLE CODE:
c
c################################################################
c
c  *****  Machine-dependent flags  *****
c
csubc: Lsubbed = is code a subroutine? (Set units below)
      Lsubbed = .false.
      if( Lsubbed )then
        isubcall = 1
      else
        isubcall = 0
      endif
c
c Whether we invoke enfs file system on Cplant:
      Lcplant = .false.
      if( Lcplant )then
        nprepath = 5
        prepath(1:nprepath) = 'enfs:'
      else
        nprepath = 0
      endif
c
c Whether we build a path, or use the path we woke up on ...
      build_path = .false.
      if( Lcplant ) build_path = .true.
c  ... and type of scheme to get directory path:
C      pathtyp = 'BUILD '
      pathtyp = 'GETCWD'
c
c Whether we check for a job name from system input:
      Ljobname  = .false.
c
c################################################################
c
c  *****  Set up timer *****
c
c itimer =1 -> system timer; =2 -> wall timer
      itimer = 1
      if( Lcplant ) itimer = 2
      call TIMINIT( itimer )
c
c  *****  Set up I/O file unit number management  *****
c
c Set the standard input, output, error unit numbers ...
      IRD = 5
      IWR = 6
      if( Lsubbed )then
        IRD = 10
        IWR = 11
      endif
      IERRFL = IWR
c  ... and initialize (null) file stuff:
      call FLINIT( IRD, IWR, IERRFL )
c
c  *****  Set up max record length size in big i/o  *****
c
      maxreclen = 1 000 000
      call BIGIOINIT( maxreclen )
c
c  *****  Set up parallel/message-passing stuff  *****
c
      nprocs = 1
      node = 0
      node0 = node
      icomm = 0
c
cmpi:start Following statements enabled for MPI, else removed
cmpi      if( .not. Lsubbed )then
cmpic       Initialize MPI system
cmpi        call UTLMPINIT( node, nprocs, node0, icomm )
cmpi      else
cmpic       Just obtain node information from preexisting MPI system
cmpi        call UTLMPINFO( node, nprocs, node0, icomm )
cmpi      endif
cmpi:end
c
c Load global MP module:
      icommlvl = 1
      call MPLOADCOMM( nprocs,node,node0,icomm, icommlvl )
c
c  *****  Set up job name (root name for output files) *****
c
      jobnm = ' '
c Set default job name (for file prefixes):
      njobnm = 4
      jobnm(1:njobnm) = 'lcao'
c
      if( Ljobname )then
c       If command level argument exists, get a job name from shell
        nargs = IARGC()
        if( nargs.gt.0 )then
c         Use argument to set job-specific directory name
          call GETARG( 1, chtmp )
          lench = INDEX( chtmp, blank ) - 1
c         For later, save this argument for a job name:
          njobnm = lench
          if( lench.gt.0 ) jobnm = chtmp
        endif
      endif
c
c  *****  Build a path to connect to correct disk directory *****
c
c  ... first, prepend to connect to correct root file system:
c
      nslash = 0
      if( nprepath.gt.0 )then
        call FLADDPATH( nprepath, prepath(1:nprepath) )
c       We will need to build a path for sure now
        build_path = .true.
c       FLADDPATH adds a closing slash to prepath, so that we will
c       need to skip an opening slash in following path builds
        nslash = 1
      endif
c
c  ... then add the full path using specified prescription:
c
      if( build_path )then
c
        if( pathtyp .eq. 'GETCWD' )then
c          ... and extend with the cwd:
c
          call GETPATH( chtmp )
          lench = INDEX( chtmp, blank )
c
          ic1 = 1
          if( Lcplant .and. lench.gt.5 )then
c           Might need to strip opening garbage for Cplant!
            if( chtmp(1:5) .ne. '/enfs' )then
c             Eliminate leading stuff before '/enfs':
              do  ic=2,lench-4
                if( chtmp(ic:ic+4) .eq. '/enfs' )then
                  ic1 = ic
                  goto 100
                endif
              enddo
            endif
          endif
  100     continue
c
          if( lench.gt.0 )then
            if( chtmp(ic1:ic1) .eq. '/' )then
c             Might need to skip opening slash if prefix prepended
              icpath = ic1 + nslash
              lenadd = lench - icpath + 1
              call FLADDPATH( lenadd, chtmp(icpath:lench) )
            else
              lenadd = lench - ic1 + 1
              call FLADDPATH( lenadd, chtmp(ic1:lench) )
            endif
          endif
c
        elseif( pathtyp .eq. 'BUILD ' )then
c         Build path a la Art Edwards:

c        ... first, complete root name:
          lench = 8
          chtmp = 'enfs/tmp'
          if( lench.gt.0 ) call FLADDPATH( lench, chtmp(1:lench) )
c
c        ... second, extend with the user name:
          call GETENV( 'USER', chtmp )
          lench = INDEX( chtmp, blank ) - 1
          if( lench.gt.0 ) call FLADDPATH( lench, chtmp(1:lench) )
c
c        ... third, if jobnm exists, set job-specific directory:
          if( njobnm.gt.0 ) call FLADDPATH( njobnm, jobnm(1:njobnm) )
c
c        ... and, fourth, if mpi, extend with node identifier:
c            This code constructs node-specific prefixes that use up
c            to four characters to identify node.
          lench = 7
          chtmp = 'node000'
          if( nprocs.gt.9999 )then
c           At most 4 chars allowed in identifier
            call STOPXERR( 'START/too many nodes' )
          elseif( nprocs.ge.1000 )then
c           Use 4 chars in identification
            call STRNUM( nstr, chnode, node, 9999 )
            chtmp(4:7) = chnode(1:4)
          elseif( nprocs.gt.1 )then
c           Default is 3 chars in identification
            call STRNUM( nstr, chnode, node, 999 )
            chtmp(5:7) = chnode(1:3)
          endif
          if( lench.gt.0 ) call FLADDPATH( lench, chtmp(1:lench) )
c
        else
c         Type of directory build is not specified
          call STOPXERR( 'STARTX: bad directory build option' )
c
c         End of directory build options
        endif
c
c       End of build_path block
      endif
c
c  ***** Take care of output file names/opens  *****
c
c  Extend job name with node number for parallel (if not node0)
c  in case different nodes share a directory, to distinguish files
c     if( nprocs.gt 1 )then
      if( node.ne.node0 )then
        call STRNUM( nstr, chnode, node,nprocs )
        jobnm(njobnm+1:njobnm+nstr) = chnode(1:nstr)
        njobnm = njobnm + nstr
      endif
      call FLSETJOB( jobnm(1:njobnm) )
c
c  Open named output file for nodes for parallel (mpi) or named jobs
c     if( node.ne.node0 )then
      if( nprocs.gt.1 .or. Ljobname .or. Lsubbed .or. IWR.ne.6 )then
        call FLNAME( 'out', jobnm )
c       We have set IWR and reserved the unit above, so use this open:
        call FLOPENU( IWR, jobnm, 'UNK', 'FOR' )
        REWIND( unit = IWR )
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STOPX
c
c
      subroutine STOPX
c---------------------------------------------------------------
c Purpose: intercept program stop, end gracefully with message
c          ***** MACHINE DEPENDENT *****
c
c Written: Peter A. Schultz, 14-June-2001 for v2.48
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c Notes:
c  This routine is to intercept ALL exits from the code.
c  It substitutes for "STOP 'stoplbl'" in such a way as to
c  allow platform-dependent final message from program, and to
c  provide a single common exit point from the code to allow
c  for any necessary final cleaning up (e.g. MPI, files) before
c  the code goes to that great electronic Valhalla (rather than
c  Purgatory, if I may mix my metaphorical milieu's)
c
      COMMON  /UTLCONFIG/ isubcall
      INTEGER  isubcall
c Our output number for errors/log messages:
      DATA IWRX / 6 /
c local declarations:
      CHARACTER*(*)  stopmsg
      CHARACTER*1    blank
      DATA           blank / ' ' /
c
c >>>> EXECUTABLE CODE:
c
c Simple stop without any message:
      call FLGETIWR( IWR )
      istop = 0
      goto 990
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STOPXERR
c
      entry STOPXERR( stopmsg )
c
c Blank or not, always want to leave note that code failed with error:
      IWR = IWRX
      if( IWR.gt.0 ) write(IWR,*) '***** ERROR: ',stopmsg
c Set stop flag to errors:
      istop = -1
c
      goto 990
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> STOPXOK
c
      entry STOPXOK( stopmsg )
c
c With normal exit, if blank, do not output message to listing file
      call FLGETIWR( IWR )
      if( stopmsg .ne. blank )then
        istop = 1
        if( IWR.gt.0 ) write(IWR,*) 'Normal exit: ',stopmsg
      else
        istop = 0
      endif
c
c Let us clean up ...
c
  990 continue
c Summarize the timer:
      call TIMFINI( IWR )
c Purge output buffer:
      call MTBUFFI( IWR )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c (clean-up and platform-dependent wrapping up should be inserted here)
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c Wrap up files:
      call FLFINI
c Wrap up mpi stuff:
cmpi: Following statement enabled for mpi, else removed
cmpi      call UTLMPFINI
c
c  ... and exit the code:
      if( isubcall .eq. 1 .and. istop .ge. 0 ) RETURN
      STOP
c
c    That's all Folks!  (really.  the very end.  code stops here.)
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETPATH
c
c
      subroutine GETPATH( jobpath )
c---------------------------------------------------------------
c     This provides a safe version of getcwd,
c     that works with a variety of compilers.
c     ***** MACHINE DEPENDENT *****
c
c     Authored by Aidan P. Thompson 8/31/01
c     character*(*) function getpath()
c
c Revision history:
c  17Dec01-PAS/2.51: Cosmetically modified for use in Quest
c---------------------------------------------------------------
c
c       Note: len, len_trim, trim are standard on f90, but not
c       on older compilers.
c
c     This routine uses getcwd, which is not a fortran standard
c     It assumes the returned value of path is in one of two forms.
c
c     a) path followed by all blanks
c
c     b) path followed by a char(0) and then garbage, possibly blanks.
c
c     To avoid overflow, we need len(path) >= length of path + 1
c
      IMPLICIT NONE
c
c Output:
      CHARACTER*(*)  jobpath
c Local:
      INTEGER  lenpath, lencwd, inull
      INTEGER  istatus,getcwd,ind0
      CHARACTER  cwdname*128, blank*1, nullch*1
      DATA  blank / ' ' /
c
c >>>> EXECUTABLE CODE:
c
      nullch = CHAR( 0 )
c Initialize path to blanks:
      lenpath = LEN( jobpath )
      do  ind0=1,lenpath
        jobpath(ind0:ind0) = blank
      enddo
c
c Get (NB: getcwd is non-std), and parse the cwd ...
      istatus = GETCWD( cwdname )
c  ... find first null and (if it exists) switch to blank:
      inull = INDEX( cwdname, nullch )
      if( inull.gt.0 ) cwdname(inull:inull) = blank
c  ... and find length of string (i.e. before first blank):
      lencwd = INDEX( cwdname, blank ) - 1
c
      if( lencwd.gt.lenpath ) call STOPXERR( 'GETPATH: path too long' )
      do  ind0=1,lencwd
        jobpath(ind0:ind0) = cwdname(ind0:ind0)
      enddo
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETTLFT
c
c
      subroutine GETTLFT( timleft )
c---------------------------------------------------------------
c Purpose: MACHINE DEPENDENT routine to return time left (secs)
c---------------------------------------------------------------
c
      IMPLICIT REAL*8  (a-h,o-z)
c
c OTHER: return big const. as routine not implemented
      timleft = 999999.9
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETTUSE
c
c
      subroutine GETTUSE( timetotl, timeuser )
c--------------------------------------------------------------
c Purpose: MACHINE DEPENDENT routine to return seconds used
c Revision history:
c  17Jul03-PAS/2.56: return both total and sys/user time
c  14Jun03-PAS/2.55: change default time array dim to 2.
c--------------------------------------------------------------
c
      IMPLICIT REAL*8  (a-h,o-z)
c
c Note on timing routine(14Jun03-PAS):
c   ntd=3 had been the historical dimension
c   ntd=3 failed on some compilers, where ntd=2 worked
c   Looked up ETIME: man and web both confirm ntd=2:
c    tarray(1) = user time
c    tarray(2) = system time
c    ETIME     = total (user+sys) time
c      ... all in R*4 seconds
c
      PARAMETER  ( ntd = 2 )
      DIMENSION  tarray(ntd)
c
ccccccccccccccccccccccccccccccccccccccccccccc
c SUN-, SGI-, DEC-UNIX: ETIME returns time used in R*4 seconds
c   (19Dec90-PAS: thanks to S.H. Lamson of GE)
      REAL*4  etime,tarray
      ttotl = ETIME( tarray )
      tuser = tarray(1)
ccccccccccccccccccccccccccccccccccccccccccccc
c IBM-UNIX: ETIME_ returns time used in R*8 seconds
C      ttotl = ETIME_( tarray )
C      tuser = tarray(1)
c
ccccccccccccccccccccccccccccccccccccccccccccc
c OTHER: just set time elapsed to zero
C      ttotl = 0.
C      tuser = 0.
ccccccccccccccccccccccccccccccccccccccccccccc
c
      timetotl = ttotl
      timeuser = tuser
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> GETWALLT
c
c
      entry GETWALLT( timwall )
c---------------------------------------------------------------
c Purpose: MACHINE DEPENDENT routine to return wall time (secs)
c Installed: 28Oct04-PAS for v2.59, to get "true" elapsed wall
c            time diagnostic for instrumentation
c---------------------------------------------------------------
c
c The SECNDS wall time routine is presumably standard, BUT
c has not been broadly tested.  If it should fail, the total
c time given by ETIME should be substituted (as above).
c
C      twall = ETIME( tarray )
      twall = SECNDS( 0.0 )
      timwall = twall
c
      RETURN
c
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MTBUFF
c
c
      subroutine MTBUFF
c--------------------------------------------------------------
c Purpose: MACHINE DEPENDENT call to flush write/print buffers
c--------------------------------------------------------------
c
c UNIX: FLUSH flushes print buffer for unit called
      call FLGETIWR( IWR )
      if( IWR.gt.0 ) call FLUSH( IWR )
c OTHER: without routine, just skip
c
      RETURN
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MTBUFFI
c
c
      entry MTBUFFI( ifile )
c--------------------------------------------------------------
c Purpose: MACHINE DEPENDENT call to flush write/print buffers
c--------------------------------------------------------------
c
c UNIX: FLUSH flushes print buffer for unit called
      if( ifile.gt.0 ) call FLUSH( ifile )
c OTHER: without routine, just skip
c
      RETURN
      END
c 
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPWRAPS
c  
c
c=======================================================================
c Provides wrappers and simple tools for the mpi code
c
c  MPSENDR4     - send an r4 message
c  MPRECVR4     - recv an r4 message
c  MPSENDR8     - send an r8 message
c  MPRECVR8     - recv an r8 message
c  MPSENDI      - send an i4 message
c  MPRECVI      - recv an i4 message
c  MPSENDTICKET - send a ticket
c  MPRECVTICKET - recv a ticket
c  MPBCASTI     - bcast an integer vector
c  MPBCAST8     - bcast an r8 vector
c  MPBCAST4     - bcast an r4 vector
c  MPREDUCI     - reduce(sum) an integer vector
c  MPREDUC8     - reduce(sum) an r8 vector
c  MPREDUC4     - reduce(sum) an r4 vector
c  MPSPLITR8    - split vector from master over processors
c  MPMERGER8    - merge vector from processors to master
c  MPSPLITCOMM  - split global communicator into local communicators
c
c=======================================================================
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBARRIER
      subroutine MPBARRIER( icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Stubbed - do nothing
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDR4
      subroutine MPSENDR4( node_to, lenmsg, r4msg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL       r4msg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVR4
      subroutine MPRECVR4( node_from, lenmsg, r4msg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL       r4msg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDR8
      subroutine MPSENDR8( node_to, lenmsg, r8msg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DIMENSION  r8msg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVR8
      subroutine MPRECVR8( node_from, lenmsg, r8msg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DIMENSION  r8msg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDI
      subroutine MPSENDI( node_to, lenmsg, imsg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DIMENSION  imsg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVI
      subroutine MPRECVI( node_from, lenmsg, imsg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DIMENSION  imsg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDTICKET
      subroutine MPSENDTICKET( itarget, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Stubbed - do nothing
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVTICKET
      subroutine MPRECVTICKET( isource, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Stubbed - do nothing
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCASTI
      subroutine MPBCASTI( node_from, lenmsg, imsg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DIMENSION  imsg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCAST8
      subroutine MPBCAST8( node_from, lenmsg, r8msg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      DIMENSION  r8msg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCAST4
      subroutine MPBCAST4( node_from, lenmsg, r4msg, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
      REAL       r4msg(*)
c Stubbed - do nothing (return message unchanged)
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUCI
c
      subroutine MPREDUCI( node_to, nvec, ivecloc, ivecred, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Message space:
      DIMENSION  ivecloc(*), ivecred(*)
c
c Stubbed - copy input vector into output
      do  i=1,nvec
        ivecred(i) = ivecloc(i)
      enddo
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUC8
c
      subroutine MPREDUC8( node_to, nvec, vecloc, vecred, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Message space:
      DIMENSION  vecloc(*), vecred(*)
c
c Stubbed - copy input vector into output
      call DCOPY( nvec, vecloc,1, vecred,1 )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUC4
c
      subroutine MPREDUC4( node_to, nvec, vecloc, vecred, icomm )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Message space:
      REAL       vecloc(*), vecred(*)
c
c Stubbed - copy input vector into output
      call SCOPY( nvec, vecloc,1, vecred,1 )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITR8
c
      subroutine MPSPLITR8( nprocl, iprocl0, iprocl,
     $     nvec,nveclocl, vec, icomm)
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Input/Output:
      DIMENSION  vec(*)
c Stubbed - do nothing
      if( nprocl .le. 1 ) RETURN
      call STOPXERR( 'stubbed mpsplitr8' )
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPMERGER8
c
      subroutine MPMERGER8( nprocl, iprocl0, iprocl,
     $     nvec,ivec1,nveclocl, vec, icomm)
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Input/Output:
      DIMENSION  vec(*)
c Stubbed - do nothing
      if( nprocl.le.1 ) RETURN
      call STOPXERR( 'stubbed mpmerger8' )
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITGCOMM
c
      subroutine MPSPLITGCOMM( mygroup )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Stubbed - do nothing
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITGCOMM
c
      subroutine MPSPLITICOMM( mygroup )
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c Stubbed - do nothing
      RETURN
      END
c####################################################################
c eigstuff stubs for serial version links
c Dev version 2.62j - distributed solver memoery in parallel
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGPROCS
c
c
      subroutine EIGPROCS( nprocs, nprow,npcol, npsolver )
c
      IMPLICIT NONE
c
      INTEGER  nprocs,nprow,npcol
      INTEGER  npsolver
c
      if( nprocs .ne. 1 ) call STOPXERR( 'serial eigprocs: nproc.ne.1' )
      npsolver = 1
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGMDIST
c
c
      subroutine EIGMDIST( matdist, norb, npsolver )
c
      IMPLICIT NONE
c
      INTEGER  matdist
      INTEGER  norb,npsolver
c
      if( npsolver.ne.1 )call STOPXERR( 'serial eigmdist: npsolv.ne.1' )
      matdist = norb*norb
c
      RETURN
      END
