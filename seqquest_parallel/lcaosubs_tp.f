c 
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPTOOLS
c  
c
c=======================================================================
c Provides wrappers and simple tools for the mpi code
c
c  MPBARRIER    - an MPI barrier
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
c  MPSPLITGCOMM  - split global communicator into local communicators
c  MPSPLITICOMM  - split global communicator into local communicators
c  MPFREECOMMK   - free MPCOMM_K communicator
c
c=======================================================================
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSYNCH
c
c
      subroutine MPBARRIER( icomm )
c---------------------------------------------------------------
c Purpose: synchronize processors
c---------------------------------------------------------------
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
      INTEGER  icomm, ierror
cmpi
      call MPI_Barrier( icomm, ierror  )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDR4
c
c
      subroutine MPSENDR4( node_to, lenmsg, r4msg, icomm )
c---------------------------------------------------------------
c Purpose: wrappers for mpi send/recv of R*4 data.
c
c Written: Peter A. Schultz, 11-February-2007, for 2.60
c          Adapted from r8 version
c---------------------------------------------------------------
c
c  node_to   - target node for outgoing message
c  node_from - originating node of incoming message
c  lenmsg    - length of message, in R*4 units
c  r4msg()   - the R*4 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      REAL       r4msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Send with a blank tag:
      msg_tag = 0
cmpi
      call MPI_SEND( r4msg, lenmsg,
     $     MPI_REAL            , node_to  , msg_tag,
     $     icomm, ierror )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVR4
c
      subroutine MPRECVR4( node_from, lenmsg, r4msg, icomm )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      REAL       r4msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c >>>> EXECUTABLE CODE:
c
c Recv with a blank tag:
      msg_tag = 0
cmpi
      call MPI_RECV( r4msg, lenmsg,
     $     MPI_REAL            , node_from, msg_tag,
     $     icomm, istatus, ierror )
c
c Find out how long the incoming message was:
cmpi
      call MPI_GET_COUNT( istatus, MPI_DOUBLE_PRECISION, lenmsg,
     $ ierror )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDR8
c
c
      subroutine MPSENDR8( node_to, lenmsg, r8msg, icomm )
c---------------------------------------------------------------
c Purpose: wrappers for mpi send/recv of R*8 data.
c
c Written: Peter A. Schultz, 21-January-2002, for v2.54
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c  14Oct04-PAS/2.59-tp: add return msg length in recv
c---------------------------------------------------------------
c
c  node_to   - target node for outgoing message
c  node_from - originating node of incoming message
c  lenmsg    - length of message, in R*8 units
c  r8msg()   - the R*8 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r8msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Send with a blank tag:
      msg_tag = 0
cmpi
      call MPI_SEND( r8msg, lenmsg,
     $     MPI_DOUBLE_PRECISION, node_to  , msg_tag,
     $     icomm, ierror )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVR8
c
      subroutine MPRECVR8( node_from, lenmsg, r8msg, icomm )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r8msg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Recv with a blank tag:
      msg_tag = 0
cmpi
      call MPI_RECV( r8msg, lenmsg,
     $     MPI_DOUBLE_PRECISION, node_from, msg_tag,
     $     icomm, istatus, ierror )
c
c Find out how long the incoming message was:
cmpi
      call MPI_GET_COUNT( istatus, MPI_DOUBLE_PRECISION, lenmsg,
     $ ierror )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDI
c
c
      subroutine MPSENDI( node_to, lenmsg, imsg, icomm )
c---------------------------------------------------------------
c Purpose: wrappers for mpi send/recv of integer data.
c---------------------------------------------------------------
c
c  node_to   - target node for outgoing message
c  node_from - originating node of incoming message
c  lenmsg    - length of message, in integer units
c  imsg()    - the integer data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  imsg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Send with a blank tag:
      msg_tag = 0
cmpi
      call MPI_SEND( imsg, lenmsg,
     $     MPI_INTEGER, node_to  , msg_tag,
     $     icomm, ierror )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVI
c
      subroutine MPRECVI( node_from, lenmsg, imsg, icomm )
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  imsg(*)
cmpi Needed for mpi:
      INTEGER    istatus(MPI_STATUS_SIZE)
c
c Recv with a blank tag:
      msg_tag = 0
cmpi
      call MPI_RECV( imsg, lenmsg,
     $     MPI_INTEGER, node_from, msg_tag,
     $     icomm, istatus, ierror )
c
c Find out how long the incoming message was:
cmpi
      call MPI_GET_COUNT( istatus, MPI_INTEGER, lenmsg,
     $ ierror )
c
      RETURN
      END
c 
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSENDTICKET
c  
c
      subroutine MPSENDTICKET( itarget, icomm )
c---------------------------------------------------------------
c Purpose: send/recv a ticket to a processor
c Written: Peter A. Schultz, 6-July-2006, for v2.60
c---------------------------------------------------------------
c     
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (lenrcpt=1)
c
      lenmsg = lenrcpt
      rcpt = 1.d0
c
      call MPSENDR8( itarget, lenmsg , rcpt, icomm )
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPRECVTICKET
c  
      subroutine MPRECVTICKET( isource, icomm )
c     
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  (lenrcpt=1)
c
      lenmsg = lenrcpt
      call MPRECVR8( isource, lenmsg , rcpt, icomm )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCASTI
c
c
      subroutine MPBCASTI( node_from, lenmsg, imsg, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for integer broadcast
c
c Written: Andrew C Pineda, 12-November-2009, for v2.63
c
c Revision history:
c  12Nov09-ACP/2.63: adapted from MPBCASTR8
c---------------------------------------------------------------
c
c  node_from - originating node of broadcast
c  lenmsg    - length of message, in integer units
c  imsg()    - the integer data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  imsg(*)
c
cmpi
      call MPI_Bcast( imsg, lenmsg, MPI_INTEGER,
     $ node_from, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPBCASTI failure' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCAST8
c
c
      subroutine MPBCAST8( node_from, lenmsg, r8msg, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r8 broadcast
c
c Written: Peter A. Schultz, 21-January-2002, for v2.59-tp
c
c Revision history:
c   6Jul06-APT/    : merged task parallel and image parallel
c---------------------------------------------------------------
c
c  node_from - originating node of broadcast
c  lenmsg    - length of message, in R*8 units
c  r8msg()   - the R*8 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r8msg(*)
c
cmpi
      call MPI_Bcast( r8msg, lenmsg, MPI_DOUBLE_PRECISION,
     $ node_from, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPBCAST8 failure' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPBCAST4
c
c
      subroutine MPBCAST4( node_from, lenmsg, r4msg, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r4 broadcast
c
c Written: Peter A. Schultz, 11-January-2007, for v2.60
c
c Revision history:
c   none
c---------------------------------------------------------------
c
c  node_from - originating node of broadcast
c  lenmsg    - length of message, in R*8 units
c  r4msg()   - the r*4 data being sent
c  icomm     - the communicator
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      REAL     r4msg
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Message space:
      DIMENSION  r4msg(*)
c
cmpi
      call MPI_Bcast( r4msg, lenmsg, MPI_REAL,
     $ node_from, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPBCAST4 failure' )
c
      RETURN
      END
c     
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUCI
c
c
      subroutine MPREDUCI( node_to, nvec, ivecloc, ivecred, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for integer reduction/sum
c
c Written: Andrew C. Pineda,  23-November-2009, for v2.63
c
c Revision history:
c   none
c---------------------------------------------------------------
c       Use standard mpi reduction as the default scheme
c
      include  'mpif.h'
      itckt = 1
      lenmsg = 1
cmpi
      call MPI_Bcast( itckt, lenmsg, MPI_INTEGER,
     &  node_to, icomm, ierr )
      call MPI_Reduce( ivecloc, ivecred, nvec, MPI_INTEGER,
     &  MPI_SUM, node_to, icomm, ierr )
      if( ierr.ne.0 ) call STOPXERR( 'MPREDUCI failure' )
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUC8
c
c
      subroutine MPREDUC8( node_to, nvec, vecloc, vecred, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r8 reduction/sum
c
c Written: Peter A. Schultz, 30-January-2007, for v2.60
c
c Revision history:
c   none
c---------------------------------------------------------------
c
c  node_to   - target of reduction
c  nvec      - length of vector, in R*8 units
c  vecloc()  - the local R*8 data being sent
c  vecred()  - the reduced vector
c  icomm     - the communicator
c
c Choices for reduction (iopt):
c   default = native mpi reduce/sum
c         1 = reduce direct to target with tickets
c         2 = simple binary tree reduction
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
      PARAMETER  ( iopt = 2 )
      PARAMETER  ( r8one = 1.d0 )
c Message space:
      DIMENSION  vecloc(*), vecred(*)
c
c >>>> EXECUTABLE CODE:
c
C      call MPNODES( nprocs )
C      call MPNODE( iproc )
c Use native rather than wrapper (do not know if k-parallel!)
cmpi
      call MPI_comm_rank( icomm, iproc, ierror )
      call MPI_comm_size( icomm, nprocs, ierror )
c
      if( nprocs .eq. 1 )then
c       Nothing to reduce, reduction is simply the local vector
        call DCOPY( nvec, vecloc,1, vecred,1 )
        RETURN
      endif
c
      if( iopt .eq. 1 )then
c       Assumimg reduction of long vec over modest procs, try to reduce
c       mpi contention by reduction directly to target using ticketing:
c
        if( iproc .eq. node_to )then
c         Target proc sends tickets, accumulates reduced vector
c
c         Start off with target's data in reduced vector:
          call DCOPY( nvec, vecloc,1, vecred,1 )
          do 100 ip=1,nprocs
            jproc = ip - 1
            if( node_to .ne. jproc )then
c             Send ticket to proc, get its data, accumulate vector
              call MPSENDTICKET( jproc, icomm )
              nvecrecv = nvec
              call MPRECVR8( jproc, nvecrecv, vecloc, icomm )
              call DAXPY( nvec, r8one, vecloc,1, vecred,1 )
            endif
  100     continue
c
        else
c         Wait for ticket, then send vector to target:
          call MPRECVTICKET( node_to, icomm )
          call MPSENDR8( node_to, nvec, vecloc, icomm )
        endif
c
      elseif( iopt .eq. 2 )then
c       Try a simple by-hand binary tree reduce, to node 0
c
c       Start off copy of local vector in reduced vector:
        call DCOPY( nvec, vecloc,1, vecred,1 )
c
        iprdel = 1
  200   continue
        if( iprdel .ge. nprocs ) goto 260
        nprset = 2*iprdel
c       Cycle over sending process:
        do 250 ipr=iprdel,nprocs-1,nprset
          iprecv = ipr - iprdel
          ipsend = ipr
          if( iproc .eq. iprecv )then
c           Send ticket to source, get its data, accumulate vector
            call MPSENDTICKET( ipsend, icomm )
            nvecrecv = nvec
            call MPRECVR8( ipsend, nvecrecv, vecloc, icomm )
            call DAXPY( nvec, r8one, vecloc,1, vecred,1 )
          elseif( iproc .eq. ipsend )then
c           Wait for ticket from target, then send vector to target:
            call MPRECVTICKET( iprecv, icomm )
            call MPSENDR8( iprecv, nvec, vecred, icomm )
c           This proc is done, get out ...
            goto 260
          endif
  250   continue
        iprdel = 2*iprdel
        goto 200
  260   continue
c
        if( node_to .ne. 0 )then
c         Send final result to right node and right location
          if( iproc .eq. 0 )then
            call MPSENDR8( node_to, nvec, vecred, icomm )
          elseif( iproc .eq. node_to )then
            nvecrecv = nvec
            call MPRECVR8( 0     , nvecrecv, vecred, icomm )
          endif
        endif
c
      else
c       Use standard mpi reduction as the default scheme
c
        itckt = 1
        lenmsg = 1
cmpi
        call MPI_Bcast( itckt, lenmsg, MPI_INTEGER,
     $   node_to, icomm, ierr )
        call MPI_Reduce( vecloc, vecred, nvec, MPI_DOUBLE_PRECISION,
     $   MPI_SUM, node_to, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Reduce vloc3' )
c
      endif
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPREDUC4
c
c
      subroutine MPREDUC4( node_to, nvec, vecloc, vecred, icomm )
c---------------------------------------------------------------
c Purpose: wrapper for r4 reduction/sum
c
c Written: Peter A. Schultz, 11-February-2007, for v2.60
c
c Revision history:
c   none
c---------------------------------------------------------------
c
c  node_to   - target of reduction
c  nvec      - length of vector, in R*8 units
c  vecloc()  - the local R*4 data being sent
c  vecred()  - the reduced vector
c  icomm     - the communicator
c
c Choices for reduction (iopt):
c   default = native mpi reduce/sum
c         1 = reduce direct to target with tickets
c         2 = simple binary tree reduction
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
      PARAMETER  ( iopt = 2 )
      REAL       r4one
      PARAMETER  ( r4one = 1.e0 )
c Message space:
      REAL       vecloc, vecred
      DIMENSION  vecloc(*), vecred(*)
c
c >>>> EXECUTABLE CODE:
c
C      call MPNODES( nprocs )
C      call MPNODE( iproc )
c Use native rather than wrapper (do not know if k-parallel!)
cmpi
      call MPI_comm_rank( icomm, iproc, ierror )
      call MPI_comm_size( icomm, nprocs, ierror )
c
      if( nprocs .eq. 1 )then
c       Nothing to reduce, reduction is simply the local vector
        call SCOPY( nvec, vecloc,1, vecred,1 )
        RETURN
      endif
c
      if( iopt .eq. 1 )then
c       Assumimg reduction of long vec over modest procs, try to reduce
c       mpi contention by reduction directly to target using ticketing:
c
        if( iproc .eq. node_to )then
c         Target proc sends tickets, accumulates reduced vector
c
c         Start off with target's data in reduced vector:
          call SCOPY( nvec, vecloc,1, vecred,1 )
          do 100 ip=1,nprocs
            jproc = ip - 1
            if( node_to .ne. jproc )then
c             Send ticket to proc, get its data, accumulate vector
              call MPSENDTICKET( jproc, icomm )
              nvecrecv = nvec
              call MPRECVR4( jproc, nvecrecv, vecloc, icomm )
              call SAXPY( nvec, r4one, vecloc,1, vecred,1 )
            endif
  100     continue
c
        else
c         Wait for ticket, then send vector to target:
          call MPRECVTICKET( node_to, icomm )
          call MPSENDR4( node_to, nvec, vecloc, icomm )
        endif
c
      elseif( iopt .eq. 2 )then
c       Try a simple by-hand binary tree reduce, to node 0
c
c       Start off copy of local vector in reduced vector:
        call SCOPY( nvec, vecloc,1, vecred,1 )
c
        iprdel = 1
  200   continue
        if( iprdel .ge. nprocs ) goto 260
        nprset = 2*iprdel
c       Cycle over sending process:
        do 250 ipr=iprdel,nprocs-1,nprset
          iprecv = ipr - iprdel
          ipsend = ipr
          if( iproc .eq. iprecv )then
c           Send ticket to source, get its data, accumulate vector
            call MPSENDTICKET( ipsend, icomm )
            nvecrecv = nvec
            call MPRECVR4( ipsend, nvecrecv, vecloc, icomm )
            call SAXPY( nvec, r4one, vecloc,1, vecred,1 )
          elseif( iproc .eq. ipsend )then
c           Wait for ticket from target, then send vector to target:
            call MPRECVTICKET( iprecv, icomm )
            call MPSENDR4( iprecv, nvec, vecred, icomm )
c           This proc is done, get out ...
            goto 260
          endif
  250   continue
        iprdel = 2*iprdel
        goto 200
  260   continue
c
        if( node_to .ne. 0 )then
c         Send final result to right node and right location
          if( iproc .eq. 0 )then
            call MPSENDR4( node_to, nvec, vecred, icomm )
          elseif( iproc .eq. node_to )then
            nvecrecv = nvec
            call MPRECVR4( 0     , nvecrecv, vecred, icomm )
          endif
        endif
c
      else
c       Use standard mpi reduction as the default scheme
c
        itckt = 1
        lenmsg = 1
cmpi
        call MPI_Bcast( itckt, lenmsg, MPI_INTEGER,
     $   node_to, icomm, ierr )
        call MPI_Reduce( vecloc, vecred, nvec, MPI_REAL,
     $   MPI_SUM, node_to, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'pvlocmat - Reduce vloc3' )
c
      endif
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITR8
c
c
      subroutine MPSPLITR8( nprocl, iprocl0, iprocl,
     $     nvec,nveclocl, vec, icomm)
c---------------------------------------------------------------
c Purpose: split over a set of processors a long vector existing
c          on the master processor
c Written: Peter A. Schultz, for task-parallel on 2.59
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
c Notes: the scheme used to split the vector to the master node
c   is the simplest possible.  The master sends pieces to each
c   node in sequence.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output:
      DIMENSION  vec(*)
c
c >>>> EXECUTABLE CODE:
c
      nveclocl = nvec
      if( nprocl .le. 1 ) RETURN
c
c Set up the split:
c  Assume nodes are contiguously sequenced, starting with
c  the source (presumably master) node (iprocl0):
      iprocmx = iprocl0 + nprocl - 1
c  Length of vector to be on first processor:
      nveclocl = nvec / nprocl
      if( nveclocl*nprocl .ne. nvec ) nveclocl = nveclocl + 1
c  Index start of local vector:
      iv1 = 1
c
c Perform the split:
c
      iv2 = iv1 + nveclocl
      nprocleft = nprocl
      nvecleft = nvecleft - nveclocl
      do 100 inode=iprocl0+1,iprocmx
        itarget = inode
        nprocleft = nprocleft - 1
        nvecleft = nvec - iv2 + 1
c       Length of vector segment to be sent:
        nvtarget = nvecleft / nprocleft
        if( nvtarget*nprocleft .ne. nvecleft ) nvtarget = nvtarget + 1
c
        if( iprocl .eq. iprocl0 )then
c         Master sends to the target proc its part of the vector
          call MPSENDR8( itarget, nvtarget, vec(iv2), icomm )
c
        elseif( iprocl .eq. itarget )then
c         Proc receives its part of vector from master
          isource = iprocl0
          nvecrecv = nvtarget
          call MPRECVR8( isource, nvecrecv, vec(iv1), icomm )
          if( nvecrecv .ne. nvtarget ) call STOPXERR( 'MPSPLITR8 err' )
          nveclocl = nvecrecv
          goto 101
        endif
        iv2 = iv2 + nvtarget
c
  100 continue
  101 continue
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPMERGER8
c
c
      subroutine MPMERGER8( nprocl, iprocl0, iprocl,
     $     nvec,ivec1,nveclocl, vec, icomm)
c---------------------------------------------------------------
c Purpose: merge onto master processor a full vector which
c          is broken up in blocks over a set of processors.
c
c Written: Peter A. Schultz, 27-October-2004, for 2.59 task-parallel
c
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
c Notes: the scheme used to merge the vector to the master node
c   is the simplest possible.  The master sends a "ticket" to
c   each node in sequence, which signals that node to send its
c   part of the vector to the master node.  The master appends
c   the incoming vector onto its existing vector, and proceeds
c   to ticket the next processor until all are done.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input/Output:
      DIMENSION  vec(*)
c
c >>>> EXECUTABLE CODE:
c
      if( nprocl.le.1 ) RETURN
      lenrcpt = 1
      rcpt = 1.d0
c
c Second, set up the merge:
c  Assume nodes are contiguously sequenced, starting with
c  the target (presumably master) node (iprocl0):
      iprocmx = iprocl0 + nprocl - 1
c  Length of vector (currently) on the local processor:
      nveclen = nveclocl
c  Index of start of existing local vector:
      iv1 = ivec1
c
c Finally, perform the merge to the root node:
c
      if( iprocl .eq. iprocl0 )then
c       Master node, collect from all nodes ...
c
        iv2 = iv1 + nveclen
        do 100 inode=iprocl0+1,iprocmx
          isource = inode
c
c         First, tell the source proc it is ok to send
          call MPSENDR8( isource, lenrcpt, rcpt, icomm )
c
c         Second, receive that proc's part of the vector
          nvecrecv = nvec - iv2 + 1
          call MPRECVR8( isource, nvecrecv, vec(iv2), icomm )
          iv2 = iv2 + nvecrecv
  100   continue
c
      else
c       Non-master, merging full vector to master ...
c
c       First, get ticket from master:
        call MPRECVR8( iprocl0, lenrcpt, rcpt, icomm )
c
c       Send local segment of vector to master:
        call MPSENDR8( iprocl0, nveclen, vec(iv1), icomm )
c
      endif
c
c    That's all, Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITGCOMM
c
c
      subroutine MPSPLITGCOMM( mygroup )
c---------------------------------------------------------------
c Purpose: generate local communicators from global communicator
c
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c >>>> EXECUTABLE CODE:
c
      call MPCOMM_G( icomm_g )
cseq:start
cseq      icomm_l = icomm_g
cseq      nprocs = 1
cseq      node = 0
cseq      node0 = 0
cseq:end
c This sets the master to the 0-node in each group:
      mykey = 0
cmpi:start
      call mpi_comm_split( icomm_g, mygroup, mykey, icomm_l, ierror )
      call mpi_comm_rank( icomm_l, node, ierror )
      call mpi_comm_size( icomm_l, nprocs, ierror )
cmpi:end
c
c Load the comm wrapper (local/image level)
      icommlvl = 2
      node0 = 0
      call MPLOADCOMM( nprocs, node, node0, icomm_l, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPSPLITICOMM
c
c
      subroutine MPSPLITICOMM( mygroup )
c---------------------------------------------------------------
c Purpose: generate local communicators from global communicator
c
c Revision history:
c   6Jul06-APT/2.60: merged task parallel and image parallel
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
cmpi Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c >>>> EXECUTABLE CODE:
c
      call MPCOMM( icomm_l )
cseq:start
cseq      icomm_k = icomm_l
cseq      nprocs = 1
cseq      node = 0
cseq      node0 = 0
cseq:end
      mykey = 0
cmpi:start
      call mpi_comm_split( icomm_l, mygroup, mykey, icomm_k, ierror )
      call mpi_comm_rank( icomm_k, node, ierror )
      call mpi_comm_size( icomm_k, nprocs, ierror )
cmpi:end
c
c Load the comm wrapper (sub-image, k-parallel level):
      icommlvl = 3
      node0 = 0
      call MPLOADCOMM( nprocs, node, node0, icomm_k, icommlvl )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MPFREECOMMK
c
c
      subroutine MPFREECOMMK()
c
      INTEGER istatus
      INCLUDE 'mpif.h'
c
      INTEGER  kparopt, icomm_k
c
      call KPFLAG_GET( kparopt )
      if( kparopt .eq. 2 )then
        call MPCOMM_K( icomm_k )
        call MPI_COMM_FREE( icomm_k, istatus )
      endif
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
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGSTUFF
c
c
c===============================================================
c Purpose: a set of routines to help handle the parallel solver stuff.
c Written; Peter A. Schultz, 9-December-2008, for 2.62 (k-parallel)
c Revision history:
c  none
c===============================================================
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGPROCS
c
c
      subroutine EIGPROCS( nprocs, nprow,npcol, npsolver )
c---------------------------------------------------------------
c Purpose: return nearly square processor grid for solver
c---------------------------------------------------------------
c
c Generate the best (most square) processor grid within nprocs
c
c nprocs = total number of processors available in this communicator
c npcol = number of columns in solver processor grid
c nprow = number of rows in solver processor grid
c npsolver = npcol*nprow = number of processors in solver grid
c
      IMPLICIT NONE
c
c Input:
      INTEGER  nprocs
c Output:
      INTEGER  npcol,nprow, npsolver
c
c Local:
      REAL     rnprocs
c
c >>>> EXECUTABLE CODE:
c
      rnprocs = DBLE( nprocs ) + 0.1
      npcol = SQRT( rnprocs )
      nprow = nprocs/npcol
      npsolver = nprow*npcol
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EIGMDIST
c
c
      subroutine EIGMDIST( matdist, norb, npsolver )
c---------------------------------------------------------------
c Purpose: return default distributed matrix size
c Revision history:
c  13Jan13-PAS/2.63a: fix a dumb mistake, very narrow failure range
c---------------------------------------------------------------
c
c norbdist = rank of matrices for which we will use a reduced
c            distributed memory allocation.
c
c norbdist is set to a somewhat larger value to ensure that
c we do not run into block-size problems, but small enough
c that using a full matrix assignment is sure to give no problems.
c
      IMPLICIT NONE
c
      INTEGER  NB
      DATA     NB / 32 /
c
c Output:
      INTEGER  matdist
c Input:
      INTEGER  norb, npsolver
c
c Local:
      INTEGER  mat, matpercpu
      INTEGER  norbdist
      DATA     norbdist / 800 /
c
c >>>> EXECUTABLE CODE
c
      mat = norb*norb
      if( norb .ge. norbdist .and. npsolver .gt. 1 )then
        matpercpu = (mat+npsolver-1) / npsolver
c       Add 32*norb to give slack for block sizes coarseness
c        and uneven distribution across processors
        matdist = matpercpu + NB*(norb+NB)
      else
c       Just give the full matrix space
        matdist = mat
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPPMEMX
c
c
      subroutine KPPMEMX( per_process, do_eigvecs,
     $ norb,nstate,  memeig )
c---------------------------------------------------------------
c Purpose: compute memory requirements (size of wmat argument) of 
c          PSCHROEDKP for wkmem using Scalapack calls.
c
c Written:  Andrew C. Pineda,  9-February-2010 for v2.62
c
c Revision history:
c  4Jan13-PAS/2.64: adapted from pschroedkpsz by acp
c---------------------------------------------------------------
c
c Return argument:
c    memeig = process value (per_process = .true.)
c           = max value over all processes (per_process = .false.)
c
      INTEGER     NB
      PARAMETER  (NB=32)
      INCLUDE  'mpif.h'
c
      INTRINSIC  MAX
c
c Input arguments
      INTEGER  norb
      LOGICAL  per_process
      LOGICAL  do_eigvecs
c
c Output arguments
      INTEGER  memeig
c
c Local variables
      integer  memeiglocal
      logical  per_process_local
c
      integer  mat, matdist
      integer  nwpreph, nwpreps, nwprep, nwsolv
      integer  nwclose, nwpsch, nwleft
      integer  iwHmat, iwOvlp, iwVecs
      integer  iwZwork, iwRwork, iwIwork, iwiwk1, iwlast, iwEnd
      integer  lwzwork,lzwork, lrwork, liwork
      integer  iwMat1, iwMat2, iwVeco
c
      integer  CONTEXT_K
      integer  nprocs_k, icomm
      integer  nprow_k, npcol_k, npsolver_k
c
      integer  ierr
c
      LOGICAL  USE_SCALAPACK_FORMULAS
      DATA    USE_SCALAPACK_FORMULAS / .true. /
C      LOGICAL  OPTIMAL_EIGENVECTORS
C      DATA    OPTIMAL_EIGENVECTORS   / .true. /
c
      LOGICAL  DEBUG_MEM
      DATA     DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DEV/KPPMEMX: do_eigvecs=',do_eigvecs
      endif
      mat = norb*norb
c  Get Local MP info
c     Image level info
      call MPCOMM( icomm )
c     K-point level info
      call MPNODES_K( nprocs_k )
c
      memeiglocal = 0
c
      if( nprocs_k .eq. 1 )then
c       Using Lapack serial eigensolver
c 
c       Only one processor assigned to the MPI context this processor
c       belongs to, set up memory to use LAPACK. 
c
        matdist = mat
        lwzwork = MAX( 4*norb, mat )
        lrwork  = 7*norb
        liwork  = 5*norb
c
c       Close serial solver memory needs
      else
c       USING SCALAPACK
c
c       Get the process grid dimensions, given the context size
c
        call EIGPROCS( nprocs_k, nprow_k, npcol_k, npsolver_k )
c
c       Determine pzhegvx workspace requirements
c
        if( USE_SCALAPACK_FORMULAS )then
          if( DEBUG_MEM )then
c           Use default mpi/blacs-avoiding naive work space:
            call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $       norb, matdist, lwzwork, lzwork, lrwork, liwork )
            lwork = lzwork
          endif
c         Get EXACT mpi/blacs-aware work space using this context:
C          CONTEXT_K = -1
          call MPCOMM_K( icomm_k )
          CONTEXT_K = icomm_k
          call BLACS_GRIDINIT( CONTEXT_K, 'R', nprow_k, npcol_k )
          call PSCHWKSZX( CONTEXT_K, npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
          call BLACS_GRIDEXIT( CONTEXT_K )
        else
c         Use default mpi/blacs-avoiding naive work space:
          call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
        endif
      endif
c
c     Initialize memory allocations within icomm_k.
c     MEM-Solver: d-C16-H, d-C16-S, d-C16-Vec, d-C16-Zwork, d-R8-Rwork)
c
cxxx264: rednudant (unused here) temp storage for solvers:
      iwHmat  = 1
      iwOvlp  = iwHmat + 2*matdist
      iwVecs  = iwOvlp + 2*matdist
      iwZwork = iwVecs + 2*matdist
      iwRwork = iwZwork + lwzwork
      iwIwork = iwRwork + lrwork
      iwiwk1  = iwIwork + liwork
      iwlast  = iwiwk1  + norb
      iwEnd   = iwlast  - 1
c
      iwMat1  = iwOvlp
      iwMat2  = MAX( iwMat1 + mat, iwVecs )
      iwVeco  = iwZwork
c
c     The solver memory allocations have the following requirements:
      nwpreph = 2*matdist + 3*mat
      nwpreps = 4*matdist + 2*mat
      nwprep = MAX( nwpreph, nwpreps )
      nwsolv = 6*matdist + lwzwork + lrwork + liwork + norb
c     nwclose = 6*matdist + 2*norb*nstate
c     because all vectors are brought back into full norb*norb ...
cxx264: nwclose (eigeenec store) needs to be dealt with (norb*nstate)
      nwclose = 6*matdist + 2*norb*norb
      nwpsch = MAX( nwprep, nwsolv, nwclose )
c
      nwleft = nwpsch - nwsolv
      if( DEBUG_MEM )then 
        write(IWR,*) 'DEV/KPPMEMX: spare solver mem=',nwleft
      endif
c
      memeiglocal = MAX( iwEnd, iwMat2+2*mat, iwVeco+2*mat )
c
c >>>>> Put it all together and report to final result
c
      per_process_local = per_process
      if( per_process_local ) then
c       Report the local memeig value
        memeig = memeiglocal
      else
c       Report the maximum value of memeig for any process in icomm
c       Need barrier to prevent a crash
        call MPBARRIER( icomm )
        icount = 1
        call MPI_ALLREDUCE( memeiglocal, memeig, icount,
     $       MPI_INTEGER, MPI_MAX, icomm,
     $       ierr )
      endif
      if( DEBUG_MEM )then
        write(IWR,'(a,4i12)') 'DEV/KPPMEMX: works,Z,R,I=',
     $                         lwzwork, lrwork, liwork
        write(IWR,'(a,3i12)') 'DEV/KPPMEMX: nwsch,memeig=',
     $                         nwpsch, memeiglocal, memeig
        write(IWR,*) 'DEV/KPPMEMX: icomm,nprocs_k=',icomm,nprocs_k
        call FLUSH( IWR )
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSOLWKSZX
c
c
      subroutine PSCHWKSZX( icontxt, nprocs, do_eigvecs,
     $ norb, matdist, lwzwork, lzwork, lrwork, liwork )
c---------------------------------------------------------------
c Purpose: compute exact memory requirements for PSCHROEDKP solvers,
c          using Scalapack (p)zhegvx calls.
c      
c Written:  Peter A. Schultz, 10-January-2013, for 2.64
c           Derived from PSCHROEDKPSZ by Andrew C. Pineda (9Feb10)
c
c Revision history:
c  none
c---------------------------------------------------------------
c
c  Work array space for
c    (1)  serial solver zhegvx
c    (2)  parallel solver pzhegvx
c
      INTEGER     NB
      PARAMETER  (NB=32)
C      INCLUDE  'mpif.h'
c
      INTRINSIC  MAX
c
c Input arguments
      INTEGER  icontxt
      INTEGER  nprocs
      LOGICAL  do_eigvecs
c
      INTEGER  norb
c
c Output arguments
      INTEGER  matdist
      INTEGER  lwzwork, lzwork, lrwork, liwork
c
c Local variables
      INTEGER  mat
      INTEGER  nprow, npcol, npsolver
c
      INTEGER  CONTEXT
      INTEGER  myrow, mycol
      INTEGER  NN, NP0, NQ0, NEIG, MQ0
      INTEGER  ANB, SQNPC, NPS, NHETRD_LWOPT, NHEGST_LWOPT
c
      LOGICAL  USE_SCALAPACK_FORMULAS
      LOGICAL  OPTIMAL_EIGENVECTORS
      DATA    USE_SCALAPACK_FORMULAS / .true. /
      DATA    OPTIMAL_EIGENVECTORS   / .true. /
c
      LOGICAL  DEBUG_MEM
      DATA     DEBUG_MEM  / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DEV/WKSZ-X: norb,nprocs=',norb, nprocs
        write(IWR,*) 'DEV/WKSZ-X: do_eigvecs=',do_eigvecs
      endif
c
      mat = norb*norb
c
      if( nprocs .eq. 1 )then
c       Using Lapack serial eigensolver
c 
c       Only one processor assigned to the MPI context this processor
c       belongs to, set up memory to use LAPACK. 
c
        matdist = mat
cxxx264        lwzwork = 2*MAX( 2*norb, mat ) ! I am sure we don't need 2mat -pas
        lwzwork =   MAX( 4*norb, mat )
        lzwork  = lwzwork / 2
        lrwork  = 7*norb
        liwork  = 5*norb
cxxx264 - need to scope a do_eigvecs-aware memory for single-proc
C        memeig = 2*mat + 2*mat + 2*norb*nstate + mat + 7*norb
c
c       Close serial solver zhegvx memory needs
      else
c       Parallel, USING SCALAPACK
c
c       Get the process grid dimensions, given the context size
        call EIGPROCS( nprocs, nprow, npcol, npsolver )
c       NOTE: following is a crude approximation to matdist, NOT exact
        call EIGMDIST( matdist, norb, npsolver )
c
        if( DEBUG_MEM )then
          write(IWR,'(a,4i12)')'WKSZX: mat,matdist=',mat,matdist
          write(IWR,'(a,4i6)')'WKSZX: nprocs,nprow,npcol,npsolver=',
     $                                nprocs,nprow,npcol,npsolver
        endif
c
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Simplified appxximate work memory requirements
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      Simplified naive values
        lwzwork = 2*matdist
        lzwork = lwzwork / 2
        lrwork = 2*matdist
        liwork = 6*norb
        if( .not. do_eigvecs )then
c         No eigenvectors, this becomes much less:
          lzwork = (NB+1) * (norb+1)
          lrwork =9*norb
        endif
c
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c    Determine EXACT memory requirements, from Scalapack docs (per ACP)
c   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      ***** Get some of the basic components of the formulas
c
        NN = MAX( norb, NB, 2 )
c       NUMROC is a SCALAPACK function to compute the amount of a 
c       distributed array that is resident on a given processor.
c       It is EXACT.
c       NP0 == number of rows assigned to root process
c           == upper bound on rows assigned to any process.
        NP0   = NUMROC( norb, NB, 0, 0, nprow )
c       NQ0 == number of columns assigned to root process
c           == upper bound on columns assigned to any process.
        NQ0   = NUMROC( norb, NB, 0, 0, npcol )
        NEIG  = norb
        MQ0   = NUMROC( MAX( NEIG,NB,2 ), NB, 0, 0, npcol )
c
        if( DEBUG_MEM )
     $  write(IWR,'(a,6i9)') 'WKSZX: norb,nb,nn,np0,nq0.mq0=',
     $                               norb,nb,nn,np0,nq0,mq0
c
c      ***** Minimum size for PZHEGVX argument WORK: lzwork
c            lzwork is the number of complex elements in wmat(iwZwork)
c
        if( do_eigvecs ) then
          lzwork = norb + NB*(NP0+MQ0+NB)
c
          if( DEBUG_MEM ) write(IWR,*) 'do_eigvecs: base lzwork=',lzwork
          if( OPTIMAL_EIGENVECTORS )then
c           For optimal performance, lzwork needs to be bigger
            if( icontxt .ge. 0 )then
c             Query existing context
              CONTEXT = icontxt
            else
c             Create trial context to get data
              call BLACS_GRIDINIT( CONTEXT, 'R', nprow, npcol )
            endif
            call BLACS_GRIDINFO( CONTEXT, nprow, npcol, myrow, mycol )
            if( DEBUG_MEM )
     $      write(IWR,*) 'optvecs: nprow,npcol=',nprow,npcol,myrow,mycol
            if( myrow .ne. -1 )then
c             OK we're in the BLACS grid/context
              ANB   = PJLAENV( CONTEXT, 3, 'PZHETTRD',
     $             'L', 0, 0, 0, 0 )
              SQNPC = INT( SQRT( DBLE( nprow * npcol ) ) )
              NPS   = MAX( NUMROC( norb, 1, 0, 0, SQNPC ), 2*ANB )
              NHETRD_LWOPT = 2*(ANB+1)*(4*NPS+2) + (NPS+4)*NPS
              NHEGST_LWOPT = (2*NP0+NQ0+NB)*NB
              lzwork = MAX( lzwork, norb+NHETRD_LWOPT, NHEGST_LWOPT )
cxxx264:      Moved gridexit outside this if-block (to balance gridinit) - A bug?
cxxx264              call BLACS_GRIDEXIT( CONTEXT )
            endif
c           Exit the BLACS, we'll recreate it on each call to the solver.
            if( icontxt .lt.0 ) call BLACS_GRIDEXIT( CONTEXT )
c
            if( DEBUG_MEM ) write(IWR,*) 'end optvec: lzwork=', lzwork
          endif
        else
c         Size of WORK array when computing eigenvalues only
          lzwork = norb + MAX( NB*(NP0+1), 3 )
          if( DEBUG_MEM ) write(IWR,*) 'eigval only: lzwork=', lzwork
        endif
c        Factor of 2 in next line  because wmat(iwZwork) is C*16, allocating into R*8
        lwzwork = 2*lzwork
c
c      ***** Minimum size of RWORK: lrwork
c
c       ICEIL(X,Y) is a SCALAPACK function that computes CEILING(DBLE(X)/DBLE(Y)).
        if( do_eigvecs )then
c         Size of RWORK array when computing eigenvalues and eigenvectors
          lrwork = 4*norb + MAX( 5*NN, NP0*MQ0 )
     $         + ICEIL( NEIG, nprow*npcol )*NN
        else
c         Size of RWORK array when computing eigenvalues only
          lrwork = 4*norb + 5*NN
        endif
cxxx264: pas: why was this 2x of lrwork in the ACP code?  Seems wrong (thinks it's c16?).
cxxx264        lrwork = 2*lrwork
c
c        ***** Minumum size of liwork.
c         Noting iwIwork is integer array, wmat is R*8, and sizeof(INTEGER) <= sizeof(REAL*8)
c
        liwork = 6*MAX( norb, nprow*npcol+1, 4 )
c       
c       Close pzhegvx parallel solver work array memory needs
      endif
c
      if( DEBUG_MEM )then
        write(IWR,'(a,3i12)') 'DEVVWKSZ-X: works Z,R,I=',
     $                         lwzwork,lrwork,liwork
        call FLUSH( IWR )
      endif
cxx:end
c
c    That's all Folks!
c
      RETURN
      END
c
c Customization of Fortran mode in Emacs. Must be in last 3000
characters of the file.
c Local Variables:
c fortran-continuation-string: "$"
c fortran-comment-line-start: "c"
c fortran-do-indent: 2
c fortran-if-indent: 2
c End:
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PEIGSOLV
c
c
      subroutine PEIGSOLV( masterlist,do_kppsolve, do_psolv,
     $ icluster, gapp,
     $ do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc, dokdefct,doblas3, ikdefct,nstbulk,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, nk, veck,wtk,
     $ eigval,eigpop,  npop, wksml,
     $ wmat )
c---------------------------------------------------------------
c Purpose: branch to complex solvers for density matrix
c
c Written: Peter A. Schultz, 4-January-2013, for 2.64 (kpp),
c          to accommodate Andrew C. Pineda kpp-solver codes
c
c Revision history:
c  none
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c Input:
      DIMENSION  masterlist(*)
      LOGICAL    do_kppsolve, do_psolv
      LOGICAL    do_eigvecs, do_eigpops
      LOGICAL    dokdefct
      LOGICAL    doblas3
      DIMENSION  veck(3,nk),wtk(nk)
      DIMENSION  elecno(2)
      DIMENSION  nstbulk(2)
c Output:
      DIMENSION  efermi(2),egap(2)
      DIMENSION  eigval(nstate*nk,nspin),eigpop(nstate*nk,nspin)
c Scratch:
      DIMENSION  wksml(norb,*)
      DIMENSION  wmat(*)
C      DIMENSION  wk(*)
c Scratch/parallel:
      INTEGER    icluster(*)
      DIMENSION  gapp(*)
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
cxxx: declare dev-version of peigsolv driver:
      if( IWR .ge. 0 )then
        write(IWR,*) 'PEIGSOLV: '//
     $               'do_eigvecs,do_eigpops,do_kppsolve,do_psolv=',
     $                do_eigvecs,do_eigpops,do_kppsolve,do_psolv
        write(IWR,*) 'DEV/PEIGSOLV: smaller PRHOIJ stripe/transpose'
      endif
c
c  Solve Schroedinger equation for the various k-vectors
      call MPNODES( nprocs )
c
      if( .not. do_kppsolve )then
c       K-parallel serial solves OR K-serial all-processor-parallel solves
c
        call PEIGSOLV0( do_psolv,
     $   icluster, gapp,
     $   do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $   idhamfl,iham0fl,iovlpfl,
     $   ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $   icloseocc, dokdefct, doblas3, ikdefct,nstbulk,
     $   etemp,edegen, elecno, efermi,egap, esumsp,
     $   norb,nstate, nk, veck,wtk,
     $   eigval,eigpop,  npop, wksml,
     $   wmat )
c -->    vmat-6s+
c
      else
c       k-parallel-parallel solves
c
c       Invokes ScaLAPACK solvers within k-point groups.
c       Should enable code to scale to 16 processors/k-point or better
c
        call PEIGSOLVKP( masterlist,
     $   icluster, gapp,
     $   do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $   idhamfl,iham0fl,iovlpfl,
     $   ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $   icloseocc, dokdefct, doblas3, ikdefct,nstbulk,
     $   etemp,edegen, elecno, efermi,egap, esumsp,
     $   norb,nstate, nk, veck,wtk,
     $   eigval,eigpop, npop, wksml,
     $   wmat )
c -->    vmat-6s+
c
      endif
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PEIGSOLV0
c
c
      subroutine PEIGSOLV0( do_psolv,
     $ icluster, gapp,
     $ do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc, dokdefct,doblas3, ikdefct,nstbulk,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, nk, veck,wtk,
     $ eigval,eigpop, npop, wksml,
     $ wmat )
c---------------------------------------------------------------
c Purpose: calling routine to solve for density matrix, given
c          hamiltonian and overlap matrices of current iteration
c
c Written: Peter A. Schultz
c
c Revision history:
c   3Jan13-PAS/2.64: cosmetics, un-drivered, and integrated
c   2May12-ACP/2.64: Added alternate entry point to support non-SCF portions of code
c   9Dec08-PAS/2.62: reduced-mem parallel-complex solve; pschroed w/nk<4
c   1May08-PAS/2.61b: patch bug - purged POP argument (edegen)
c  25Jun07-PAS/2.60: merge tp and serial
c   7Jul06-PAS/2.60: k-parallel dmat generation added
c  14Aug05-PAS/2.59: master-only matrix files, read multi-k Hams
c  15Feb05-PAS/2.59: options to impose 0K occs, and defect sampling
c  25Sep04-PAS/tp0.1: clean-up task parallel routines
c  30Apr03-DMC     : Added a parallel version of RHOIJ
c  22Apr03-DMC     : Added a second version of a parallel SCHROED to take
c                    advantage when nk > 1.
c   7Apr03-DMC     : Made parallel by modifying SCHROED to be parallel.
c                    Does not appear to be sufficient work in POP or RHOIJ
c                    to warrant changing at this time.
c   5Oct01-PAS/2.50: remove backspace, and compress restart
c  20Jul01-PAS/2.49: spin-polarized dft, fpop made scratch, dmat
c  17May01-PAS/2.47: extract listing output (eigenspectrum)
c  23Aug99-PAS/2.38: install zhegvx/selected eigvec solver
c  10May99-PAS/2.35: eigenfunction grid density, pass npop
c   9Apr99-PAS/2.34: complex lapack eigensolver installed
c  31Mar99-PAS/2.33: clean out green fcn stuff, tighten mem use
c   9Feb99-PAS/2.31: pass write unit into routine
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  14Jul98-PAS/2.23: pass option for un-chol'd d-Ham into schroed
c  10Oct97-PAS/2.21: change convention for fpop
c  13Sep93-PAS: rhoofr extracted out of routine, calling arguments
c               rationalized, dimensions fixed, renamed "eigsolv"
c---------------------------------------------------------------
c
c  On input, idhamfl contains full input delta-hamiltonian.
c  To obtain full hamiltonian, delta-H needs to be combined with
c  H0 (computed in setup and saved in iham0fl; overlap in iovlpfl),
c  done in routine "schroed", which also solves the Schroedinger
c  eqn. Routine "pop" populates the resulting levels and finds the
c  fermi energy, then "rhoij" takes the eigenfcns and occupations
c  and builds and saves to disk the density matrix for each k-vector.
c     
c  When this routine is done:
c   (1) eigenvectors are saved to "ivecfl"
c   (2) density matrix is saved to "idmat"
c   (3) energy-weighted density matrix is saved to "iemat"
c   (4) eigenvalue spectrum is returned in "eigval()" (not written)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
c Input:
      LOGICAL    do_psolv
      LOGICAL    do_eigvecs, do_eigpops
      LOGICAL    dokdefct
      LOGICAL    doblas3
      DIMENSION  veck(3,nk),wtk(nk)
      DIMENSION  elecno(2)
      DIMENSION  nstbulk(2)
c Output:
      DIMENSION  efermi(2),egap(2)
      DIMENSION  eigval(nstate*nk,nspin),eigpop(nstate*nk,nspin)
c Scratch:
      DIMENSION  wksml(norb,*)
      DIMENSION  wmat(*)
c Scratch/parallel:
      INTEGER    icluster(*)
      DIMENSION  gapp(*)
c
c Local:
      DATA       itimlog / 5 /
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
c  Solve Schroedinger equation for the various k-vectors
c
      mat = norb*norb
      nmat = nk*mat
c
      call FLGETIWR( IWR )
cxxx: declare dev-version of peigsolv0:
      if( IWR .ge. 0 )
     $write(IWR,*) 'PEIGSOLV0: do_psolv,do_eigvecs,do_eigpops=',
     $                         do_psolv,do_eigvecs,do_eigpops
c
c Partition space within wmat() for schroed:
      iw1 = 1
      iww = iw1 + mat
      nwsch = norb*(4*norb + 2*nstate)
      iwx = iww + nwsch
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      if( iproc .eq. masterp )then
c       Go to start of idhamfl, to get ready for dHam reads
        REWIND( unit=idhamfl )
c       Go to start of ivecfl, to get ready for eigvec writes
        REWIND( unit=ivecfl )
      endif
c
      if( itimlog.gt.0 ) call TIMER( 'eigsolv/begin      ' )
c
C     if( .not. do_psolv )then
      if( nprocs.eq.1 .or.
     $    (  .not. do_psolv .and.
     $      ( norb.lt.800 .or. nk.gt.4 .or. nprocs.lt.(4*nk) )  )      
     $                      )then
c
c       Record path through solver: single-proc solves
        if( IWR .gt. 0 ) write(IWR,*) 'PEIGSOLV: single-proc solves'
c
        do  ispin=1,nspin
c
c         Parallel over k-points - single-proc solves
          call PSCHROED2( do_eigvecs,
     $     idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $     norb,nstate,nk, eigval(1,ispin), wksml,    wksml(1,9), scut,
c -->                                       rwork-7ns seval(n)-s
     $     wmat(iw1), wmat(iww), wmat(iwx) )
c -->      vmat-s     wmat-6ks   iwork-6ns
c
          if( itimlog.gt.1 ) call TIMER( 'eigsolv/schroedkp  ' )
        enddo
c
      else
c       All-proc parallel solvers
c
c       Record path through solver: parallel
        if( IWR.ge.0 ) write(IWR,*) 'PEIGSOLV: all-parallel solves'
c
        call PSCHROED( do_eigvecs,
     $   idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $   norb,nstate,nk, eigval, wksml,
c -->                            wksml-ns
     $   nspin, icluster, gapp,
     $   wmat )
c -->    wmat-s
c
        if( itimlog.gt.1 ) call TIMER( 'eigsolv/pschroed   ' )
c
      endif
c
      if( .not. do_eigpops )then
c       Skip populations and density matrix construction
        RETURN
      endif
c
c  Evaluate density matrix in diagonal (eigfcn) representation
c
      nstk = nstate*nk
c
      npop = 0
      esumsp = zero
      occfull = occlvl(nspin)
c
      do  ispin=1,nspin
        elnumbr = elecno(ispin)
c
        if( ( icloseocc .gt. 0 ) .or. dokdefct )then
c         Use closed-shell occupations, or defect sampling
          call POPCLOSE( icloseocc,dokdefct,ikdefct,
     $     efermi(ispin),egap(ispin), nstbulk(ispin),
     $     etemp,edegen, elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),numpop,
     $     wmat )
c -->      wlvl-s(nstate)
c
        else
c         Use standard metallic occupations
          mspin = 1
          call POP( mspin, popnumbr,
     $     efermi(ispin),egap(ispin),
     $     etemp,        elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),numpop,
     $     wmat(1), wmat(1+nstk),wmat(1+2*nstk),
     $     wmat(1+3*nstk), wmat(1+4*nstk) )
c -->      elvl wlvl klvl nlvl flvl - all scratch (nk*state)
        endif
c
        npop = MAX( npop,numpop )
        esumsp = esumsp + esum1e
      enddo
c
      if( itimlog.gt.0 ) call TIMER( 'eigsolv/pop        ' )
c
      if( .not. do_eigvecs )then
c       Without eigvecs, skip density matrix construction
        RETURN
      endif
c
c  Evaluate density matrix in orbital basis
c
      mat = norb*norb
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idmatfl )
        REWIND( unit=iematfl )
        if( nspin.eq. 2 )then
          REWIND( unit=idmatsfl )
          REWIND( unit=iematsfl )
        endif
      endif
      iw1 = 1
      iw2 = iw1 + mat
      iw3 = iw2 + mat
c
      idmfile = idmatfl
      iemfile = iematfl
      do  ispin=1,nspin
c
        if( nk .ge. nprocs )then
c         K-parallel, single-processor dmat(k) construction
c
          call PRHOIJ2( ivecfl, idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s(n,nst)
c
        else
c         All-processor parallel: distribute each k-point
c
c         Set up sizes of stripe-optimized vectors
          call RHOSTRIP( jstripe, norb, nprocs )
          iw2 = iw1 + jstripe*norb
          iw3 = MAX( iw1 + 2*norb*npop, iw2 + jstripe*norb )
c
          call PRHOIJ( ivecfl,idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s
c ------>  2*norb*npop       ->  2*norb*npop
c ------>  stripe     stripe <-  2*norb*npop
c ------>  stripe     mat
c         Total prhoij mem = max(A,B)
c            A= norb*( max( 2*npop, 2*jstripe ) + 2*npop )
c            B= jstripe*norb + norb^2
c
        endif
c       Switch to spin-dn:
        idmfile = idmatsfl
        iemfile = iematsfl
c
      enddo
      if( itimlog.gt.0 ) call TIMER( 'eigsolv/rhoij      ' )
c
c    That's all Folks!
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
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PERLSOLV
c
c
      subroutine PERLSOLV( icluster, gapp, nlr,nlc,
     $ idosolv,scut, nspin,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc, doblas3,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, wtk, eigval,eigpop,fpop,npop, iwork,ifail, lwork8,
     $ vmat,  tempio )
c---------------------------------------------------------------
c Purpose: solve real (double precision) Schroedinger equation
c
c Written: Dave Raczkowski, Peter A. Schultz, and Dean P. McCullough
c
c Revision history:
c   3Jan13-PAS/2.64: cosmetics, un-drivered, and integrated
c   2May12-ACP/2.64: Added flags for non-scf, and distribute >2GB matrix
c  20Dec07-PAS/2.60d: close BLACS INIT with EXIT in PSOLVSIZ
c  30Aug07-RPM/2.60: infrostructure to reduce parallel memory
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files, pass scut
c  15Feb05-PAS/2.59: add option to impose closed shell occs
c  27Oct04-PAS/tp0.1: clean-up/debug task-parallel routines
c   5Jul03-DMC: Need eigsolv to be norb*nk.  However in wkmem, was only
c               setting up nstate*nk.
c  10Feb03-DMC: Modified extensively to use ScaLAPACK and PBLAS routines.
c   8Mar02-PAS/2.52: clean unused/const d0; closed-spin bug
c   5Oct01-PAS/2.50: remove backspace, and compress restart
c  20Jul01-PAS/2.49: spin-polarized dft, compact dmat
c  21Jun01-PAS/2.48: replace STOPs
c  17May01-PAS/2.47: extract listing output (eigenspectrum)
c  20Dec00-PAS/2.46: patch for small basis problems
c  20Aug99-PAS/2.38: bugfix in call to dsygvx; full dsygvx
c                    implemented, doblas3 passed
c  10May99-PAS/2.35: eigenfunction grid density, pass npop
c   2Apr99-PAS/2.33: tighten big memory usage
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  11Feb99-PAS/2.31: adapted for archival code
c---------------------------------------------------------------
c
c Notes:
c  Call parameter idosolv was used to control which version of the
c  LAPACK routine DSYGV ws used.  Here, we will always use PDSYGVX.
c  The parameter will be used only to control whether all, or only
c  nstate eigenvalues/eigenvectors will be requested.
c  DMC 12/Feb/03
c
c  The size of iwork in the calling routine needs to be increased to
c  at least 6*norb (vs. 5*norb).
c
c  The doblas3 flag is ignored in solvers.
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c     Upper limit on size of array that PDGEMR2D can distribute.
c     MAXBLOCKR= 2**31/sizeof(DOUBLE PRECISION) = 2**31/8.
      PARAMETER  ( maxblockr = 268 435 456 )
c
ctp: Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c
c Input:
      DIMENSION  wtk(1), elecno(nspin)
      DIMENSION  vmat(nlr,nlc,*)
      LOGICAL    doblas3
c Output:
      DIMENSION  efermi(nspin),egap(nspin)
      DIMENSION  eigval(nstate,nspin),eigpop(nstate,nspin)
c Scratch:
      DIMENSION  fpop(norb)
      DIMENSION  iwork(6*norb),ifail(norb)
c Scratch/parallel:
      INTEGER    icluster(*)
      DOUBLE PRECISION  gapp(*)
      DIMENSION  tempio(norb*norb)
c
ctp: parallel local
c Local/parallel
      PARAMETER  ( NB = 32 )
      INTEGER    ICONTXT
      INTEGER    DESCA(9), DESCWORK(9)
      DATA       itimlog / 3 /
      LOGICAL    do_eigvecs, do_eigpops
      DATA       do_eigvecs, do_eigpops / .true., .true. /
      LOGICAL    DEBUG_SOLV
      DATA       DEBUG_SOLV / .true. /
c Local:
      LOGICAL    closeshl
      DIMENSION  npops(2)
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
      CHARACTER  jobz*1
      LOGICAL    pos,neg
      DATA  zero,half,one,two / 0.d0,0.5d0,1.d0,2.d0 /
      DATA  elcerr,elconv / 1.d-6,1.d-11 /
c     for constructing number of occupied levels
      DATA  occnil / 1.d-9 /
c
c >>>> EXECUTABLE CODE:
c
      if( itimlog.gt.0 ) call TIMER( 'erlsolv/begin  ' )
      call FLGETIWR( IWR )
cxxx: declare dev-version of perlsolv
      if( IWR .ge. 0 ) write(IWR,*) 'PERLSOLV: kpp-copmatible version'
c
      if( do_eigvecs )then
c       Compute eigenvalues and eigenvectors
        jobz = 'V'
      else
c       Compute eigenvalues only
        jobz = 'N'
      endif
c
      mat = norb*norb
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c Initialize ScaLAPACK Grid
c
c     Create a close to square grid of processors for solver
      rnprocs = DBLE( nprocs ) + 0.1
      npcol = SQRT( rnprocs )
      nprow = nprocs / npcol
      npsolver = nprow*npcol
c
c  Create BLACS context using system context, 
c  i.e. the local MP communicator
c
      ICONTXT = icomm
      call BLACS_GRIDINIT( ICONTXT, 'R', nprow, npcol )
      call BLACS_GRIDINFO( ICONTXT, nprow, npcol, myrow, mycol )
c
      if( myrow .ne. -1 )then
c       Initialize matrix descriptions, for procs in solver context
c
c       Set up the descriptor for a distributed matrix
        call DESCINIT( DESCA, norb,norb, NB,NB, 0,0,
     $   ICONTXT, nlr, INFO )
        if( info.ne.0 ) call STOPXERR( 'perlsolv: INIT for DESCA' )
c        nlocr = NUMROC(norb,NB,myrow,0,nprow)
c        nlocc = NUMROC(norb,NB,mycol,0,npcol)
c       This call to PSOLVSIZES is not uesed (yet).
        call PSOLVSIZES( ICONTXT, norb,nstate,nlocr,nlocc,mylwork )
c
c       Set up the descriptor for full local matrix
        call DESCINIT( DESCWORK, norb,norb, norb,norb, 0,0,
     $   ICONTXT, norb, INFO )
        if( info.ne.0 ) call STOPXERR( 'perlsolv: INIT for DESCWORK' )
      endif
c
c Set files to correct locations
c
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idhamfl )
      endif
c
      do 1000 ispin=1,nspin
        if( myrow .eq. -1 ) goto 800
c
c       Read in array to work space, and  distribute to all processors.
c       If array is too big to fit on one processor, this can be done in
c       several steps.
c
c       If nspin = 1, this read not done in the original serial code.
c       However, I (DMC) do not  know where it comes from, and the fill
c       might not be distributed across the processors, thus reload it.
c       Found I also need to rewind it.
c
c       Retrieve all of delta-Ham onto master ...
        if( iproc .eq. masterp )then
          call READBIG( idhamfl, mat, tempio )
        endif
c        ... and distribute d-Ham onto all processors.
        if( mat .le. maxblockr )then
c          Matrix small enough to distribute from full local storage
           call PDGEMR2D( norb,norb,
     $          tempio, 1,1, DESCWORK,
     $          vmat(1,1,1), 1,1, DESCA,  ICONTXT )
        else
c         Distribute the matrix in (column) stripes that are smaller than 2GB each.
c         Should be good as long the distributed matrix is no larger than 2GB on any processor.
          nblock = maxblockr / norb
          nleft = norb - nblock
          noffset = 1
          moffset = 1
  151     continue
          if( noffset.le.norb )then
            call PDGEMR2D( norb,nblock,
     $           tempio(moffset), 1,1, DESCWORK,
     $           vmat(1,1,1), 1, noffset, DESCA,  ICONTXT )
            noffset = noffset + nblock
            moffset = moffset + norb*nblock
            if( nleft.ge.nblock ) then
c             Have enough to do another full block.
              nleft = nleft - nblock
            else
c             Have less than a full block to send.
              nblock = nleft
              nleft = 0
            endif
            goto 151
          endif
        endif
c
c       Retrieve iter-independent Hamiltonian ...
        if( iproc .eq. masterp )then
          REWIND( unit=iham0fl )
          call READBIG( iham0fl, mat, tempio )
        endif
c        ... and distribute over processors:
        if( mat .le. maxblockr )then
          call PDGEMR2D( norb,norb,
     $         tempio, 1,1, DESCWORK,
     $         vmat(1,1,2), 1,1, DESCA,  ICONTXT )
        else
c         Distribute the matrix in (column) stripes that are smaller than 2GB each.
c         Should be good as long the distributed matrix is no larger than 2GB on any processor.
          nblock = maxblockr / norb
          nleft = norb - nblock
          noffset = 1
          moffset = 1
  152     continue
          if( noffset.le.norb )then
            call PDGEMR2D( norb,nblock,
     $           tempio(moffset), 1,1, DESCWORK,
     $           vmat(1,1,2), 1, noffset, DESCA,  ICONTXT )
            noffset = noffset + nblock
            moffset = moffset + norb*nblock
            if( nleft.ge.nblock )then
c             Have enough to do another full block.
              nleft = nleft - nblock
            else
c             Have less than a full block to send.
              nblock = nleft
              nleft = 0
            endif
            goto 152
          endif
        endif
c
c       Retrieve overlap matrix ...
        if( iproc .eq. masterp )then
          REWIND( unit=iovlpfl )
          call READBIG( iovlpfl , mat, tempio )
        endif
c        ... and distribute over processors:
        if( mat .le. maxblockr )then
           call PDGEMR2D( norb,norb,
     $          tempio, 1,1, DESCWORK,
     $          vmat(1,1,3), 1,1, DESCA,  ICONTXT )
        else
c         Distribute the matrix in (column) stripes that are smaller than 2GB each.
c         Should be good as long the distributed matrix is no larger than 2GB on any processor.
          nblock = maxblockr / norb
          nleft = norb - nblock
          noffset = 1
          moffset = 1
  153     continue
          if( noffset.le.norb )then
            call PDGEMR2D( norb,nblock,
     $           tempio(moffset), 1,1, DESCWORK,
     $           vmat(1,1,3), 1, noffset, DESCA,  ICONTXT )
            noffset = noffset + nblock
            moffset = moffset + norb*nblock
            if( nleft.ge.nblock )then
c             Have enough to do another full block.
              nleft = nleft - nblock
            else
c             Have less than a full block to send.
              nblock = nleft
              nleft = 0
            endif
            goto 153
          endif
        endif
c
c       Combine iter-independent Ham and del-Ham to get total Ham
c
        do  i=1,norb
c          do  ii=1,norb
c            vmat(ii,i,2) = vmat(ii,i,2) + vmat(ii,i,1)
c          enddo
          call PDAXPY( norb, one, vmat(1,1,1), 1,i, DESCA, 1,
     $                            vmat(1,1,2), 1,i, DESCA, 1 )
        enddo
c
crpm  Do we really need this now that we're using PDAXPY? Guess it can't hurt
c        call BLACS_BARRIER( ICONTXT, 'All' )
c
c  Total hamiltonian is now in vmat(,,2), overlap in vmat(,,3)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                    Solve Schroedinger Equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c lwork is workspace available within vmat(,,4), and normally you
c would set lwork=mat(dist) and that would be more than ample.
c But for small problems (e.g. H-atom), we crash as the routine wants
c more than the small mat it gets.  AAAArrrrrgggghhhh!
c Hence I pass in the full free space available ...
c
        lwork = lwork8
c
        if( idosolv .eq.1 )then
          jstate = norb
        elseif( idosolv .eq. 2 )then
          jstate = nstate
        else
          call STOPXERR( 'perlsolv: invalid solver option idosolv' )
        endif
c
        rmone = 2*PDLAMCH( ICONTXT, 'S' )
c
        call PDSYGVX( 1,jobz, 'I','L',
     $   norb,  vmat(1,1,2),1,1,DESCA,  vmat(1,1,3),1,1,DESCA,
     $   zero,zero,1,jstate,  rmone, nstout, nvout,
     $   fpop           ,rmone,
     $   vmat(1,1,1),1,1,DESCA,  vmat(1,1,4),
c -->    Z                       WORK
     $   lwork,iwork,6*norb,ifail,icluster,gapp,info )
c
        if( info .ne. 0 )then
          call STOPXERR( 'perlsolv: pdsygvx error' )
        endif
        call DCOPY( nstate, fpop,1, eigval(1,ispin),1 )
c
c       If no eigvecs, we are done with this spin-Ham:
        if( .not. do_eigvecs ) goto 1000
c
c       Save the the eigenvectors in vmat(,,1) to ivecfl
c
c       Bring back eigenvectors to master processor
        if( mat .le. maxblockr )then
          call PDGEMR2D( norb,norb,
     $         vmat(1,1,1), 1,1, DESCA,
     $         tempio, 1,1, DESCWORK,  ICONTXT )
         else
c          Distribute the matrix in (column) stripes that are smaller than 2GB each.
c          Should be good as long the distributed matrix is no larger than 2GB on any processor.
           nblock = maxblockr / norb
           nleft = norb - nblock
           noffset = 1
           moffset = 1
  154      continue
           if( noffset.le.norb )then
             call PDGEMR2D( norb,nblock,
     $            vmat(1,1,1), 1,1, DESCA, 
     $            tempio(moffset), 1, noffset, DESCWORK,  ICONTXT )
             noffset = noffset + nblock
             moffset = moffset + norb*nblock
             if( nleft.ge.nblock )then
c              Have enough to do another full block.
               nleft = nleft - nblock
             else
c              Have less than a full block to send.
               nblock = nleft
               nleft = 0
             endif
             goto 154
           endif
        endif
c
c       Processors not in ScaLAPACK solver context return here:
  800   continue
Cc
Cc       Record eigenvectors and values
Cc
Cc       Broadcast eigenvectors from master to all processors
C        call MPI_Bcast( vmat(1,1,4), norb*nstate, MPI_DOUBLE_PRECISION,
C     $   0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast eigvec' )
Cc
c       Save eigenvectors to disk for later
        if( do_eigvecs .and. ( iproc .eq. masterp ) )then
          call WRITBIG( ivecfl, norb*nstate, tempio )
        endif
c
c       Close loop over spins:
 1000 continue
c
c  If some processors not used in the ScaLAPACK context ...
      if( nprocs .ne. npsolver )then
c        ... also need to broadcast eigvals to those not having them
        call MPI_Bcast( eigval(1,1),nstate*nspin, MPI_DOUBLE_PRECISION,
     $   0, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast eigval' )
      endif
ctp:
      if( itimlog.gt.1 ) call TIMER( 'erlsolv/schroed' )
c
      if( .not. do_eigpops )then
c       Skip populations and density matrix
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   Construct level populations, i.e. dmat in diagonal representation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      npop = 0
      esumsp = zero
      occfull = occlvl(nspin)
c
c     Set closed shell flag if to populate using 0K
      closeshl = .false.
      if( icloseocc .eq. 2 ) closeshl = .true.
c
      do 2000 ispin=1,nspin
        elnumbr = elecno(ispin)
c
        call POP1K( closeshl,
     $   fermilvl,gap,ehomo,elumo,
     $   etemp,edegen, elnumbr,occfull,
     $   nstate,numpop, wtk(1),
     $   eigval(1,ispin),eigpop(1,ispin), fpop )
c
        efermi(ispin) = fermilvl
        egap(ispin) = gap
        npops(ispin) = numpop
        npop = MAX( npop, numpop )
c
        esum1e = zero
        do  i=1,numpop
          esum1e = esum1e + eigpop(i,ispin)*eigval(i,ispin)
        enddo
        esumsp = esumsp + wtk(1)*esum1e
 2000 continue
c
      if( itimlog.gt.1 ) call TIMER( 'erlsolv/pop    ' )
c
      if( .not. do_eigvecs )then
c       Without eigvecs, skip density matrix construction
        RETURN
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       Construct density-matrix, e-matrix over basis functions
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if( nspin .eq. 2 )then
c       Need to recover eigvecs ... they have been stepped on
        if( iproc .eq. masterp )then
          REWIND( unit=ivecfl )
        endif
      endif
      if( iproc .eq. masterp )then
        REWIND( unit=idmatfl )
        REWIND( unit=iematfl )
        if( nspin .eq. 2 )then
c         Rewind the spin-dn mats files
          REWIND( unit=idmatsfl )
          REWIND( unit=iematsfl )
        endif
      endif
      idmfile = idmatfl
      iemfile = iematfl
c
      do 4000 ispin=1,nspin
        if( myrow .eq. -1 ) goto 2995
c
        numpop = npops(ispin)
        if( nspin .eq. 2 )then
c         Retrieve spin-dependent eigenvectors (in mem if non-spin)
          if( iproc .eq. masterp )then
            call READBIG( ivecfl, norb*nstate, tempio )
          endif
c          ... and distribute over processors
          if( mat .le. maxblockr )then
c           Matrix is small enough to be distributed monolithically via BLACS ...
            call PDGEMR2D( norb,norb,
     $           tempio, 1,1, DESCWORK,
     $           vmat(1,1,1), 1,1, DESCA,  ICONTXT )
          else
c           Matrix is larger than BLACS can handle in one chunk ...
c           Distribute matrix in (column) stripes smaller than 2GB each.
c           Should be good as long the distributed matrix <2GB on any processor.
            nblock = maxblockr / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 1
 2155       continue
            if( noffset.le.norb )then
              call PDGEMR2D( norb,nblock,
     $             tempio(moffset), 1,1, DESCWORK,
     $             vmat(1,1,1), 1, noffset, DESCA,  ICONTXT )
              noffset = noffset + nblock
              moffset = moffset + norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 2155
            endif
          endif
        endif
c
        if( numpop .lt. 1 )then
c         No occupied levels: dmat=emat=0, and skip construction
crpm      Testing whether we can use PDLASET instead of MKZERO here:
c          call MKZERO( mat, vmat(1,1,3) )
c          call MKZERO( mat, vmat(1,1,4) )
          call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,3),1,1,DESCA )
          call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,4),1,1,DESCA )
          goto 2990
        endif
c
c       Scale level occupations with appropriate sample weights
        wt = wtk(1)
        do  n=1,numpop
          fpop(n) = wt*eigpop(n,ispin)
        enddo
c
c >>>>> Closed shell, use Raczkowski blas code:
c
c       Make a copy of the eigenvectors ...
        do  i=1,norb
           call PDCOPY( norb,vmat(1,1,1), 1,i, DESCA, 1,
     $                       vmat(1,1,2), 1,i, DESCA, 1 )
        enddo
c
c        ni = (numpop-1) / (NB*npcol)
c        ... and scale copy vmat(:,i,2) by occupations fpop(i)
c        do  i=0,ni
c          ki = i*NB
c          kii = (i*npcol + mycol)*NB
c          nj = MIN( NB, numpop-kii )
c          if( nj .ge. 0 )then
c            do  j=1,nj
c              do  jj=1,norb
c                vmat(jj,ki+j,2) = vmat(jj,ki+j,2)*fpop(kii+j)
c              enddo
c            enddo
c          endif
c        enddo
crpm    Replace the explicit math with a call to PDSCAL:
        do  i=1,numpop
          call PDSCAL( norb, fpop(i),vmat(1,1,2),1,i,DESCA, 1 )
        enddo
c
c       To make dscal happy (not return a floating point error)
crpm    Testing whether PDLASET will work here:
c        call MKZERO( mat, vmat(1,1,3) )
c        call MKZERO( mat, vmat(1,1,4) )
        call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,3),1,1,DESCA )
        call PDLASET( 'A',norb,norb,zero,zero,vmat(1,1,4),1,1,DESCA )
c
c       Make density matrix:
        call PDGEMM( 'n','t', norb,norb,numpop,one,
     $   vmat(1,1,2),1,1,DESCA, vmat(1,1,1),1,1,DESCA,
     $   zero, vmat(1,1,3),1,1,DESCA )
c
c       Make e-matrix:
c
c        ... scale pop-wt'd eigenvectors with eigenvalues
c        do  i=0,ni
c          ki = i*NB
c          kii = (i*npcol + mycol)*NB
c          nj = MIN( NB, numpop-kii )
c          if( nj .ge. 0 )then
c            do  j=1,nj
c              do  jj=1,norb
c                vmat(jj,ki+j,2) = vmat(jj,ki+j,2)*eigval(kii+j,ispin)
c              enddo
c            enddo
c          endif
c        enddo
crpm    Replace the explicit math with a call to PDSCAL:
        do  i=1,numpop
          call PDSCAL( norb, eigval(i,ispin),vmat(1,1,2),1,i,DESCA, 1 )
        enddo
c
        call PDGEMM( 'n','t', norb,norb,numpop,one,
     $   vmat(1,1,2),1,1,DESCA, vmat(1,1,1),1,1,DESCA,
     $   zero, vmat(1,1,4),1,1,DESCA )
c
c       Scale dmat and emat by weight function:
        wt = wtk(1) 
crpm        call DSCAL( 2*mat, wt, vmat(1,1,3),1 )
        do  i=1,norb
          call PDSCAL( norb,wt,vmat(1,1,3),1,i,DESCA, 1 )
          call PDSCAL( norb,wt,vmat(1,1,4),1,i,DESCA, 1 )
        enddo
c
c     Zero out the imaginary portion (upper triangle) of the two matrices
crpm  Using the PDLASET routine to do this. This is a slightly non-standard
c     way to call the routine, but it insures that the upper triangle
c     without the diagonal is zeroed. The default call also zeros the diagonal,
c     so for this call, you basically give arguments for the first off-diag
c     element, and then zero N-1 elements in the triangle.
        call PDLASET( 'U',norb-1,norb-1,zero,zero,
     $       vmat(1,1,3),1,2,DESCA )
        call PDLASET( 'U',norb-1,norb-1,zero,zero,
     $       vmat(1,1,4),1,2,DESCA )
c        nic = (norb-1) / (NB*npcol)
c        nir = (norb-1) / (NB*nprow)
c        do  ic=0,nic
c          do  ir=0,nir
c            jc = ic*npcol + mycol
c            jr = ir*nprow + myrow
c            nc = MIN( NB, norb - jc*NB )
c            nr = MIN( NB, norb - jr*NB )
c            kr = ir*NB
c            kc = ic*NB
c            if( jc .gt. jr )then
c              do  i=1,nr
c                do  j=1,nc
c                  vmat(i+kr,j+kc,3) = zero
c                  vmat(i+kr,j+kc,4) = zero
c                enddo
c              enddo
c            elseif( jc. eq. jr )then
c              do  i=1,nr-1
c                do  j=i+1,nc
c                  vmat(i+kr,j+kc,3) = zero
c                  vmat(i+kr,j+kc,4) = zero
c                enddo
c              enddo
c            endif
c          enddo
c        enddo
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c               Write out density-matrix, e-matrix:
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       If we skipped constructing zero-matrices, return here:
 2990   continue
c
c       Non-context processors return here:
 2995   continue
        if( itimlog.gt.2 ) call TIMER( 'erlsolv/dmat compute' )
c
c       Save density matrix dmat to idmatfl:
c
        if( myrow .ne. -1 )then
c         Reassemble full dmat on master from distributed dmat
          if( mat .le. maxblockr )then
            call PDGEMR2D( norb,norb,
     $           vmat(1,1,3), 1,1, DESCA,
     $           tempio, 1,1, DESCWORK,  ICONTXT )
          else
c           Should be good as long distributed matrix <2GB on any processor.
            nblock = maxblockr/norb
            nleft = norb - nblock
            noffset = 1
            moffset = 1
 3156       continue
            if( noffset.le.norb )then
              call PDGEMR2D( norb,nblock,
     $             vmat(1,1,3), 1, noffset, DESCA,
     $             tempio(moffset), 1,1, DESCWORK,  ICONTXT )
              noffset = noffset + nblock
              moffset = moffset + norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 3156
            endif
          endif
        endif
        if( itimlog.gt.3 ) call TIMER( 'erlsolv/dmat collect' )
c
Cc       Broadcast dmat from master to all processors
C        call MPI_Bcast( vmat(1,1,1), mat, MPI_DOUBLE_PRECISION,
C     $   0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast dmat' )
C        if( itimlog.gt.3 ) call TIMER( 'erlsolv Bcast dmat' )
c
        if( iproc .eq. masterp )then
c         Only master saves dmat to idmatfl
          call WRITBIG( idmfile, mat, tempio )
          if( itimlog.gt.3 ) call TIMER( 'erlsolv/dmat write' )
        endif
c
c       Save e-matrix emat to iematfl:
c
        if( myrow .ne. -1 )then
c         Reassemble full emat on master from distributed emat
          if( mat .le. maxblockr )then
            call PDGEMR2D( norb,norb,
     $           vmat(1,1,4), 1,1, DESCA,
     $           tempio, 1,1, DESCWORK,  ICONTXT )
          else
c           Should be good as long distributed matrix <2GB on any processor.
            nblock = maxblockr / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 1
 3157       continue
            if( noffset.le.norb )then
              call PDGEMR2D( norb,nblock,
     $             vmat(1,1,4), 1, noffset, DESCA,
     $             tempio(moffset), 1,1, DESCWORK,  ICONTXT )
              noffset = noffset + nblock
              moffset = moffset + norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 3157
            endif
          endif
        endif
        if( itimlog.gt.3 ) call TIMER( 'erlsolv/emat collect' )
c
Cc       Broadcast emat from master to all processors
C        call MPI_Bcast( vmat(1,1,1), mat, MPI_DOUBLE_PRECISION,
C     $   0, icomm, ierr )
C        if( ierr.ne.0 ) call STOPXERR( 'perlsolv: Bcast emat' )
C         if( itimlog.gt.3 ) call TIMER( 'erlsolv/emat Bcast' )
c
        if( iproc .eq. masterp )then
          call WRITBIG( iemfile, mat, tempio )
          if( itimlog.gt.3 ) call TIMER( 'erlsolv/emat write' )
        endif
c
c       Switch to spin-dn:
        idmfile = idmatsfl
        iemfile = iematsfl
c
c       Close loop over spins:
 4000 continue
c
cxxx: should GRIDEXIT be used to clear context for out-of-context procs, too?
      if( myrow .ne. -1 ) call BLACS_GRIDEXIT( ICONTXT )
c
      if( itimlog.gt.1 ) call TIMER( 'erlsolv/rhoij  ' )
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> DISTRMAT
c
c
      subroutine DISTRMAT( norb, mat, maxblockr,
     $ DESCWORK,DESCA,ICONTXT,
     $ rmatfull,rmatdist )
c---------------------------------------------------------------
c Purpose: distribute matrix over processes, with >2GB handling
c
c **** this is just draft, not tested or invoked yet *****
c
c Revision history:
c   2Jan12-PAS/2.64: extracted from ACP code in persolv
c---------------------------------------------------------------
c
      IMPLICIT NONE
c
c Input:
      INTEGER  norb, mat, maxblockr
      INTEGER  ICONTXT
      INTEGER  DESCA(9), DESCWORK(9)
      REAL*8   rmatfull(*)
c Output:
      REAL*8   rmatdist(*)
c Local declarations:
      INTEGER  nblock, nleft,noffset,moffset
c
c >>>> EXECUTABLE CODE:
c
crpm_tp          call PDGEMR2D( norb,norb, vmat(1,1,4), 1,1,
c
      if( mat .le. maxblockr )then
c       Matrix is small enough to be distributed monolithically via BLACS ...
c
        call PDGEMR2D( norb,norb,
     $                 rmatfull, 1,1, DESCWORK,
     $                 rmatdist, 1,1, DESCA, ICONTXT )
c
      else
c       Distribute the matrix in (column) stripes that are smaller than 2GB each.
c       Should be good as long the distributed matrix is no larger than 2GB on any processor.
c
        nblock = maxblockr / norb
        nleft = norb - nblock
        noffset = 1
        moffset = 1
  151   continue
        if( noffset.le.norb )then
          call PDGEMR2D( norb,nblock,
     $                   rmatfull(moffset), 1,1,       DESCWORK,
     $                   rmatdist,          1,noffset, DESCA, ICONTXT )
          noffset = noffset + nblock
          moffset = moffset + norb*nblock
          if( nleft.ge.nblock ) then
c           Have enough to do another full block.
            nleft = nleft - nblock
          else
c           Have less than a full block to send.
            nblock = nleft
            nleft = 0
          endif
          goto 151
        endif
      endif
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSOLVSIZES
c
c
      subroutine PSOLVSIZES( ICONTXT, norb,neig,nlocrow,nloccol,lwork )
c---------------------------------------------------------------
c Purpose: determine solver needs for serial solver (dsygvx)
c Written: Richard P. Muller
c Revision history:
c  17Jan13-PAS/2.64: tidied up a little during kpp solver work
c---------------------------------------------------------------
c
c On input:
c    ICONTXT:  <0 indicates need to create and probe context for info
c            0,>0 indicates probe existing context (ICONTXT) for info
c Note:
c     The output here is only the _minimum_ required to solve,
c     and not the amount needed to make this optimal - 17Jan13
c
      IMPLICIT NONE

crpm  This is only here for the MPI_ALLREDUCE command. If we don't need
c     the command, the include statement can go:
      INCLUDE  'mpif.h'
c
      INTEGER     NB
      PARAMETER  (NB=32)
c     NB      Blocking size (=32)
c
c Functions (these are scalapack functions)
      INTEGER  NUMROC,ICEIL
c
c Input arguments
      INTEGER  ICONTXT
      INTEGER  norb,neig
c     norb    Number of orbitals
c     neig    Number of eigenvalues to compute
c
c Output arguments
      INTEGER  nlocrow,nloccol,lwork
c     nlocrow Number of rows required for local storage
c     nloccol Number of columns required for local storage
c     lwork   Size needed for the secratch array WORK (minimum)
c
c Local variables
      INTEGER  nn,np0,mq0,nprow,npcol,myrow,mycol,nenb,
     $         nprocs,iproc,masterp,icomm,npsolver,JCONTXT
c
      DOUBLE PRECISION  rnprocs
      LOGICAL  init_blacs
c
c >>>> EXECTUABLE CODE:
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c  Create a close to square grid of processors for solver
      rnprocs = DBLE( nprocs ) + 0.1
      npcol = SQRT( rnprocs )
      nprow = nprocs / npcol
      npsolver = nprow*npcol
c
c  Create the BLACS context if we need to:
      if( ICONTXT .lt. 0 )then
c       We do not have an existing context ... create one
        JCONTXT = icomm
        call BLACS_GRIDINIT( JCONTXT, 'R', nprow, npcol )
        init_blacs = .true.
      else
c       We have an existing context, use it ...
        JCONTXT = ICONTXT
        init_blacs = .false.
      endif
c  Get information about BLACS processor grid:
      call BLACS_GRIDINFO( JCONTXT, nprow, npcol, myrow, mycol )
c
      nn = MAX( norb,NB )
      nenb = MAX( neig,NB )
      np0 = NUMROC( nn,NB,0,0,nprow )
      mq0 = NUMROC( nenb,NB,0,0,npcol )
c
      nlocrow = NUMROC( norb,NB,myrow,0,nprow )
      nloccol = NUMROC( norb,NB,mycol,0,npcol )
c
      lwork = 5*norb + MAX( 5*nn,np0*mq0 + 2*NB*NB)
     $     + ICEIL( neig,nprow*npcol )*nn
c
c     nlocrow/col have to be at least 1, since we're using them
c     as the leading dimension of matrices
      nlocrow = MAX( nlocrow, 1 )
      nloccol = MAX( nloccol, 1 )
c
      if( init_blacs )then
c       Close up this BLACS Context if we created it
        if( myrow .ne. -1 ) call BLACS_GRIDEXIT( JCONTXT )
      endif
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
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOSTRIP
c
c
      subroutine RHOSTRIP( jstripe, norb, nprocs )
c---------------------------------------------------------------
c Purpose: width of mat-stripe (larger or "equitable" widths)
c Written: Peter A. Schultz, 31Dec19, for 2.68b
c---------------------------------------------------------------
      IMPLICIT NONE
      INTEGER  norb,nprocs,    jstripe
c
      jstripe = (norb+nprocs-1) / nprocs
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOIJ
c
c
      subroutine PRHOIJ( ivecfl,idmatfl,iematfl, norb,nstate,nk,
     $ eigval,eigpop,fpop,npop, wtk,
     $ dmat,emat, eigvec )
c---------------------------------------------------------------
c Purpose: construct density matrix in orbital basis, parallel
c
c Written: Peter A. Schultz,  9-April-1999, for v2.34
c
c Revision history:
c  29Mar20-PAS/2.68e: READB1 to read partial ivecfl(npop<nstate)
c  31Dec19-PAS/2.68b: reduced-mem eigvec(npop)
c  11May07-PAS/2.60: merged tp and serial
c  12Jul06-PAS/2.60: stripe-distributed dmat generation
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  30Jul05-PAS/2.59: master-only matrix files
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  DMC: task parallel
c  18Jan02-PAS/2.51b: replace code that led to segmentation flt.
c  20Jul01-PAS/2.49: spin-polarized dft, fpop made scratch, dmat
c  10May99-PAS/2.35: eigenfunction grid density
c   9Apr99-PAS/2.34: complex lapack eigensolver implemented
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  10Oct97-PAS/2.21: change fpop/wtk convention
c  12Aug94-PAS/2.11: cleanup/purge complex/clean args
c---------------------------------------------------------------
c
c mem/cpu optimization notes:
c  MUST use transpose despite (2*norb*npop) mem cost
c    - stride in dm-gen loops leads to 4-5 X cost increase
c  pre-weight eigvec(n;) wtih sqrt(wtk*pop) in transpose not helpful
c    - eliminates 1 (second-order) op in dm-gen, but only negligible gain
c    - code risk, might want neg wtk in future for some unforeseen purpose
c Possible future ideas:
c   skip emat in scf, only do emat after scf (complicated, extra solve)
c   block or block-cyclic dm-mat gen (complicated code)
c Now, this mem cost below core solvers, cpu cost is marginal, leave alone
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
      DIMENSION  eigval(nstate,nk),eigpop(nstate,nk),wtk(nk)
c Scratch arrays:
      DIMENSION  fpop(nstate)
c
c     *** NB: dmat used as 2*norb*npop   scratch below, in addition
c         to norb^2 output.  Currently, the idle emat is contiguous
c         to dmat, so this is safe, even for npop  =norb.
      DIMENSION  dmat(norb,*),emat(norb,norb)
      DIMENSION  eigvec(2,npop  ,norb)
c
c Local declarations:
      DATA       lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      mat = norb*norb
c
c  Get Local MP info
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
      if( lstout.gt.0 .and. nprocs.gt.1 .and. iproc.eq.masterp )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DMAT: all-parallel stripe/merge dmat(k)'
        write(IWR,*) 'PRHOIJ(06): compact stripe+transpose'
      endif
c
      j1 = 1
      jn = norb
      jlen = norb
c Take care of arithmetic for stripe mgmt, using (1):
c   (1) most equitable, rather than
c   (2) blocked (i.e. BLACS compatible), or (3) 0-favored
c  All assign mem allocations of wider stripe (everyone has same larger size)
c  Thinner work allocations (shorter width by 1) will go to first processors
c  First, determine width of (wider) stripe in matrix
C      jstripe = (norb+nprocs-1) / nprocs ! replaced with call (so
C      peigsolv can use
      call RHOSTRIP( jstripe, norb, nprocs )
c  Allocate memory locations (not used - idea for combined mem later)
C      j0dmat = 0
C      j0emat = j0dmat + jstripe
C      j0vec  = j0emat + jstripe
c
      jlen = jstripe - 1
      nshort = ( jstripe * nprocs ) - norb
      j1 = 1 + iproc*jlen
      if( iproc.ge.nshort )then
c       Last nlong procs get bigger stripe (by 1), assume 0:(nprocs-1) id
        jlen = jlen + 1
        j1 = j1 + iproc - nshort
      endif
      jn = j1 - 1 + jlen
c
c    Bloch vector loop:
c
      do 1000 k=1,nk
c
c       Retrieve (only occupied npop, not all nstate) eigenvectors ...
c       NB: dmat used as 2*norb*npop   scratch space here
        if( iproc .eq. masterp )then
c         A read of the partial vector, npop out of nstate saved
          call READB1 ( ivecfl, 2*norb*nstate, 2*norb*npop, dmat )
        endif
        call MPBCAST8( masterp, 2*norb*npop  , dmat  , icomm )
c
c       Every processor works on every k-point
c
c       Transpose (only occupied!) eigvecs for use below:
        call TREIGVEC( norb,npop  , dmat,eigvec )
        numpop = npop
c
c       Scale fpop by k-point sampling weight for use in dmat loops
        wt = wtk(k)
        do  n=1,numpop
          fpop(n) = wt*eigpop(n,k)
        enddo 
c       
c       Clear out dmat, emat stripes
        nmatstr = jlen*norb
        call MKZERO( nmatstr, dmat )
        call MKZERO( nmatstr, emat )
c       
        jdm = 0
        jem = 0
        do 400 j=j1,jn
          jdm = jdm + 1
          jem = jem + 1
c
c         Do imaginary part first ...
          do 210 i=1,j-1
            dmi = zero
            emi = zero
            do 200 n=1,numpop 
              occnk = fpop(n)
              ti = occnk*( eigvec(2,n,i)*eigvec(1,n,j)
     $                   - eigvec(1,n,i)*eigvec(2,n,j) )
              dmi = dmi + ti 
              emi = emi + ti*eigval(n,k) 
  200       continue
            dmat(i,jdm) = dmi 
            emat(i,jem) = emi
  210     continue
c           
c          ... and then switch to real part:
          do 310 i=j,norb
            dmr = zero
            emr = zero
            do 300 n=1,numpop
              occnk = fpop(n)
              tr = occnk*( eigvec(1,n,i)*eigvec(1,n,j)
     $                   + eigvec(2,n,i)*eigvec(2,n,j) )
              dmr = dmr + tr
              emr = emr + tr*eigval(n,k)
  300       continue
            dmat(i,jdm) = dmr
            emat(i,jem) = emr
  310     continue
c
c         Close loop over j/norb
  400   continue
c
        if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat compute' )
c
c >>>>> Store dmat and emat away
c
c Here, dmat&emat in distributed stripes and we are done with eigvecs()
c This version, we merge to masterp and write to disk to save
c
c       Must do emat() first, as it is in the way of the full dmat
        if( nprocs.gt.1 )then
c         Merge emat to masterp
          nvec = mat
          nmatstr = jlen*norb
c         Now stripe frontloaded locally, so  ...
c          ivec1 = (j1-1)*norb + 1
          ivec1 = 1
c
          call MPMERGER8( nprocs, masterp,iproc,
     $                    nvec,ivec1,nmatstr, emat,
     $                    icomm )
          if( lstout.gt.3 ) call TIMER( 'prhoij/emat merged' )
        endif
c
c       Save full emat to iematfl for use in the force section
        if( iproc .eq. masterp )then
          call WRITBIG( iematfl, mat, emat )
          if( lstout.gt.3 ) call TIMER( 'prhoij/emat writ' )
        endif
c       And now done with the emat(), so its space is free, switch to dmat ...
c
        if( nprocs.gt.1 )then
c         Merge dmat to masterp
          nvec = mat
          nmatstr = jlen*norb
c         Now stripe frontloaded locally, so  ...
c          ivec1 = (j1-1)*norb + 1
          ivec1 = 1
c
          call MPMERGER8( nprocs, masterp,iproc,
     $                    nvec,ivec1,nmatstr, dmat,
     $                    icomm )
          if( lstout.gt.3 ) call TIMER( 'prhoij/dmat merged' )
        endif
c
c       Save full dmat to idmatfl for use later
        if( iproc .eq. masterp )then
          call WRITBIG( idmatfl, mat, dmat )
          if( lstout.gt.3 ) call TIMER( 'prhoij/dmat writ' )
        endif
c
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> RHOIJ2
c
c
      subroutine PRHOIJ2( ivecfl,idmatfl,iematfl, norb,nstate,nk,
     $ eigval,eigpop,fpop,npop, wtk,
     $ dmat,emat, eigvec )
c---------------------------------------------------------------
c Purpose: construct density matrix in orbital basis, k-parallel
c
c Written: Peter A. Schultz,  9-April-1999, for v2.34
c
c Revision history:
c  10Jul06-PAS/2.60: k-parallel version prhoij2 written
c  30Jul05-PAS/2.59: master-only matrix files
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  DMC: task parallel
c  18Jan02-PAS/2.51b: replace code that led to segmentation flt.
c  20Jul01-PAS/2.49: spin-polarized dft, fpop made scratch, dmat
c  10May99-PAS/2.35: eigenfunction grid density
c   9Apr99-PAS/2.34: complex lapack eigensolver implemented
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c  10Oct97-PAS/2.21: change fpop/wtk convention
c  12Aug94-PAS/2.11: cleanup/purge complex/clean args
c---------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
      DIMENSION  eigval(nstate,nk),eigpop(nstate,nk),wtk(nk)
c Scratch arrays:
      DIMENSION  fpop(nstate)
c     *** NB: dmat used as 2*norb*npop   scratch below, in addition
c         to norb^2 output.  Currently, the idle emat is contiguous
c         to dmat, so this is safe, even for npop  =norb.
      DIMENSION  dmat(norb,*),emat(norb,norb)
      DIMENSION  eigvec(2,npop  ,norb)
c
c Local declarations:
      DATA       lstout / 1 /
c
c >>>> EXECUTABLE CODE:
c
      mat = norb*norb
cdev:
      call FLGETIWR( IWR )
      write(IWR,*) 'DEV/PRHOIJ2-pij03: npop+pull up stripe-arith'
ctp:
c  Get local MP info:
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
      if( lstout.gt.0 .and. nprocs.gt.1 .and. iproc.eq.masterp )then
        call FLGETIWR( IWR )
        write(IWR,*) 'DMAT: k-parallel/idle-master version'
      endif
c
c    Bloch vector loop:
c
      do 1000 kblk=1,nk, nprocs
c       We will start with master(=0) doing first k-pt in set:
        kproc0 = -1
c
        kfirst = kblk
        klast = kblk + nprocs - 1
        if( klast.gt.nk )then
c         More procs than k-pts. Truncate loop
          klast = nk
c         Make masterp=0 an idle proc, so as to hide some i/o
          kproc0 = 0
        endif
c
c       >>>> Retrieve eigenvectors and send to proper worker proc
c
        kproc = kproc0
        kdmat = 0
        do 100 k=kfirst,klast
          kproc = kproc + 1
          if( kproc.eq.nprocs ) kproc = 0
c
          nvecs = 2*norb*npop
          if( iproc.eq.masterp )then
c           Retrieve (only occupied, not all nstate) eigvecs for this k-pt
            call READBIG( ivecfl, nvecs, eigvec )
            if( iproc.eq.kproc )then
c              ... and keep on master to generate dmat(k)
              call DCOPY( nvecs, eigvec,1, dmat,1 )
              kdmat = k
            else
c              ... or send to worker to generate dmat(k)
              lenmsg = nvecs
              call MPSENDR8( kproc, lenmsg, eigvec, icomm )
            endif
          elseif( iproc.eq.kproc )then
c            ... receive eigvec on worker
            lenmsg = nvecs
            call MPRECVR8( masterp, lenmsg, dmat  , icomm )
            kdmat = k
          endif
c
  100   continue
c
        if( kdmat.le.0 ) goto 500
c
c       >>>> Compute density matrix dmat, energy-weighted, emat
c
        if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat assignd' )
c       We are working on dmat(k):
        k = kdmat
c
c       Transpose (only occupied) eigvecs for use below:
        call TREIGVEC( norb,npop  , dmat,eigvec )
        numpop = npop
c
        do  n=1,numpop
          fpop(n) = eigpop(n,k)
C          do  i=1,norb
C            eigvec(1,n,i) = dmat(1,i,n)
C            eigvec(2,n,i) = dmat(2,i,n)
C          enddo
        enddo
c
c       Clear out dmat, emat:
        call MKZERO( mat, dmat )
        call MKZERO( mat, emat )
c
c       Scale fpop by k-point sampling weight for use below.
        wt = wtk(k)
        do  n=1,numpop
          fpop(n) = wt*fpop(n)
        enddo
c
        do 400 j=1,norb
c
c         Do imaginary part first ...
          do 210 i=1,j-1
            dmi = zero
            emi = zero
            do 200 n=1,numpop
              occnk = fpop(n)
              ti = occnk*( eigvec(2,n,i)*eigvec(1,n,j)
     $                   - eigvec(1,n,i)*eigvec(2,n,j) )
              dmi = dmi + ti
              emi = emi + ti*eigval(n,k)
  200       continue
            dmat(i,j) = dmi
            emat(i,j) = emi
  210     continue
c
c          ... and then switch to real part:
          do 310 i=j,norb
            dmr = zero
            emr = zero
            do 300 n=1,numpop
              occnk = fpop(n)
              tr = occnk*( eigvec(1,n,i)*eigvec(1,n,j)
     $                   + eigvec(2,n,i)*eigvec(2,n,j) )
              dmr = dmr + tr
              emr = emr + tr*eigval(n,k)
  300       continue
            dmat(i,j) = dmr
            emat(i,j) = emr
  310     continue
c
c         Close loop over j/norb
  400   continue
        if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat compute' )
c
c       Drop here if no dmat to work on ...
  500   continue
c
c       >>>> Put dmat and emat into storage for use later
c
c       Currently, all sent to masterp and written to disk from there
        kproc = kproc0
        do 600 k=kfirst,klast
          kproc = kproc + 1
          if( kproc.eq.nprocs ) kproc = 0
c
          if( iproc .eq. masterp )then
c           Get dmats and emats from workers, and store
c           Use ticketing to keep from overwhelming master/commun.
c
c           Save dmat to idmatfl for use later
            if( iproc.ne.kproc )then
c             First, send kproc a ticket ... ok to send
              call MPSENDTICKET( kproc, icomm )
              lenmsg = mat
              call MPRECVR8( kproc, lenmsg, dmat, icomm )
            endif
            call WRITBIG( idmatfl, mat, dmat )
c
c           Save emat to iematfl for use in the force section
            if( iproc.ne.kproc )then
c             First, send kproc another ticket ... ok to send
              call MPSENDTICKET( kproc, icomm )
              lenmsg = mat
              call MPRECVR8( kproc, lenmsg, emat, icomm )
            endif
            call WRITBIG( iematfl, mat, emat )
c
            if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat writes' )
          elseif( iproc.eq.kproc )then
c           Send dmat, emat to master for storage
c           First, wait for ticket from master:
            call MPRECVTICKET( masterp, icomm )
c           Send dmat(k) to master:
            call MPSENDR8( masterp, mat, dmat, icomm )
c           Get another ticket from master:
            call MPRECVTICKET( masterp, icomm )
c           Send emat(k) to master:
            call MPSENDR8( masterp, mat, emat, icomm )
            if( lstout.gt.3 ) call TIMER( 'prhoij2/edmat sent' )
          endif
  600   continue
c
c       End loop over blocks of k-points
 1000 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSCHROED
c
c
      subroutine PSCHROED( do_eigvecs,
     $ idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $ norb,nstate,nk, eigval, wksml,
     $ nspin, icluster, gapp,
     $ wmat )
c---------------------------------------------------------------
c Purpose: solve Schroedinger equation given input ham and ovlp,
c          single solve, distributed over processors
c
c Written: Peter A. Schultz,  8-April-1999, for v2.34
c
c Revision history:
c  22Aug17-PAS/2.68: fix idle-proc (solver context) bug
c  10Jan13-PAS/2.64: cosmetics, un-drivered, and integrated
c   2May12-ACP/2.64: Added alternate entry point to support non-SCF portions of code
c   2May12-ACP/2.64: Workaround for 2GB limit in ScaLAPACK
c   9Dec08-PAS/2.62: reduce mem complex-parallel solve
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  18Aug05-PAS/2.59: master-only matrix files, dHam from disk
c  27Oct04-PAS/tp0.1: clean-up task-parallel routines
c   7Apr03-DMC: Parallelized by moving loop on nspin into routine,
c               and using ScaLAPACK to solve eigen problem.
c  24Jul01-PAS/2.49: spin-polarized dft
c  21Jun01-PAS/2.48: replace STOPs
c  17May01-PAS/2.47: extract writes
c  23Aug99-PAS/2.38: install zhegvx expert driver eigensolve
c  10May99-PAS/2.35: eigvals cut from ivecfl output
c   9Apr99-PAS/2.34: complex lapack eigensolver installed
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   5Sep98-PAS/2.25: print out spectral limits
c  14Jul98-PAS/2.23: option to use unchol'd delta-Ham
c  12Aug94-PAS/2.11: cleanup
c---------------------------------------------------------------
c
c Notes:
c  The call parameter idosolv was used to control which version of the
c  LAPACK routine ZHEGV was used.  Here, we always use PZHEGVX.
c  parameter will be used only to control whether the number of
c  eigenvectors/values returned is all (norb) or just nstate.
c  DMC 8Apr03
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
c     MAXBLOCK= 2**31/sizeof(DOUBLE COMPLEX) = 2**31/16.
      PARAMETER ( maxblockz = 67 108 864 )
c
ctp: Pull in the mpi stuff
      INCLUDE  'mpif.h'
c
      LOGICAL    do_eigvecs
c Output:
      DIMENSION  eigval(nstate,nk,nspin)
c Scratch:
      DIMENSION  wksml(*)
      DIMENSION  wmat(*)
c Scratch/parallel:
      INTEGER    icluster(*)
      DIMENSION  gapp(*)
c
c Local declarations:
      CHARACTER  jobz*1
      DATA       iprint / 1 /
ctp: Local for ScaLAPACK
      PARAMETER  (NB = 32)
      INTEGER    CONTEXT
      INTEGER    DESCA(9), DESCWORK(9)
c
      LOGICAL    USE_SCALAPACK_FORMULAS
      DATA       USE_SCALAPACK_FORMULAS / .true. /
      LOGICAL    DEBUG_MEM
      DATA       DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      if( iprint.gt.0 .or. DEBUG_MEM )then
        call FLGETIWR( IWR )
        write(IWR,*) 'PSCHROED: parallel distributed solver'
      endif
c
c     Do eigenvectors with jobz='V' (jobz="N" turns off eigenvectors):
      if( do_eigvecs )then
c       Compute eigenvalues and eigenvectors
        jobz = 'V'
      else
c       Compute eigenvalues only
        jobz = 'N'
      endif
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c  Initialize ScaLAPACK Grid
c
c  Define processor grid dimensions of largest nearly square rectangular array
c  of processors within nprocs. npsolver is the number of processors in the grid
c  (npsolver.le.nprocs).
c
      call EIGPROCS( nprocs, nprow,npcol, npsolver )
      if( iprint.gt.0 ) write(IWR,'(a,4i8)')
     $ ' SOLVER grid:  procs, rowsXcols=solver =',
     $      nprocs,  nprow,npcol, npsolver
c
c  Create BLACS context using system context, 
c  i.e. the local MP communicator
      CONTEXT = icomm
      call BLACS_GRIDINIT( CONTEXT, 'R', nprow, npcol )
      call BLACS_GRIDINFO( CONTEXT, nprow, npcol, myrow, mycol )
c
      if( myrow .ne. -1 )then
c       Initialize matrix descriptions, for procs in solver context
c
c       Set up the descriptor for a distributed matrix
c       First, need to limit lead dimension in local space
c       nblkrow = ceil(norb/NB) = # of blocks in a row
        nblkrow = (norb+NB-1)/NB
c       npblkrow = ceil(nblkrow/nprow)
        npblkrow = (nblkrow+nprow-1)/nprow
        LLDA = npblkrow*NB
        call DESCINIT( DESCA, norb,norb, NB,NB, 0,0,
     $   CONTEXT, LLDA, INFO )
        if( INFO.ne.0 ) call STOPXERR( 'pschroed - DESCINIT/DESCA' )
c
c       Set up the descriptor for full local matrix
c
c       call DESCINIT( DESC, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, INFO)
c       MB=M, NB=N means only the processor (RSRC,CSRC) has any data.
        call DESCINIT( DESCWORK, norb,norb, norb,norb, 0,0,
     $   CONTEXT, norb, INFO )
        if( INFO.ne.0 ) call STOPXERR( 'pschroed - DESCINIT/DESCWORK' )
      else
c       This processor not in solver context, and leaves
        goto 2000
      endif
      if( iprint.gt.2 )then
        write(IWR,*) 'DESC= dtype ctxt m n mb nb csrc rsrc lld'
        write(IWR,'(a,9i8)') '  DESCA=',DESCA
        write(IWR,'(a,9i8)') '  DWORK=',DESCWORK
      endif
c
c  Initialize memory allocations:
c
      mat = norb*norb
      call EIGMDIST( matdist, norb, npsolver )
c
c Allocate "simple" memory sizes for solver work arrays (using crude assumptions).
c Work (Zwork) is complex, r8 length (in wmat) and c16 length (in solver) differ 2x
      lwzwork = mat
      lzwork = lwzwork / 2
C      lrwork = mat
      lrwork = matdist
      liwork = 6*norb
c
      if( DEBUG_MEM )then
        write(IWR,'(a,4i10)') '  PSCH: old native works Z,R,I=',
     $                         lwzwork,lrwork,liwork
c       Use default mpi/blacs-avoiding naive work space:
        call PSCHWKSZ0(            npsolver  , do_eigvecs,
     $   norb, matdist, lwzwork, lzwork, lrwork, liwork )
        lwork = lzwork
        write(IWR,'(a,4i10)') '  PSCH: wksz0  works Z,R,I=',
     $                         lwzwork,lrwork,liwork
      endif
c
      if( USE_SCALAPACK_FORMULAS )then
        if( DEBUG_MEM .or. iprint.gt.2 )
     $  write(IWR,*) '  PSCH: WKSZ/SCAL FORMULAS, CONTEXT=',CONTEXT
c       Get EXACT mpi/blacs-aware work space using this context:
        call PSCHWKSZX( CONTEXT  , npsolver  , do_eigvecs,
     $   norb, matdist, lwzwork, lzwork, lrwork, liwork )
        lwork = lzwork
        if( DEBUG_MEM .or. iprint.gt.2 )
     $  write(IWR,'(a,4i10)') '  PSCH: wkszx  works Z,R,I=',
     $                         lwzwork,lrwork,liwork
      else
c       Use default mpi/blacs-avoiding naive work space:
        call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $   norb, matdist, lwzwork, lzwork, lrwork, liwork )
        lwork = lzwork
      endif
c
      if( .not. do_eigvecs )then
c       No eigenvectors, this becomes much less:
        lzwork = (NB+1) * (norb+1)
        lrwork =9*norb
      endif
c
c The solver memory allocations have the following requirements:
      nwpreph = 2*matdist + 3*mat
      nwpreps = 4*matdist + 2*mat
      nwprep = MAX( nwpreph, nwpreps )
c     Using a naive allocation for solver work arrays
      nwsolv = 6*matdist + lwzwork + lrwork + liwork + norb
C      nwclose = 6*matdist + 2*norb*nstate
c      because vectors are brought back into full norb*norb ...
      nwclose = 6*matdist + 2*norb*norb
      nwpsch = MAX( nwprep, nwsolv, nwclose )
c
      nwleft = nwpsch - nwsolv
      if( DEBUG_MEM )then
        write(IWR,*) '  PSCHROED: spare solver mem=',nwleft
      endif
c
c MEM-Solver: d-C16-H, d-C16-S, d-C16-Vec, d-C16-Zwork, d-R8-Rwork)
c
      iwHmat  = 1
      iwOvlp  = iwHmat + 2*matdist
      iwVecs  = iwOvlp + 2*matdist
      iwZwork = iwVecs + 2*matdist
      iwRwork = iwZwork + lwzwork
      iwIwork = iwRwork + lrwork
      iwiwk1  = iwIwork + liwork
      iwlast  = iwiwk1  + norb
      iwEnd   = iwlast  - 1
c
      if( nwsolv .ne. iwEnd )then
c       Our memory arithmetic has failed this idiot-check
        write(*,*) '***** ERROR/PSCHROED: nwsolv,iwsolv=',nwsolv,iwEnd
        call STOPXERR( 'Code ERROR: nwsolv mem arith in pschroed' )
      endif
c
c MEM-Prep: (d-C16-H, d-C16-S=full-R8-mat, full-C16-mat)
      iwMat1 = iwOvlp
c     Put full-C16-mat after full-R8-mat AND after dist-S
      iwMat2 = MAX( iwMat1 + mat, iwVecs )
c MEM-After: (d-C16-H, d-C16-S, d-C16-Vec, full-C16-Vec)
c     Put full output vectors after dist-vectors:
      iwVeco = iwZwork
c
c Set files to correct locations
c
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idhamfl )
      endif
c
      do 100 ispin=1,nspin
c
        if( iproc .eq. masterp )then
          REWIND( unit=iham0fl )
          REWIND( unit=iovlpfl )
        endif
c
c       Bloch vector loop:
c
        do 1000 k=1,nk
c
          if( idosolv .eq. 1 )then
            jstate = norb
          elseif( idosolv .eq. 2 )then
            jstate = nstate
          else
            call STOPXERR( 'pschroed - invalid solver option' )
          endif
c
c         Form H0+dH=full Hamiltonian(spin,k), and distribute
c
          if( iproc .eq. masterp )then
c           Retrieve reference Hamiltonian H0 (packed) from iham0fl
            call READBIG( iham0fl, mat, wmat(iwMat2) )
c           Retrieve delta-Hamiltonian dH (packed) from idhamfl
            call READBIG( idhamfl, mat, wmat(iwMat1) )
c           Combine H0 and dH to get full Hamiltonian (packed)
            do  ij=1,mat
              wmat(iwMat1-1+ij) = wmat(iwMat1-1+ij) + wmat(iwMat2-1+ij)
            enddo
c           Load Ham (packed real) into lower triangle of complex matrix
            call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
c-->                      n    real*8 nxn    complex*16 nxn (lower triangular)
          endif
c         Distribute complex Ham over processors
          if( norb*norb .le. maxblockz )then
c           The matrix is small enough that pzgemr2d can swallow it whole.    
            call PZGEMR2D( norb,norb,
c                From full matrix on masterp ...
     $           wmat(iwMat2), 1,1, DESCWORK,
c                 ... into distributed matrix:
     $           wmat(iwHmat), 1,1, DESCA,  CONTEXT )
          else
c           Distribute the matrix in (column) stripes that are smaller than 2GB each.
c           Should be good as long the distributed matrix is no larger than 2GB on any processor.
            nblock = maxblockz / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 0
  151       continue
            if( noffset.le.norb )then
              call PZGEMR2D( norb, nblock,
     $             wmat(iwMat2+moffset), 1,1, DESCWORK,
     $             wmat(iwHmat), 1, noffset, DESCA,  CONTEXT )
              noffset = noffset + nblock
              moffset = moffset + 2*norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 151
            endif
          endif
c
c         Retrieve and distribute overlap matrix
c
          if( iproc .eq. masterp )then
c           Retrieve overlap matrix (packed) from iovlpfl
            call READBIG( iovlpfl , mat, wmat(iwMat1) )
c           Load overlap (packed) into lower triangle complex matrix
            call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
          endif
c         Distribute complex overlap over processors
          if( norb*norb .le. maxblockz )then
c           The matrix is small enough that pzgemr2d can swallow it whole.    
            call PZGEMR2D( norb,norb,
     $           wmat(iwMat2), 1,1, DESCWORK,
     $           wmat(iwOvlp), 1,1, DESCA,  CONTEXT )
          else
c           Distribute the matrix in (column) stripes that are smaller than 2GB each.
c           Should be good as long the distributed matrix is no larger than 2GB on any processor.
            nblock = maxblockz / norb
            nleft = norb - nblock
            noffset = 1
            moffset = 0
  152       continue
            if( noffset.le.norb )then
              call PZGEMR2D( norb, nblock,
     $             wmat(iwMat2+moffset), 1,1, DESCWORK,
     $             wmat(iwOvlp), 1, noffset, DESCA,  CONTEXT )
              noffset = noffset + nblock
              moffset = moffset + 2*norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 152
            endif
          endif
c
c         Ham is in wmat(iwHmat), overlap in wmat(iwOvlp)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Solve Schroedinger equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          rmone = 2*PDLAMCH( CONTEXT, 'S' )
c
          call PZHEGVX( 1,jobz, 'I', 'L',
     $     norb, wmat(iwHmat),1,1,DESCA, wmat(iwOvlp),1,1,DESCA,
c -->            H(2,n,n)-i              S(2,n,n)-i
     $     -100d0,-20d0, 1,jstate, rmone, nstout,nvout,
     $     wksml,
c <--      eigval(1,k,ispin),
     $     rmone, wmat(iwVecs),1,1,DESCA,
     $     wmat(iwZwork),lzwork, wmat(iwRwork),lrwork,
     $     wmat(iwIwork),liwork, wmat(iwiwk1),
     $     icluster,gapp, INFO )
c  -->     int(6n)               int(n)
c
c         Check if something went wrong in parallel solver
          if( INFO .ne. 0 )then
            call FLGETIERR( IERRFL )
            write(IERRFL,*) 'PSCHROED/PZHEGVX error, INFO=',INFO
            call STOPXERR( 'PSCHROED/PZHEGVX solver error' )
          endif
c
c         Put eigenvalues into the right place
          call DCOPY( nstate, wksml,1, eigval(1,k,ispin),1 )
c
c         If no eigvecs, we are done with this Ham:
          if( .not. do_eigvecs ) goto 1000
c
c         Collect (only nstate!) eigenvectors onto master
          if( norb*nstate .le. maxblockz )then
c           The matrix is small enough to retrieve whole.
            call PZGEMR2D( norb, nstate,
     $           wmat(iwVecs), 1,1, DESCA,
     $           wmat(iwVeco), 1,1, DESCWORK,  CONTEXT )
          else
c           Retrieve the matrix in (column) stripes that are smaller than 2GB each.
c           Should be good as long the distributed matrix is no larger than 2GB on any processor.
            nblock = maxblockz / norb
            nleft = nstate - nblock
            noffset = 1
            moffset = 0
 1502       if( noffset.le.nstate )then
              call PZGEMR2D( norb, nblock,
     $             wmat(iwVecs+moffset), 1,1, DESCA,
     $             wmat(iwVeco), 1, noffset, DESCWORK,  CONTEXT )
              noffset = noffset + nblock
              moffset = moffset + 2*norb*nblock
              if( nleft.ge.nblock )then
c               Have enough to do another full block.
                nleft = nleft - nblock
              else
c               Have less than a full block to send.
                nblock = nleft
                nleft = 0
              endif
              goto 1502
            endif
          endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c              Record eigenvalues and eigenvectors
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c         Save eigenvectors to ivecfl for later
          if( do_eigvecs .and. ( iproc .eq. masterp ) )then
            call WRITBIG( ivecfl, 2*norb*nstate, wmat(iwVeco) )
          endif
c
 1000   continue
  100 continue
c
c Non-solver procs return here:
c
 2000 continue
c
c If some processors not used in the ScaLAPACK solver context ...
      if( nprocs .ne. npsolver )then
c        ... also need to broadcast eigvals to those not having them
        call MPI_Bcast( eigval(1,1,1), nstate*nk*nspin,
     $   MPI_DOUBLE_PRECISION,
     $   masterp, icomm, ierr )
        if( ierr.ne.0 ) call STOPXERR( 'eigsolv - Bcast eigval' )
      endif
c
      if( myrow .ne. -1 ) call BLACS_GRIDEXIT( CONTEXT )
      call MPBARRIER( icomm )
c      call MPI_ABORT(icomm, ierrors, istatus)
c
c    That's all Folks!
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
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSCHROED2
c
c
      subroutine PSCHROED2( do_eigvecs,
     $ idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $ norb,nstate,nk, eigval, rwork, seval, scut,
     $ vmat, wmat, iwork )
c---------------------------------------------------------------
c Purpose: solve Schroedinger equation given input ham and ovlp,
c          with k-resolved distribution over processors
c          using LAPACK on each processor.
c
c Written: Peter A. Schultz,  8-April-1999, for v2.34
c
c Revision history:
c   2May12-ACP/2.64: Added flag to skip eigvecs to support non-SCF
c  25Jun07-PAS/2.60: merge tp and serial
c   6Jul06-APT/2.60: merged task parallel and image parallel
c  18Aug05-PAS/2.59: master-only matrix files; read dHams
c  25Sep04-PAS/tp0.1: clean-up task-parallel routines
c  22Apr03-DMC: This version of a parallel schroed routine is designed
c               to take advantage of nk>1, by parallelizing on k.
c  24Jul01-PAS/2.49: spin-polarized dft
c  21Jun01-PAS/2.48: replace STOPs
c  17May01-PAS/2.47: extract writes
c  23Aug99-PAS/2.38: install zhegvx expert driver eigensolve
c  10May99-PAS/2.35: eigvals cut from ivecfl output
c   9Apr99-PAS/2.34: complex lapack eigensolver installed
c  15Feb99-PAS/2.31: routine to intercept perhaps too-big i/o
c  10Dec98-PAS/2.29: promote source to explicit double precision
c   5Sep98-PAS/2.25: print out spectral limits
c  14Jul98-PAS/2.23: option to use unchol'd delta-Ham
c  12Aug94-PAS/2.11: cleanup
c---------------------------------------------------------------
c
c When this routine is done
c   (1) it writes out eigenvectors to ivecfl
c   (2) it returns eigenvalues in "eigval()"
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
ctp: Pull in the mpi stuff:
      INCLUDE  'mpif.h'
c MPI information
      INTEGER status(MPI_STATUS_SIZE)
c
c Input:
      LOGICAL    do_eigvecs
c Output: eigenvalue spectrum
      DIMENSION  eigval(nstate,nk)
c Scratch arrays:
      DIMENSION  rwork(*)
      DIMENSION  seval(*)
      DIMENSION  iwork(norb,*)
      DIMENSION  vmat(norb*norb)
      DIMENSION  wmat(norb*norb,*)
c
c Local declarations:
      CHARACTER  jobz*1,uplo*1
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
cxxx: declare dev-version of pschroedkp:
      if( IWR.ge.0 ) write(IWR,*)
     $ 'PSCHROED2: k-parallel single-processor solver, nk=',nk
      mat = norb*norb
c
c Do eigenvectors with jobz='V' (jobz="N" turns off eigenvectors):
      if( do_eigvecs )then
c       Compute eigenvalues and eigenvectors
        jobz = 'V'
      else
c       Compute eigenvalues only
        jobz = 'N'
      endif
      uplo = 'L'
c
c  Get Local MP info
c
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
      if( iproc .eq. masterp )then
        REWIND( unit=iham0fl )
        REWIND( unit=iovlpfl )
      endif
c
c    Bloch vector loop:
c
      do 1000 kblk=1,nk, nprocs
        do 1010 kk=1,nprocs
          k = kk + kblk - 1
          if( k .gt. nk ) goto 1011
c
c         This sets it up so that masterp goes *last*
          kproc = MOD( kk, nprocs )
c
          if( iproc .eq. masterp )then
c
c           Retrieve delta-Hamiltonian from idhamfl
            call READBIG( idhamfl, mat, vmat )
c
c           Retrieve reference Hamiltonian from iham0fl
            call READBIG( iham0fl, mat, wmat(1,1) )
c
c           Add del-Ham to ref-Ham to get total Ham
            do  ij=1,mat
              wmat(ij,1) = wmat(ij,1) + vmat(ij)
            enddo
c
c           Retrieve overlap matrix from iovlpfl
            call READBIG( iovlpfl , mat, wmat(1,2) )
c
            if( kproc .ne. masterp )then
c             Send H and S together to node that will do this work
              lenmsg = 2*mat
              call MPSENDR8( kproc, lenmsg, wmat(1,1), icomm )
C              call MPI_SEND( wmat(1,1), 2*mat, MPI_DOUBLE_PRECISION,
C     $         kproc,k, icomm, ierr )
C              if( ierr.ne.0 ) call STOPXERR( 'pschroed2 - Send err' )
            endif
c
          elseif( iproc .eq. kproc )then
c
            lenmsg = 2*mat
            call MPRECVR8( masterp, lenmsg, wmat(1,1), icomm )
C            call MPI_RECV( wmat(1,1), lenmsg,  MPI_DOUBLE_PRECISION,
C     $       masterp,k, icomm, status, ierr )
C              if( ierr.ne.0 ) call STOPXERR( 'pschroed2 - Send err' )
c
          endif
c
c         Only the node(s) with work will continue
          if( iproc .ne. kproc ) goto 1010
c
c         Load overlap into lower triangle full complex matrix wmat(3:4)
          call C16LOW( norb, wmat(1,2), wmat(1,3) )
c
c         Load Ham into lower triangle of full complex matrix wmat(1:2)
          call DCOPY( mat, wmat(1,1),1, vmat,1 )
          call C16LOW( norb, vmat     , wmat(1,1) )
c
c         Need 3*norb-2 for r*8 workspace "rwork"
c         lwork is (complex) length for workspace "l"
          lwork = (mat) / 2
c
c         Ham is in wmat(,1:2), overlap in wmat(,3:4)
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                  Solve Schroedinger equation
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
C         We will always use the expert driver for parallel.
c         We use idosolv only to control how many eigenvalues/vectors
c         to solve for.
c
          if( idosolv .eq. 1 )then
            jstate = norb
          elseif( idosolv .eq. 2 )then
            jstate = nstate
          else
            call STOPXERR( 'pschroed2 - invalid solver option' )
          endif
c
c          ***** Selected eigenvalue Hermitian eigenproblem *****
c
          call ZHEGVX( 1,jobz, 'I', 'L',
     $     norb, wmat(1,1),norb,         wmat(1,3),norb,
c -->            H(2,n,n)-i              S(2,n,n)-i
     $     -100d0,-20d0, 1,jstate, 2*DLAMCH('S'), nstout,
     $     eigval(1,k),
     $     wmat(1,5),norb,    vmat     ,lwork, rwork,
c -->      eigvec(2,n,nstate) work(2,lwork)    (7n)
     $     iwork(1,2),iwork, INFO )
c -->      int(5n)    int(n)
c
          if( INFO .ne. 0 )then
            call STOPXERR( 'pschroed2 - diagonalization error' )
          endif
c
          if( do_eigvecs )then
c           Put eigenvectors in the right place ...
            call DCOPY( 2*norb*nstate, wmat(1,5),1, wmat(1,1),1 )
          endif
c
 1010   continue
 1011   continue
c
c       Now collect results of work at master, and store from there
c
        do 1020 kk=1,nprocs
          k = kk + kblk - 1
          if( k .gt. nk ) goto 1021
c
          kproc = MOD( kk, nprocs )
c
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c          Send eigenvalues to everyone from kproc that owns this k-pt
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          call MPBCAST8( kproc, nstate, eigval(1,k), icomm )
c
c         Without eigvecs, skip their retrieval and storage
          if( .not. do_eigvecs ) goto 1020
c
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                          Record eigenvectors
c         * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
          nvec = 2*norb*nstate
c
c         Transfer the eigenvectors to masterp for storage
          if( iproc.eq.kproc .and. iproc.ne.masterp )then
c           Send eigvecs to masterp (safe - Bcast above shields)
            call MPSENDR8( masterp, nvec, wmat(1,1), icomm )
c
          elseif( iproc.eq.masterp .and. kproc.ne.masterp )then
c           Receive eigvecs from kproc (safe - Bcast above shields)
            call MPRECVR8( kproc, nvec, wmat(1,1), icomm )
c
          elseif( iproc.eq.masterp .and. kproc.eq.masterp )then
c           Receives stomped on masterp eigvecs, must refresh
c           (because masterp does last k-pt, not first)
            call DCOPY( nvec, wmat(1,5),1, wmat(1,1),1 )
          endif
c
c         Save eigenfunction coefficients to ivecfl
          if( iproc .eq. masterp )then
            call WRITBIG( ivecfl, 2*norb*nstate, wmat(1,1) )
          endif
c
 1020   continue
 1021   continue
c
 1000 continue
C      call MPBARRIER( icomm )
C      call MPI_ABORT( icomm, ierrors, istatus )
c
c    That's all Folks!
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
c
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PEIGSOLVKP
c
c
      subroutine PEIGSOLVKP( masterlist,
     $ icluster, gapp,
     $ do_eigvecs,do_eigpops, idosolv,scut, nspin,
     $ idhamfl,iham0fl,iovlpfl,
     $ ivecfl, idmatfl,idmatsfl, iematfl,iematsfl,
     $ icloseocc,dokdefct,doblas3, ikdefct,nstbulk,
     $ etemp,edegen, elecno, efermi,egap, esumsp,
     $ norb,nstate, nk, veck,wtk, eigval,eigpop,npop, wksml,
     $ wmat )
c---------------------------------------------------------------
c Purpose: k-parallel-parallel solves for complex density matrices,
c          hamiltonian and overlap matrices of current iteration
c
c Written: Andrew C. Pineda, 16-May-2012, for 2.64
c          Adapted (forked) from  peigsolv
c
c Revision history:
c   4Dec13-PAS/2.64: cosmetics, un-drivered, integrated
c---------------------------------------------------------------
c
c  On input, idhamfl contains full input delta-hamiltonian.
c  To obtain full hamiltonian, delta-H needs to be combined with
c  H0 (computed in setup and saved in iham0fl; overlap in iovlpfl),
c  done in routine "schroed", which also solves the Schroedinger
c  eqn. Routine "pop" populates the resulting levels and finds the
c  fermi energy, then "rhoij" takes the eigenfcns and occupations
c  and builds and saves to disk the density matrix for each k-vector.
c     
c  When this routine is done:
c   (1) eigenvectors are saved to "ivecfl" (if do_eigvecs)
c   (2) density matrix is saved to "idmat" (and if do_eigpops)
c   (3) energy-weighted density matrix is saved to "iemat"
c   (4) eigenvalue spectrum is returned in "eigval()" (not written)
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c
      PARAMETER  ( zero = 0.d0 )
c
c Input:
      LOGICAL    do_eigvecs, do_eigpops
      LOGICAL    dokdefct
      LOGICAL    doblas3
      DIMENSION  veck(3,nk),wtk(nk)
      DIMENSION  elecno(2)
      DIMENSION  nstbulk(2)
c Output:
      DIMENSION  efermi(2),egap(2)
      DIMENSION  eigval(nstate*nk,nspin),eigpop(nstate*nk,nspin)
c Scratch:
      DIMENSION  wksml(norb,*)
      DIMENSION  wmat(*)
c Scratch/parallel:
      INTEGER    icluster(*)
c ACP/space for list of kp-masters used in pschroedkp
      INTEGER    masterlist(*)
      DIMENSION  gapp(*)
      DATA       itimlog / 5 /
c
c Local:
      DIMENSION  occlvl(2)
      DATA       occlvl / 2.d0,1.d0 /
c
c >>>> EXECUTABLE CODE:
c
cdev
      call FLGETIWR( IWR )
      if( IWR.gt.0 )then
        write(IWR,*) 'PEIGSOLVKP: k-parallel-parallel solver'
        write(IWR,*) 'KPP/DEV: smaller prhoij(06) stripe/transpose'
      endif
c  Solve Schroedinger equation for the various k-vectors
c
      mat = norb*norb
      nmat = nk*mat
c
c Partition space within wmat() for schroed:
      iw1 = 1
      iww = iw1 + mat
      nwsch = norb*(4*norb + 2*nstate)
      iwx = iww + nwsch
c
      call MPNODES( nproc )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      if( iproc .eq. masterp )then
c       Go to start of idhamfl, to get ready for dHam reads
        REWIND( unit=idhamfl )
c       Go to start of ivecfl, to get ready for eigvec writes
        REWIND( unit=ivecfl )
      endif
c
      if( itimlog.gt.0 ) call TIMER( 'eigsolvkp/begin      ' )
c
c      if( nproc.eq.1 .or. nk.gt.1 )then
      if( nproc.eq.1) then
c$$$      if( nproc.eq.1 .or. nk.gt.4 .or. norb.lt.800 .or.
c$$$     $    nproc.lt.(4*nk) )then
c       Parallel over k-points - single-proc solves
c
        do  ispin=1,nspin
c
          call PSCHROED2( do_eigvecs,
     $     idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $     norb,nstate,nk, eigval(1,ispin), wksml,    wksml(1,9),
c -->                                       rwork-7ns seval(n)-s
     $     scut,
     $     wmat(iw1), wmat(iww), wmat(iwx) )
c -->      vmat-s     wmat-6ks   iwork-6ns
c
          if( itimlog.gt.1 ) call TIMER( 'eigsolvkp/psch2 k-par' )
        enddo
c
      else
c       Parallel solvers
c
        call PSCHROEDKP( do_eigvecs,
     $   idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $   norb,nstate,nk, eigval, wksml,
c -->                            wksml-ns
     $   nspin, icluster, gapp,
     $   wmat, masterlist )
c -->    wmat-s
c
        if( itimlog.gt.1 ) call TIMER( 'eigsolvkp/pschkp scalap' )
      endif
c
      if( .not. do_eigpops )then
c       Skip populations and density matrix construction
        RETURN
      endif
c
c  Evaluate density matrix in diagonal (eigfcn) representation
c
      nstk = nstate*nk
c
      npop = 0
      esumsp = zero
      occfull = occlvl(nspin)
c
      do  ispin=1,nspin
        elnumbr = elecno(ispin)
c
        if( ( icloseocc .gt. 0 ) .or. dokdefct )then
c         Use closed-shell occupations, or defect sampling
          call POPCLOSE( icloseocc,dokdefct,ikdefct,
     $     efermi(ispin),egap(ispin), nstbulk(ispin),
     $     etemp,edegen, elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),numpop,
     $     wmat )
c -->      wlvl-s(nstate)
c
        else
c         Use standard metallic occupations
          mspin = 1
          call POP( mspin, popnumbr,
     $     efermi(ispin),egap(ispin),
     $     etemp,        elnumbr,esum1e, occfull,
     $     nstate,nk, wtk, eigval(1,ispin),eigpop(1,ispin),numpop,
     $     wmat(1), wmat(1+nstk),wmat(1+2*nstk),
     $     wmat(1+3*nstk), wmat(1+4*nstk) )
c -->      elvl wlvl klvl nlvl flvl - all scratch (nk*state)
        endif
c
        npop = MAX( npop,numpop )
        esumsp = esumsp + esum1e
      enddo
c
      if( itimlog.gt.0 ) call TIMER( 'eigsolvkp/pop        ' )
c
      if( .not. do_eigvecs )then
c       Without eigvecs, skip density matrix construction
        RETURN
      endif
c
c  Evaluate density matrix in orbital basis
c
      mat = norb*norb
c
      call MPNODES( nproc )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idmatfl )
        REWIND( unit=iematfl )
        if( nspin.eq. 2 )then
          REWIND( unit=idmatsfl )
          REWIND( unit=iematsfl )
        endif
      endif
      iw1 = 1
      iw2 = iw1 + mat
      iw3 = iw2 + mat
c
      idmfile = idmatfl
      iemfile = iematfl
      do  ispin=1,nspin
c
        if( nk .ge. nproc )then
c         Then distribute full k-points over processors
c
          call PRHOIJ2( ivecfl, idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s(n,nst)
c ------> 2*norb*nstate  mat     2*norb*nstate
c
        else
c         Distribute every k-point over all processors
c
c         Set up sizes of stripe-optimized vectors
          call RHOSTRIP( jstripe, norb, nproc )
          iw2 = iw1 + jstripe*norb
          iw3 = MAX( iw1 + 2*norb*npop, iw2 + jstripe*norb )
c
          call PRHOIJ( ivecfl,idmfile,iemfile,
     $     norb,nstate,nk,
     $     eigval(1,ispin),eigpop(1,ispin),wksml,npop, wtk,
     $     wmat(iw1), wmat(iw2), wmat(iw3) )
c -->      dmat-(2)s  emat-s     eigvec-2s
c ------>  2*norb*npop       ->  2*norb*npop
c ------>  stripe     stripe <-  2*norb*npop
c ------>  stripe     mat
c
        endif
c       Switch to spin-dn:
        idmfile = idmatsfl
        iemfile = iematsfl
c
      enddo
      if( itimlog.gt.0 ) call TIMER( 'eigsolvkp/rhoij      ' )
c
c    That's all Folks!
c
      RETURN
      END

c
c
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> KPMASTERLIST
c
c
      subroutine KPMASTERLIST( nprocs, iproc, master, icomm,
     $     nprocs_k, iproc_k, master_k, icomm_k,
     $     masterlist, imasters )
c
c  Build a list of kpmasters and compute the number of kpmasters.
c
c  This is needed for the PSCHROEDKP k-parallel solver.
c
      INTEGER  nprocs, iproc, master, icomm
      INTEGER  nprocs_k, iproc_k, master_k, icomm_k
      INTEGER  masterlist, imasters
      INTEGER  jproc
      INTEGER  istate
      DIMENSION masterlist(nprocs)
c    
      imasters = 0
      lenmsg = 1
      do  jproc=0,nprocs-1
        if( iproc .eq. jproc )then
          if( iproc_k .eq. 0 )then
            istate = 0
          else
            istate = 1
          endif
          if( iproc .ne. master )then
            call MPSENDI( master, lenmsg, istate, icomm )
          endif
        endif
        if( iproc .eq. master ) then
          if( jproc .ne. master )then
            call MPRECVI( jproc,  lenmsg, istate, icomm )
          endif
          if( istate.eq.0 )then
            imasters = imasters + 1
            masterlist(imasters) = jproc
          endif
        endif
      enddo
c
c Distribute the list of kp-masters to all processes.      
c
      lenmsg = 1
      call MPBCASTI( master, lenmsg, imasters, icomm )
      lenmsg = imasters
      call MPBCASTI( master, lenmsg, masterlist, icomm )
      do jproc=imasters+1,nprocs, 1
        masterlist(jproc) = -1
      enddo
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
c     
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PSCHROEDKP
c
c
      subroutine PSCHROEDKP( do_eigvecs,
     $     idhamfl,iham0fl,iovlpfl,ivecfl, idosolv,
     $     norb,nstate,nk, eigval, wksml,
     $     nspin, icluster, gapp,
     $     wmat, masterlist )
c---------------------------------------------------------------
c Purpose: solve Schroedinger equation given input ham and ovlp,
c     fully k-parallel-parallel solve. 
c     
c Written: Andrew C. Pineda,  9-February-2010 for v2.64,
c          (adapted from "pschroed" routine)
c     
c Revision history:
c  16May12-ACP/2.64: Removed ENTRY point.
c  2May12-ACP/2.64: Added alternate entry point to support non-SCF portions of code
c  17Apr12-ACP/2.64: Work around ScaLAPACK 2GB limit
c---------------------------------------------------------------
c     
c  Notes:
c     Each k-point can have multiple processes associated with it.
c     Is so, SCALAPACK is used instead of LAPACK.
c     The master process, with access to all the data on disk, 
c     distributes ham and ovlp for eachk-point to a k-point master.
c     Each k-master leads a parallel solvewithin its MPI context,
c     BLACS context derived from the k-point communicator). 
c
c     The call parameter idosolv was used to control which version of the
c     LAPACK routine ZHEGV was used.  Here, we always use PZHEGVX.
c     parameter will be used only to control whether the number of
c     eigenvectors/values returned is all (norb) or just nstate.
c     Note: this means that SCF treatment of linear dependence
c       (idsolv=3) is not possible with this routine as it stands.
c     DMC 8Apr03
c
      IMPLICIT DOUBLE PRECISION  (a-h,o-z)
c     MAXBLOCK= 2**31/sizeof(DOUBLE COMPLEX) = 2**31/16.
      PARAMETER  ( maxblockz = 67 108 864 )
c
ctp: Pull in the mpi stuff
      INCLUDE  'mpif.h'
c
c Input:
      LOGICAL    do_eigvecs
c Output:
      DIMENSION  eigval(nstate,nk,nspin)
c Scratch:
      DIMENSION  wksml(*)
      DIMENSION  wmat(*)
c Scratch/parallel:
      INTEGER    icluster(*)
      DIMENSION  gapp(*)
      INTEGER    masterlist(*)
c
c Local declarations:
      CHARACTER  jobz*1
ctp: Local for ScaLAPACK
      PARAMETER  (NB = 32)
      INTEGER    ANB, SQNPC
      INTEGER    CONTEXT_K
      INTEGER    DESCA(9), DESCWORK(9)
      INTEGER    procindex
c
      LOGICAL    USE_SCALAPACK_FORMULAS
      LOGICAL    OPTIMAL_EIGENVECTORS
      DATA  USE_SCALAPACK_FORMULAS / .true. /
      DATA  OPTIMAL_EIGENVECTORS   / .true. /
      LOGICAL    DEBUG_MEM
      DATA       DEBUG_MEM / .false. /
c
c >>>> EXECUTABLE CODE:
c
      call FLGETIWR( IWR )
cxxx: declare dev-version of pschroedkp:
      if( IWR.ge.0 ) write(IWR,*)
     $ 'PSCHROEDKP: k-parallel-parellel solves (opt mem)'
c
      mat = norb*norb
c
c     Get Local MP info
c
c     Image level
      call MPNODES( nprocs )
      call MPNODE( iproc )
      call MPNODE0( masterp )
      call MPCOMM( icomm )
c
c     K-parallel level
      call MPNODES_K( nprocs_k )
      call MPNODE_K( iproc_k )
      call MPNODE0_K( master_k )
      call MPCOMM_K( icomm_k )
c
c     Build a list of the ranks of the kpmasters within the icomm
c     communicator and get a count of kpmasters. The returned list of
c     ranks is in increasing order.
c
      call KPMASTERLIST( nprocs, iproc, masterp, icomm,
     $     nprocs_k, iproc_k, master_k, icomm_k,
     $     masterlist, num_masters )
c
      if( IWR.ge.0 )then
        write(IWR,9131) 'number of kpp masterss=',num_masters
         if( DEBUG_MEM ) write(IWR,9132)
     $     ' kpp masterlist ranks=',(masterlist(i),i=1, num_masters)
      endif
 9131 format( 4x,a,1x,i4)
 9132 format( 4x,a,(16i5))
c
c     Broadcast rank in icomm communicator of kp-master to processes in
c     icomm_k and store in my_kpmaster
c
      lenmsg = 1
      my_kpmaster = iproc
      call MPBCASTI( master_k, lenmsg, my_kpmaster, icomm_k )
c
      if( nprocs_k.eq.1 )then
c       SINGLE PROCESSOR SOLVE - Use LAPACK
c
        matdist = mat
C        lwzwork  = 2*MAX( 2*norb, mat )
        lwzwork  = MAX( 4*norb, mat )
        lzwork  = lwzwork / 2
        lwork  = lzwork
        lrwork  = 7*norb
        liwork  = 5*norb
c
c       Cloase single-processor solve work sizes
      else
c       PRALLEL SOLVER - Use SCALAPACK
c
c      Initialize ScaLAPACK Grid for k-point to be handled by this processor.
c      Define processor grid dimensions of largest nearly square
c       rectangular array of processors within icomm_k.
c       npsolver_k is the number of processors in the grid (npsolver_k.le.nprocs_k).
c
        call EIGPROCS( nprocs_k, nprow_k, npcol_k, npsolver_k )
        if( IWR.gt.0 ) write(IWR,'(a,4i6)') ' procs, rowsxcols=solver:',
     $                nprocs_k,  nprow_k,npcol_k, npsolver_k
c
c     Create BLACS context using system context from 
c     k-point communicator.
c
        CONTEXT_K = icomm_k
        call BLACS_GRIDINIT( CONTEXT_K, 'R', nprow_k, npcol_k )
        call BLACS_GRIDINFO( CONTEXT_K,
     $       nprow_k, npcol_k,
     $       myrow_k, mycol_k )
c
c     myrow_k = -1 if iproc > npsolver_k?
c
c     Initialize matrix descriptions, for procs in solver contexts
c
        if( myrow_k .ne. -1 )then
c         Process in solver context, set up matrix descriptors
c
c         First, set up descriptor for distributed matrix
c
c         Need to limit lead dimension in local space
c          nblkrow = ceiling(DBLE(norb)/DBLE(NB)) = # of blocks in a row
c          npblkrow = ceiling(DBLE(nblkrow)/DBLE(nprow_k))
c          Lead dimension LLDA is row bloacks times blocking factor
          nblkrow = (norb+NB-1) / NB
          npblkrow = (nblkrow+nprow_k-1) / nprow_k
          LLDA = npblkrow*NB
          call DESCINIT( DESCA,  norb,norb, NB,NB,  0, 0,
     $         CONTEXT_K, LLDA, INFO )
          if( INFO.ne.0 )
     $      call STOPXERR( 'pschroedkp - DESCINIT/DESCA' )
c
c        Set up the descriptor for full local matrix.
c
c        call DESCINIT( DESC, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, INFO)
c        MB=M, NB=N means only the processor (RSRC,CSRC)=(0,0) has any data.
c
          call DESCINIT( DESCWORK,  norb,norb, norb,norb,  0, 0,
     $         CONTEXT_K, norb, INFO )
          if( INFO.ne.0 )
     $      call STOPXERR(' pschroedkp - DESCINIT/DESCWORK' )
        else
c         This processor is not in solver context, and leaves
          goto 990
        endif
c
        if( USE_SCALAPACK_FORMULAS )then
c         Get EXACT mpi/blacs-aware work space using this context:
          call PSCHWKSZX( CONTEXT_K, npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
        else
c         Use default mpi/blacs-avoiding naive work space:
          call PSCHWKSZ0(            npsolver_k, do_eigvecs,
     $     norb, matdist, lwzwork, lzwork, lrwork, liwork )
          lwork = lzwork
        endif
c
        if( DEBUG_MEM ) write(IWR,'(a,3i12)')'DEV/SCHKPP: works Z,R,I=',
     $                                        lwzwork,lrwork,liwork
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                Do the memory allocations for parallel solves
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     MEM-Solver: d-C16-H, d-C16-S, d-C16-Vec, d-C16-lwork, d-R8-rwork)
c
c     The solver memory allocations have the following requirements:
        nwpreph = 2*matdist + 3*mat
        nwpreps = 4*matdist + 2*mat
        nwprep  = MAX( nwpreph, nwpreps )
        nwsolv  = 6*matdist + lwzwork + lrwork + liwork + norb
c     nwclose = 6*matdist + 2*norb*nstate
c     because vectors are brought back into full norb*norb ...
        nwclose = 6*matdist + 2*norb*norb
        nwpsch  = MAX( nwprep, nwsolv, nwclose )
c
        nwleft = nwpsch - nwsolv
        if( DEBUG_MEM )then
          write(IWR,*) 'PSCHROEDKP: spare solver mem=',nwleft
        endif
c
c     end of  if (nprocs_k.gt.1) then
      endif
c
      iwHmat  = 1
      iwOvlp  = iwHmat + 2*matdist
      iwVecs  = iwOvlp + 2*matdist
      iwZwork = iwVecs + 2*matdist
c
      iwRwork = iwZwork + lwzwork
      iwIwork = iwRwork + lrwork
      iwiwk1  = iwIwork + liwork
      iwlast  = iwwk1   + norb
      iwEnd   = iwlast  - 1
c
c     MEM-Prep: (d-C16-H, d-C16-S=full-R8-mat, full-C16-mat)
c
      iwMat1 = iwOvlp
c
c     Put full-C16-mat after full-R8-mat AND after dist-S
c
      iwMat2 = MAX( iwMat1 + mat, iwVecs )
c
c     MEM-After: (d-C16-H, d-C16-S, d-C16-Vec, full-C16-Vec)
c     Put full output vectors after dist-vectors:
c
      iwVeco = iwZwork
c
c     Set files to correct locations
c
      if( iproc .eq. masterp )then
        REWIND( unit=ivecfl )
        REWIND( unit=idhamfl )
      endif
c
      do 100 ispin=1,nspin
c
        if( iproc .eq. masterp )then
          REWIND( unit=iham0fl )
          REWIND( unit=iovlpfl )
        endif
c
        if( do_eigvecs )then
c         Compute eigenvalues and eigenvectors
          jobz = 'V'
        else
c         Compute eigenvalues only
          jobz = 'N'
        endif
c
c     Bloch vector loop: individual k point problems are distributed to
C     blocks of processes controled by the k-point masters.
c
        do 1000 kblk = 1, nk, num_masters
          my_k = -1
          do 1010 kk = 1, num_masters
            k = kk + kblk - 1  
c     k out of range, everyone skips to end of loop if they are not the image master
c     which must remain to do I/O.
            if( k .gt. nk )then
              if( (my_k.eq.-1) .and. (iproc.eq.masterp) )then
                my_k = 0   
              endif
              goto 1011
            endif
c     Since masterlist is in increasing order, we ensure that masterp does the last k-point
            procindex = 1 + MOD( kk, num_masters )
c     Select the k-point master for this k.
            kmasterpr = masterlist(procindex)
            if( my_kpmaster .eq. kmasterpr )then
c             Set parameters for the group of processes working on this k-point.
              my_k = k
              if( idosolv .eq. 1 )then
                jstate = norb
              elseif( idosolv .eq. 2 )then
                jstate = nstate
              else
                call STOPXERR( 'PSCHROEDKP: invalid solver option' )
              endif
            endif
c
c     Form H0+dH=full Hamiltonian(spin,k), and distribute it to the k-point master
c
            if( iproc .eq. masterp )then
c             Retrieve reference Hamiltonian H0 (packed) from iham0fl
              call READBIG( iham0fl, mat, wmat(iwMat2) )
              if( iproc.ne.kmasterpr )then
                lenmsg = mat
                call MPSENDR8( kmasterpr,lenmsg,wmat(iwMat2),icomm )
              endif
c             Retrieve delta-Hamiltonian dH (packed) from idhamfl
              call READBIG( idhamfl, mat, wmat(iwMat1) )
              if( iproc.ne.kmasterpr )then
                lenmsg = mat
                call MPSENDR8( kmasterpr,lenmsg,wmat(iwMat1),icomm )
              endif
            endif
c
c     Receive full Hamiltonian on the k-point master at wmat(iwMat2).
            
            if( iproc .eq. kmasterpr )then
              if( iproc .ne. masterp )then
                lenmsg = mat
                call MPRECVR8( masterp,lenmsg,wmat(iwMat2),icomm )
                lenmsg = mat
                call MPRECVR8( masterp,lenmsg,wmat(iwMat1),icomm )
              endif
c             Combine H0 and dH to get full Hamiltonian (packed)
              do  ij=1,mat
                wmat(iwMat1-1+ij) =
     $               wmat(iwMat1-1+ij) + wmat(iwMat2-1+ij)
              end do
c             Load Ham (packed real) into lower triangle of complex matrix
              call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
c----------------------------> n   real*8 nxn     complex*16 nxn
c---------------------------->                  (lower triangular)
c
c             Close if(iproc.eq.kmasterpr) 
            endif
c
c     Distribute complex Ham over the subset of processes belonging to CONTEXT_K using SCALAPACK.
c     If I have done this right this is the group of processes headed up by kpmasterpr.
c
c     iwMat2 --> full matrix on kmasterpr
c     iwHmat --> distributed matrix within icomm_k
c
            if( my_kpmaster .eq. kmasterpr )then
              if( nprocs_k.eq.1 )then
c               Copy iwMat2 to iwHmat
                call DCOPY( 2*mat, wmat(iwMat2), 1,
     $               wmat(iwHmat), 1)
              else
c     Distribute the matrix if iproc is part of the correct CONTEXT_K/icomm_k.
c     All processes in this CONTEXT_K/icomm_k must participate in this call.
                if( norb*norb.le.maxblockz )then
c                 Matrix small enough that pzgemr2d can swallow it whole.    
                  call PZGEMR2D( norb,norb,
     $                 wmat(iwMat2), 1,1, DESCWORK,
     $                 wmat(iwHmat), 1,1, DESCA,
     $                 CONTEXT_K )
                else
c                 Distribute the matrix in (column) stripes smaller than 2GB each.
c                 Good if distributed matrix is no larger than 2GB on any processor.
                  nblock = maxblockz / norb
                  nleft = norb - nblock
                  noffset = 1
                  moffset = 1
  151             continue
                  if( noffset.le.norb )then
                    call PZGEMR2D( norb,nblock,
     $                   wmat(iwMat2-1+moffset), 1, 1, DESCWORK,
     $                   wmat(iwHmat), 1, noffset, DESCA,
     $                   CONTEXT_K )
                    noffset = noffset + nblock
                    moffset = moffset + 2*norb*nblock
                    if( nleft.ge.nblock )then
c                     Have enough to do another full block.
                      nleft = nleft - nblock
                    else
c                     Have less than a full block to send.
                      nblock = nleft
                      nleft = 0
                    endif
                    goto 151
                  endif
                endif
              endif
            endif
c
c     Retrieve overlap matrix on masterp and send it to k-point master
c
            if( iproc .eq. masterp )then
c             Retrieve overlap matrix (packed) from iovlpfl
              call READBIG( iovlpfl , mat, wmat(iwMat1) )
              if( iproc .ne. kmasterpr )then
                lenmsg = mat
                call MPSENDR8( kmasterpr,lenmsg,wmat(iwMat1),icomm )
              endif
            endif
c
c     Receive complex overlap on k-point master
c
            if( iproc .eq. kmasterpr )then
              if( iproc .ne. masterp )then
                lenmsg = mat
                call MPRECVR8( masterp,lenmsg,wmat(iwMat1),icomm )
              endif
c             Load overlap (packed) into lower triangle complex matrix
              call C16LOW( norb, wmat(iwMat1), wmat(iwMat2) )
c----------------------------n   real*8 nxn    complex*16 nxn
c----------------------------(lower triangular)
            endif
c
c     Distribute complex overlap over processes controlled by k-point master
c
            if( my_kpmaster .eq. kmasterpr )then
              if( nprocs_k.eq.1 )then
c               Copy complex n X n matrix into place.
                call DCOPY( 2*mat, wmat(iwMat2), 1,
     $               wmat(iwOvlp), 1 )
              else
c     Distribute the matrix if iproc is part of the correct CONTEXT_K/icomm_k.
c     All processes in this CONTEXT_K/icomm_k must participate in this call.
                if( norb*norb .le. maxblockz )then
c                 Matrix is small enough that pzgemr2d can swallow it whole.                      
                  call PZGEMR2D( norb,norb,
     $                 wmat(iwMat2), 1,1, DESCWORK,
     $                 wmat(iwOvlp), 1,1, DESCA,
     $                 CONTEXT_K )
                else
c                 Distribute matrix in (column) stripes smaller than 2GB each.
c                 Good if distributed matrix no larger than 2GB on any processor.
                  nblock = maxblockz / norb
                  nleft = norb - nblock
                  noffset = 1
                  moffset = 1
  152             continue
                  if( noffset.le.norb )then
c     We handle the offset to the full matrix here so that ScaLAPACK doesn't have to.
c     NB: This foreshadows what we want to do to distribute the matrix in stripes.
                    call PZGEMR2D( norb,nblock,
     $                   wmat(iwMat2-1+moffset), 1,1, DESCWORK,
     $                   wmat(iwOvlp), 1,noffset, DESCA,
     $                   CONTEXT_K )
                    noffset = noffset + nblock
                    moffset = moffset + 2*norb*nblock
                    if( nleft.ge.nblock )then
c                     Have enough to do another full block.
                      nleft = nleft - nblock
                    else
c                     Have less than a full block to send.
                      nblock = nleft
                      nleft = 0
                    endif
                    goto 152
                  endif
                endif
              endif
            endif
c     
c     End of 1st kk loop (Eigenproblem distribution)
c     
 1010     continue
 1011   continue
c
c     Ham is in wmat(iwHmat), overlap in wmat(iwOvlp)
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Solve Schroedinger equation
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Processors in a context do the solve for their assigned k-point
c
c     If I'm part of a CONTEXT_K that is processing a k-point (my_k>=0),
c     call the SCALAPACK solver. Note that there are up to num_masters
c     parallel solves happening here as there are num_masters CONTEXT_Ks. 
c
        if( my_k .ge. 1 )then
c
          if( nprocs_k.eq.1 )then
c           Only one processor assigned to this k-point, use LAPACK
c
            rmone = 2*DLAMCH('S')
c
            call ZHEGVX( 1,   jobz, 'I',  'L', norb,
     $       wmat(iwHmat), norb,  wmat(iwOvlp),norb,
c -->        A             LDA    B            LDB
     $       -100d0, -20d0,    1, jstate,  rmone,
c -->        VL      VU c      IL IU       ABSTOL
     $       nstout,   wksml,
c -->        M         W[N]
     $       wmat(iwVecs), norb,
c -->        Z             LDZ  [Z(LDZ, max(1,M))]
     $       wmat(iwZwork), lwork,
c -->        WORK           LWORK    [WORK(LWORK)]
     $       wmat(iwRwork),  wmat(iwIwork),  wmat(iwiwk1),
c -->        RWORK[7*N]      IWORK[5*N]      IFAIL[N]
     $       INFO )
c
c           If something failed, go do postmortem
            if( INFO .ne. 0 ) goto 1310
c
c           Close serial LAPACK solve
          else
c           Parallel Scalapack solver
c
            if( myrow_k .ne. -1 ) then
c
              rmone = 2*PDLAMCH( CONTEXT_K, 'S' )
c
              call PZHEGVX( 1, jobz, 'I', 'L', norb,
     $             wmat(iwHmat), 1,1,DESCA, 
     $             wmat(iwOvlp), 1,1,DESCA,
     $             -100d0,-20d0, 1,jstate, rmone,
     $             nstout, nvout, wksml, rmone, 
c <--                             eigval(1,k,ispin),
     $             wmat(iwVecs), 1,1,DESCA,
     $             wmat(iwZwork), lwork,
     $             wmat(iwRwork), lrwork,
     $             wmat(iwIwork), liwork,
     $             wmat(iwiwk1),
     $             icluster,  gapp, INFO )
c
c             If something failed, go do postmortem
              if( INFO .ne. 0 ) goto 1320
c
c             Close parallel solver branch
            endif
          endif
        endif
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                      Record eigenvalues and eigenvectors
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Retrieve results from kpmasters
        if( my_k .ge. 0 )then
          do 2010 kk = 1, num_masters
            k = kk + kblk - 1
c           If k out of range, everyone skips to end of loop
            if( k .gt. nk )then
              goto 2011
            endif
            procindex = 1 + MOD( kk,num_masters )
            kmasterpr = masterlist(procindex)
c
c     Put eigenvalues into the right place on masterp.
c
            if( iproc .eq. kmasterpr )then
              if( iproc .eq. masterp )then
                call DCOPY( nstate, wksml, 1,
     $               eigval(1,k,ispin), 1 )
              else
                lenmsg = nstate
                call MPSENDR8( masterp,lenmsg,wksml,icomm )
              endif
            endif
            if( iproc .eq. masterp )then
              if( iproc .ne. kmasterpr )then
                lenmsg = nstate
                call MPRECVR8( kmasterpr,lenmsg,
     $               eigval(1,k,ispin),icomm )
              endif
            endif
c
            if( do_eigvecs )then
c     Optionally retrieve eigenvectors for current k-point from wmat(iwVecs) and
c     put them in wmat(iwVeco) on the k-point master (kmasterpr).
c     We only retrieve the vectors for nstate states (length norb) we need.
c     There should be at least nstate states in wmat(iwVecs), more if
c     jstate=norb.
              if( my_kpmaster .eq. kmasterpr )then
                if( nprocs_k.eq.1 )then
                  call DCOPY( 2*norb*nstate, wmat(iwVecs), 1,
     $                 wmat(iwVeco), 1 )
                else
c                 Matrix small enough that PZGEMR2D won't choke reading it all at once.
                  if( norb*nstate.le.maxblockz )then
                    call PZGEMR2D( norb, nstate,
     $                   wmat(iwVecs), 1, 1, DESCA,
     $                   wmat(iwVeco), 1, 1, DESCWORK,
     $                   CONTEXT_K )
                  else
c                   Retrieve matrix in (column) stripes smaller than 2GB each.
c                   Good if distributed matrix is no larger than 2GB on any processor.
                    nblock = maxblockz / norb
                    nleft = nstate - nblock
                    noffset = 1
                    moffset = 1
  153               continue
                    if( noffset.le.nstate )then
                      call PZGEMR2D( norb, nblock,
     $                     wmat(iwVecs), 1, noffset, DESCA,
     $                     wmat(iwVeco-1+moffset), 1, 1, DESCWORK,
     $                     CONTEXT_K )
                      noffset = noffset + nblock
                      moffset = moffset + 2*norb*nblock
                      if( nleft.ge.nblock )then
c                       Have enough to do another full block.
                        nleft = nleft - nblock
                      else
c                       Have less than a full block to send.
                        nblock = nleft
                        nleft = 0
                      endif
                      goto 153
                    endif
                  endif
                endif
              endif
c
c     Send eigenvectors from kpmaster for this k-point to masterp
c     if kpmaster is not the same as masterp.
c
              if( iproc .eq. kmasterpr )then
                if( iproc .ne. masterp )then
c                 Send only the nstate (complex) vectors we need.
c                 lenmsg is 2*norb*nstate since the wmat data is complex.
                  lenmsg = 2*norb*nstate
                  call MPSENDR8( masterp,lenmsg,wmat(iwVeco),icomm )
                endif
              endif
c
c             On masterp, save the eigenvectors to ivecfl.
c
              if( iproc .eq. masterp )then
c               Receive nstate vectors we need from kmasterpr if needed. 
c               The factor of 2 is because the data is complex.
                if( iproc .ne. kmasterpr )then
                  lenmsg = 2*norb*nstate
                  call MPRECVR8( kmasterpr,lenmsg,wmat(iwVeco),
     $                 icomm )
                endif
c               wmat(iwVeco) on masterp is filled either by PZGEMR2D or MPRECVR8
c               nstate <= jstate so we should always have data.
                isize = 2*norb*nstate
                call WRITBIG( ivecfl, isize, wmat(iwVeco) )
              endif
            endif
c
c           End of 2nd kk loop
 2010     continue
 2011     continue
        endif
c       Close kblk loop
 1000 continue
c     Close ispin loop
  100 continue
c
c Non-solver processors re-enter here:
  990 continue
c
c  Need to broadcast eigvals to everyone
c
      lenmsg = nstate*nk*nspin
      call MPBCAST8( masterp, lenmsg, eigval(1,1,1), icomm )
c
c  Everyone exits their BLACS context
c
      if( nprocs_k.gt.1 )then
        if( myrow_k .ne. -1 ) call BLACS_GRIDEXIT( CONTEXT_K )
      endif
c
      call MPBARRIER( icomm )
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                           NORMAL EXIT
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      RETURN
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                     SOLVER ERROR POSTMORTEMS
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c           ***** Lapack/ZHEGVX error postmortem *****
 1310 continue
      call FLGETIERR( IERRFL )
      write(IERRFL,*) '***** ERROR: PSCHROEDKP/ZHEGVX, INFO=',INFO
      if( INFO .lt. 0 )then
        write(IERRFL,*) -INFO,'th argument of ZHEGVX had  illegal value'
        if( INFO .eq. -20 ) write(IERRFL,*) 'lwork = ', lwork
      else
        if( INFO .le. norb )then
          write(IERRFL,*)  INFO,' eigenvectors failed to converge'
        else
          write(IERRFL,*)
     $     'factorization of overlap matrix could not be completed '//
     $     'no eigenvalues or eigenvectors computed'
        endif
      endif
      call FLUSH( IERRFL )
c
      call STOPXERR( 'PSCHROEDKP/ZHEGVX solver error' )
      goto 1399
c
c           ***** Scalapack/PZHEGVX error postmortem *****
 1320 continue
      call FLGETIERR( IERRFL )
      write(IERRFL,*) '***** ERROR: PSCHROEDKP/PZHEGVX, INFO=',INFO
c
      if(MOD( INFO,2 ) .ne. 0 )then
        write(IERRFL,*) 'One or more eigenvectors failed to converge'
        call print_ifail_array( IERRFL, wmat(iwiwk1), norb )
      endif
      if( MOD( INFO/2 ,2 ) .ne. 0 )then
        write(IERRFL,*) 'One or more eigenvectors '//
     $   'could not be reorthogonalized due to insufficient workspace'
        write(IERRFL,*) 'Here is the list of clusters:'
c       The cluster list is a zero terminated list
        ikl = 1
        write(IERRFL,*) 'Cluster #    Starts    Ends   Size    Gap'
 1321   write(IERRFL,*)  ikl, icluster(2*ikl-1),icluster(2*ikl),
     $    icluster(2*ikl)-icluster(2*ikl-1)+1, gapp(ikl)
        if( icluster(2*ikl+1).ne.0 )then
          ikl = ikl + 1
          goto 1321
        endif
        write(IERRFL,*)'Here are the workspace sizes:'
        write(IERRFL,*)'lwork,lrwork,liwork=',lwork,lrwork,liwork
      endif
      if( MOD( INFO/4, 2 ) .ne. 0 )
     $ write(IERRFL,*) 'Space limit prevented PZHEGVX from ' //
     $  'computing all eigenvectors between VL(-100) and VU(-20).'
      if( MOD( INFO/8, 2 ) .ne. 0 )
     $ write(IERRFL,*) 'Scalapack PZSTEBZ failed to compute eigenvalues'
      if( mod( INFO/16, 2 ) .ne. 0 )then
        write(IERRFL,*) 'Overlap matrix is not positive definite'
        write(IERRFL,*) 'norb=',norb
        call print_ifail_array( IERRFL, wmat(iwiwk1), 1 )
      endif
      if( INFO.eq.-28 )
     $  write(IERRFL,*) 'LWORK= ', lwork,', but needs=',wmat(iwZwork)
      if( INFO.eq.-30 )
     $  write(IERRFL,*) 'LRWORK=',lrwork,', but needs=',wmat(iwRwork)
      if( INFO.eq.-32 )
     $  write(IERRFL,*) 'LRWORK=',liwork,', but needs=',wmat(iwIwork)
c
      call STOPXERR( 'PSCHROEDKP/PZHEGVX solver error' )
c
 1399 continue
c
c    That's all Folks!
c
      RETURN
      END
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> print_ifail_array
c
      subroutine print_ifail_array(iunit, iarray, n)
      INTEGER  iunit
      INTEGER  n
      INTEGER  iarray(n)
      INTEGER  j
      do  j=1,n
        write(iunit,*)'ifail(',j,')=', iarray(j)
      enddo
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
