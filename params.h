! ####&
! namelen is the length of filenames (in characters)
! maxlay is the maximum number of layers allowed in the model
! maxtr is the maximum number of traces allowed
! maxseg: maximum # of segments (should be 3*maxlay for 1st-order multiples
! maxph: maximum number of phases per trace
! buffsize is the max. line length assumed for reading files.
! mk is the maximum number of layers handled by mineos
      integer namelen, maxlay,maxtr,maxseg,maxph
      integer buffsize,milay,malay,ndatadmax
      integer mk
      parameter (namelen=40, milay = 5, malay=80)
      parameter (maxlay=300, maxtr=13, maxseg=45)
      parameter (maxph=40000,buffsize=120)
      parameter (ndatadmax=100)
      parameter (mk=350)

      integer invertype
      parameter (invertype=1)! 1:Only P 2: P&S


! Units for reading and writing
      integer iounit1,iounit2
      parameter (iounit1=1,iounit2=2)

! pi: duh. ztol: tolerance considered to be equiv. to zero.
      real pi,ztol
      parameter (pi=3.141592653589793,ztol=1.e-7)

! nsamp is the number of samples per trace.
      integer maxsamp
      parameter (maxsamp=5000)

! P anisotropy ratio (w.r.t S)
      real S2P_aniso
      parameter (S2P_aniso=1.)

      real widthh, shifttp,shiftts
      parameter (widthh=0.1,shifttp=20,shiftts=20)

! Rotation to output: 0 is EW/Z/NS, 1 is T/Z/R, 2 is SV/SH/P
      integer out_rott
      parameter (out_rott=1)

! Thickness minimum for one layer
      real thickmin
      parameter (thickmin=2)


! For De-homo
       integer cut
       real dH_min, DH_max
       parameter  (cut = 4)
        parameter (dH_min = 0)
       parameter (dH_max = 350)

! Scaling Parameters
       real vpvs
       parameter  (vpvs = 1.73)

! Earth radius in meters
       real rearth
       parameter  (rearth = 6371000)

! maximum distance between mode and period
       real maxdis
       parameter (maxdis=30.)

! max acceptable raylquo
        real maxrq
        parameter (maxrq=1000.)
