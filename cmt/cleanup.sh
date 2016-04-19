if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/cern.ch/sw/contrib/CMT/v1r20p20090520; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=Rich1MDCS -version=v1r0 -path=/afs/cern.ch/user/a/anandi/cmtuser/Brunel_v49r2p1/Rich $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

