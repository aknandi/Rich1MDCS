# echo "Setting Rich1MDCS v1r0 in /afs/cern.ch/user/a/anandi/cmtuser/Brunel_v49r2p1/Rich"

if test "${CMTROOT}" = ""; then
  CMTROOT=/afs/cern.ch/sw/contrib/CMT/v1r20p20090520; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh

tempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then tempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt setup -sh -pack=Rich1MDCS -version=v1r0 -path=/afs/cern.ch/user/a/anandi/cmtuser/Brunel_v49r2p1/Rich  -no_cleanup $* >${tempfile}; . ${tempfile}
/bin/rm -f ${tempfile}

