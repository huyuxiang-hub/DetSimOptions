# echo "cleanup DetSimOptions v0 in /junofs/users/huyuxiang/juno_centos7/offline/Simulation/DetSimV2"

if test "${CMTROOT}" = ""; then
  CMTROOT=/cvmfs/juno.ihep.ac.cn/centos7_amd64_gcc830/Pre-Release/J20v2r0-branch/ExternalLibs/CMT/v1r26; export CMTROOT
fi
. ${CMTROOT}/mgr/setup.sh
cmtDetSimOptionstempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if test ! $? = 0 ; then cmtDetSimOptionstempfile=/tmp/cmt.$$; fi
${CMTROOT}/mgr/cmt cleanup -sh -pack=DetSimOptions -version=v0 -path=/junofs/users/huyuxiang/juno_centos7/offline/Simulation/DetSimV2  $* >${cmtDetSimOptionstempfile}
if test $? != 0 ; then
  echo >&2 "${CMTROOT}/mgr/cmt cleanup -sh -pack=DetSimOptions -version=v0 -path=/junofs/users/huyuxiang/juno_centos7/offline/Simulation/DetSimV2  $* >${cmtDetSimOptionstempfile}"
  cmtcleanupstatus=2
  /bin/rm -f ${cmtDetSimOptionstempfile}
  unset cmtDetSimOptionstempfile
  return $cmtcleanupstatus
fi
cmtcleanupstatus=0
. ${cmtDetSimOptionstempfile}
if test $? != 0 ; then
  cmtcleanupstatus=2
fi
/bin/rm -f ${cmtDetSimOptionstempfile}
unset cmtDetSimOptionstempfile
return $cmtcleanupstatus

